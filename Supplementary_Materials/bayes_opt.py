# Import libraries
import numpy as np
import scipy as sp
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("darkgrid"); sns.set_context("talk")
import GPy
import GPyOpt   # This will do the Bayesian optimization
import os
import sdf_helper as sh
import shutil
import time
import subprocess
import sys

# Functions
def loss_func(coeffs,parameters,PICEnergy,Delta=1):
    # List Constants
    q = 1.6e-19 #Elementary charge (C)
    Z = 1           #Atomic number
    mi = 1.67e-27   #Proton mass (kg)
    lamb = 1e-6   #Wavelength of light (m)
    n_crit = 1e21/((lamb/1e-6)**2)  #if lambda is in um, this gives answer in cm^-3
    epsilon0 = 8.854e-12 #In F/m
    me = 9.11*10**-31   #Electorn mass (kg)
    c = 3*10**8     #Speed of light (m/s)]
    JtoMeV = 6.242e12   #Converts J to MeV

    tau = parameters[0]*(10**-15)
    intensity = 10**parameters[1]
    targ_dens = parameters[2]*n_crit*1e6    #In m^-3
    D = parameters[3]*10**-6
    Lg = parameters[4]*10**-6
    tsci = (10**int(np.format_float_scientific(tau)[-3:]))

    multiplier = [1,tsci/(10**int(np.format_float_scientific(intensity)[-3:])),tsci/(10**int(np.format_float_scientific(targ_dens)[-3:])),tsci/(10**int(np.format_float_scientific(D)[-3:])),tsci/(10**int(np.format_float_scientific(Lg)[-3:]))]

    acc_time = coeffs[:,0]*tau*multiplier[0] + coeffs[:,1]*intensity*multiplier[1] + coeffs[:,2]*targ_dens*multiplier[2] + coeffs[:,3]*D*multiplier[3] + coeffs[:,4]*Lg*multiplier[4]

    TeHot = JtoMeV*me*c**2 * (np.sqrt(1+(intensity*(lamb/1e-6)**2)/(1.37e18))-1)
    wp = np.sqrt((Z * q**2 * targ_dens)/(mi*epsilon0))  #Ion plasma frequency in SI
    ta = np.float64((wp*acc_time)/(np.sqrt(2*np.e)))   #Ion acceleration time    

    Mora_E = 2*TeHot*(np.log(ta + np.sqrt(ta**2 + 1)))**2
    
    if np.abs(PICEnergy - Mora_E) <= Delta:
        return 0.5 * (PICEnergy - Mora_E)**2
    else:
        return Delta * (np.abs(PICEnergy - Mora_E) - 0.5 * Delta)
    
def make_dir(params):
    new_dir = f'EPOCH_SIM_t{params[0]}_I{params[1]}_p{params[2]}_D{params[3]}_Lg{params[4]}'
    os.system(f'mkdir {new_dir}')
    return new_dir

def write_job(title):
    f = open('job.sh', 'w')
    f.write(
       f"""#!/bin/sh
# File: Job.sh
# Description: Template for running a job with slurm
# Submit job from command line: sbatch <path_to_Job.sh>
# Check the status of your jobs: squeue -u $USER
# Cancel a job: scancel <PID>  (PID is found with squeue)

# ----------------------------------------------------------------------------
# SLURM COMMAND       DESCRIPTION
# ----------------------------------------------------------------------------
# --job-name          Job name
# --time              Walltime (max time limit) in hours:minutes:seconds
# --nodes             Number of nodes to run on
# --nodelist          Which nodes to run on (if irrelevant, omit this line)
# --ntasks-per-node   Number of processes to use per node
# --mem-per-cpu       Minimum memory required per allocated CPU - [K|M|G|T]
# --mail-type         Types of events that send an email notification
# --mail-user         The email to send notifications to

#SBATCH --job-name={title}
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tamminga.2@osu.edu

# ----------------------------------------------------------------------------
# MY PARAMETER        DESCRIPTION
# ----------------------------------------------------------------------------
# NPROC               Number of processors to use per node (--ntasks-per-node)
# EXEC                Path to the executable
# INPUT               Path to the input file (if N/A, set to empty string)
# GETTOP              Path to the gettop script
# USETOP              Whether or not to run gettop (monitors top metrics)
# TOPFREQ             How often to sample top metrics (in seconds)

NPROC=10
EXEC="$HOME/epoch-4.18.0-2022-12-14/epoch1d/bin/epoch1d"
INPUT="input.deck"

# load modules
module load intel
module load openmpi/4.0.1

# output simple commands
set -x

# change to submission directory
cd $SLURM_SUBMIT_DIR

# output hostname and nodes to files
echo $SLURM_SUBMIT_HOST > hostname.txt 
echo $SLURM_JOB_NODELIST > nodesfile.txt 

# run simulation -- the "--mca..openib" portion disables infiniband in
# order to suppress a warning in the output
echo $SLURM_SUBMIT_DIR | mpirun --mca btl ^openib -np $NPROC $EXEC $INPUT
""")
    f.close()
    return

def write_deck(params):
    f = open('input.deck', 'w')
    f.write(
        f"""begin:constant
  lambda0 = 1 * micron
  tau = {params[0]} * femto
  tau0 = tau / (2*sqrt(loge(2)))
  omega = 2 * pi * c / lambda0
  n_crit = critical(omega)
  thickness = {params[3]} * micron
  Lg = {params[4]} * micron
  l_intensity = 1.0e{params[1]}
  
  #This laser parameter setup follows example_decks/laser_focus.deck 
  # These two set the beam focus
  #w_0 = 15 * micron # 15 micron FWHM	
  #x_spot = 16 * micron # Distance from x_min to spot

  # These are the parameters calculated for driving the laser
  # These should not need to be modified
  #x_R = pi * w_0^2 / lambda0 # Rayleigh range
  #RC = x_spot * (1.0 + (x_R/x_spot)^2) # Radius of curvature on x_min
  #w_bnd = w_0 * sqrt( 1.0 + (x_spot/x_R)^2) # Spot size at x_min
  #gouy = atan(x_spot/x_R) # Gouy phase shift at x_min

  ###
  
  den_targ = {params[2]} * n_crit #epoch is in meters cubed
  den_min = 0
  den_max = thickness
  preplasma_bound = den_min + Lg * log10(0.001)
end:constant

begin:control
  # Final time of simulation
  t_end = 5 * pico

  # Size of domain
  x_min = -100 * micron
  x_max = 150 * micron

  # Number of steps
  nx = (x_max - x_min) / (0.01 * lambda0)

  stdout_frequency = 10 
  smooth_currents = T
  field_ionisation = T #tells the system we are tracking ionisation
  physics_table_location = /home/tamminga.2/epoch-4.18.0-2022-12-14/epoch1d/src/physics_packages/TABLES
end:control

begin:boundaries
  bc_x_min = simple_laser
  bc_x_max = open
end:boundaries

begin:species
  name = deuteron
  charge = 1.0
  mass = 3670.483 #i subtracted 1 electron
  density = if((x gt preplasma_bound) and (x lt den_min),den_targ*exp((x-den_min)/Lg),0)
  density = if((x gt den_min) and (x lt den_max),den_targ,density(deuteron))
  denisty = if((x eq den_min),den_targ,density(deuteron))
  temp_ev = 100
  npart_per_cell = 100
  ionisation_electron_species = (electron)
end:species

begin:species
  name = electron
  charge = -1.0
  mass = 1.0
  temp_ev = 100
  density = density(deuteron)
  npart_per_cell = 100
end:species

begin:laser
  boundary = x_min
  intensity_w_cm2 = l_intensity
  lambda = lambda0
  t_end = 4*tau
  t_profile = gauss(time, 2*tau0, tau0)
  phase = pi / 2
end:laser

begin:output_global
  force_final_to_be_restartable = T
end:output_global

begin:output
  name = Bayes
  dt_snapshot = 20 * femto

  particles = always 
  particle_weight = always
  # Properties on grid
  grid = always
  ex = always
  ey = always
  ez = always
  bx = never
  by = never
  bz = never
  jx = never
  jy = never
  jz = never
  px = always
  py = always
  pz = always
  ekbar = never
  charge_density = never
  number_density = always + species
  temperature = always + species
  total_energy_sum = always + species
  mass = always
  charge = always
  particle_energy = always
end:output""")
    f.close()
    return

def max_energy():
    dir = 'SDF'
    files = sorted(os.listdir(dir))
    fig = plt.figure()
    JtoMeV = 6.242e18/1e6
    max_E=10 #This is in MeV
    bins=max_E*10 #10 bins per MeV, could change
    data = sh.getdata(dir + '/' + files[-1], verbose=False)
    Energy_p =  data.Particles_Ek_deuteron.data*JtoMeV
    w_p = data.Particles_Weight_deuteron.data
    ax2=fig.add_subplot(122)
    H_p = ax2.hist(Energy_p,weights=w_p, bins=bins, range=(0,max_E))
    x_bins = H_p[1][0:-1]
    E = x_bins[np.where(H_p[0])[0].max()]
    return E

def Mora_Energy(coeffs,parameters):
    # List Constants
    q = 1.6e-19 #Elementary charge (C)
    Z = 1           #Atomic number
    mi = 1.67e-27   #Proton mass (kg)
    lamb = 1e-6   #Wavelength of light (m)
    n_crit = 1e21/((lamb/1e-6)**2)  #if lambda is in um, this gives answer in cm^-3
    epsilon0 = 8.854e-12 #In F/m
    me = 9.11*10**-31   #Electorn mass (kg)
    c = 3*10**8     #Speed of light (m/s)]
    JtoMeV = 6.242e12   #Converts J to MeV

    tau = parameters[0]*(10**-15)
    intensity = 10**parameters[1]
    targ_dens = parameters[2]*n_crit*1e6    #In m^-3
    D = parameters[3]*10**-6
    Lg = parameters[4]*10**-6
    tsci = (10**int(np.format_float_scientific(tau)[-3:]))

    multiplier = [1,tsci/(10**int(np.format_float_scientific(intensity)[-3:])),tsci/(10**int(np.format_float_scientific(targ_dens)[-3:])),tsci/(10**int(np.format_float_scientific(D)[-3:])),tsci/(10**int(np.format_float_scientific(Lg)[-3:]))]

    acc_time = coeffs[:,0]*tau*multiplier[0] + coeffs[:,1]*intensity*multiplier[1] + coeffs[:,2]*targ_dens*multiplier[2] + coeffs[:,3]*D*multiplier[3] + coeffs[:,4]*Lg*multiplier[4]

    TeHot = JtoMeV*me*c**2 * (np.sqrt(1+(intensity*(lamb/1e-6)**2)/(1.37e18))-1)
    wp = np.sqrt((Z * q**2 * targ_dens)/(mi*epsilon0))  #Ion plasma frequency in SI
    ta = np.float64((wp*acc_time)/(np.sqrt(2*np.e)))   #Ion acceleration time    

    Mora_E = 2*TeHot*(np.log(ta + np.sqrt(ta**2 + 1)))**2

    return Mora_E
    
def TestFinish():
    result = subprocess.run('/home/tamminga.2/running.py', capture_output=True, text=True)
    length = len(str(result.stdout))
    return length

# Parameters to be varied are:
#   - tau: beam temporal FWHM (fs)
#   - l_intensity: laser intensity I0 (W cm^-2)
#   - den_targ: target number density (nc)
#   - thickness: target thickness (um)
#   - Lg: preplasma scale length (um)

# We are looking at the output of Ei, or ion energy

# Variable array order is tau, l_intensity (exponent only), den_targ, thickness, Lg - picked from Djordjevic paper
par_min = [20,17,80,5,0]
par_max = [500,21,120,25,10]
mora_coeffs = [0,1.3,0,0,0,0,0,0,0] # Standard coeffs. in the original Mora model
mora_min = -10
mora_max = 10
base_len = TestFinish()   #Sets the base length of the query response to if our code is still running. Should be 85


# Set up Bayesian parameters
bounds = [{'name':'tau','type':'continuous','domain':(mora_min,mora_max)},
          {'name':'I0','type':'continuous','domain':(mora_min,mora_max)},
          {'name':'p','type':'continuous','domain':(mora_min,mora_max)},
          {'name':'D','type':'continuous','domain':(mora_min,mora_max)},
          {'name':'Lg','type':'continuous','domain':(mora_min,mora_max)}]
my_acquisition_type = 'EI'
my_model_type = 'GP'

# Loop through multiple pic simulations to prove robustness
for i in range(10):
  print(f'Starting {i}')
  # Randomly pick our starting values
  tau = int(np.random.uniform(low=par_min[0], high=par_max[0]))
  l_intensity = int(np.random.uniform(low=par_min[1], high=par_max[1]))
  den_targ = 100 #This was not varied because the Djordjevic paper mentioned a lack of dependence on this parameter
  #int(np.random.uniform(low=par_min[2], high=par_max[2]))
  thickness = int(np.random.uniform(low=par_min[3], high=par_max[3]))
  Lg = int(np.random.uniform(low=par_min[4], high=par_max[4]))
  params = [tau,l_intensity,den_targ,thickness,Lg]

  # Make new directory and write our input deck and job.sh
  folder = make_dir(params)
  path = os.getcwd()
  dir = path + '/' + folder
  os.chdir(dir)
  write_job(folder)
  write_deck(params)

  # Run PIC Code
  os.system('sbatch job.sh')

  #ENTER WAIT FUNCTION
  done = True
  while(done):
      curr_len = TestFinish()
      if curr_len == base_len:   #The program is done
          done = False
      else:
          time.sleep(15)

  # Get PIC Values
  shutil.copy('/home/tamminga.2/file_organizer.py', os.getcwd())
  os.system('python3 file_organizer.py')
  PIC_Energy = max_energy()
  plt.close()

  # Define function for Bayesian Optimization
  def BoptFunc(coeffs):
      return loss_func(coeffs,parameters=params,PICEnergy=PIC_Energy,Delta=1)
    
  myBopt = GPyOpt.methods.BayesianOptimization(f=BoptFunc,
                                               domain=bounds,
                                               acquisition_type=my_acquisition_type,
                                               exact_feval=True,
                                               model_type=my_model_type)
  # Enter the loop that will
  max_iter=1
  max_time=100
  eps=1.e-8

  num_iter=50
  for i in range(num_iter):
      myBopt.run_optimization(max_iter,max_time,eps)

  reshapecoeff = np.reshape(myBopt.x_opt, newshape=(1,len(myBopt.x_opt)))
  MoraE = Mora_Energy(reshapecoeff,params)

  f = open('RESULTS.txt', 'w')
  f.write('Coefficients: ' + str(myBopt.x_opt) + '\r') #Save Parameters
  f.write('Loss: ' + str(myBopt.fx_opt) + '\r') #Save loss
  f.write('Mora Energy: ' + str(MoraE) + '\r') #Save Mora Energy
  f.write('PIC Energy: ' + str(PIC_Energy) + '\r') #Save PIC Energy
  f.close()

  # Change to original directory
  os.chdir(path)
  print(f'Done with iteration {i}')