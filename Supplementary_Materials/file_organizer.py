# Import libraries
import os
import shutil

# Set directory and file prefix lists
dir = 'SDF'

# Make new directories
os.mkdir(dir)

# Get files in current directory
path = os.getcwd()
files = [i for i in os.listdir(path)if os.path.isfile(os.path.join(path,i))]

# Move files from current directory to their new directoreis
mv_files = [i for i in files if '.sdf' in i]
for k in range(len(mv_files)):
    new_path = path + '/' + dir + '/' + mv_files[k]
    shutil.move(mv_files[k], new_path)