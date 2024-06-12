# BayesStatsProject
*By Nathaniel Tamminga*

## Description
For my final stats project, I decided to use Bayesian Optimization through GPyOpt to optimize a function that estimates the maximum energy of ions accelerated in laser-plasma interactoin (specifically TNSA). This model is known as the Mora model. My project is an extension of the work performed by Djordjevic in these two papers (https://doi.org/10.1088/1361-6587/ac172a and https://doi.org/10.1063/5.0045449). In these papers, Djordjevic uses 1D PIC simulations to and nerual networks to analyze the Mora model and make it more accuate.

The Mora model has a term called the acceleration time which is dependent on the laser pulse width $\tau$. Djordjevic suggests that this acceleration time is linearly dependent on the laser pulse width ($\tau$), laser intensity ($I_0$), target thickness (D), and prepulse scale length (Lg). In this project, I adjusted the acceleration time equation to fit suggestions from the Djordjevic papers.

$t_a = \alpha_1*\tau + \alpha_2*I_0 + \alpha_3*n_0 + \alpha_4*D + \alpha_5*Lg$

In the Djordjevic papers, he also varied $n_0$, but decided it had no notable effect on the ion energy. I included it for completeness. 

I randomly selected values in the ranges used by Djordjevic, then ran a 1D PIC sim to find the maximum ion energy accelerated from the target. I then compared this to the Mora model with the same loss function used by Djordjevic. I ran a Bayesian Optimization on this function, with the $\alpha$ parameters being set by the Bayesian Optimization to minimize the loss function, making the Mora model accurately predict the outcome of the PIC sims. I ran this code for 11 different sims with randomly selected conditions to test the robustness of the predicted improvements made to the Mora model.

The code used for this can be found in the bayes_opt.ipynb file.

## Results
Below are the results of my 11 runs.

Input Parameters
|Run    |$\tau$     |$I_0$              |$n_0$      |D          |Lg         |
|-------|-----------|-------------------|-----------|-----------|-----------|
|1      |50 fs      |$10^{20} W cm^{-2}$|$100 n_c$  |8 $\mu m$  |8 $\mu m$  |
|2      |88 fs      |$10^{19} W cm^{-2}$|$100 n_c$  |12 $\mu m$ |5 $\mu m$  |
|3      |103 fs     |$10^{18} W cm^{-2}$|$100 n_c$  |22 $\mu m$ |7 $\mu m$  |
|4      |139 fs     |$10^{17} W cm^{-2}$|$100 n_c$  |9 $\mu m$  |9 $\mu m$  |
|5      |166 fs     |$10^{17} W cm^{-2}$|$100 n_c$  |21 $\mu m$ |5 $\mu m$  |
|6      |183 fs     |$10^{18} W cm^{-2}$|$100 n_c$  |10 $\mu m$ |1 $\mu m$  |
|7      |313 fs     |$10^{17} W cm^{-2}$|$100 n_c$  |21 $\mu m$ |4 $\mu m$  |
|8      |352 fs     |$10^{17} W cm^{-2}$|$100 n_c$  |6 $\mu m$  |7 $\mu m$  |
|9      |364 fs     |$10^{18} W cm^{-2}$|$100 n_c$  |8 $\mu m$  |3 $\mu m$  |
|10     |456 fs     |$10^{17} W cm^{-2}$|$100 n_c$  |10 $\mu m$ |5 $\mu m$  |
|11     |487 fs     |$10^{20} W cm^{-2}$|$100 n_c$  |5 $\mu m$  |4 $\mu m$  |

Note: $n_0$ was unvaried due to the analysis by Djordjevic. It is also given in terms of the plasma critical density, $n_c$.

Optimization Results:
|Run    |$\alpha_1$     |$\alpha_2$     |$\alpha_3$     |$\alpha_4$     |$\alpha_5$     |Loss   |Mora Energy    |PIC Energy     |
|-------|---------------|---------------|---------------|---------------|---------------|-------|---------------|---------------|
|1      |-10.0          |-4.6534804     |7.64248169     |-4.22944659    |10.0           |1.018  |11.218 MeV     |9.7 MeV        |
|2      |0.06210955     |-6.33299868    |0.33894258     |-3.58003688    |1.70506284     |0.288  |4.441 MeV      |5.2 MeV        |
|3      |1.99713051     |-9.7275129     |0.35871941     |-3.41899509    |2.10753574     |0.002  |0.436 MeV      |0.5 MeV        |
|4      |9.54353805     |-8.25479491    |-6.06038722    |4.80694701     |-4.68524272    |0.0002 |0.020 MeV      |0.0 MeV        |
|5      |3.68739488     |7.04902376     |8.28054182     |7.67713727     |-7.50161491    |0.0006 |0.0356 MeV     |0.0 MeV        |
|6      |1.86909854     |-2.61244158    |4.58296491     |-0.86815546    |3.55097377     |0.024  |1.620 MeV      |1.4 MeV        |
|7      |-6.73127417    |-3.0693978     |7.51636318     |1.22456461     |3.53573298     |0.0008 |0.0598 MeV     |0.1 MeV        |
|8      |6.10474638     |-2.61296928    |-0.74514379    |-3.53476367    |0.44624222     |0.0002 |0.020 MeV      |0.0 MeV        |
|9      |6.02975377     |-7.4031575     |0.3387455      |-2.57702812    |1.84323673     |8.42e-5|1.287 MeV      |1.3 MeV        |
|10     |-5.25743831    |0.89534154     |8.06854283     |3.9441878      |-4.91798095    |9.38e-5|0.114 MeV      |0.1 MeV        |
|11     |3.42680065     |2.21027944     |1.73885014     |-6.96945185    |3.55784892     |8.28   |1.12 MeV       |9.9 MeV        |

## Conclusions
My conclusions are that I see no decernable trends in the coefficients that would lead to a robust extension to the Mora model. The Bayesian optimization worked good for the most part, particularly with runs 3,9, and 10 and with the exception of run 11.

More files or figures can be presented upon request. The result text files and python scripts are attatched in a .zip folder. This work was done using the EPOCH PIC code on the Ohio State Unity Cluster.