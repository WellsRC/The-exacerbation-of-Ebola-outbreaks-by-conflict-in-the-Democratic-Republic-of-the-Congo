# The exacerbation of Ebola outbreaks by conflict in the Democratic Republic of the Congo
Chad R. Wells, Abhishek Pandey, Martial L. Ndeffo Mbah, Bernard-A. Gaüzère, Denis Malvy, Burton H. Singer, Alison P. Galvani

Copyright 2019, Chad Wells et al. All rights reserved. Released under the GNU AGPL 3.

The MATLAB code provided here will run the fitting and analysis for the modelling portion of the manuscript.

# Data
DataNorthKivu.mat - Provides the data used to fit the model

# Fitting
BM - Runs a primary search of the parameter space and calculates the log-likelihood

BMFocus - Runs a more focused search around the more likely parameter sets to improve the log-likelihood and obtain additional unique parameter sets in the sampling process for the analysis

FitEpiSimNK - Computes the negative log-likelihood of the model fit in case fmincon or other minimization algorithms want to be applied

# Model
EpiSimDVG - The system of ODEs. It returns the different latent classes, infectious classes, and cumulative incidence

EpiSimDVGAlt - The system of ODEs. It returns the different latent classes, infectious classes, cumulative incidence, cumulative recovery/rmeoval, and cumulative latent infections

EpiSimNKV - Runs the Euler method with step-size 0.1 to solve the system of ODEs (EpiSimDVG). It returns the different latent classes, infectious classes, and cumulative incidence over the specified time

EpiSimDCVGREff- Runs the Euler method with step-size 0.1 to solve the system of ODEs (EpiSimDVGAlt). This functino was used in the analysis portion. It returns It returns the different latent classes, infectious classes, cumulative incidence, cumulative recovery/rmeoval, cumulative latent infections, R Effective and the effects of conflict applied in the model.

TDCD - Computes the value of the effect of conflict at time t. Contains the information about the type, location, and time of attack. As well, contains the weights for the impact of the attack based on the location.

TriDist - Computes the impact of the attacks included in the model relative to time t. There is an increase prior to the attack and exponetial decay after the attack

SatF- The function returns the value of the specified saturation functions

# Analysis
ValRunNKL- Samples the parameter sets, then calculates output from the model. These results are saved to a ModelFit.mat

Figure2- Produces a figure based on the output from ValRunNKL

