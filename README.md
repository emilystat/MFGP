# MFGP-code
The repository contains sample code as a realization of Multivariate Fused Gaussian Process (MFGP).

Both folders contain examples using MFGP. In folder "2var_non_stationary", simulated data containing two independent variables are generated from non-stationary and asymmetric process using conditional approach.
In folder "3var_stationary", data containing three variables are simulated from stationary process.

In both folders, "Y_simu.mat" contains the response variable generated using "generate_data.m". 
"main.m" is the main folder to do the job.
Functions in "EM_MCAR_2Vars.m" and "EM_MCAR_3Vars.m" are the MFGP algorithms for 2 variables and 3 variables respectively.
"crps.m" is the metric function applied in the paper to conduct comparisons.
"fminsearchbnd.m" contains the optimization function used in the algorithm.
"Predict_missings_v2.m" and "Predict_missings_v3.m" are used for prediction purpose.
"variogram.m" is implemented to compute variogram.

All the other files are copied directely from SuiteSparse software (https://people.engr.tamu.edu/davis/suitesparse.html) to conduct Cholesky decomposition.
