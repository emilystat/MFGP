The repository contains sample code as a realization of Multivariate Fused Gaussian Process (MFGP).

Both folders contain examples using MFGP. In folder "2var_non_stationary", simulated data containing two independent variables are generated from non-stationary and asymmetric process using conditional approach.

In folder "3var_stationary", data containing three independent variables are simulated from stationary process.

In both folders, "Y_simu.mat" contains the response variable generated using "generate_data.m". 

"main.m" is the main folder to do the job.

Functions in "EM_MCAR_2Vars.m" and "EM_MCAR_3Vars.m" are the MFGP algorithms for 2 variables and 3 variables respectively.

"crps.m" is the metric function applied in the paper to conduct comparisons. ref(1)

"fminsearchbnd.m" contains the optimization function used in the algorithm. ref(2)

"Predict_missings_v2.m" and "Predict_missings_v3.m" are used for prediction purpose.

"variogram.m" is implemented to compute variogram. ref(3)

All the other files are copied directely from SuiteSparse software ref(4) to conduct Cholesky decomposition.




Reference:

(1) Durga Lal Shrestha (2021). Continuous rank probability score (https://www.mathworks.com/matlabcentral/fileexchange/47807-continuous-rank-probability-score), MATLAB Central File Exchange. Retrieved July 3, 2021

(2) John D'Errico (2021). fminsearchbnd, fminsearchcon (https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon), MATLAB Central File Exchange. Retrieved July 3, 2021.

(3) Wolfgang Schwanghart (2021). Experimental (Semi-) Variogram (https://www.mathworks.com/matlabcentral/fileexchange/20355-experimental-semi-variogram), MATLAB Central File Exchange. Retrieved July 3, 2021.

(4) SUITESPARSE : A SUITE OF SPARSE MATRIX SOFTWARE (https://people.engr.tamu.edu/davis/suitesparse.html)
