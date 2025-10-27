The files include R sampling code and MATLAB code for WAIC decomposition.The R sampling code consists of Stan code and main code; running the main code will generate sampling data.Use the MATLAB decomposition code to perform decomposition calculations on the example data.

1. rt.csv ， resp.csv: Example data
2. m.R:  Uses rstan for MCMC sampling, including the joint model and conditional models.
3. WAIC_calculate_total.m:  The main MATLAB file for WAIC calculation
4. JointpdfMatrix.m ,  Var_Cov.m:  Auxiliary functions for WAIC calculation
5. m.stan，m_rt_give_resp.stan，m_resp_give_rt.stan：Stan code for joint and conditional sampling
6. samples.csv:  Joint model sampling data of the example data using Stan