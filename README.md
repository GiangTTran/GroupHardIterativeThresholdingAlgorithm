# GroupHardIterativeThresholdingAlgorithm
Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward, 2018.

Reference: "Learning Dynamical Systems and Bifurcation via Group Sparsity", https://arxiv.org/abs/1709.01558
1. group_hard_iterative_thresholding.m

      Solve  
 
           min_{c1,...,cm} \|D1 * c1 - b1\|_2^2 + ... \|Dm * cm - bm\|_2^2 + gamma * \|[c1,..., cm]\|_{2,0}  

      where \|A\|_{2,0} = number of non-zero rows of the matrix A.   

2. dictionary1d.m

      Dictionary matrix for 1D data: include all monomials of degree up to p
      
3. dictionary3.m

      Dictionary matrix for 3D data: include all multivariable monomials of degree up to p

4. time_derivative.m

      Velocity approximation from data using 1st/2nd-order approximation of time derivative (excluding the initial and/or the final time data)

5. time_derivative1.m

      Velocity approximation from data using 1st/2nd-order approximation of time derivative (including the initial and the final time data)

6. Test_logistic_100experiments.m

      Test our group hard iterative thresholding algorithm on recovering the logistic equation
           xdot = alpha * x * (1-x)
      where data comes from two sources corresponding to alpha = 0.05 and alpha = 0.23. 
       
      Run the experiment 100 times and report the probability of recovering correct terms in the governing equations. Compare       with the results using l^0-model that treats each data separately.
7. Test_Lorenz_5sets.m

      Test our group hard iterative thresholding algorihm on recovering the well-known Lorenz system where data comes from           five sources corresponding to five bifurcation parameters.
       
      Output the recovered coefficient matrix and compute the relative error.
8. Test_Lorenz_5sets_100experiments.m

      Test our group hard iterative thresholding algorihm on recovering the well-known Lorenz system where data comes from           five sources corresponding to five bifurcation parameters.
       
      Run the experiment 100 times and report the probability of recovering correct terms in the governing equations as well         as the rellative error ranges for all five data sets.

      
