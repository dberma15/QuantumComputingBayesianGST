# QuantumComputingBayesianGST

The following code is written in R and stan. Code that implements this is written in Matlab.

## Purpose 
Frequentist statistics requires Monte Carlo simulations in order to get error bars on estimates. For large and computationally difficult problems, this might require a cluster. However, this can be avoided with Bayesian statistics, which provides error bars on estimates. Rather than use MLE to perform Gate Set Tomogrpahy, this code uses Bayesian Statistics. This code is initially run using Matlab code (which cannot be posted). The Matlab code runs the R code, which implements Rstan. The results of Rstan are then fed into Matlab to complete the analysis. 


## Input to R

The input to R is two files temp.data and temp.init. Both are generated form the program, but data.init is more associated with initial values and data.temp is associated with constants or parameters.

### Data.Temp

measurementResults: a structure of a vector with dimensions (number of rows for each gate) x (number of gates * number of columns for each gate). It contains the number of counts observed and should equivalent to lining up the gates. 
       
    If input into R as: 
    structure(c(10, 5, 5, 0, 5, 0, 5, 5, 5, 5, 0, 5, 0, 5, 5, 10, 5, 0, 5, 5, 0, 5, 0, 10, 5, 5, 5, 5, 5, 10, 5, 5, 5, 5, 0, 5, 5, 5, 5, 5, 0, 0, 5, 0, 5, 5, 10, 5, 0, 5, 5, 10, 5, 10, 5, 5, 5, 5, 10, 5, 10, 5, 5, 0), Dim = c(4,16))
    It should be stored in R as:
    10  5  5  0  5  0  5  5  5  5  0  5  0  5  5 10
     5  0  5  5  5  5  5 10  5  5  0  5  5 10  5  5
     5  5  0  5  0  0  5  5  0  5  5 10  5  5 10  5
     0  5  5  5 10  5  5  5  5  5  0  5 10  5  5  0
     
trials: An integer containing the number of trials performed for each ijk measurement

Fnumber: An integer representing the number of F gates

chiDimensions: An integer representing the number of rows (or columns) in the chi representation of each matrix in operator form.

rhoDimensions: An integer representing the number of rows (or columns) in the initial (or final) state in oeprator form.

chiDimensionsSquared: ChiDimensions squared

chiDimensionsDoubled: ChiDimensions doubled

chiDimensionsDoubledSquared: ChiDimensions doubled, then squared

rhoDimensionsDoubled: rhoDimensions squared

rhoVectorDimensions: The number of dimensions in the superoperator form of the initial/final state. This twice rhoDimensionsDoubled

rhoDimensionsDoubledSquared: rhoDimensions doubled, then squared

numberOfGates: number of gates

F_map: A maping of the SPAM gates to the gates that are to be estimated.

A_new: The transformation that changes the vectorized superoperator G into the process matrix Chi, which is formed using the pauli matricies and the Kronecker product: A_new[n,]=vector(transpose(T) _ j x P_i) where x is the tensor product. vec(chi)=A_new * Vec(G). See equaitons 25 and 26. 

PauliDimensions: Number of pauli matrices

ceq_real: An array that is (number of gates x rows per gate x columns per gate) and contains real parts of the trace of the product of the i * j * k Pauli Matrices. This is constant so it is easy to calculate outside of Stan and import it in. 

ceq_real: An array that is (number of gates x rows per gate x columns per gate) and contains imaginary parts of the trace of the product of the i * j * k Pauli Matrices. This is constant so it is easy to calculate outside of Stan and import it in. 

constraint_real: Imports the physicality constraints on the real parts of the process matrices. See equation 28 for origin. 

constraint_imag: Imports the physicality constraints on the imaginary parts of the process matrices. See equation 28 for origin. 

alpha: should be c(1,1). Acts as the parameters generating a uniform Dirichlet distribution prior. 

## Data.Init

chi_Init: initialized values from GST of the process matrices in block form. It is an array with dimensions (number of gates) x (2 * rows in chi) x (2 * columns in chi)
      
      Note: Block form is when a matrix A is rewritten as
            Re(A) -Im(A)
            Im(A)  Re(A)
chi_Init_Real: the real part of the process matrix. It is an array with dimensions (number of gates) x (rows in chi) x (columns in chi)

chi_Init_Imag: the imaginary part of the process matrix. It is an array with dimensions (number of gates) x (rows in chi) x (columns in chi)

rho_init: initialized values from GST of the initial state in block form. It is an array with dimensions (2 * rows in rho) x (2 * columns in rho)

rho_init_Real: the real part of the GST estimate of the initial state. It is an array with dimensions (rows in rho) x (columns in rho)

rho_init_Imag: the imaginary part of the GST estimate of the initial state. It is an array with dimensions (rows in rho) x (columns in rho)

rho_init_chol_real: the cholesky decompositions of the real part of the GST estimate of the initial state. It is an array with dimensions (rows in rho) x (columns in rho)

rho_init_chol_imag: the cholesky decompositions of the imaginary part of the GST estimate of the initial state. It is an array with dimensions (rows in rho) x (columns in rho)

E_init: initialized values from GST of the final state in block form. It is an array with dimensions (2 * rows in rho) x (2 * columns in rho)

E_init_Real: the real part of the GST estimate of the final state. It is an array with dimensions (rows in rho) x (columns in rho)

E_init_Imag: the imaginary part of the GST estimate of the final state. It is an array with dimensions (rows in rho) x (columns in rho)

E_init_chol_real: the cholesky decompositions of the real part of the GST estimate of the final state. It is an array with dimensions (rows in rho) x (columns in rho)

E_init_chol_imag: the cholesky decompositions of the imaginary part of the GST estimate of the final state. It is an array with dimensions (rows in rho) x (columns in rho)

pi1_init: Initial values of a simplex so that the trace of rho can sum to 1

pi2_init: Initial values of a simplex so that the trace of E can sum to 1

sds1_init: square root of pi1

sds2_init: square root of pi2
