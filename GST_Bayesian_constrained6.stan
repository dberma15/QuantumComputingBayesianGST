data{
  int trials;    #integer indicating where the positive data starts
  int Fnumber;    #number of F gates involved
  int chiDimensionsDoubled; #dimensions of a gate
  int chiDimensionsDoubledSquared; #dimensions of a gate
  int chiDimensionsSquared; #dimensions of a gate
  #int measurements[trials];
  int numberOfGates; #number of gates
  int rhoDimensionsDoubled; #Dimensions of rho and E
  int rhoDimensionsDoubledSquared; #Dimensions of rho and E squared
  int rhoDimensions; #number of dimensions in rho and E
  int chiDimensions; #number of dimensions in the gate
  int measurementResults[numberOfGates,chiDimensionsSquared]; #1s and 0s indicating subject's choice. Rawdata equals 1 means the subject choose option B.
  real measurementCounts[numberOfGates,1,chiDimensionsSquared];
  matrix[chiDimensions,chiDimensions] Sigma; #Previous parameter for prior distribution. Not needed anymore
  matrix[rhoDimensions,rhoDimensions] Sigma2; #Previous parameter for prior distribution. Not needed anymore
  int F1[Fnumber]; #Number of gate components for the SPAM gates
  int F_map[Fnumber,1]; #Vector of mapping of gates to SPAM gates. Only 1 gate per SPAM gate.
  int rhoVectorDimensions; #size of the vector form of rho
  real<lower=0> nu; #parameter for prior distribution
  vector<lower=0>[rhoDimensions] alpha; #parameter for prior distribution
  int<lower=1> d; #value for fidelity calculation.
  matrix[chiDimensionsDoubled,chiDimensionsDoubled]Gates_Actual[numberOfGates]; #The Actual gates , used for trace distance
  int zero;
  matrix[chiDimensions, chiDimensions] chi_LGST_Real[numberOfGates];
  matrix[chiDimensions, chiDimensions] chi_LGST_Imag[numberOfGates];
 
  int chiDimensionsSquaredDoubled;
  matrix[chiDimensionsSquaredDoubled,chiDimensionsSquaredDoubled] A;
  real ceq_real[numberOfGates,chiDimensions,chiDimensions] ;
  real ceq_imag[numberOfGates,chiDimensions,chiDimensions];
  matrix[chiDimensions,chiDimensions] constraint_real;
  matrix[chiDimensions,chiDimensions] constraint_imag;
  int PauliDimensions;
  
  real constraint_real_std;
  real constraint_imag_std;
  real GmeasuredImag_std;
  real Gmeasured_std;
  }
 
#parameters of interest
parameters{

  #They must be specified as such so that they can be described as positive semi-definite
  cholesky_factor_corr[rhoDimensions] rho_chol_Real; #correlation cholesky decomposition matrix for the real part of rho.
  cholesky_factor_corr[rhoDimensions] rho_chol_Imag; #correlation cholesky decomposition matrix for the imaginary part of rho
  simplex[rhoDimensions] pi1; #Simplex for creating a trace of 1 for the real part of rho
 
  cholesky_factor_corr[rhoDimensions] E_chol_Real; #correlation cholesky decomposition matrix for the real part of E
  cholesky_factor_corr[rhoDimensions] E_chol_Imag; #correlation cholesky decomposition matrix for the imaginary part of E
  simplex[rhoDimensions] pi2; #Simplex for creating a trace of 1 for the real part of E

  #Array of matrices. The length of the array is the number of gates that are being estimated.
  matrix<lower=-1.05,upper=1.05>[chiDimensions,chiDimensions] chi_Real[numberOfGates];
  matrix<lower=-1.05,upper=1.05>[chiDimensions,chiDimensions] chi_Imag[numberOfGates];
}

transformed parameters{
    #matrix array of the gates that are going to be estimated
    matrix<lower=-1.05,upper=1.05>[chiDimensionsDoubled, chiDimensionsDoubled] chi[numberOfGates];
   
    #matrix array of the initial state to be estimated, the real part of it and the imaginary part of it.
    matrix[rhoDimensionsDoubled, rhoDimensionsDoubled] rho;
    matrix[rhoDimensions, rhoDimensions] rho_Real;
    matrix[rhoDimensions, rhoDimensions] rho_Imag;
   
    #matrix array of the initial state to be estimated, the real part of it and the imaginary part of it.
    matrix[rhoDimensionsDoubled, rhoDimensionsDoubled] E;
    matrix[rhoDimensions, rhoDimensions] E_Real;
    matrix[rhoDimensions, rhoDimensions] E_Imag;
    matrix<lower=0,upper=1>[chiDimensions,chiDimensions] Gmeasured[numberOfGates];

    #Parameter to ensure that I-E is positive semidefinite.
    #matrix[rhoDimensions,rhoDimensions] E_Real_One_minus;
    #cov_matrix[rhoDimensions] E_Real_One_minus;

    vector[rhoDimensions] sds1;
    vector[rhoDimensions] sds2;
   
    matrix[rhoDimensionsDoubled,rhoDimensions] rho_temp;
    vector[rhoVectorDimensions] rhoVector;

    matrix[rhoDimensionsDoubled,rhoDimensions] E_temp;
    row_vector[rhoVectorDimensions] EVector;
   
    matrix[chiDimensionsDoubled,chiDimensionsDoubled] chi_Process[numberOfGates];
    matrix[chiDimensions,chiDimensions] chi_Process_real[numberOfGates];
    matrix[chiDimensions,chiDimensions] chi_Process_imag[numberOfGates];
    #cov_matrix[chiDimensions] chi_Process_real[numberOfGates];
    #cov_matrix[chiDimensions] chi_Process_imag[numberOfGates];
    matrix[chiDimensionsDoubled,chiDimensions] chi_k_temp;
    vector[chiDimensionsSquaredDoubled] chi_k_Vector;
    vector[chiDimensionsSquaredDoubled] tempProcess_k;

    matrix[chiDimensionsDoubled, chiDimensionsDoubled] chiReverse[numberOfGates];
    matrix<lower=-.05,upper=.05>[chiDimensions, chiDimensions] GmeasuredImag[numberOfGates];
	matrix[numberOfGates,PauliDimensions] const_real_sum;
    matrix[numberOfGates,PauliDimensions] const_imag_sum;
	
    #takes the square roots of each element in pi1 and pi2 so that sds2*t(sds2) is a square matrix with pi2 on the diagonal and the same for pi1
    for (i in 1:rhoDimensions) {
        sds2[i] <- sqrt(pi2[i]);
        sds1[i] <- sqrt(pi1[i]);
    }   
   
   
    rho_Real <- multiply_lower_tri_self_transpose(diag_pre_multiply(sds1, rho_chol_Real) );      #Finds E_Real by multiplying t(rho_chol_Real)*sds1*t(sds1)*rho_chol_Real
    rho_Imag <- (rho_chol_Imag-(rho_chol_Imag)');     #Creates a hollow matrix with zeros on the diagonal
    rho<-append_row(append_col(rho_Real,multiply(-1,rho_Imag)),append_col(rho_Imag,rho_Real));     #creates a block matrix


    E_Real <- multiply_lower_tri_self_transpose(diag_pre_multiply(sds2, E_chol_Real) );     #Finds E_Real by multiplying t(E_chol_Real)*sds2*t(sds2)*E_chol_Real
    E_Imag <- (E_chol_Imag-(E_chol_Imag)');     #Creates a hollow matrix with zeros on the diagonal
    E<-append_row(append_col(E_Real,multiply(-1,E_Imag)),append_col(E_Imag,E_Real));     #creates a block matrix

    #E_Real_One_minus<-(diag_matrix(rep_vector(1, rows(E_Real)))-(E_Real)); #ensures I-E is positive semidefinite

    #Creates a block matrix of [Re(chi),-Im(chi); Im(chi), Re(chi)] for form the new chi matrix
    for (k in 1:numberOfGates){
        chi[k]<-append_row(append_col(chi_Real[k],multiply(-1,chi_Imag[k])),append_col(chi_Imag[k],chi_Real[k]));
    }
   
    for (k in 1:numberOfGates){
        chiReverse[k]<-append_row(append_col(multiply(-1,chi_Imag[k]),chi_Real[k]),append_col(chi_Real[k],chi_Imag[k]));
    }

    #Takes the first column of blocks for reformatting into a vector
    for (j in 1:rhoDimensions){
        for (i in 1:rhoDimensions){
            rho_temp[i,j]<-rho[i,j];
            E_temp[i,j]<-E[i,j];
        }
        for (i in (rhoDimensions+1):(rhoDimensionsDoubled)){
            rho_temp[i,j]<-rho[i,j];
            E_temp[i,j]<-E[i,j];
        }
    }
    #reformats into a vector
    rhoVector<-to_vector(rho_temp');
    EVector<-to_row_vector(E_temp');
   
    for (k in 1:numberOfGates){
        for (j in 1:Fnumber){
            for (i in 1:Fnumber){
                Gmeasured[k][i,j]<-EVector*(chi[F_map[j,1]]*(chi[k]*(chi[F_map[i,1]]*rhoVector)));
                GmeasuredImag[k][i,j]<-EVector*(chiReverse[j]*(chiReverse[k]*(chiReverse[i]*rhoVector)));
            }
        }   
    }

    for (k in 1:numberOfGates){
        for (j in 1:chiDimensions){
            for (i in 1:chiDimensions){
                chi_k_temp[i,j]<-chi[k][i,j];
            }
            for (i in (chiDimensions+1):(chiDimensionsDoubled)){
                chi_k_temp[i,j]<-chi[k][i,j];
            }
        }
        chi_k_Vector<-to_vector(chi_k_temp');
        tempProcess_k<-inverse(A)*chi_k_Vector;
        #count<-1;
        for (i in 1:chiDimensions){
            for (j in 1:chiDimensions){
                chi_Process_real[k][i,j]<-tempProcess_k[chiDimensions*(i-1)+j];
                chi_Process_imag[k][i,j]<-tempProcess_k[chiDimensions*(i-1)+j+chiDimensionsSquared];
                #count<-count+1;
            }
        }
        chi_Process[k]<-append_row(append_col(chi_Process_real[k],multiply(-1,chi_Process_imag[k])),append_col(chi_Process_imag[k],chi_Process_real[k]));
    }
	
	for (k in 1:numberOfGates){
        for (rr in 1:PauliDimensions){
            const_real_sum[k,rr]<-zero;
            const_imag_sum[k,rr]<-zero;
        }
    }


    for (k in 1:numberOfGates){
        for (rr in 1:PauliDimensions){
            for(mm in 1:chiDimensions){
                for(nn in 1:chiDimensions){
                    const_real_sum[rr,k]<-const_real_sum[rr,k]+(chi_Process_real[k][mm,nn]*ceq_real[rr,mm,nn])-chi_Process_imag[k][mm,nn]*ceq_imag[rr,mm,nn];
                    const_imag_sum[rr,k]<-const_imag_sum[rr,k]+(chi_Process_imag[k][mm,nn]*ceq_real[rr,mm,nn])+chi_Process_real[k][mm,nn]*ceq_imag[rr,mm,nn];
                }
            }
        }
    }

}
model {
    int count;
    real measurementRates;

    matrix[chiDimensionsDoubled,chiDimensionsDoubled] FI[numberOfGates];
    matrix[chiDimensionsDoubled,chiDimensionsDoubled] FJ[numberOfGates];
    matrix[chiDimensionsDoubled,chiDimensionsDoubled] FIimag[numberOfGates];
    matrix[chiDimensionsDoubled,chiDimensionsDoubled] FJimag[numberOfGates];

     

   
    rho_chol_Real ~ lkj_corr_cholesky(nu); #Uniform distribution prior for cholesky matrix
    rho_chol_Imag ~ lkj_corr_cholesky(nu); #Uniform distribution prior for cholesky matrix
    pi1~dirichlet(alpha); #uniform distribution prior for simplex
   
    E_chol_Real ~ lkj_corr_cholesky(nu); #Uniform distribution prior for cholesky matrix
    E_chol_Imag ~ lkj_corr_cholesky(nu); #Uniform distribution prior for cholesky matrix
    pi2~dirichlet(alpha); #uniform distribution prior for simplex


    #uniform prior on chi_real and chi_imag using uniform distribution from -1 to 1
    for (k in 1:numberOfGates){
        for (i in 1:chiDimensions){
            for (j in 1:chiDimensions){
                chi_Real[k][i,j] ~ uniform(-1.05,1.05);#cauchy(chi_LGST_Real[k][i,j],.05); #
                chi_Imag[k][i,j] ~ uniform(-1.05,1.05);#cauchy(chi_LGST_Imag[k][i,j],.05); #
			}
        }
    }


   
    #editable prior
    for (k in 1:numberOfGates){
        for (rr in 1:PauliDimensions){
            constraint_real[rr,k]~normal(const_real_sum[rr,k],constraint_real_std);
            constraint_imag[rr,k]~normal(const_imag_sum[rr,k],constraint_imag_std);
        }
    }       
   
    #Uses the measurements as a bernoulli distribuion of the Gmeasured value.
    #editable prior
    for (k in 1:numberOfGates){
        count<-1;
        for (i in 1:Fnumber){
            for (j in 1:Fnumber){
                zero~normal(GmeasuredImag[k][i,j],GmeasuredImag_std);
                measurementCounts[k,1,count]~normal(Gmeasured[k][i,j],Gmeasured_std);
                count<-count+1;
            }
        }   
    }
}