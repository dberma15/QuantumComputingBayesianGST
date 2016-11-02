
# Load the rstan and other package
library(Rcpp)
library(inline)
library(rstan)
library(MASS)
library(ggmcmc)
library(coda)
library(Matrix)
library(R.matlab)
library(matrixStats)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#*******************************Loaded Data**************************************
#********************************************************************************
# The first step is to set the working directory, e.g.: 
workingDirectory="C:\\Users\\daniel\\Documents\\QuantumComputingIndependentStudy\\gst_matlab"
setwd(workingDirectory)
#Loads the necessary data and filename
source('temp.data.R')
source('temp.inits.R')
outputfilename=readLines('outputfilename.txt',n=1)




#Sets the model you are working with and translates and compiles it once. 
translatedModel<-stanc(file="GST_bayesian_constrained6.stan",model_name='CPTbehavioral')
compiledModel<-stan_model(stanc_ret=translatedModel,verbose=FALSE)
#Creates a list for F_map to import into RStan
F_map2<-list()
for (i in 1:Fnumber){
  F_map2[[i]]<-F_map[i]
}


#ensures that pi2 and pi1 are proper simplexes by splitting the differences of extra values.
if (sum(pi1_init)>1){
  pi1_init<-pi1_init-(sum(pi1_init)-1)/2
}
if (sum(pi2_init)>1){
  pi2_init<-pi2_init-(sum(pi2_init)-1)/2
}
if (sum(pi1_init)<1){
  pi1_init<-pi1_init+(1-sum(pi1_init))/2
}
if (sum(pi2_init)<1){
  pi2_init<-pi2_init+(1-sum(pi2_init))/2
}

measurementCounts<-array(0,c(numberOfGates,1,chiDimensionsSquared))
for (i in 1:numberOfGates){
  for (j in 1:chiDimensionsSquared){
    # measurementCounts[i,1,j]<-sum(measurementResults[i,,j])
    measurementCounts[i,1,j]<-measurementResults[i,j]
  }
}
#measurementCounts<-round(measurementCounts,digits=1)
chiDimensionsSquaredDoubled<-2*chiDimensionsSquared
measurementCounts<-measurementCounts/trials
zero<-0
chi_Init_Imag<-pmin(chi_Init_Imag,.99999)
chi_Init_Real<-pmin(chi_Init_Real,.99999)
chi_Init<-pmin(chi_Init,.99999)
constraint_real_std<-.1#.2;
constraint_imag_std<-.1#.2;
GmeasuredImag_std<-.05#.1;
Gmeasured_std<-.05#.01;

#normally this would be done with GenerateKroneckerBasis function
expData = list('trials'=trials,
               'Fnumber'=Fnumber,
               'chiDimensions'=chiDimensions,
               'rhoDimensions'=rhoDimensions,
               'chiDimensionsSquared'=chiDimensionsSquared,
               'chiDimensionsDoubled'=chiDimensionsDoubled,
               'chiDimensionsDoubledSquared'=chiDimensionsDoubledSquared,
               'rhoDimensionsDoubled'=rhoDimensionsDoubled,
               'rhoVectorDimensions'=rhoVectorDimensions,
               'rhoDimensionsDoubledSquared'=rhoDimensionsDoubledSquared,
               'measurementResults'=measurementResults,
               'numberOfGates'=numberOfGates,
               'F_map'=F_map2,
               'nu'=1,
               'Gates_Actual'=Gates_Actual,
               'alpha'=alpha,
               'zero'=zero,
               'chi_LGST_Real'=chi_Init_Real,
               'chi_LGST_Imag'=chi_Init_Imag, 
               'A'=A_new,
               'PauliDimensions'=PauliDimensions,
               'ceq_real'=ceq_real,
               'ceq_imag'=ceq_imag,
               'constraint_real'=constraint_real,
               'constraint_imag'=constraint_imag,
               'chiDimensionsSquaredDoubled'=chiDimensionsSquaredDoubled,
               'measurementCounts'=measurementCounts,
		   'constraint_real_std'=constraint_real_std,
  		   'constraint_imag_std'=constraint_imag_std,
  	 	   'GmeasuredImag_std'=GmeasuredImag_std,
               'Gmeasured_std'=Gmeasured_std)


inits = function(){
  return(
    list(chi=chi_Init,
         chi_Real=chi_Init_Real,
         chi_Imag=chi_Init_Imag,               
         
         rho=rho_init,
         rho_Real=rho_init_real,
         rho_Imag=rho_init_imag,
         rho_chol_Real=rho_init_chol_real,
         rho_chol_Imag=rho_init_chol_imag,
         
         E=E_init,
         E_Real=E_init_real,
         E_Imag=E_init_imag,
         E_chol_Real=E_init_chol_real,
         E_chol_Imag=E_init_chol_imag,
         
         sds1=sds1_init,
         sds2=sds2_init,
         
         pi1=pi1_init,
         pi2=pi2_init,
         traceDist=traceDist_Init)
  )
}

#Bayesian estimation
#hierarchical= sampling(compiledModel,data=expData, init=inits,iter=30000, chains=3,thin=5,warmup=1000, pars=c('chi', 'rho','E', 'const_real_sum','const_imag_sum'))#,
hierarchical= sampling(compiledModel,data=expData, init=inits,iter=40000, chains=10,thin=5,warmup=1000, pars=c('chi', 'rho','E', 'const_real_sum','const_imag_sum'))#,
print(hierarchical, digits=6)
#Calculation of means and covariance matrices
mcmc_hierarchical<-mcmc.list(lapply(1:ncol(hierarchical), function(x) mcmc(as.array(hierarchical)[,x,]))) #for checking of iteration values change
# traceplot(mcmc_hierarchical)
stan_trace(hierarchical,pars="chi[1,1,1]")
stan_dens(hierarchical,pars="chi[1,1,1]")
stan_trace(hierarchical,pars="chi[2,1,1]")
stan_dens(hierarchical,pars="chi[2,1,1]")
stan_trace(hierarchical,pars="chi[2,2,2]")
stan_dens(hierarchical,pars="chi[2,2,2]")
covarmatrix<-cov(as.matrix(hierarchical))
relevantCols<-grep('rho\\[|E\\[|chi\\[|traceDist',colnames(covarmatrix))
relevantRows<-grep('rho\\[|E\\[|chi\\[|traceDist',rownames(covarmatrix))
relevantCovar<-diag(covarmatrix[relevantRows,relevantCols])

#Calculates the means
parameterMeans<-colMedians(as.matrix(hierarchical))
names(parameterMeans)<-names(colMeans(as.matrix(hierarchical)))

chiLoc<-grep('chi\\[',names(parameterMeans))
chi<-array(parameterMeans[chiLoc],c(numberOfGates,chiDimensionsDoubled,chiDimensionsDoubled))

ELoc<-grep('E\\[',names(parameterMeans))
E<-matrix(parameterMeans[ELoc],rhoDimensionsDoubled,rhoDimensionsDoubled)

rhoLoc<-grep('rho\\[',names(parameterMeans))
rho<-matrix(parameterMeans[rhoLoc],rhoDimensionsDoubled,rhoDimensionsDoubled,byrow=TRUE)

traceDistLoc<-grep('traceDist\\[',names(parameterMeans))
traceDist_measured<-matrix(parameterMeans[traceDistLoc],numberOfGates,1,byrow=TRUE)

#reformats the estimates into proper matrices
chi_mat<-array(0, dim = c(numberOfGates,chiDimensions,chiDimensions))
for (j in 1:numberOfGates){
  chi_mat[j,,]<-(chi[j,1:chiDimensions,1:chiDimensions]+
                   1i*chi[j,(chiDimensions+1):chiDimensionsDoubled,1:chiDimensions]) 
}
E_mat<-t(matrix(t(E[1:rhoDimensions,1:rhoDimensions]),rhoDimensions,rhoDimensions))+
  1i*t(matrix(t(E[(1+rhoDimensions):rhoDimensionsDoubled,1:rhoDimensions]),rhoDimensions,rhoDimensions))
rho_mat<-matrix(t(rho[1:rhoDimensions,1:rhoDimensions]),rhoDimensions,rhoDimensions)+
  1i*matrix(t(rho[(1+rhoDimensions):rhoDimensionsDoubled,1:rhoDimensions]),rhoDimensions,rhoDimensions)

#computes the standard deviations and formats them into matrices
chiCovar<-array(relevantCovar[grep('chi',names(relevantCovar))],dim = c(numberOfGates,chiDimensionsDoubled,chiDimensionsDoubled))
chi_Std<-sqrt(chiCovar)

rhoCovar<-array(relevantCovar[grep('rho',names(relevantCovar))],dim = c(rhoDimensionsDoubled,rhoDimensionsDoubled))
rhoStd<-sqrt(rhoCovar)

ECovar<-array(relevantCovar[grep('E',names(relevantCovar))],dim = c(rhoDimensionsDoubled,rhoDimensionsDoubled))
EStd<-sqrt(ECovar)

traceDistCovar<-array(relevantCovar[grep('traceDist',names(relevantCovar))],dim = c(numberOfGates,1))
traceDistStd<-sqrt(traceDistCovar)


EStd<-t(matrix(t(EStd[1:rhoDimensions,1:rhoDimensions]),rhoDimensions,rhoDimensions))+
  1i*t(matrix(t(EStd[(1+rhoDimensions):rhoDimensionsDoubled,1:rhoDimensions]),rhoDimensions,rhoDimensions))
rhoStd<-matrix(t(rhoStd[1:rhoDimensions,1:rhoDimensions]),rhoDimensions,rhoDimensions)+
  1i*matrix(t(rhoStd[(1+rhoDimensions):rhoDimensionsDoubled,1:rhoDimensions]),rhoDimensions,rhoDimensions)

chiStd<-array(0, dim = c(numberOfGates,chiDimensions,chiDimensions))
for (j in 1:numberOfGates){
  chiStd[j,,]<-(chi_Std[j,1:chiDimensions,1:chiDimensions]+
                  1i*chi_Std[j,(chiDimensions+1):chiDimensionsDoubled,1:chiDimensions]) 
}

hierarchical_raw<-rstan::extract(hierarchical)


rho_actual<-matrix(c(1,0,0,0),2,2)
E_actual<-rho_actual

estimatedmeasurements<-array(0,c(numberOfGates,numberOfGates,numberOfGates))
for (j in 1:numberOfGates){
  for (k in 1:numberOfGates){
    for (i in 1:numberOfGates){
      estimatedmeasurements[k,i,j]<-matrix(E_mat,1,4)%*%(chi_mat[j,,]%*%chi_mat[k,,]%*%chi_mat[i,,])%*%matrix(rho_mat,4,1)
    }
  }
}
#Saves the matrices and strandard deviations and exports them
estimates = list('chi_mat'=chi_mat,
                 'rho_mat'=rho_mat,
                 'E_mat'=E_mat,
                 'chiStd'=chiStd,
                 'rhoStd'=rhoStd,
                 'EStd'=EStd,
                 'raw_chi'=hierarchical_raw$chi,
                 'traceDist_measured'=traceDist_measured,
                 'traceDistStd'=traceDistStd)
outputfilename=paste(workingDirectory,outputfilename,sep='\\')
outputfilename

writeMat(con=outputfilename,x=estimates)
