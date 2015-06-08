
#' Initialize the parameters for the DLM used for a specific ration.
#' 
#' @param tMax Maximum growth period in weeks.
#' @param k2 Value of \eqn{k_2} in the \eqn{F} matrix.
#' @param FF The \eqn{F} matrix.
#' @param GG The \eqn{G} matrix.
#' @param V The \eqn{V} matrix.
#' @param W The \eqn{W} matrix.
#' @param m0 Prior latent mean \eqn{m_0}.
#' @param C0 Prior latent covariance matrix \eqn{C_0}.
#' @return A list containing all the parameters
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
iniDLM<-function(tMax=15,
                 k2=1.549,
                 FF=matrix(data=c(1,0,0.044,k2),ncol=2), 
                 GG=matrix(c(1,0,1,1),nrow=2), 
                 V=matrix(c(0.066,0.027,0.027,0.012),ncol=2),#matrix(data=c(0.066,0.00013,0.00013,0.0019),ncol=2),#matrix(c(0.25,0.013,0.013,0.73),ncol=2), 
                 W=matrix(c(2.1,-0.124,-0.124,0.112),ncol=2),#matrix(c(0,0,0,0.112),ncol=2),#matrix(data=c(0,0,0,0.0066),ncol=2),#matrix(c(2.194,1.167,1.167,1.601),ncol=2), 
                 m0=matrix(c(26.49,6),nrow=2), 
                 C0=matrix(c(4.26,0.32,0.32,0.53),ncol=2)
){
  model<-list(t=t)
  model$k2<-k2
  model$FF<-t(FF)    # note FF in dlm package is actually F' in W&H book!
  model$GG<-GG
  model$V<-V
  model$W<-W
  model$m0<-m0
  model$C0<-C0
  return(model)
}


####################################################################


#' Estimate values of \eqn{k_{1,t}} in \eqn{F_t} of the DLM and calculate the
#' covariance matrices to be used when have to calculate the transition
#' probabilities.
#' 
#' Estimate using formula \eqn{k_{1,t} = k_1*a_t^(-0.25)} where \eqn{a_t} is the prior at time \eqn{t}. 
#' 
#' @param iniDLM The DLM created using \link{iniDLM}.
#' @param Y A matrix with 2 observation rows (weight and intake). One column for each week. Used to estimate
#'   \eqn{k_{1,t}}. Could be mean values observed in the herd.
#' @param k1 Value of \eqn{k_1} used in the estimation.
#' @return A list containing all the parameters. Element L contains the covariance matrix for transition probabilities. 
#'   Element JFF is a (2x2) matrix indicating where \eqn{F_t'} depends on \eqn{t} (i.e. a 1 is put where \eqn{k_t}) is
#'   included). Element X contains the \eqn{k_t} values. See the vignette (\code{vignette("dlm")}) in the \pkg{dlm}
#'   package for further info about the specification.
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
buildDLM<-function(iniDLM,Y,k1=0.044) {
  tMax<-dim(Y)[2]
  k<-rep(NA,tMax)
  m<-vector("list", tMax)    # means of posterior - preallocate an empty list of length tMax
  C<-vector("list", tMax)    # covariance matrices of posterior 
  L<-vector("list", tMax)    # covariance matrices of transition probabilities
  FF<-t(iniDLM$FF)           # F def used in W&H book
  for(t in 1:tMax){
    #Prior
    if(t==1){
      a<-iniDLM$GG %*% iniDLM$m0
      R<-iniDLM$GG %*% iniDLM$C0 %*% t(iniDLM$GG) + iniDLM$W
    }
    else{
      a<-iniDLM$GG %*% m[[t-1]] 
      R <-iniDLM$GG %*% C[[t-1]]  %*% t(iniDLM$GG) + iniDLM$W
    }
    # One step forcast
    k[t]<-7*k1/(a[1]^0.25)
    FF[1,2]<-k[t]
    f<-t(FF) %*% a
    Q<-t(FF) %*% R %*% FF + iniDLM$V
    #Posterior (we see Y_t here)
    A<-R %*% FF %*% solve(Q)
    m[[t]] <-a + A %*% (Y[,t]-f)   
    C[[t]] <-R  - A %*% Q %*% t(A) 
    #Compute the covariance matrix for transition probabilities
    L[[t]] <-R-C[[t]] 
  }
  #L[[1]]<-L[[2]]
  
  #iniDLM$m<-m
  iniDLM$C<-C
  #iniDLM$k<-k
  iniDLM$L<-L
  iniDLM$JFF<-matrix(c(0,1,0,0), nrow=2) # F' time dependent in entry (2,1) (dlm package specification way)
  iniDLM$X<-matrix(k,ncol=1) # value of k_t (dlm package specification way)
  return(iniDLM)
}


####################################################################


#' Initialize the parameters for the DGLM 
#' 
#' @param tMax Maximum growth period in weeks.
#' @param nf Sample size to estimate the sample variance in the pen (we suppose this value is constant and the sample is taken from the herd)
#' @param W System variance (we suppose this value is 0 in the DGLM)
#' @param alpha0 Initial shape parameter at time t_0
#' @param beta0 Initial rate parameter at time t_0
#' @param c0 coefficient of G_t in the dglm model (Based on the simulation test, this value is 1) 
#' 
#' @return A list containing all the parameters
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
iniDGLM<-function(tMax=15, nf=35, W=0, alpha0=130, beta0=17, c0=1) {
  model<-list(tMax=tMax)  
  model$nf=35
  model$W=0 #? I consider this number as zero  based on the Anders paper and talking with him. 
  model$alpha0=130
  model$beta0=17
  model$co=1
  return(model)
}

####################################################################
