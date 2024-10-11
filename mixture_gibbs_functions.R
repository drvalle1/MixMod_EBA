#Samples z

update.z=function(nmat,Nmat,N.minus.n,pi1,theta,ngroup,nloc,z){
    #pre-calculate some useful quantities
    log.theta=log(theta)
    log.pi1=log(pi1)
    log.1.minus.pi1=log(1-pi1)
    
    #calculate lprob
    lprob=array(NA,dim=c(nloc,nheight,ngroup))
    for (j in 1:nheight){
      for (k in 1:ngroup){
        lprob[,j,k]=nmat[,j]*log.pi1[k,j]+
                    N.minus.n[,j]*log.1.minus.pi1[k,j]
      }
    }
    
    #calculate max for each loc and height
    soma=apply(lprob,c(1,3),sum)
    max1=matrix(apply(soma,1,max),nloc,ngroup)
    soma1=soma-max1
    prob=soma1/rowSums(soma1)
    
    #sample z
    z=rep(NA,nloc)
    for (i in 1:nloc){
      lprob=rep(NA,ngroup)
      for (k in 1:ngroup) lprob[k]=sum(nmat[i,]*log.pi1[k,]+N.minus.n[i,]*log.1.minus.pi1[k,])+log.theta[k]
      max1=max(lprob)
      lprob1=lprob-max1
      prob1=exp(lprob1)
      prob2=prob1/sum(prob1)
      tmp=rmultinom(1,size=1,prob2)
      z[i]=which(tmp==1)
    }
    z  
}

#--------------------------------------------
#Samples theta and v parameters

update.theta=function(ntk,ngroup,gamma1,ntransect){

    #sample v from a beta distribution
    n.greater.k=t(apply(ntk[,ngroup:1],1,cumsum))[,ngroup:1]
    
    #get theta from v1 using the stick-breaking equation 
    theta=v1=matrix(NA,ntransect,ngroup)
    tmp=rep(1,ntransect)
    for (i in 1:(ngroup-1)){
        v1[,i]=rbeta(ntransect,ntk[,i],n.greater.k[,i+1]+gamma1)
        theta[,i]=v1[,i]*tmp
        tmp=tmp*(1-v1[,i])
    }
    theta[,ngroup]=tmp    

    #output results
    theta
}

#--------------------------------------------
#Samples theta

update.pi1=function(z,nmat,N.minus.n,ngroup,nheight){
  pi1=matrix(NA,ngroup,nheight)
  for (k in 1:ngroup){
    cond=z==k
    soma=sum(cond)
    if (soma>1){
      nmat1=colSums(nmat[z==k,])
      N.minus.n1=colSums(N.minus.n[z==k,])
    }  
    if (soma==1){
      nmat1=nmat[z==k,]
      N.minus.n1=N.minus.n[z==k,]
    }  
    if (soma==0){
      nmat1=N.minus.n1=rep(0,nheight)      
    } 
    pi1[k,]=rbeta(nheight,nmat1+1,N.minus.n1+1)
  }
  pi1
}

#--------------------------------------------
#Calculate llk

calc.llk=function(z,nmat,Nmat,ngroup,nheight,pi1,nloc){
  llk=matrix(NA,nloc,nheight)
  for (k in 1:ngroup){
    cond=z==k
    for (j in 1:nheight){
      llk[cond,j]=dbinom(nmat[cond,j],size=Nmat[cond,j],prob=pi1[k,j],log=T)  
    }
  }
  sum(llk)
}