mixture.gibbs=function(nmat,Nmat,ngroup,ntransect,trans.id,
                       ngibbs,burnin,gamma1){
    
    #useful pre-calculated quantities 
    nloc=nrow(nmat)
    nheight=ncol(nmat)
    N.minus.n=Nmat-nmat

    #initial parameter values
    z=sample(1:ngroup,size=nloc,replace=T)
    theta=matrix(1/ngroup,ntransect,ngroup)
    tmp=runif(ngroup*nheight)
    pi1=matrix(tmp,ngroup,nheight)

    #to store results from gibbs sampler
    store.pi1=matrix(NA,ngibbs,nheight*ngroup)
    store.theta=matrix(NA,ngibbs,ngroup*ntransect)
    # store.z=matrix(NA,ngibbs,nloc) #too big
    store.logl=rep(NA,ngibbs)
    
    #run gibbs sampler
    norder=50
    max.llk=-Inf

    for (i in 1:ngibbs){
        print(i)
        
        #sample group allocation vector z
        ltheta=log(theta)
        lpi1=log(pi1)
        l1.minus.pi1=log(1-pi1)
        #z=true.z$z
        z=samplez(ltheta=ltheta,
                  nmat=nmat,
                  Nminusn=N.minus.n,
                  lpi1=lpi1,
                  l1minuspi1=l1.minus.pi1,
                  ngroup=ngroup,
                  nloc=nloc,
                  nheight=nheight,
                  randu=runif(nloc),
                  TransID=trans.id-1)

        #re-order groups if necessary
        if (i%%norder==0 & i<burnin){
            med=apply(theta,2,mean)
            order1=order(med,decreasing=T)
            theta=theta[,order1]
            pi1=pi1[order1,]
            znew=rep(NA,nloc)
            for (j in 1:ngroup){
                cond=z==order1[j]
                znew[cond]=j
            }
            z=znew
        }
        
        #sample pi1
        pi1=update.pi1(z=z,
                       nmat=nmat,
                       N.minus.n=N.minus.n,
                       ngroup=ngroup,
                       nheight=nheight)

        #calculate the number of locations in each transect 
        #assigned to each group
        ntk=calc_ntk(z=z-1,
                     ngroup=ngroup,nloc=nloc,ntransect=ntransect,
                     TransID=trans.id-1)
        
        #sample theta
        theta=update.theta(ntk=ntk,ngroup=ngroup,gamma1=gamma1,
                           ntransect=ntransect)

        #get loglikelihood
        logl=calc.llk(z=z,nmat=nmat,Nmat=Nmat,
                      ngroup=ngroup,nheight=nheight,
                      pi1=pi1,nloc=nloc)
        
        #store if MLE
        if (logl>max.llk & i>burnin){
          max.llk=logl
          MLE.pi1=pi1
          MLE.theta=theta
          MLE.z=z
          MLE.iter=i
        }
        
        #store results
        store.logl[i]=logl
        store.pi1[i,]=pi1
        store.theta[i,]=theta
    }
    
    #output MCMC results after discarding the burn-in phase
    seq1=burnin:ngibbs
    list(pi1=store.pi1[seq1,],
         theta=store.theta[seq1,],
         logl=store.logl,
         MLE.pi1=MLE.pi1,
         MLE.theta=MLE.theta,
         MLE.z=MLE.z,
         MLE.iter=MLE.iter,
         last.z=z)
}