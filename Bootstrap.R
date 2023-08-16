###Bootstrap parametrico estimación del CI###
#data, datos observados
#intial_values valores iniciales del modelo
#nboot, numero de interacciones a realizar
#model tipo de modelo a emplear

SIT.boot<-function(data,initial_values,nboot,model=c("constant","time.dependent")){

  if(model=="constant"){
    npar=2
    inits<-c(-1,-1)
    names(inits)<-c("beta_1","r")
  }  else{
    npar=3
    inits<-c(-1,-1,-1)
    names(inits)<-c("beta_0","beta_2","r")
  }
ndat<-nrow(data)
boot.par<-matrix(NA,nrow = nboot,ncol=npar)
boot.pred.I<-matrix(NA,nrow = nboot,ncol=ndat)
boot.pred.T1<-matrix(NA,nrow = nboot,ncol=ndat)

for (i in 1:nboot) {
  r.infected<-rpois(ndat,data[,"I"])
  r.removed<-rpois(ndat,data[,"T1"])
  
  obs.r<-cbind(r.infected,r.removed)
  colnames(obs.r)<-c("I","T1")
  
  
  model.optim<-SIT.optim(init.pars = inits,initial_values = valores_iniciales,data=obs.r,timepoints = timepoints,model = model,test=T)
  
  boot.par[i,]<-model.optim$Param.hat
  boot.pred.I[i,]<-model.optim$Predicted[,"I"]
  boot.pred.T1[i,]<-model.optim$Predicted[,"T1"]
}
#cuantiles para el intervalo de confianza 
par.CI<-apply(boot.par, 2, quantile,probs=c(0.025,0.5,0.975))
par.CI.I<-apply(boot.pred.I, 2, quantile,probs=c(0.025,0.5,0.975))
par.CI.T1<-apply(boot.pred.T1, 2, quantile,probs=c(0.025,0.5,0.975))
results<-list(parameter=par.CI,infected=par.CI.I,tratament=par.CI.T1)
return(results)
}


