#SIT optim
#Description
#donde S, son todos los que pueden adquirir la enfermedad
#I, individuos enfermos - los que desarrollaron sintomas 
#T1, individuos que ingresan a tratamiento
require(deSolve)
library(deSolve)
SIT.optim<-function(init.pars,data,initial_values,timepoints,model=c("constant","time.dependent")
                    ,test=FALSE,length.test=10,...){

  sit.model<-switch(model,
                    constant=sit_mod.c,
                    time.dependent=seit_model)
  
 #Optimización de modelos 
   if(test){
    optim.list<-list()
    rmse.vec<-rep(NA,length=length.test)
    predict.list<-list()
    inits.sampl<-seq(-0.01,-10,length=length.test)
                     for(i in 1:length.test){
      if(model=="constant"){
        init.pars<-rep(inits.sampl[i],2)                  #modelo beta constante
        names(init.pars)<-c("beta_1","r")
      }else{init.pars<-rep(inits.sampl[i],3)
      names(init.pars)<-c("beta_0","beta_2","r")
      }
      optim.list[[i]]<-optim(par=init.pars,fn=seit.negll,method="Nelder-Mead"
                             ,datos=data,initial_values=initial_values
                             ,timepoints=timepoints,model=model)
      ith.hats<-exp(optim.list[[i]]$par)
      names(ith.hats)<-names(init.pars)
      ith.out<-lsoda(initial_values,timepoints,sit.model,ith.hats)
      I.rmse<-sqrt(mean((ith.out[,"I"]-data[,"I"])^2))
      T1.rmse<-sqrt(mean((ith.out[,"T1"]-data[,"T1"])^2))
      rmse.vec[i]<-I.rmse+T1.rmse
      predict.list[[i]]<-ith.out
    }
    min.index<-which.min(rmse.vec)
    tmp.optim<-optim.list[[min.index]]
    predicted<-predict.list[[min.index]]
    params.hats<-exp(tmp.optim$par)
    names(params.hats)<-names(init.pars)
    loglik<--(tmp.optim$value)
    BIC.mod<--2*loglik+log(nrow(data))*length(params.hats)
    rmse<-rmse.vec[min.index]
  }else{
    
 #Optimizacion modelo beta dependiente de tiempo  
    
   tmp.optim<-optim(par=init.pars,fn=seit.negll,method="Nelder-Mead"
                   ,datos=data,initial_values=initial_values,timepoints=timepoints, model= model)
    params.hats<-exp(tmp.optim$par)
    names(params.hats)<-names(init.pars)
    loglik<--(tmp.optim$value)
    BIC.mod<--2*loglik+log(nrow(data))*length(params.hats) 
    predicted<-lsoda(initial_values,timepoints,sit.model,params.hats)
    rmse<-sqrt(mean((predicted[,"I"]-data[,"I"])^2))+sqrt(mean((predicted[,"T1"]-data[,"T1"])^2))
  }
  res<-list(Param.hat=params.hats,loglik=loglik,BIC=BIC.mod,Predicted=predicted,RMSE=rmse)
return(res)
  }

#SIT model - constant
#modelo SIT TB
seit_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  I = state_values [2]        # infectious
  T1 = state_values [3]        # recovered
  N = S+I+T1
  
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      beta_1 =  beta_0*(1+beta_2*cos((pi*current_timepoint)/52))
      dS = (-beta_1 * S * (I/N))+lambda*(N)-mu*S
      dI = ((beta_1 *S*(I/N))-(r*I) -(mu*I))
      dT1 = (r*I)-(mu*T1)
      
      # combine results
      results = c (dS,dI, dT1)
      list (results)
    }
  )
}

#beta constante
sit_mod.c = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  I = state_values [2]        # infectious
  T1 = state_values [3]        # recovered
  N = S+I+T1
  
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      dS = (-beta_1 * S * (I/N))+lambda*(N)-mu*S
      dI = ((beta_1 *S*(I/N))-(r*I) -(mu*I))
      dT1 = (r*I)-(mu*T1)
      
      # combine results
      results = c (dS,dI, dT1)
      list (results)
    }
  )
}

# beta 1 en funci?n del tiempo - Time dependent
beta_1.fun<-function(t,beta_0,beta_2){
  beta_0*(1+beta_2*cos((pi*t)/52))
  
}

# Verosimilitud negativa
seit.negll<-function(par,datos,initial_values,timepoints,model=c("constant","time.dependent")){
  sit.model<-switch(model,
                    constant=sit_mod.c,
                    time.dependent=seit_model)
  mu = (1+(4.2/100000))^(1/52)-1
  lambda=(1+(15.8/100000))^(1/52)-1
  parameter_list<-c(exp(par),mu=mu,lambda=lambda)
  inits.mod<-initial_values[c("S","I","T1")]
  output1 = lsodar (inits.mod, timepoints, sit.model, parameter_list)
  I.predict<-output1[-1,"I"]
  T.predict<-output1[-1,"T1"]
  Neg.test<-sum(I.predict<0)+sum(T.predict<0)
  if(Neg.test>0){
    negll<-.Machine$double.xmax
  }else{
  Infectados<-datos[,1]
  Tratamiento<-datos[,2]
  I.loglik<-dpois(Infectados,I.predict,log=TRUE)
  T.loglik<-dpois(Tratamiento,T.predict,log=TRUE)
  negll<--(sum(I.loglik)+sum(T.loglik))
  if(is.nan(negll)|is.infinite(negll)){negll<-.Machine$double.xmax}}
  return(negll)
  
}




