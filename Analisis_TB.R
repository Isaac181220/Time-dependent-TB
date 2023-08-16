############################
#Total
#Cargar previamente los datos (Data_TB_atlantico) y funciones (SIT_functions,bootstrap)
############################
#obs.all corresponde a los datos observados
#best.molde.list lista de modelos que mejor se ajustan
#best.model.boot lista de modelos e intervalos de confianza para los parametros utilizando
#técnicas Bootstrap
require(deSolve)
source("SIT_functions_2.R")
source("Bootstrap_12-08-2022")

enfermos_T<-casos.df$Enfermos[53:415]
tratamiento_T<-casos.df$Tratamiento[53:415]
year.T<-as.numeric(substring(casos.df$Semana_ID[53:415],1,4))
obs.all=cbind(I=enfermos_T,T1=tratamiento_T)
sampling.year<-c(2015:2021)
best.model.list<-list()
best.model.boot<-list()

#mu.vec es el vector de parametros para las dinamicas vitales, nacimientos 

mu.vec<-c((1+(4.2/100000))^(1/52)-1,(1+(4/100000))^(1/52)-1,
          (1+(4/100000))^(1/52)-1, (1+(5.1/100000))^(1/52)-1,
          (1+(4.7/100000))^(1/52)-1,(1+(6.63/100000))^(1/52)-1,
          (1+(5.48/100000))^(1/52)-1)

#lambda.vec es el vector de parametros para las dinamicas vitales, defunciones 
lambda.vec<-c((1+(15.8/100000))^(1/52)-1,(1+(16/100000))^(1/52)-1,
              (1+(14.9/100000))^(1/52)-1,(1+(16.2/100000))^(1/52)-1,
              (1+(16/100000))^(1/52)-1,(1+(16.2/100000))^(1/52)-1,
              (1+(16.2/100000))^(1/52)-1)

#n.vec vector de población total 
n.vec<-c(1235506,1259491,1275681,1329198,1329198,1447878,1472209)

#suscep.vec vector de individuos susceptibles para cada año de estudio
suscep.vec<-c(1235381,1259044,1274977,1328193,1393746,1446230,1257643)

infected.pred<-(NULL)
tratament.pred<-(NULL)
boot=T

for (i in 1:7) {
  
  ith.year<-sampling.year[i]
  
  #Modelo - Total 
  
  obs.year_i=obs.all[year.T == ith.year,]
  
  
  #Valores iniciales ajustar a la primera semana
  #epidemiologica que exista en los datos
  N = n.vec[i]  # Total population   
  W = suscep.vec[i]        # susceptible hosts
  I = obs.year_i[1,1]           # infectious hosts
  T1 = obs.year_i[1,2]           # treatment hosts
  
  
  # reco hosts
  # Tasa mortalidad para colombia anual = 4.7*10^5
  # Tasa natalidad anual para colombia = 16*10^5
  # Para convertir tasa anual a tasa semanal
  # (1+Tasa Anual)^(1/#Semanas Anho)-1
  mu <-mu.vec[i]
  lambda<-lambda.vec[i]
  
  #Valores iniciales
  valores_iniciales = c (S = W, I, T1 )
  
  #Periodo
  timepoints = seq (2,52, by=1)
  if(i==7){
    timepoints=seq(2,51,by=1)
  }
  
  # Optimizaci?n del modelo
  inits<-c(-1,-1)
  names(inits)<-c("beta_1","r")
  model.constant<-SIT.optim(init.pars = inits,initial_values = valores_iniciales,data=obs.year_i[timepoints,],timepoints = timepoints,model = "constant",test=T)
  
  inits<-c(-1,-1,-1)
  names(inits)<-c("beta_0","beta_2","r")
  model.time.dep<-SIT.optim(init.pars = inits,initial_values = valores_iniciales,data=obs.year_i[timepoints,],timepoints = timepoints,model = "time.dependent",test=T)
  
  BIC.vec<-c(model.constant$BIC,model.time.dep$BIC)
  best.model<-which.min(BIC.vec)
  if (best.model==1){
    best.model.list[[i]]<-model.constant
    model.name<-"constant"
  }else{
    best.model.list[[i]]<-model.time.dep
    model.name<-"time.dependent"
  }
  if(boot){
    best.model.boot[[i]]<-SIT.boot(data = best.model.list[[i]]$Predicted,
                                   initial_values = valores_iniciales, nboot = 200,
                                   model = model.name)
  }
    
  print(ith.year)
  infected.pred<-c(infected.pred,c(I,best.model.list[[i]]$Predicted[,"I"]))
  tratament.pred<-c(tratament.pred,c(T1,best.model.list[[i]]$Predicted[,"T1"]))
}

n.par<-sum(sapply(best.model.list,function(x)length(x$Param.hat)))
log.lik<-sum(sapply(best.model.list, function(x)x$loglik))
n.all<-2*nrow(obs.all)
BIC.all<-log(n.all)*n.par-2*log.lik


#Modelo - Total 
years.vec<-2015:2021
tiff("Fig1.tif",width = 7.5,height = 8.75,units = "in",res = 300,type = "cairo",family = "times",compression = "lzw")
par(mfrow=c(7,2),mar=c(3,3,1,1),mgp=c(1.75,0.5,0),tcl=-0.3)
for (i in 1:7) {
  ith.year<-sampling.year[i]
  my.ylab<-paste("Infected")
  my.xlab<-paste("Weeks")
  plot(obs.all[year.T == ith.year,1],ylab="Infected",xlab="Epidemiologic Weeks",pch=19,bty="l")
  points(best.model.boot[[i]]$infected[2,],type = "l",col="red")
  points(best.model.boot[[i]]$infected[1,],type ="l", col="blue",lty=2)
  points(best.model.boot[[i]]$infected[3,],type ="l", col="blue",lty=2)
  mtext(years.vec[i],adj=0)
  
  plot(obs.all[year.T == ith.year,2],ylab="Treatment",xlab="Epidemiologic Week",pch=19,bty="l")
  points(best.model.boot[[i]]$tratament[2,],type = "l",col="red")
  points(best.model.boot[[i]]$tratament[1,],type ="l", col="blue",lty=2)
  points(best.model.boot[[i]]$tratament[3,],type ="l", col="blue",lty=2)
  
}
dev.off()






#Modelo conjunto todos los años
par(mfrow=c(1,2))
plot(obs.baq_T[,1],ylab="Infectados 2015-2021",xlab="Semanas 2015-2021",pch=19,bty="l")
points(infected.pred,type = "l",col="red")
plot(obs.baq_T[,2],ylab="Tratamiento 2015-2021",xlab="Semana 2015-2021",pch=19,bty="l")
points(tratament.pred,type = "l",col="red")
