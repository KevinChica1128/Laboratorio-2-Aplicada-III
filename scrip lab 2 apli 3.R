# Kevin García - 1533173
# Alejandro Vargas - 1525953
# Laboratorio 2 - Aplicada III
# Datos Importaciones

install.packages("far")
install.packages("ade4")
#base de datos
importaciones=read.table("clipboard",header=T,row.names =1)
View(importaciones)
#creacion de matrices con datos faltantes
N=length(importaciones$Colombia)*length(importaciones)#total de observaciones
N_1=N*0.05#numero de datos faltantes al 5%
N_2=N*0.1#numero de datos faltantes al 10%
N_3=N*0.15#numero de datos faltantes al 15%
N_4=N*0.2#numero de datos faltantes al 20%
#funcion (a)
remplazardatos=function(a,datos){
  mustcol1=sample(c(1:6),a,replace = T)
  mustrow1=sample(c(1:20),a,replace = T)
  matr=datos
  for(i in 1:a){
    matr[mustrow1[i],mustcol1[i]]=NA
  }
  
  return(matr)  
}
#matriz al 5%
importaciones_1=remplazardatos(N_1,importaciones)
#matriz al 10%
importaciones_2=remplazardatos(N_2,importaciones)
#matriz al 15%
importaciones_3=remplazardatos(N_3,importaciones)
#matriz al 20%
importaciones_4=remplazardatos(N_4,importaciones)

#****Funciones manuales para ACP 
# Función para sd (1/n)
sd2 <- function (x) {
  
  sqrt(sum((x - mean(x))^2) / (length(x)))
  
} 
#funcion para el calculo de la matriz Z
matrZ=function(datos,n,p){
  datos
  #Estandarización
  Z=datos
  for(i in 1:p){
    Z[,i]=(Z[,i]-mean(Z[,i]))/sd2(Z[,i])
  }
  
  return(Z)
}
#funcion para el calculo del ACP
#ACP_manual la entrada l recibe los vaores 1,2 y 3 
#para obtener 1:(correlacion) 2:(valores y vectores) 3:(componentes)
#la entrada tipo recibe valores 0:(VARIABES) 1:(INDIVIDUOS)
ACP_manual=function(datos,n,p,l,tipo){
  Z=matrZ(datos,n,p)
  if(tipo==1){
    R=(t(Z)%*%as.matrix(Z))*(1/n)
    dv=eigen(R)
    CP=as.matrix(Z)%*%dv$vectors
  }
  if(tipo==0){
    N=diag(sqrt(1/n),n,n)
    R=N%*%as.matrix(Z)%*%t(Z)%*%N
    dv=eigen(R)
    CP=t(Z)%*%N%*%dv$vectors
  }
  if(l==1){
    return(R)
  }
  if(l==2){
    return(dv)
  }
  if(l==3){
    return(CP)
  }
  
}
#****ACP MANUAL DATOS COMPLETOS
#ACP (VARIABLES)
ACP_manual(importaciones,20,6,2,0) #valores y vectores propios
ACP_manual(importaciones,20,6,3,0) #componentes
#ACP (INDIVIDUOS)
ACP_manual(importaciones,20,6,2,1) #valores y vectores propios
ACP_manual(importaciones,20,6,3,1) #componentes
#****ACP MANUAL DATOS COMPLETOS (NIPALS)
#La entrada c recibe 0:(componentes principales) y 1:(vectores propios)
#la entrada tipo recibe valores 0:(VARIABES) 1:(INDIVIDUOS)
ACP_NIPALS=function(datos,n,p,c,tipo){
  X=matrZ(datos,n,p)
  if(tipo==1){
    X0=X
    T=matrix(NA,n,p)
    P=matrix(NA,p,p)
  }
  if(tipo==0){
    N=diag(sqrt(1/n),n,n)
    X0=t(X)%*%N
    T=matrix(NA,p,p)
    P=matrix(NA,n,p)
  }
  for(h in 1:p){
    t1=as.matrix(X0[,1])
    for(i in 1:100){
      P11=(t(X0)%*%t1)/(as.numeric(t(t1)%*%t1))
      nP11=as.numeric(t(P11)%*%P11)
      P1=1/sqrt(nP11)*P11
      t1=as.matrix(X0)%*%P1
    }
    T[,h]=t1
    P[,h]=P1
    X1=X0-t1%*%t(P1)
    X0=X1
  }
  if(c==1){
    return(T)
  }
  if(c==0){
    return(P)
  }
}
#**ACP-NIPALS(INDIVIDUOS) MANUAL
P=ACP_NIPALS(importaciones,20,6,1,1) #componentes principales
T=ACP_NIPALS(importaciones,20,6,0,1) #vectores propios
valp=diag(t(P)%*%as.matrix(P)*(1/20)) #valores propios
#**ACP-NIPALS(VARIABLES) MANUAL
PV=ACP_NIPALS(importaciones,20,6,1,0) #componentes principales
TV=ACP_NIPALS(importaciones,20,6,0,0) #vectores propios
valpV=diag(t(PV)%*%as.matrix(PV)*(1/20)) #valores propios
#****ACP NO MANUAL (ADE4)
library(ade4)
#**ACP(INDIVIDUOS)
acpnipals=dudi.pca(importaciones,scannf=F,nf=6)
acpnipals$li #componentes primeras 2 dim
acpnipals$c1 #vectores propios primeras 2 dim
acpnipals$eig #valores propios
#**ACP(VARIABLES)
#?????



#ACP MANUAL DATOS FALTANTES (NIPALS)
ACP_NIPALSNA=function(datos,n,p){
  library("far")
  njm=colMeans(datos,na.rm=TRUE); njs=colSums(datos,na.rm=TRUE)
  nj=njs/njm
  X=scale(datos)*sqrt(nj/(nj-1))
  X0=X
  P=matrix(NA,p,p); T=matrix(NA,n,p)
  P1i=matrix(NA,p,1)
  for(h in 1:p){
    
    t1=X0[,1]
    
    for(i in 1:100){
      
      for(j in 1:p){
        
        j1 =na.omit(cbind(X0[,j],t1))
        P1i[j]=sum(j1[,1]*j1[,2])/sum(j1[,2]^2)
      }
      P[,h]=P1i
      Portn=orthonormalization(P[,1:h]); P1=Portn[,h]
      for(i in 1:n){
        
        i1=na.omit(cbind(X0[i,],P1))
        t1[i]=sum(i1[,1]*i1[,2])/sum(i1[,2]^2)
      }
      T[,h]=t1
      Tortg=orthonormalization(T[,1:h],norm=FALSE); t1=Tortg[,h]
    } # end i
    P[,h]=P1
    T[,h]=t1
    X1=X0-t1%*%t(P1); X0=X1
  } # end h
  
  L=diag(t(T)%*%T)/n
  r.nipNA=list(T,P,L); return(r.nipNA)
} 
#**ACP martriz al 5%
ACP_NIPALSNA(importaciones_1,20,6)
#**ACP martriz al 10%
ACP_NIPALSNA(importaciones_2,20,6)
#**ACP martriz al 15%
ACP_NIPALSNA(importaciones_3,20,6)
#**ACP martriz al 20%
ACP_NIPALSNA(importaciones_4,20,6)
#****ACP NO MANUAL (ADE4)
library("ade4")
#**ACP martriz al 5%
acpnipals5NA<-nipals(importaciones_1,nf=6)
acpnipals5NA$li #Coordenadas de los individuos
acpnipals5NA$co #Coordenadas de las variables
acpnipals5NA$tab #Matriz con los valores faltantes estimados
acpnipals5NA$eig #Valores propios
acpnipals5NA$c1 #Vectores propios
par(mfrow=c(1,2))
s.label(acpnipals$li,boxes=FALSE)
s.label(acpnipals5NA$li,boxes=FALSE) ## Nube de individuos (matriz con NA)
s.corcircle(acpnipals5NA$co) ## Observe lo que pasa con el circulo de correlacciones
#**ACP martriz al 10%
acpnipals10NA<-nipals(importaciones_2,6)
acpnipals10NA$li #Coordenadas de los individuos
acpnipals10NA$co #Coordenadas de las variables
acpnipals10NA$tab #Matriz con los valores faltantes estimados
acpnipals10NA$eig #Valores propios
acpnipals10NA$c1 #Vectores propios
par(mfrow=c(1,2))
s.label(acpnipals$li,boxes=FALSE)
s.label(acpnipals10NA$li,boxes=FALSE) ## Nube de individuos (matriz con NA)
s.corcircle(acpnipals10NA$co) ## Observe lo que pasa con el circulo de correlacciones
#**ACP martriz al 15%
acpnipals15NA<-nipals(importaciones_3,6)
acpnipals15NA$li #Coordenadas de los individuos
acpnipals15NA$co #Coordenadas de las variables
acpnipals15NA$tab #Matriz con los valores faltantes estimados
acpnipals15NA$eig #Valores propios
acpnipals15NA$c1 #Vectores propios
par(mfrow=c(1,2))
s.label(acpnipals$li,boxes=FALSE)
s.label(acpnipals15NA$li,boxes=FALSE) ## Nube de individuos (matriz con NA)
s.corcircle(acpnipals15NA$co) ## Observe lo que pasa con el circulo de correlacciones
#**ACP martriz al 20%
acpnipals20NA<-nipals(importaciones_4,6)
acpnipals20NA$li #Coordenadas de los individuos
acpnipals20NA$co #Coordenadas de las variables
acpnipals20NA$tab #Matriz con los valores faltantes estimados
acpnipals20NA$eig #Valores propios
acpnipals20NA$c1 #Vectores propios
par(mfrow=c(1,2))
s.label(acpnipals$li,boxes=FALSE)
s.label(acpnipals20NA$li,boxes=FALSE) ## Nube de individuos (matriz con NA)
s.corcircle(acpnipals20NA$co) ## Observe lo que pasa con el circulo de correlacciones

#--------------------------------------------------------#
#ACP EM 
#FUNCIÓN PARA REEMPLAZAR DATOS FALTANTES POR VALORES INICIALES
reempNA=function(datos,n,p){
  X0=datos
  for (i in 1:p) {
    for (j in 1:n) {
      if(is.na(X0[j,i])=='TRUE'){
        X0[j,i]=mean(X0[,i],na.rm = 'TRUE')
    }
  }
  }
  return(X0)
}

#Función para estandarizar la matriz sin tener en cuenta NAS y con la desviación sobre n
library('plyr')
matrZNA=function(datos,n,p){
  datos
  #Estandarización
  Z=datos
  for(i in 1:p){
    Z[,i]=(Z[,i]-mean(Z[,i],na.rm = TRUE))/((sqrt((n-count(Z[,i],vars = 'NA')$freq -1)/(n-count(Z[,i],vars = 'NA')$freq)))*sd(Z[,i],na.rm = TRUE))
      }
  return(Z)
}

#Función algoritmo ACP EM
ACPEM=function(datos,n,p){
  datos
  L=0
  Z=matrZNA(datos,n,p)
  Z0=reempNA(Z,n,p)
  for (L in 0:500) {
    VL=ACP_manual(Z0,n,p,2,1)$vectors  #Matriz de vectores propios
    CL=ACP_manual(Z0,n,p,3,1)          #Componentes principales
    Z0=CL%*%VL                         #Reconstitucion de la matriz estandarizada
    }
  result=list(VL,CL,Z0,Z)
  return(result)
}


#**ACP EM matriz al 5%
ACPEM(importaciones_1,20,6)
#**ACP EM matriz al 10%
ACPEM(importaciones_2,20,6)
#**ACP EM matriz al 15%
ACPEM(importaciones_3,20,6)
#**ACP EM matriz al 20%
ACPEM(importaciones_4,20,6)

mean(importaciones_1[,1],na.rm = TRUE)
sqrt((20-count(importaciones_1[,1],vars = 'NA')$freq -1)/(20-count(importaciones_1[,1],vars = 'NA')$freq))*sd(importaciones_1[,1],na.rm = TRUE)
sqrt((20-count(importaciones_1[,1],vars = 'NA')$freq)/(20-count(importaciones_1[,1],vars = 'NA')$freq -1))*sd2(importaciones[,1])

#ACPEM (ACP iterativo) con función de R 
library(missMDA)


#PCAEM Matriz del 5% de NA
## Validación cruzada para elegir el número de ejes..
?estim_ncpPCA
nb5 <- estim_ncpPCA(importaciones_1,ncp.max=6,method = 'EM')
## Imputación
imp5NA<-imputePCA(data.frame(importaciones_1),ncp=2)
imp5NA$completeObs
?imputePCA
## ACP Impute ade4
pca.em5 <- dudi.pca(imp5NA$completeObs,nf=3,scannf = FALSE)
pca.em5$li #Componentes principales individuos
pca.em5$co #Componentes principales variables
pca.em5$c1 #Vectores propios
pca.em5$eig #Valores propios
pca.em5$tab #Matriz modificada

#PCAEM Matriz del 10% de NA
## Validación cruzada para elegir el número de ejes..
?estim_ncpPCA
nb10 <- estim_ncpPCA(importaciones_2,ncp.max=6)
## Imputación
imp10NA<-imputePCA(data.frame(importaciones_2),ncp=2)
imp10NA$completeObs
?imputePCA
## ACP Impute ade4
pca.em10 <- dudi.pca(imp10NA$completeObs,nf=3,scannf = FALSE)
pca.em10$li #Componentes principales individuos
pca.em10$co #Componentes principales variables
pca.em10$c1 #Vectores propios
pca.em10$eig #Valores propios
pca.em10$tab #Matriz modificada

#PCAEM Matriz del 15% de NA
## Validación cruzada para elegir el número de ejes..
?estim_ncpPCA
nb15 <- estim_ncpPCA(importaciones_3,ncp.max=6,method = 'EM')
## Imputación
imp15NA<-imputePCA(data.frame(importaciones_3),ncp=2)
imp15NA$completeObs
?imputePCA
## ACP Impute ade4
pca.em15 <- dudi.pca(imp15NA$completeObs,nf=3,scannf = FALSE)
pca.em15$li #Componentes principales individuos
pca.em15$co #Componentes principales variables
pca.em15$c1 #Vectores propios
pca.em15$eig #Valores propios
pca.em15$tab #Matriz modificada

#PCAEM Matriz del 20% de NA
## Validación cruzada para elegir el número de ejes..
?estim_ncpPCA
nb20 <- estim_ncpPCA(importaciones_4,ncp.max=6,method = 'EM')
## Imputación
imp20NA<-imputePCA(data.frame(importaciones_4),ncp=2)
imp20NA$completeObs
?imputePCA
## ACP Impute ade4
pca.em20 <- dudi.pca(imp20NA$completeObs,nf=3,scannf = FALSE)
pca.em20$li #Componentes principales individuos
pca.em20$co #Componentes principales variables
pca.em20$c1 #Vectores propios
pca.em20$eig #Valores propios
pca.em20$tab #Matriz modificada


#Comparación poder descriptivo 2 primeros ejes:
sum(acpnipals$eig[1:2])/sum(acpnipals$eig) #Completa
sum(acpnipals5NA$eig[1:2])/sum(acpnipals5NA$eig) #Nipals 5%
sum(acpnipals10NA$eig[1:2])/sum(acpnipals10NA$eig) #Nipals 10%
sum(acpnipals15NA$eig[1:2])/sum(acpnipals15NA$eig) #Nipals 15%
sum(acpnipals20NA$eig[1:2])/sum(acpnipals20NA$eig) #Nipals 20%
sum(pca.em5$eig[1:2])/sum(pca.em5$eig) #ACP-EM 5%
sum(pca.em10$eig[1:2])/sum(pca.em10$eig) #ACP-EM 10%
sum(pca.em15$eig[1:2])/sum(pca.em15$eig) #ACP-EM 15%
sum(pca.em20$eig[1:2])/sum(pca.em20$eig) #ACP-EM 20%

#Correlacion componente 1 completos vs faltantes
abs(cor(acpnipals$li[1],acpnipals5NA$li[,1])) #5% nipals
abs(cor(acpnipals$li[1],acpnipals10NA$li[,1])) #10% nipals
abs(cor(acpnipals$li[1],acpnipals15NA$li[,1])) #15% nipals
abs(cor(acpnipals$li[1],acpnipals20NA$li[,1])) #20% nipals
abs(cor(acpnipals$li[1],pca.em5$li[,1])) #5% ACP-EM
abs(cor(acpnipals$li[1],pca.em10$li[,1])) #10% ACP-EM
abs(cor(acpnipals$li[1],pca.em15$li[,1])) #15% ACP-EM
abs(cor(acpnipals$li[1],pca.em20$li[,1])) #20% ACP-EM

#Comparación nube de individuos
x11()
s.label(acpnipals$li,boxes=FALSE)
x11()
par(mfrow=c(1,2))
s.label(acpnipals20NA$li,boxes=FALSE) ## Nube de individuos Nipals 20% NA 
s.label(pca.em20$li,boxes=FALSE) ## Nube de individuos ACP-EM 20% NA

#Comparación nube de variables
rbind()
x11()
s.corcircle(-acpnipals$co)

circulonipals20NA<-cbind(acpnipals20NA$co[,1],-acpnipals20NA$co[,2])
circuloacpem20NA<-cbind(-pca.em20$co[,1],-pca.em20$co[,2])
rownames(circuloacpem20NA)<-c("Colombia","Brasil","Chile","Argentina","Ecuador","Peru")
x11()
par(mfrow=c(1,2))
s.corcircle(circulonipals20NA) ## Nube de variables Nipals 20% NA
s.corcircle(circuloacpem20NA) ## Nube de variables ACP-EM 20% NA

#Ortogonalidad en las componentes
t(acpnipals$li[,1])%*%acpnipals$li[,2] #ACP completo
t(acpnipals5NA$li[,1])%*%acpnipals5NA$li[,2] #Nipals 5%
t(acpnipals10NA$li[,1])%*%acpnipals10NA$li[,2] #Nipals 10%
t(acpnipals15NA$li[,1])%*%acpnipals15NA$li[,2] #Nipals 15%
t(acpnipals20NA$li[,1])%*%acpnipals20NA$li[,2] #Nipals 20%
t(pca.em5$li[,1])%*%pca.em5$li[,2] #ACP-EM 5%
t(pca.em10$li[,1])%*%pca.em10$li[,2] #ACP-EM 10%
t(pca.em15$li[,1])%*%pca.em15$li[,2] #ACP-EM 15%
t(pca.em20$li[,1])%*%pca.em20$li[,2] #ACP-EM 20%

#Ortonormalidad en los vectores propios
#Norma
sqrt(acpnipals$c1[,1]%*%acpnipals$c1[,1]) #Vector propio 1 completa
sqrt(acpnipals$c1[,2]%*%acpnipals$c1[,2]) #Vector propio 2 completa
sqrt(acpnipals5NA$c1[,1]%*%acpnipals5NA$c1[,1]) #Vector propio 1 nipals 5%
sqrt(acpnipals5NA$c1[,2]%*%acpnipals5NA$c1[,2]) #Vector propio 2 nipals 5%
sqrt(acpnipals10NA$c1[,1]%*%acpnipals10NA$c1[,1]) #Vector propio 1 nipals 10%
sqrt(acpnipals10NA$c1[,2]%*%acpnipals10NA$c1[,2]) #Vector propio 2 nipals 10%
sqrt(acpnipals15NA$c1[,1]%*%acpnipals15NA$c1[,1]) #Vector propio 1 nipals 15%
sqrt(acpnipals15NA$c1[,2]%*%acpnipals15NA$c1[,2]) #Vector propio 2 nipals 15%
sqrt(acpnipals20NA$c1[,1]%*%acpnipals20NA$c1[,1]) #Vector propio 1 nipals 20%
sqrt(acpnipals20NA$c1[,2]%*%acpnipals20NA$c1[,2]) #Vector propio 2 nipals 20%
sqrt(pca.em5$c1[,1]%*%pca.em5$c1[,1]) #Vector propio 1 ACP-EM 5%
sqrt(pca.em5$c1[,2]%*%pca.em5$c1[,2]) #Vector propio 2 ACP-EM 5%
sqrt(pca.em10$c1[,1]%*%pca.em10$c1[,1]) #Vector propio 1 ACP-EM 10%
sqrt(pca.em10$c1[,2]%*%pca.em10$c1[,2]) #Vector propio 2 ACP-EM 10%
sqrt(pca.em15$c1[,1]%*%pca.em15$c1[,1]) #Vector propio 1 ACP-EM 15%
sqrt(pca.em15$c1[,2]%*%pca.em15$c1[,2]) #Vector propio 2 ACP-EM 15%
sqrt(pca.em20$c1[,1]%*%pca.em20$c1[,1]) #Vector propio 1 ACP-EM 20%
sqrt(pca.em20$c1[,2]%*%pca.em20$c1[,2]) #Vector propio 2 ACP-EM 20%
#ortogonalidad
t(acpnipals$c1[,1])%*%acpnipals$c1[,2] #ACP completo
t(acpnipals5NA$c1[,1])%*%acpnipals5NA$c1[,2] #Nipals 5%
t(acpnipals10NA$c1[,1])%*%acpnipals10NA$c1[,2] #Nipals 10%
t(acpnipals15NA$c1[,1])%*%acpnipals15NA$c1[,2] #Nipals 15%
t(acpnipals20NA$c1[,1])%*%acpnipals20NA$c1[,2] #Nipals 20%
t(pca.em5$c1[,1])%*%pca.em5$c1[,2] #ACP-EM 5%
t(pca.em10$c1[,1])%*%pca.em10$c1[,2] #ACP-EM 10%
t(pca.em15$c1[,1])%*%pca.em15$c1[,2] #ACP-EM 15%
t(pca.em20$c1[,1])%*%pca.em20$c1[,2] #ACP-EM 20%

#Imputación
matrZ(importaciones,20,6) #Matriz original estandarizada
acpnipals20NA$tab
pca.em20$tab

dudi.pca
nipals
