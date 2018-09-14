# Kevin García - 1533173
# Alejandro Vargas - 1525953
# Laboratorio 2 - Aplicada III
# Datos Importaciones

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
valpV=diag(t(P)%*%as.matrix(P)*(1/20)) #valores propios
#****ACP NO MANUAL (ADE4)
library(ade4)
#**ACP(INDIVIDUOS)
acpnipals=dudi.pca(importaciones,scannf=F,nf=2)
acpnipals$li #componentes primeras 2 dim
acpnipals$c1 #vectores propios primeras 2 dim
acpnipals$eig #valores propios
#**ACP(VARIABLES)
#?????