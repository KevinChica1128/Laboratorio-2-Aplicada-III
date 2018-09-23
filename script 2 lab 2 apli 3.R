# Kevin García - 1533173
# Alejandro Vargas - 1525953
# Laboratorio 2 - Aplicada III
# Datos Importaciones

library(ade4) ## Contiene NIPALS
library(missMDA) ## Contiene ACP-EM (Acp iterativo)
library(FactoMineR)
library(factoextra)
#Datos
importaciones=read.table("clipboard",header=T,row.names =1)
#ACP datos completos
acp1 <- dudi.pca(importaciones,scannf=FALSE,nf=9)
biplot(acp1)
fviz_pca(acp1)
acp1$eig ## Valores Propios
acp1$li ## Componentes (Coordenadas en Rp)
acp1$l1 ## vectores propios v
acp1$co ## Coordenadas en Rn
acp1$c1 ## vectores propios u 

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

#FUNCIÓN NIPALS DATOS COMPLETOS
fnipals <- function(Xi)
{
  a <- qr(Xi)$rank	# rango de Xi
  
  p <- ncol(Xi); n <- nrow(Xi)
  Xo <- sqrt(n/(n-1))*scale(Xi)
  
  T <- matrix(1,n,a) # a compals 
  P <- matrix(1,p,a)
  
  for(h in 1:a)
  {
    t1 <- as.matrix(Xo[,1])
    
    for(e in 1:100)
    {
      P1i <- (t(Xo)%*%t1)/sum(t1^2)
      nP1i <- sqrt(sum(P1i^2))
      P1 <- P1i/nP1i		# vect unitario
      t1 <- Xo%*%P1
    }
    
    T[,h] <- t1
    P[,h] <- P1
    X1 <- Xo - t1%*%t(P1)	# deflacta
    Xo <- X1
  }
  
  L <- diag(t(T)%*%T)/n
  
  r.nip <- list(T,P,L)
  return(r.nip)
  
}	# End, Nipals datos completos de rango a.


nipals2 <- fnipals(datacom1)

T <- nipals2[[1]]  ## Componentes Principales

cbind(T[,1:2],acp1$li[,1:2]) ## Ok

P <- nipals2[[2]]  ## Vectores propios en Rp (U)

cbind(P[,1:2],acp1$c1[,1:2]) ## Ok

L <- nipals2[[3]]  ## Valores propios

rbind(L,acp1$eig)  ## ok
