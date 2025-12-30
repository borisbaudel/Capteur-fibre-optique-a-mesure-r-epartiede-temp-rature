# Photodetector_Response_Deconvolution_with_Tikhonov_regularisaton_v1.R

#
# 
# 
#
#
#
#
##########

rm(list=ls())

# Paramètres

f_3dB<- 400e6 # Bande passante à -3dB : tau*2*pi*f_3dB=1
tau <- 1/(2*pi*f_3dB)/1e-9# Constante de temps système du premier ordre (ns)

# Noyau de convolution : Réponse impulsionnelle photodiode
hPhD <- Vectorize(function(t) if (t<0) return(0) else 1/tau*exp(-t/tau),"t")

Kernel <- function(x,y) h*hPhD(x-y) # Noyau de la convolution 

N <- 501 # Nombre de points
h <- 1/(N-1) # pas
tmin <- 0; tmax <- 10
T <- seq(tmin,tmax,len=N)

# Opérateur discrétisé : Matrice A[i,j] 
A <- outer(T, T, Kernel)# dim(A) : matrice N x N

# Signal optique (qu'il faudra retrouver )
# Impulsion gaussienne
fopt <- Vectorize(function(t,p) p[1]*dnorm(t,p[2],p[3]) ,"t")#+ fopt2(t)

p0<-c(1, 2, 0.5)# c(amplitude, position centrale, écart-type)
curve(hPhD,tmin,tmax,col="red");grid() # réponse impulsionnelle
curve(fopt(x,p=p0),add=TRUE,tmin,tmax) # Impulsion optique

f<- fopt(T,p=p0) # Signal optique 

# Signal du photodétecteur
g0<- convolve(hPhD(T),rev(f),type="open")/N # ? convolve 
g <- g0[1:N]# Signal mesuré length(g0) =1001

lines(T,g*5, col="blue")

svdA <- svd(A)

# Matrices of the SVD decomposition
D <- diag(svdA$d)	# svdA$d : singular values # object.size(D)
U <- svdA$u		# svdY$u : U matrix	# object.size(U)
V <- svdA$v		# svdY$v : V matrix

# plot of singular values of the SVD decomposition,
# in order to appreciate the significant singular values
#
#windows()
#plot(svdA$d[svdA$d>svdA$d[1]/1e6],main=paste("j=",j),xlim=c(0,40),log="y")

windows()
plot(svdA$d,main=paste("j="),log="y",ylim=c(min(svdA$d)/1000,max(svdA$d)))
noise <- 0.2e-3
error <- rnorm(length(g),0,noise)
points(abs(t(U)%*%error), col="red")
abline(h=sd(abs(t(U)%*%error)))
alpha<- 25*sd(abs(t(U)%*%error))
abline(h=alpha); locator(1)$y
dev.off()

Ylim <- c(-0.1,1.1)
# Filtre
Psi <-  diag(svdA$d/(svdA$d^2+alpha^2))	

TitrePlot <- bquote(paste("Paramètre de régularisation : ", alpha==.(alpha)))
plot(T,f,type="l",main=TitrePlot,ylim=Ylim) #; lines(T,g, col="blue")

error <- rnorm(length(g),0,noise) # bruit gaussien
g1 <- g + error
fa <- V%*%Psi%*%t(U)%*%g1# Signal reconstitué

lines(T,g1*10, col="blue")
lines(T,fa, col="red")





