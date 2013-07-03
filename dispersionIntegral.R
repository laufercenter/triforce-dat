#     rotation schema
#
#      y                     y
#      |                     |
#      |                     |
#      |______x        z_____|
#     /                     /
#    /                     /
#   z                     x
#     /\                    /\
# |___|                 |___|
#   +                     +




library(multicore)


args <- commandArgs(TRUE)
species <- as.numeric(args[1])




DERIVATIVE_LEVEL=0
STORE_PHI=0

#the resolution of our grid
grid=8
res_lambda = grid
res_psi = grid
res_PHI = grid
res_d = grid
res_eps = grid
res_sig = grid

#limits for the parameters
max_lambda = pi/2
max_psi = pi
max_PHI = pi
max_d = 8
max_eps = 16.736
max_sig = 7.0


res_theta=100000
max_theta=pi

dimensions=6
fd=0.0025


ex = matrix(c(1,0,0),3,1)
ey = matrix(c(0,1,0),3,1)
ez = matrix(c(0,0,1),3,1)

THRESHOLD_NUMERICAL = 0.0001

#normal vector to the plane connecting the origin, the integration origin and the interace center of the circular region
nOrigin = matrix(c(0,0,1),3,1)




# realisticEps = c(0.000000, 0.339967, 0.339967, 0.305240, 0.311815, 0.106908, 0.264953, 0.247135, 0.229317, 0.211499, 0.259964, 0.251055, 0.242146, 0.195998, 0.418722, 0.440104, 0.332840, 0.890899, 0.141225, 0.325000, 0.295992, 0.315061, 0.306647, 0.300001, 0.374177, 0.356359, 0.339967, 0.473602, 0.526699, 0.604920, 0.202590, 0.195998)
# realisticSig = 10*c(0.00000000, 0.35982400, 0.45773000, 1.92376000, 0.25522400, 0.06568880, 0.06568880, 0.06568880, 0.06568880, 0.06568880, 0.06276000, 0.06276000, 0.06276000, 0.06568880, 1.67360000, 0.41840000, 0.01158970, 0.41840000, 3.74342000, 0.71128000, 0.87864000, 0.63638600, 0.88031400, 0.71128000, 0.83680000, 1.04600000, 0.35982400, 0.00137235, 0.00071128, 0.00033723, 0.07656720, 0.05230000)
# realisticClass=c(0,1,1,1,0,2,2,2,2,2,2,2,2,2,0,0,0,0,0,3,4,4,4,4,5,6,0,0,0,0,0,0)


realisticEps = c(0.339967, 0.264953,  0.325000, 0.295992, 0.374177, 0.356359)
realisticSig = 10*c(0.35982400, 0.06568880, 0.71128000, 0.87864000, 0.83680000, 1.04600000)
realisticClass=c(0,1,2,3,4,5)
classProb=c(30.6, 50, 8.4, 9.0, 0.1, 1.8)





LJ <- function(x,eps,sig){
	(4*eps*sig^12)/(x^12) - (4*eps*sig^6)/(x^6)
	
}


LJField<-function(ts, v, eps0, sig0, d, eps1, sig1){
	#assume that v is normalised
	res=c()
	atm = c(-d,0,0)
	
	for(i in 1:length(ts)){
		t=ts[i]
		x=v*t
		f=0
		s0 = veclen(x)
		s1 = veclen(atm-x)
		f = LJ(s0,eps0,sig0) + LJ(s1,eps1,sig1)
		res=c(res,f)
	}
	
	res
}


probeLJFieldSigBuffered<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	wd=30
	res=c()
	for(i in 1:length(thetas)){
		theta=thetas[i]
		
		i_theta = floor((theta)/thetaStep)
		p=i_theta*thetaStep
		l=(theta-p)/thetaStep
		
		sig=(1-l)*sigBuffer[i_theta+1]+l*sigBuffer[i_theta+2]
		
		res=c(res,sig)
	}
	res

}


probeLJFieldEpsBuffered<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	wd=30
	res=c()
	for(i in 1:length(thetas)){
		theta=thetas[i]
		
		i_theta = floor((theta)/thetaStep)
		p=i_theta*thetaStep
		l=(theta-p)/thetaStep
		
		eps=(1-l)*epsBuffer[i_theta+1]+l*epsBuffer[i_theta+2]
		
		res=c(res,eps)
	}
	res

}


probeLJFieldSig<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	wd=20
	res=c()
	for(i in 1:length(thetas)){
		theta=thetas[i]
		v=c(cos(theta), cos(phi)*sin(theta), sin(theta)*sin(phi))
		m=optimize(LJField,interval=c(0,wd), v=v, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$minimum
		res=c(res,m)
	}
	res

}

probeLJFieldEps<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	res=c()
	ms=probeLJFieldSig(thetas,phi,eps0,sig0,d,eps1,sig1)
	for(i in 1:length(ms)){
		theta=thetas[i]
		v=c(cos(theta), cos(phi)*sin(theta), sin(theta)*sin(phi))
		m=ms[i]
		e=LJField(m,v,eps0,sig0,d,eps1,sig1)
		res=c(res,e)
	}
	res
}


probeLJField<-function(theta,phi,eps0,sig0,d,eps1,sig1){
	m=probeLJFieldSig(theta,phi,eps0,sig0,d,eps1,sig1)
	v=c(cos(theta), cos(phi)*sin(theta), sin(theta)*sin(phi))
	e=LJField(m,v,eps0,sig0,d,eps1,sig1)
	c(e,m)
}

thetaConcave<-function(phi, psi, lambda){
	s=cos(phi)^2*sin(psi)^2*(-cos(lambda)^2 + cos(psi)^2 + cos(phi)^2*sin(psi)^2)
	if(s<0){
		#print(c("perturbation error: ",s))
		s=0
	}
	a=(cos(lambda)*cos(psi) + sqrt(s)) / (cos(psi)^2 + cos(phi)^2*sin(psi)^2)
	if(a > 1){
		#print(c("acos perturbation error: ",a))
		a=1
	}
	if(a < -1){
		#print(c("acos perturbation error: ",a))
		a=-1
	}
	acos(a)
}


thetaConvex<-function(phi, psi, lambda){
	s=cos(phi)^2*sin(psi)^2*(-cos(lambda)^2 + cos(psi)^2 + cos(phi)^2*sin(psi)^2)
	if(s<0){
		#print(c("perturbation error: ",s))
		s=0
	}
	a=(cos(lambda)*cos(psi) - sqrt(s)) / (cos(psi)^2 + cos(phi)^2*sin(psi)^2)
	if(a > 1){
		#print(c("acos perturbation error: ",a))
		a=1
	}
	if(a < -1){
		#print(c("acos perturbation error: ",a))
		a=-1
	}
	acos(a)
}

targetEps<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	sin(thetas)*probeLJFieldEpsBuffered(thetas,phi,eps0,sig0,d,eps1,sig1)
}

targetSig<-function(thetas,phi,eps0,sig0,d,eps1,sig1){
	sin(thetas)*probeLJFieldSigBuffered(thetas,phi,eps0,sig0,d,eps1,sig1)
}

integrateSigThetaConcave<-function(phis,psi,lambda,eps0,sig0,d,eps1,sig1){
	res=c()
	for(i in 1:length(phis)){
		phi=phis[i]
		theta=thetaConcave(phi,psi,lambda)
		integral=integrate(targetSig, lower=0, upper=theta, phi=phi,  eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val
		res=c(res, integral)
	}
	res
}

integrateSigThetaConvex<-function(phis,psi,lambda,eps0,sig0,d,eps1,sig1){
	res=c()
	for(i in 1:length(phis)){
		phi=phis[i]
		theta=thetaConvex(phi,psi,lambda)
		integral=integrate(targetSig, lower=0, upper=theta, phi=phi,  eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val
		res=c(res, integral)
	}
	res
}

integrateEpsThetaConcave<-function(phis,psi,lambda,eps0,sig0,d,eps1,sig1){
	res=c()
	for(i in 1:length(phis)){
		phi=phis[i]
		theta=thetaConcave(phi,psi,lambda)
		integral=integrate(targetEps, lower=0, upper=theta, phi=phi,  eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val
		res=c(res, integral)
	}
	res
}

integrateEpsThetaConvex<-function(phis,psi,lambda,eps0,sig0,d,eps1,sig1){
	res=c()
	for(i in 1:length(phis)){
		phi=phis[i]
		theta=thetaConvex(phi,psi,lambda)
		integral=integrate(targetEps, lower=0, upper=theta, phi=phi,  eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val
		res=c(res, integral)
	}
	res
}

saveTable <-function(filename, dimensions, headerPHI, headerPsi, headerLambda, headerD, headerEps, headerSig, tbl){
	DERIVATIVE_LEVEL=0
	STORE_PHI=0
	
	headerPHIdat = array(0,dim=c(4,length(headerPHI)))
	headerPsidat = array(0,dim=c(4,length(headerPsi)))
	headerLambdadat = array(0,dim=c(4,length(headerLambda)))
	headerDdat = array(0,dim=c(4,length(headerD)))
	headerEpsdat = array(0,dim=c(4,length(headerEps)))
	headerSigdat = array(0,dim=c(4,length(headerSig)))
	
	tbldat = array(0,dim=c(4,length(tbl)))
	for(i in 1:length(headerPHI)){
		headerPHIdat[,i] = float2CharArray(headerPHI[i])
	}
	for(i in 1:length(headerPsi)){
		headerPsidat[,i] = float2CharArray(headerPsi[i])
	}
	for(i in 1:length(headerLambda)){
		headerLambdadat[,i] = float2CharArray(headerLambda[i])
	}
	for(i in 1:length(headerD)){
		headerDdat[,i] = float2CharArray(headerD[i])
	}
	for(i in 1:length(headerEps)){
		headerEpsdat[,i] = float2CharArray(headerEps[i])
	}
	for(i in 1:length(headerSig)){
		headerSigdat[,i] = float2CharArray(headerSig[i])
	}
	for(i in 1:length(tbl)){
		tbldat[,i] = float2CharArray(tbl[i])
	}

	write(file=filename,DERIVATIVE_LEVEL, append=FALSE)
	write(file=filename,STORE_PHI, append=TRUE)
	write(file=filename,dimensions,append=TRUE)
	write(file=filename,c(1,res_PHI,1,res_psi,1,res_lambda, 1,res_d, 1,res_eps, 1,res_sig),ncol=12,append=TRUE)
	write(file=filename,headerPHIdat,ncol=4,append=TRUE)
	write(file=filename,headerPsidat,ncol=4,append=TRUE)
	write(file=filename,headerLambdadat,ncol=4,append=TRUE)
	write(file=filename,headerDdat,ncol=4,append=TRUE)
	write(file=filename,headerEpsdat,ncol=4,append=TRUE)
	write(file=filename,headerSigdat,ncol=4,append=TRUE)
	write(file=filename,tbldat,ncol=4,append=TRUE)
}




#returns a rotation matrix around the Z axis for a certain angle theta
rotz <- function(theta){
	matrix(c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1),3,3)
}

#returns a rotation matrix around the X axis for a certain angle theta
rotx <- function(theta){
	matrix(c(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta)),3,3)
}

#calculates the three dimensional cross product between two vectors a and b
crossproduct <-function(a,b){
	c = matrix(0,3,1)
	c[1] = a[2]*b[3] - a[3]*b[2]
	c[2] = a[3]*b[1] - a[1]*b[3]
	c[3] = a[1]*b[2] - a[2]*b[1]
	c
}

#matrix-vector multiplication and reformatting of the output into a proper vector
mvMultiply <- function(A,b){
	matrix(c(A %*% b),3,1)
}

dot <- function(a,b){
	c(t(a) %*% b)
}

#normalization of a vector
normalise <- function(a){
	a / sqrt(c(t(a) %*% a))
}

#gives angle between two normal vectors
normalAngle <-function(a,b){
	acos(c(t(a) %*% b))
}

veclen <- function(a){
	sqrt(c(t(a) %*% a))
}

angle <-function(a,b){
	acos(c(t(a) %*% b) / (veclen(a)*veclen(b)))
}

nangle <-function(a,b){
	acos(c(t(a) %*% b))
}

#converts a PHI angle into a phi angle
PHI2phi <- function(PHI, psi, lambda){

	

	#print(c("PHI2phi",PHI))
	if(lambda==0)	res=0
	else if(psi==lambda && PHI==0) res=pi/2
	else if(PHI==0) res=0
	else if(isWithinNumericalLimits(psi+lambda,pi) && isWithinNumericalLimits(PHI,pi)) res=pi
	else{
		border = mvMultiply(rotz(-lambda), ex)
		v = mvMultiply(rotz(-psi), ex)
		a = border-ex
		aRot = mvMultiply(rotx(PHI), a)
	
		ip = ex + aRot
	
		nIntersection = crossproduct(v,ip)
		nIntersection = normalise(nIntersection)
		phiIntersection = normalAngle(nIntersection,nOrigin)
		res = phiIntersection
	}



	#res = acos((-cos(PHI) * cos(psi) * sin(lambda)+cos(lambda) * sin(psi))/(sqrt(cos(psi)^2 * sin(lambda)^2-2 * cos(PHI) * cos(lambda) * cos(psi) * sin(lambda) * sin(psi)+(cos(lambda)^2+sin(PHI)^2 * sin(lambda)^2) * sin(psi)^2)))


	res


}



phiLimitPlain <-function(psi,lambda){
	if(lambda==0) limit=0
	else if(psi==0) limit=pi/2
	else if(psi<=lambda) limit=pi/2
	else if(psi+lambda>=pi) limit=pi/2
	else {
		if(psi>pi/2) psi=pi-psi

		x = (cos(lambda)^2-cos(psi)^2) * csc(psi)^2
		if(x<0){
			print(c("INVALID LIMIT, REPLACING ",x," WITH 0 ",lambda,psi,cos(lambda)^2,cos(psi)^2,csc(psi)^2))
			x=0
		}
		limit=acos(sqrt(x))
	}
	limit
}


PHILimitPlain <-function(psi,lambda){
	reverse=FALSE

	if(psi+lambda>pi){
		limit = acos(tan(psi)*cot(lambda))
	}
	else
	if(psi+lambda==pi){
		limit = pi
	}
	else{

		if(psi>pi/2){
			reverse=TRUE
			psi = pi-psi
		}

		if(lambda==0) limit=pi/2
		else
		if(psi==0) limit=0
		else
		if(psi<lambda) limit=pi/2
		else
		if(psi==lambda) limit=0
		else{
			limit = acos(cot(psi)*tan(lambda))
		}

		if(reverse){
			limit = pi-limit
		}
	}


	limit
}




isWithinNumericalLimits <- function(x,l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) r = TRUE
	else r = FALSE
	
	r
}






csc <-function(x){
	1/sin(x)
}

cot <- function(x){
	1/tan(x)
}





isAtPositiveDiscontinuity <-function(i_psi,i_lambda){
	b=FALSE
	psi = headerPsi[i_psi+1]
	lambda = headerLambda[i_lambda+1]


	if(i_psi+1 < res_psi){
		vpsi = headerPsi[i_psi+2]
		vlambda = headerLambda[i_lambda+1]

		if(psi<lambda){
			if(vpsi >= vlambda) b = TRUE
		}
		else if(psi+lambda<pi){
			if(vpsi+vlambda>=pi) b=TRUE
		}
	}

	if(i_lambda+1 < res_lambda){
		vpsi = headerPsi[i_psi+1]
		vlambda = headerLambda[i_lambda+2]

		if(psi<lambda){
			if(vpsi >= vlambda) b = TRUE
		}
		else if(psi+lambda<pi){
			if(vpsi+vlambda>=pi) b=TRUE
		}
	}

	b

}




isAtNegativeDiscontinuity <-function(i_psi,i_lambda){
	b=FALSE
	psi = headerPsi[i_psi+1]
	lambda = headerLambda[i_lambda+1]


	if(i_psi > 0){
		vpsi = headerPsi[i_psi]
		vlambda = headerLambda[i_lambda+1]

		if(psi+lambda>=pi){
			if(vpsi+vlambda<pi) b=TRUE
		}
		else if(psi>=lambda){
			if(vpsi < vlambda) b = TRUE
		}
	}

	if(i_lambda > 0){
		vpsi = headerPsi[i_psi+1]
		vlambda = headerLambda[i_lambda]

		if(psi+lambda>=pi){
			if(vpsi+vlambda<pi) b=TRUE
		}
		else if(psi>=lambda){
			if(vpsi < vlambda) b = TRUE
		}
	}

	b

}

isWithinLimitsConcave <-function(PHI, psi, lambda){
	valid=TRUE
	
	if(PHI<0 || psi<0 || lambda<0 || PHI>max_PHI || psi > max_psi || lambda > max_lambda) valid = FALSE	
	
	
	valid
}








integralSigConcave <-function(PHI, psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#print(c("icx",phi,PHI,psi,lambda))
	#calculate the maximal phi value to which we can integrate
	#phiLimit = PHI2phi(pi/2,psi,lambda)

	#print(c("enterconcave",PHI,psi,lambda))

	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<lambda){


			phiLimit = phiLimitPlain(psi,lambda)
		
			if(PHI==0) phi=pi


			if(phi>=pi/2){
				A = A + abs(integrate(integrateSigThetaConcave,lower=phi,upper=pi,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
			}
			else{
				A = A + abs(integrate(integrateSigThetaConcave,lower=pi/2,upper=pi,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
				A = A + abs(integrate(integrateSigThetaConvex,lower=phi,upper=pi/2,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
			}

	
	}
	else{

		phiLimit = phiLimitPlain(psi,lambda)
		PHILimit = PHILimitPlain(psi,lambda)
		#print(c("row1",phi,phiLimit,PHILimit))

		if(psi+lambda >= pi){
			if(PHI >= PHILimit || isWithinNumericalLimits(PHI,PHILimit)){
				A = A - abs(integrate(integrateSigThetaConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(integrateSigThetaConcave,lower=0,upper=ip,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
		}
		else{


			if(PHI > PHILimit){
				A = A + abs(integrate(integrateSigThetaConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(integrateSigThetaConcave,lower=0,upper=ip,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
		}

	}
	
	A

}






integralEpsConcave <-function(PHI, psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#print(c("icx",phi,PHI,psi,lambda))
	#calculate the maximal phi value to which we can integrate
	#phiLimit = PHI2phi(pi/2,psi,lambda)

	#print(c("enterconcave",PHI,psi,lambda))

	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<lambda){


			phiLimit = phiLimitPlain(psi,lambda)
		
			if(PHI==0) phi=pi


			if(phi>=pi/2){
				A = A + abs(integrate(integrateEpsThetaConcave,lower=phi,upper=pi,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
			}
			else{
				A = A + abs(integrate(integrateEpsThetaConcave,lower=pi/2,upper=pi,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
				A = A + abs(integrate(integrateEpsThetaConvex,lower=phi,upper=pi/2,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1)$val)
			}

	
	}
	else{

		phiLimit = phiLimitPlain(psi,lambda)
		PHILimit = PHILimitPlain(psi,lambda)
		#print(c("row1",phi,phiLimit,PHILimit))

		if(psi+lambda >= pi){
			if(PHI >= PHILimit || isWithinNumericalLimits(PHI,PHILimit)){
				A = A - abs(integrate(integrateEpsThetaConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(integrateEpsThetaConcave,lower=0,upper=ip,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
		}
		else{


			if(PHI > PHILimit){
				A = A + abs(integrate(integrateEpsThetaConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(integrateEpsThetaConcave,lower=0,upper=ip,psi=psi,lambda=lambda, eps0=eps0, sig0=sig0, d=d, eps1=eps1, sig1=sig1, stop.on.error=FALSE)$val)
		}

	}
	
	A

}








gradient <- function(integrator, PHI, psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i){

	switch(modePHI, 
		central={pfd_PHI=fd; nfd_PHI= -fd; cfd_PHI=0},
		forward={pfd_PHI=2*fd; nfd_PHI=0; cfd_PHI=fd},
		backward={pfd_PHI=0; nfd_PHI= -2*fd; cfd_PHI= -fd},
		invalid={pfd_PHI=0; nfd_PHI= 0; cfd_PHI= 0},
	)
	switch(modepsi, 
		central={pfd_psi=fd; nfd_psi= -fd; cfd_psi=0},
		forward={pfd_psi=2*fd; nfd_psi=0; cfd_psi=fd},
		backward={pfd_psi=0; nfd_psi= -2*fd; cfd_psi= -fd},
		invalid={pfd_psi=0; nfd_psi= 0; cfd_psi= 0},
	)
	switch(modelambda, 
		central={pfd_lambda=fd; nfd_lambda= -fd; cfd_lambda=0},
		forward={pfd_lambda=2*fd; nfd_lambda=0; cfd_lambda=fd},
		backward={pfd_lambda=0; nfd_lambda= -2*fd; cfd_lambda= -fd},
		invalid={pfd_lambda=0; nfd_lambda= 0; cfd_lambda= 0},
	)
	


	res = array(0,dim=c(3))
	res[1] = (integrator(PHI+pfd_PHI, psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI+nfd_PHI, psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	res[2] = (integrator(PHI, psi+pfd_psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI, psi+nfd_psi, lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	res[3] = (integrator(PHI, psi, lambda+pfd_lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI, psi, lambda+nfd_lambda, eps0, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	
	res[4] = (integrator(PHI, psi, lambda, eps0+fd, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI, psi, lambda, eps0-fd, sig0, d, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	res[5] = (integrator(PHI, psi, lambda, eps0, sig0+fd, d, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI, psi, lambda, eps0, sig0-fd, d, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	res[6] = (integrator(PHI, psi, lambda, eps0, sig0, d+fd, eps1, sig1, modePHI, modepsi, modelambda, i) - integrator(PHI, psi, lambda, eps0, sig0, d-fd, eps1, sig1, modePHI, modepsi, modelambda, i))/(2*fd)
	res
	
}




hessian <- function(integrator, PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	res = array(0,dim=c(3,3))

	switch(modePHI, 
		central={pfd_PHI=fd; nfd_PHI= -fd; cfd_PHI=0},
		forward={pfd_PHI=2*fd; nfd_PHI=0; cfd_PHI=fd},
		backward={pfd_PHI=0; nfd_PHI= -2*fd; cfd_PHI= -fd},
		invalid={pfd_PHI=0; nfd_PHI= 0; cfd_PHI= 0},
	)
	switch(modepsi, 
		central={pfd_psi=fd; nfd_psi= -fd; cfd_psi=0},
		forward={pfd_psi=2*fd; nfd_psi=0; cfd_psi=fd},
		backward={pfd_psi=0; nfd_psi= -2*fd; cfd_psi= -fd},
		invalid={pfd_psi=0; nfd_psi= 0; cfd_psi= 0},
	)
	switch(modelambda, 
		central={pfd_lambda=fd; nfd_lambda= -fd; cfd_lambda=0},
		forward={pfd_lambda=2*fd; nfd_lambda=0; cfd_lambda=fd},
		backward={pfd_lambda=0; nfd_lambda= -2*fd; cfd_lambda= -fd},
		invalid={pfd_lambda=0; nfd_lambda= 0; cfd_lambda= 0},
	)


	#f PHI PHI
	res[1,1] = (integrator(PHI+pfd_PHI, psi, lambda, modePHI, modepsi, modelambda, i) - 2*integrator(PHI+cfd_PHI, psi, lambda, modePHI, modepsi, modelambda, i) + integrator(PHI+nfd_PHI, psi, lambda, modePHI, modepsi, modelambda, i))/fd^2
	#f PHI psi
	res[1,2] = (integrator(PHI+pfd_PHI, psi+pfd_psi, lambda, modePHI, modepsi, modelambda, i) - integrator(PHI+pfd_PHI, psi+nfd_psi, lambda, modePHI, modepsi, modelambda, i) - integrator(PHI+nfd_PHI, psi+pfd_psi, lambda, modePHI, modepsi, modelambda, i) + integrator(PHI+nfd_PHI, psi+nfd_psi, lambda, modePHI, modepsi, modelambda, i))/(4*fd^2)
	#f PHI lambda
	res[1,3] = (integrator(PHI+pfd_PHI, psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) - integrator(PHI+pfd_PHI, psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i) - integrator(PHI+nfd_PHI, psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) + integrator(PHI+nfd_PHI, psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i))/(4*fd^2)

	#f psi PHI
	res[2,1] = res[1,2]
	#f psi psi
	res[2,2] = (integrator(PHI, psi+pfd_psi, lambda, modePHI, modepsi, modelambda, i) - 2*integrator(PHI, psi+cfd_psi, lambda, modePHI, modepsi, modelambda, i) + integrator(PHI, psi+nfd_psi, lambda, modePHI, modepsi, modelambda, i))/fd^2
	#f psi lambda
	res[2,3] = (integrator(PHI, psi+pfd_psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) - integrator(PHI, psi+pfd_psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i) - integrator(PHI, psi+nfd_psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) + integrator(PHI, psi+nfd_psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i))/(4*fd^2)


	#f lambda  PHI
	res[3,1] = res[1,3]
	#f lambda psi
	res[3,2] = res[2,3]
	#f psi psi
	res[3,3] = (integrator(PHI, psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) - 2*integrator(PHI, psi, lambda+cfd_lambda, modePHI, modepsi, modelambda, i) + integrator(PHI, psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i))/fd^2

	res
}







calculateDispersionField<-function(i_PHI, i_psi, i_lambda, i_eps, i_sig, i_d, speciesEps, speciesSig){

	d = max_d * i_d/(res_d-1)
	eps = max_eps * i_eps/(res_eps-1)
	sig = max_sig * i_sig/(res_sig-1)
	lambda = max_lambda*i_lambda/(res_lambda-1)
	#lambda==0 will not occur in reality. this point would have extreme derivatives, so we start with a circle of very small radius instead
	lambda = max(0.01,lambda)
	psi = (i_psi * (max_psi)/(res_psi-1))
	PHI = max_PHI * i_PHI/(res_PHI-1)
	
	res=c()


	if(isWithinLimitsConcave(PHI,psi,lambda)){
		deps = integralEpsConcave(PHI,psi,lambda,eps,sig,d,speciesEps,speciesSig)
		dsig = integralSigConcave(PHI,psi,lambda,eps,sig,d,speciesEps,speciesSig)
	}
	else{
		deps = NaN
		dsig = NaN
	}
	
	res=c(res,deps)
	res=c(res,dsig)
	
	modePHI="central"
	modepsi="central"
	modelambda="central"
	moded="central"
	modeeps="central"
	modesig="central"

	
	if(i_PHI==0) modePHIConcave="forward"
	else if(i_PHI==(res_PHI-1)) modePHIConcave="backward"
	else modePHIConcave="central"
	if(i_psi==0) modepsiConcave="forward"
	else if(i_psi==(res_psi-1)) modepsiConcave="backward"
	else modepsiConcave="central"
	if(i_lambda==0) modelambdaConcave="forward"
	else if(i_lambda==(res_lambda-1)) modelambdaConcave="backward"
	else modelambdaConcave="central"
	
	derivableConcave=TRUE
	derivableConvex=TRUE


	if(!isWithinLimitsConcave(PHI+2*fd,psi,lambda) && !isWithinLimitsConcave(PHI-2*fd,psi,lambda)) modePHIConcave="invalid"
	else if(modePHIConcave!="backward" && !isWithinLimitsConcave(PHI+2*fd,psi,lambda)) modePHIConcave="backward"
	else if(modePHIConcave!="forward" && !isWithinLimitsConcave(PHI-2*fd,psi,lambda)) modePHIConcave="forward"


	if(modepsiConcave!="backward" && isAtPositiveDiscontinuity(i_psi,i_lambda)) modepsiConcave="backward"
	else if(modepsiConcave!="forward" && isAtNegativeDiscontinuity(i_psi,i_lambda)) modepsiConcave="forward"

	if(modelambdaConcave!="backward" && isAtPositiveDiscontinuity(i_psi,i_lambda)) modelambdaConcave="backward"
	else if(modelambdaConcave!="forward" && isAtNegativeDiscontinuity(i_psi,i_lambda)) modelambdaConcave="forward"


	if(!isWithinLimitsConcave(PHI,psi+2*fd,lambda) && !isWithinLimitsConcave(PHI,psi-2*fd,lambda)) modepsiConcave="invalid"
	else if(modepsiConcave!="backward" && !isWithinLimitsConcave(PHI,psi+2*fd,lambda)) modepsiConcave="backward"
	else if(modepsiConcave!="forward" && !isWithinLimitsConcave(PHI,psi-2*fd,lambda)) modepsiConcave="forward"

	if(!isWithinLimitsConcave(PHI,psi,lambda+2*fd) && !isWithinLimitsConcave(PHI,psi,lambda-2*fd)) modelambdaConcave="invalid"
	else if(modelambdaConcave!="backward" && !isWithinLimitsConcave(PHI,psi,lambda+2*fd)) modelambdaConcave="backward"
	else if(modelambdaConcave!="forward" && !isWithinLimitsConcave(PHI,psi,lambda-2*fd)) modelambdaConcave="forward"


	if(is.nan(dsig) || is.nan(deps) || !derivableConcave){
		geps = NaN
		gsig = NaN
	}
	else{
		geps = gradient(integralEpsConcave,PHI,psi,lambda,eps,sig,d,speciesEps,speciesSig, modePHIConcave, modepsiConcave, modelambdaConcave)
		gsig = gradient(integralSigConcave,PHI,psi,lambda,eps,sig,d,speciesEps,speciesSig, modePHIConcave, modepsiConcave, modelambdaConcave)
	}

	res=c(res,geps)
	res=c(res,gsig)
	
	res
	
}















#empty grid data

dataEps=array(0,dim=c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig))
dataSig=array(0,dim=c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig))

gradientsEps = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig))
gradientsSig = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig))

headerLambda=array(0,dim=c(res_lambda))
headerPsi=array(0,dim=c(res_psi))
headerPHI=array(0,dim=c(res_PHI))
headerD=array(0,dim=res_d)
headerEps=array(0,dim=res_eps)
headerSig=array(0,dim=res_sig)

headerTheta=array(0,dim=res_theta)

epsBuffer=array(0,dim=c(res_theta))
sigBuffer=array(0,dim=c(res_theta))




#first create the header
for(i_lambda in 0:(res_lambda-1)){
	lambda = max_lambda*i_lambda/(res_lambda-1)
	headerLambda[i_lambda+1]=lambda
}

#iterate over psi angles
for(i_psi in 0:(res_psi-1)){
	psi = (i_psi * (max_psi)/(res_psi-1))
	headerPsi[i_psi+1]=psi
}

for(i_PHI in 0:(res_PHI-1)){
	PHI = max_PHI * i_PHI/(res_PHI-1)
	headerPHI[i_PHI+1]=PHI
}

for(i_d in 0:(res_d-1)){
	d = max_d*i_d/(res_d-1)
	headerD[i_d+1]=d
}

for(i_sig in 0:(res_sig-1)){
	sig = max_sig*i_sig/(res_sig-1)
	headerSig[i_sig+1]=sig
}
for(i_eps in 0:(res_eps-1)){
	eps = max_eps*i_eps/(res_eps-1)
	headerEps[i_eps+1]=eps

}

for(i_theta in 0:(res_theta-1)){
	theta = max_theta*i_theta/(res_theta-1)
	headerTheta[i_theta+1]=theta

}
thetaStep=headerTheta[2]



speciesEps = realisticEps[species+1]
speciesSig = realisticSig[species+1]


starttime = proc.time()[3]

#iterate over lambda angles
for(i_d in 0:(res_d-1)){
	d = max_d * i_d/(res_d-1)


	for(i_eps in 0:(res_eps-1)){
		eps = max_eps * i_eps/(res_eps-1)
		
		for(i_sig in 0:(res_sig-1)){
			sig = max_sig * i_sig/(res_sig-1)
			
			if(i_sig==0 && i_eps==0 && i_d==0) print(paste(round(100*i_d/(res_d)),"% completed",sep=""))
			else{
				maxtasks=res_d*res_eps*res_sig
				cmptasks=i_d*(res_eps*res_sig) + i_eps*(res_sig) + i_sig
				print(paste(round(100*cmptasks/maxtasks,digits=2),"% completed. eta: ",round((maxtasks*t/cmptasks - t)/(60*60), digits=2),"h",sep=""))
			}
			
			
			#create par buffer
			res = mclapply(headerTheta, probeLJField, phi=0, eps0=eps, sig0=sig, d=d, eps1=speciesEps, sig1=speciesSig)
			for(i in 1:res_theta){
				epsBuffer[i]=res[[i]][1]
				sigBuffer[i]=res[[i]][2]
			}
			
			
			
			
			for(i_lambda in 0:(res_lambda-1)){

				#iterate over psi angles
				for(i_psi in 0:(res_psi-1)){
					
					res = mclapply(0:(res_PHI-1), calculateDispersionField, i_psi=i_psi, i_lambda=i_lambda, i_eps=i_eps, i_sig=i_sig, i_d=i_d, speciesEps=speciesEps, speciesSig=speciesEps)
					
					for(i_PHI in 0:(res_PHI-1)){
					
						dataEps[i_PHI+1, i_psi+1, i_lambda+1, i_d+1, i_eps+1, i_sig+1] = res[[i_PHI+1]][1]
						dataSig[i_PHI+1, i_psi+1, i_lambda+1, i_d+1, i_eps+1, i_sig+1] = res[[i_PHI+1]][2]
						gradientsEps[,i_PHI+1, i_psi+1, i_lambda+1, i_d+1, i_eps+1, i_sig+1] = res[[i_PHI+1]][3:8]
						gradientsSig[,i_PHI+1, i_psi+1, i_lambda+1, i_d+1, i_eps+1, i_sig+1] = res[[i_PHI+1]][9:14]
					}
					
				}

			}
			
			stoptime=proc.time()[3]
			t=(stoptime-starttime)			
			
		}
	}
	
	
	
}
print(paste(100,"% completed",sep=""))



print("saving data")
#now, save to file

source("floatconversion.R")


saveTable(paste("dispersionFieldEps",species,".csv",sep=""),c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig),headerPHI, headerPsi, headerLambda, headerD, headerEps, headerSig, dataEps)
saveTable(paste("dispersionFieldSig",species,".csv",sep=""),c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig),headerPHI, headerPsi, headerLambda, headerD, headerEps, headerSig, dataSig)

print("saving derivatives and hessians")
for(i in 1:dimensions){
	saveTable(paste("dispersionFieldEps",species,"_",(i-1),".csv", sep=""),c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig),headerPHI, headerPsi, headerLambda, headerD, headerEps, headerSig, gradientsEps[i,,,,,,])
	saveTable(paste("dispersionFieldSig",species,"_",(i-1),".csv", sep=""),c(res_PHI,res_psi,res_lambda,res_d,res_eps,res_sig),headerPHI, headerPsi, headerLambda, headerD, headerEps, headerSig, gradientsSig[i,,,,,,])
}


