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


FAILSAFE=TRUE
DODERIVATIVES=FALSE


#the resolution of our grid
res_lambda = 12
res_psi = 12
res_PHI = 12

#limits for the parameters
max_lambda = pi/2
max_psi = pi
max_PHI = pi

dimensions=3
fd=0.025


ex = matrix(c(1,0,0),3,1)
ey = matrix(c(0,1,0),3,1)
ez = matrix(c(0,0,1),3,1)

THRESHOLD_NUMERICAL = 0.0001

#normal vector to the plane connecting the origin, the integration origin and the interace center of the circular region
nOrigin = matrix(c(0,0,1),3,1)


saveTable <-function(filename, dimensions, headerPHI, headerPsi, headerLambda, tbl, tblgradients, tblhessians){

	headerPHIdat = array(0,dim=c(8,length(headerPHI)))
	headerPsidat = array(0,dim=c(8,length(headerPsi)))
	headerLambdadat = array(0,dim=c(8,length(headerLambda)))
	tbldat = array(0,dim=c(8,length(tbl)))
	tblgradientsdat = array(0,dim=c(8,length(tblgradients)))
	tblhessiansdat = array(0,dim=c(8,length(tblhessians)))
	for(i in 1:length(headerPHI)){
		headerPHIdat[,i] = double2CharArray(headerPHI[i])
	}
	for(i in 1:length(headerPsi)){
		headerPsidat[,i] = double2CharArray(headerPsi[i])
	}
	for(i in 1:length(headerLambda)){
		headerLambdadat[,i] = double2CharArray(headerLambda[i])
	}
	for(i in 1:length(tbl)){
		tbldat[,i] = double2CharArray(tbl[i])
	}
	for(i in 1:length(tblgradients)){
		tblgradientsdat[,i] = double2CharArray(tblgradients[i])
	}
	for(i in 1:length(tblhessians)){
		tblhessiansdat[,i] = double2CharArray(tblhessians[i])
	}

	write(file=filename,3,append=FALSE)
	write(file=filename,c(1,res_PHI,1,res_psi,1,res_lambda),ncol=3,append=TRUE)
	write(file=filename,headerPHIdat,ncol=8,append=TRUE)
	write(file=filename,headerPsidat,ncol=8,append=TRUE)
	write(file=filename,headerLambdadat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
	write(file=filename,tblgradientsdat,ncol=8,append=TRUE)
	write(file=filename,tblhessiansdat,ncol=8,append=TRUE)
}


saveTableRaw <-function(filename, dimensions, header, tbl){

	headerdat = array(0,dim=c(length(header),8))
	tbldat = array(0,dim=c(length(tbl),8))
	write(file=filename,c(length(dimensions),nrow(header)),append=FALSE)
	write(file=filename,dimensions,append=TRUE)
	for(i in 1:length(header)){
		headerdat[i,] = header[i]
	}
	for(i in 1:length(tbl)){
		tbldat[i,] = tbl[i]
	}

	write(file=filename,headerdat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
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


PHI2HalfSphereAngle <- function(PHI, psi, lambda){

	border = mvMultiply(rotz(-lambda), ex)
	v = mvMultiply(rotz(-psi), ex)
	a = border-ex
	aRot = mvMultiply(rotx(PHI), a)
	
	ip = ex + aRot
	ipv = mvMultiply(rotz(psi),ip)

	res = normalAngle(normalise(ipv),ex)
	


	res


}



tmp=c(0,0)


phiLimitPlain <-function(psi,lambda){
	if(lambda==0) limit=0
	else if(psi<=lambda) limit=pi/2
	else if(psi+lambda>=pi) limit=pi/2
	else {
		if(psi>pi/2) psi=pi-psi

		x = (cos(lambda)^2-cos(psi)^2) * csc(psi)^2
		if(x<0){
			print(c("INVALID LIMIT, REPLACING ",x," WITH 0"))
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

arrayPHILimitPlain <-function(detail){
	a = array(0,dim=c(detail,detail))

	for(i in 1:detail){
		psi = (i-1)*pi/(detail-1)

		for(j in 1:detail){
		lambda = (j-1)*(pi/2)/(detail-1)

		a[i,j] = PHILimitPlain(psi,lambda)
		}
	}

	a

}

arrayphiLimitPlain <-function(detail){
	a = array(0,dim=c(detail,detail))

	for(i in 1:detail){
		psi = (i-1)*pi/(detail-1)

		for(j in 1:detail){
		lambda = (j-1)*(pi/2)/(detail-1)

		a[i,j] = phiLimitPlain(psi,lambda)
		}
	}

	a

}


maxPHI <- function(psi, lambda){
	#tmp <<- c(psi,lambda)
	hsboundary=pi

	if(lambda==0) res=0
	else
	if(psi<lambda) res = pi/2
	else
	if(psi-lambda>=pi/2) res=0
	else
	if(psi+lambda > pi/2){
		if(psi-lambda >=pi/2) res=0
		else res = halfSpherePHI(psi,lambda)
	}
	else res = acos(cot(psi)*tan(lambda))

	res
}

isWithinNumericalLimits <- function(x,l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) r = TRUE
	else r = FALSE
	
	r
}


halfSpherePHI <- function(psi,lambda)
{
	if(abs(psi-pi)<=THRESHOLD_NUMERICAL && abs(lambda-pi/2)<=THRESHOLD_NUMERICAL) res=pi/2
	else
	if(abs(pi/2+lambda - psi) <= THRESHOLD_NUMERICAL) res = 0
	else
	if(lambda==0){
		res=0
	}
	else{
		rho = acos(cos(lambda)*csc(psi))
		if(is.nan(rho)){
		stop(c("HALF SPHERE NAN EVALUATION",psi,lambda))
	}
	
		v = mvMultiply(rotx(rho),ey)
		n = mvMultiply(rotz(psi),ex) * cos(lambda)
		t0 = normalise(v-n)
		v2 = mvMultiply(rotz(psi-lambda),ex)
		t1 = normalise(v2-n)
		res = normalAngle(t0,t1)
	}
	
	if(is.nan(res)){
		stop(c("HALF SPHERE NAN EVALUATION",psi,lambda))
	}

	res
}





csc <-function(x){
	1/sin(x)
}

cot <- function(x){
	1/tan(x)
}
#gives theta values dependent on phi, psi and lambda angles for a convex spherical arc
arcConvex <- function(phi, psi, lambda){
	if(lambda==0){
		v = array(0,length(phi))
	}
	else{
		s = cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		for(i in 1:length(phi))
		if(s[i]<0){
			print(c("perturbation error: s=",s[i]))
			s[i] = 0
		}
		v = 1+(-cos(lambda) * cos(psi)+sqrt(s))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	}

	#print(c("fphi",phi))
	#print(c("fpsi",psi))
	#print(c("flambda",lambda))
	#print(c("fv",v))

	
	#for(i in 1:length(phi))
	#	v[i] = min(v[i],pi/2)

	v
}

#gives theta values dependent on phi, psi and lambda angles for a concave spherical arc
arcConcave <- function(phi, psi, lambda){
	if(lambda==0){
		v = array(0,length(phi))
	}
	else{
		s = cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		for(i in 1:length(phi))
		if(s[i]<0){
			print(c("perturbation error: s=",s[i]))
			s[i] = 0
		}
		v = 1-(cos(lambda) * cos(psi)+sqrt(s))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	}

	#for(i in 1:length(phi))
	#	v[i] = min(v[i],pi/2)


	v
}



phiLimitAbs <-function(psi,lambda){
	#print(c(psi,lambda))
	#if(abs(cos(psi)) >= abs(cos(lambda))) limit=pi/2
	#else limit=acos(sqrt((cos(lambda)^2-cos(psi)^2) * csc(psi)^2))
	#print(c("phiLimitAbs",maxPHI(psi,lambda)))

	limit = PHI2phi(maxPHI(psi,lambda),psi,lambda)

	limit
}


isAtPositiveConcaveDiscontinuity <-function(i_psi,i_lambda){
	b=FALSE
	if(i_psi+1 < res_psi){
		if(isWithinNumericalLimits(headerPsiConcave[i_psi+2]+headerLambdaConcave[i_lambda+1],pi)) b=TRUE
	}
	if(i_lambda+1 < res_lambda){
		if(isWithinNumericalLimits(headerPsiConcave[i_psi+1]+headerLambdaConcave[i_lambda+2],pi)) b=TRUE
	}
	b
}

isAtNegativeConcaveDiscontinuity <-function(i_psi,i_lambda){
	b=FALSE
	if(isWithinNumericalLimits(headerPsiConcave[i_psi+1]+headerLambdaConcave[i_lambda+1],pi)) b=TRUE
	if(i_psi > 0){
		if(isWithinNumericalLimits(headerPsiConcave[i_psi]+headerLambdaConcave[i_lambda+1],pi)) b=TRUE
	}
	if(i_lambda > 0){
		if(isWithinNumericalLimits(headerPsiConcave[i_psi+1]+headerLambdaConcave[i_lambda],pi)) b=TRUE
	}
	b
}


isWithinLimitsConvex <-function(PHI, psi, lambda){
	valid = TRUE

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	#if(PHI2HalfSphereAngle(PHI,psi,lambda)-THRESHOLD_NUMERICAL > pi/2){
	#	valid=FALSE
	#}

	if(PHI<0 || psi<0 || lambda<0 || PHI>max_PHI || psi > max_psi || lambda > max_lambda) valid = FALSE	
	
	valid
}

isWithinLimitsConcave <-function(PHI, psi, lambda){
	valid=TRUE

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	#if(psi<=lambda){
	#	valid=FALSE
	
	#}
	#else{
		#if(PHI2HalfSphereAngle(PHI,psi,lambda)-THRESHOLD_NUMERICAL > pi/2){
		#	valid=FALSE
		#}
#
#	}
	
	if(PHI<0 || psi<0 || lambda<0 || PHI>max_PHI || psi > max_psi || lambda > max_lambda) valid = FALSE	
	
	
	valid
}






integralConvex <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#print(c("icx",phi,PHI,psi,lambda))
	#calculate the maximal phi value to which we can integrate
	#phiLimit = PHI2phi(pi/2,psi,lambda)


	#print(c("enterconvex",PHI,psi,lambda))
	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<=lambda){
		#if(FAILSAFE && PHI2HalfSphereAngle(PHI,psi,lambda)-THRESHOLD_NUMERICAL > pi/2){
		#	A = NaN
		#}
		#else{


			phiLimit = phiLimitPlain(psi,lambda)
		
			if(PHI==0) phi=pi

			#print(c("row0",phi,phiLimit))


			#calculations for convex arcs
			#if(phi<=pi/2){
			#	A = A + abs(integrate(arcConvex,lower=0,upper=phi,psi=psi,lambda=lambda)$val)
			#}
			#else{
			#	A = A + abs(integrate(arcConvex,lower=0,upper=pi/2,psi=psi,lambda=lambda)$val)
			#	A = A + abs(integrate(arcConcave,lower=pi/2,upper=pi-phi,psi=psi,lambda=lambda)$val)
			#}




			if(phi>=pi/2){
				A = A + abs(integrate(arcConcave,lower=phi,upper=pi,psi=psi,lambda=lambda)$val)
			}
			else{
				A = A + abs(integrate(arcConcave,lower=pi/2,upper=pi,psi=psi,lambda=lambda)$val)
				A = A + abs(integrate(arcConvex,lower=phi,upper=pi/2,psi=psi,lambda=lambda)$val)
			}







		#}
		
	}
	else if(psi+lambda>=pi){
		A = NaN
	}
	else{



		#there are no convex regions above pi/2 (they are handled by the mirror-side)

		#we reject every PHI that is "on the other side of the half-sphere"
		#print(c("icvx",psi,lambda,PHI,halfSpherePHI(psi,lambda)))
		#if(psi-lambda>=pi/2 || (psi+lambda>pi/2 && PHI > halfSpherePHI(psi,lambda))){
		#	A= -1
		#}
		#else{


		#if(FAILSAFE && PHI2HalfSphereAngle(PHI,psi,lambda)-THRESHOLD_NUMERICAL > pi/2){
		#	A = NaN
		#}
		#else{
		
			phiLimit = phiLimitPlain(psi,lambda)
			PHILimit = PHILimitPlain(psi,lambda)

			#print(c("row1",phi,phiLimit,PHILimit))

			#for the convex case
			#if(PHI < PHILimit){
			#	A = A - abs(integrate(arcConcave,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
			#	#print(abs(integrate(arcConcave,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, stop.on.error=FALSE)$val))
			#	ip=phiLimit
			#}
			#else ip=phi

			
			#A = A + abs(integrate(arcConvex,lower=0,upper=ip,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
			#print(abs(integrate(arcConvex,lower=0,upper=ip,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)


			if(PHI > PHILimit){
				A = A + abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)

		#}
	}
	A
}



integralConcave <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
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
				A = A + abs(integrate(arcConcave,lower=phi,upper=pi,psi=psi,lambda=lambda)$val)
			}
			else{
				A = A + abs(integrate(arcConcave,lower=pi/2,upper=pi,psi=psi,lambda=lambda)$val)
				A = A + abs(integrate(arcConvex,lower=phi,upper=pi/2,psi=psi,lambda=lambda)$val)
			}

	
	}
	else{

		phiLimit = phiLimitPlain(psi,lambda)
		PHILimit = PHILimitPlain(psi,lambda)
		#print(c("row1",phi,phiLimit,PHILimit))

		if(psi+lambda >= pi){
			if(PHI >= PHILimit || isWithinNumericalLimits(PHI,PHILimit)){
				A = A - abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
		}
		else{


			if(PHI > PHILimit){
				A = A + abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A - abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda, stop.on.error=FALSE)$val)
		}

	}
	
	A

}


#integralConcave <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
#	#convert PHI to phi
#	phi = PHI2phi(PHI,psi,lambda)
#	#calculate the maximal phi value to which we can integrate
#	
#	
#
#	A=0
#
#	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
#	if(psi<=lambda){
#		phiLimit = phiLimitAbs(psi,lambda)
#		#this is empty on purpose, calculating these values makes only theoretically sense, we will not use them
#	}
#	else{
#		if(psi-lambda>=pi/2 || (psi+lambda>pi/2 && PHI > halfSpherePHI(psi,lambda))){
#			A= NaN
#		}
#		else{
#
#
#			phiLimit = phiLimitAbs(psi,lambda)
#			PHILimit = maxPHI(psi,lambda)
#			#for the concave case
#			if(PHI > PHILimit){
#				A = A - abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
#				ip=phiLimit
#			}
#			else ip=phi
#			
#			A = A + abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda)$val)
#		}
#	}
#	A
#		
#
#}


gradientIntegralConvex <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	g = gradient(integralConvex,PHI,psi,lambda, modePHI, modepsi, modelambda, i);
	g[i]
 
}
gradientIntegralConcave <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	g = gradient(integralConcave,PHI,psi,lambda, modePHI, modepsi, modelambda, i);
	g[i]
 
}





gradient <- function(integrator, PHI, psi, lambda, modePHI, modepsi, modelambda, i){

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
	res[1] = (integrator(PHI+pfd_PHI, psi, lambda, modePHI, modepsi, modelambda, i) - integrator(PHI+nfd_PHI, psi, lambda, modePHI, modepsi, modelambda, i))/(2*fd)
	res[2] = (integrator(PHI, psi+pfd_psi, lambda, modePHI, modepsi, modelambda, i) - integrator(PHI, psi+nfd_psi, lambda, modePHI, modepsi, modelambda, i))/(2*fd)
	res[3] = (integrator(PHI, psi, lambda+pfd_lambda, modePHI, modepsi, modelambda, i) - integrator(PHI, psi, lambda+nfd_lambda, modePHI, modepsi, modelambda, i))/(2*fd)
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













#empty grid data
dataConvex = array(0,dim=c(res_PHI,res_psi,res_lambda))
dataConcave = array(0,dim=c(res_PHI,res_psi,res_lambda))



gradientsConvex = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda))
hessiansConvex = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))
gradientsGradientsConvex = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))
hessiansGradientsConvex = array(0,dim=c(dimensions, dimensions, dimensions, res_PHI,res_psi,res_lambda))

gradientsConcave = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda))
hessiansConcave = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))
gradientsGradientsConcave = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))
hessiansGradientsConcave = array(0,dim=c(dimensions, dimensions, dimensions, res_PHI,res_psi,res_lambda))

headerLambdaConvex=array(0,dim=c(res_lambda))
headerPsiConvex=array(0,dim=c(res_psi))
headerPHIConvex=array(0,dim=c(res_PHI))

headerLambdaConcave=array(0,dim=c(res_lambda))
headerPsiConcave=array(0,dim=c(res_psi))
headerPHIConcave=array(0,dim=c(res_PHI))


#first create the header
for(i_lambda in 0:(res_lambda-1)){
	lambda = max_lambda*i_lambda/(res_lambda-1)

	headerLambdaConvex[i_lambda+1]=lambda
	headerLambdaConcave[i_lambda+1]=lambda

	#iterate over psi angles
	for(i_psi in 0:(res_psi-1)){

		psi_convex = (i_psi * (max_psi)/(res_psi-1))
		psi_concave = (i_psi * (max_psi)/(res_psi-1))

		headerPsiConvex[i_psi+1]=psi_convex
		headerPsiConcave[i_psi+1]=psi_concave

		for(i_PHI in 0:(res_PHI-1)){
			PHI_convex = max_PHI * i_PHI/(res_PHI-1)
			PHI_concave = max_PHI * i_PHI/(res_PHI-1)

			headerPHIConvex[i_PHI+1]=PHI_convex
			headerPHIConcave[i_PHI+1]=PHI_concave
		}
	}
}


#iterate over lambda angles
for(i_lambda in 0:(res_lambda-1)){
	print(c(round(100*i_lambda/(res_lambda)),"% completed"))

	lambda = max_lambda*i_lambda/(res_lambda-1)


	#lambda==0 will not occur in reality. this point would have extreme derivatives, so we start with a circle of very small radius instead
	lambda = max(0.01,lambda)



	#iterate over psi angles
	for(i_psi in 0:(res_psi-1)){

		psi_convex = (i_psi * (max_psi)/(res_psi-1))
		psi_concave = (i_psi * (max_psi)/(res_psi-1))





		for(i_PHI in 0:(res_PHI-1)){
			PHI_convex = max_PHI * i_PHI/(res_PHI-1)
			PHI_concave = max_PHI * i_PHI/(res_PHI-1)



			#print(c("convex: ",i_PHI,PHI_convex,i_psi,psi_convex,i_lambda,lambda))
			#print(c("concave: ",i_PHI,PHI_concave,i_psi,psi_concave,i_lambda,lambda))

			#print("dataConvex")
			if(isWithinLimitsConvex(PHI_convex,psi_convex,lambda))
				dconv = integralConvex(PHI_convex,psi_convex,lambda)
			else dconv = NaN

			#print("dataConcave")
			#if(isWithinLimitsConcave(PHI_concave,psi_concave,lambda))
				#dconc = integralConcave(PHI_concave,psi_concave,lambda)
			#else dconc = NaN

			dataConvex[i_PHI+1, i_psi+1, i_lambda+1] = dconv
			dataConcave[i_PHI+1, i_psi+1, i_lambda+1] = dconc
			modePHI="central"
			modepsi="central"
			modelambda="central"



			if(i_PHI==0) modePHIConvex="forward"
			else if(i_PHI==(res_PHI-1)) modePHIConvex="backward"
			else modePHIConvex="central"
			if(i_psi==0) modepsiConvex="forward"
			else if(i_psi==(res_psi-1)) modepsiConvex="backward"
			else modepsiConvex="central"
			if(i_lambda==0) modelambdaConvex="forward"
			else if(i_lambda==(res_lambda-1)) modelambdaConvex="backward"
			else modelambdaConvex="central"
				
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

			#if(!isWithinLimitsConvex(PHI_convex+2*fd,psi_convex,lambda) && !isWithinLimitsConvex(PHI_convex-2*fd,psi_convex,lambda)) modePHIConvex="invalid"
			#else if(modePHIConvex!="backward" && !isWithinLimitsConvex(PHI_convex+2*fd,psi_convex,lambda)) modePHIConvex="backward"
			#else if(modePHIConvex!="forward" && !isWithinLimitsConvex(PHI_convex-2*fd,psi_convex,lambda)) modePHIConvex="forward"

			#if(!isWithinLimitsConvex(PHI_convex,psi_convex+2*fd,lambda) && !isWithinLimitsConvex(PHI_convex,psi_convex-2*fd,lambda)) modepsiConvex="invalid"
			#else if(modepsiConvex!="backward" && !isWithinLimitsConvex(PHI_convex,psi_convex+2*fd,lambda)) modepsiConvex="backward"
			#else if(modepsiConvex!="forward" && !isWithinLimitsConvex(PHI_convex,psi_convex-2*fd,lambda)) modepsiConvex="forward"

			#if(!isWithinLimitsConvex(PHI_convex,psi_convex,lambda+2*fd) && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda-2*fd)) modelambdaConvex="invalid"
			#else if(modelambdaConvex!="backward" && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda+2*fd)) modelambdaConvex="backward"
			#else if(modelambdaConvex!="forward" && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda-2*fd)) modelambdaConvex="forward"


			if(!isWithinLimitsConcave(PHI_concave+2*fd,psi_concave,lambda) && !isWithinLimitsConcave(PHI_concave-2*fd,psi_concave,lambda)) modePHIConcave="invalid"
			else if(modePHIConcave!="backward" && !isWithinLimitsConcave(PHI_concave+2*fd,psi_concave,lambda)) modePHIConcave="backward"
			else if(modePHIConcave!="forward" && !isWithinLimitsConcave(PHI_concave-2*fd,psi_concave,lambda)) modePHIConcave="forward"


			if(modepsiConcave!="backward" && isAtPositiveConcaveDiscontinuity(i_psi,i_lambda)) modepsiConcave="backward"
			else if(modepsiConcave!="forward" && isAtNegativeConcaveDiscontinuity(i_psi,i_lambda)) modepsiConcave="forward"

			if(modelambdaConcave!="backward" && isAtPositiveConcaveDiscontinuity(i_psi,i_lambda)) modelambdaConcave="backward"
			else if(modelambdaConcave!="forward" && isAtNegativeConcaveDiscontinuity(i_psi,i_lambda)) modelambdaConcave="forward"


			if(!isWithinLimitsConcave(PHI_concave,psi_concave+2*fd,lambda) && !isWithinLimitsConcave(PHI_concave,psi_concave-2*fd,lambda)) modepsiConcave="invalid"
			else if(modepsiConcave!="backward" && !isWithinLimitsConcave(PHI_concave,psi_concave+2*fd,lambda)) modepsiConcave="backward"
			else if(modepsiConcave!="forward" && !isWithinLimitsConcave(PHI_concave,psi_concave-2*fd,lambda)) modepsiConcave="forward"

			if(!isWithinLimitsConcave(PHI_concave,psi_concave,lambda+2*fd) && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda-2*fd)) modelambdaConcave="invalid"
			else if(modelambdaConcave!="backward" && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda+2*fd)) modelambdaConcave="backward"
			else if(modelambdaConcave!="forward" && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda-2*fd)) modelambdaConcave="forward"


			#print(c(i_PHI,i_psi,i_lambda))
			#print(c(modePHIConcave, modepsiConcave, modelambdaConcave))

			#print("gradientsConvex")

			#if(is.nan(dconv) || !derivableConvex) gconv = NaN
			#else gconv = gradient(integralConvex,PHI_convex,psi_convex,lambda, modePHIConvex, modepsiConvex, modelambdaConvex)
			#gradientsConvex[,i_PHI+1, i_psi+1, i_lambda+1] = gconv

			#print("gradientsConcave")

			if(is.nan(dconc) || !derivableConcave) gconc = NaN
			else gconc = gradient(integralConcave,PHI_concave,psi_concave,lambda, modePHIConcave, modepsiConcave, modelambdaConcave)
			gradientsConcave[,i_PHI+1, i_psi+1, i_lambda+1] = gconc

			#print("hessiansConvex")

			#if(is.nan(dconv) || !derivableConvex) hconv = NaN
			#else hconv = hessian(integralConvex, PHI_convex, psi_convex, lambda, modePHIConvex, modepsiConvex, modelambdaConvex)
			#hessiansConvex[,,i_PHI+1, i_psi+1, i_lambda+1] = hconv

			#print("hessiansConcave")

			if(is.nan(dconc) || !derivableConcave) hconc = NaN
			else hconc = hessian(integralConcave, PHI_concave, psi_concave, lambda, modePHIConcave, modepsiConcave, modelambdaConcave)
			hessiansConcave[,,i_PHI+1, i_psi+1, i_lambda+1] = hconc

			#print("done")

		}
	}
}
print(c(100,"% completed"))


if(!DODERIVATIVES) print("we will not generate tables for derivatives")
if(DODERIVATIVES) print("generating tables for both derivatives")

if(DODERIVATIVES){
#build tables for the components of the gradient
#iterate over lambda angles
for(i_lambda in 0:(res_lambda-1)){
	print(c(round(100*i_lambda/(res_lambda)),"% completed"))
	lambda = max_lambda*i_lambda/(res_lambda-1)



	#iterate over psi angles
	for(i_psi in 0:(res_psi-1)){

		psi_convex = (i_psi * (max_psi)/(res_psi-1))
		psi_concave = (i_psi * (max_psi)/(res_psi-1))




		for(i_PHI in 0:(res_PHI-1)){
			PHI_convex = max_PHI * i_PHI/(res_PHI-1)
			PHI_concave = max_PHI * i_PHI/(res_PHI-1)



			#print(c("convex: ",PHI_convex,psi_convex,lambda))
			#print(c("concave: ",PHI_concave,psi_concave,lambda))

			dconv = dataConvex[i_PHI+1, i_psi+1, i_lambda+1]
			dconc = dataConcave[i_PHI+1, i_psi+1, i_lambda+1]

			modePHI="central"
			modepsi="central"
			modelambda="central"



			if(i_PHI==0) modePHIConvex="forward"
			else if(i_PHI==(res_PHI-1)) modePHIConvex="backward"
			else modePHIConvex="central"
			if(i_psi==0) modepsiConvex="forward"
			else if(i_psi==(res_psi-1)) modepsiConvex="backward"
			else modepsiConvex="central"
			if(i_lambda==0) modelambdaConvex="forward"
			else if(i_lambda==(res_lambda-1)) modelambdaConvex="backward"
			else modelambdaConvex="central"
				
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

			#if(!isWithinLimitsConvex(PHI_convex+2*fd,psi_convex,lambda) && !isWithinLimitsConvex(PHI_convex-2*fd,psi_convex,lambda)) modePHIConvex="invalid"
			#else if(modePHIConvex!="backward" && !isWithinLimitsConvex(PHI_convex+2*fd,psi_convex,lambda)) modePHIConvex="backward"
			#else if(modePHIConvex!="forward" && !isWithinLimitsConvex(PHI_convex-2*fd,psi_convex,lambda)) modePHIConvex="forward"

			#if(!isWithinLimitsConvex(PHI_convex,psi_convex+2*fd,lambda) && !isWithinLimitsConvex(PHI_convex,psi_convex-2*fd,lambda)) modepsiConvex="invalid"
			#else if(modepsiConvex!="backward" && !isWithinLimitsConvex(PHI_convex,psi_convex+2*fd,lambda)) modepsiConvex="backward"
			#else if(modepsiConvex!="forward" && !isWithinLimitsConvex(PHI_convex,psi_convex-2*fd,lambda)) modepsiConvex="forward"

			#if(!isWithinLimitsConvex(PHI_convex,psi_convex,lambda+2*fd) && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda-2*fd)) modelambdaConvex="invalid"
			#else if(modelambdaConvex!="backward" && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda+2*fd)) modelambdaConvex="backward"
			#else if(modelambdaConvex!="forward" && !isWithinLimitsConvex(PHI_convex,psi_convex,lambda-2*fd)) modelambdaConvex="forward"


			if(!isWithinLimitsConcave(PHI_concave+2*fd,psi_concave,lambda) && !isWithinLimitsConcave(PHI_concave-2*fd,psi_concave,lambda)) modePHIConcave="invalid"
			else if(modePHIConcave!="backward" && !isWithinLimitsConcave(PHI_concave+2*fd,psi_concave,lambda)) modePHIConcave="backward"
			else if(modePHIConcave!="forward" && !isWithinLimitsConcave(PHI_concave-2*fd,psi_concave,lambda)) modePHIConcave="forward"

			if(!isWithinLimitsConcave(PHI_concave,psi_concave+2*fd,lambda) && !isWithinLimitsConcave(PHI_concave,psi_concave-2*fd,lambda)) modepsiConcave="invalid"
			else if(modepsiConcave!="backward" && !isWithinLimitsConcave(PHI_concave,psi_concave+2*fd,lambda)) modepsiConcave="backward"
			else if(modepsiConcave!="forward" && !isWithinLimitsConcave(PHI_concave,psi_concave-2*fd,lambda)) modepsiConcave="forward"

			if(!isWithinLimitsConcave(PHI_concave,psi_concave,lambda+2*fd) && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda-2*fd)) modelambdaConcave="invalid"
			else if(modelambdaConcave!="backward" && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda+2*fd)) modelambdaConcave="backward"
			else if(modelambdaConcave!="forward" && !isWithinLimitsConcave(PHI_concave,psi_concave,lambda-2*fd)) modelambdaConcave="forward"





			#for(i_c in 1:3){
		#		if(is.nan(dconv)) gconv = NaN
		#		else gconv = gradient(gradientIntegralConvex,PHI_convex,psi_convex,lambda, modePHIConvex, modepsiConvex, modelambdaConvex,i_c)
		#		gradientsGradientsConvex[i_c,,i_PHI+1, i_psi+1, i_lambda+1] = gconv
		#	}
			for(i_c in 1:3){
				if(is.nan(dconc)) gconc = NaN
				else gconc = gradient(gradientIntegralConcave,PHI_concave,psi_concave,lambda, modePHIConcave, modepsiConcave, modelambdaConcave,i_c)
				gradientsGradientsConcave[i_c,,i_PHI+1, i_psi+1, i_lambda+1] = gconc
			}

		#	for(i_c in 1:3){
		#		if(is.nan(dconv)) hconv = NaN
		#		else gconv = hessian(gradientIntegralConvex,PHI_convex,psi_convex,lambda, modePHIConvex, modepsiConvex, modelambdaConvex,i_c)
		#		hessiansGradientsConvex[i_c,,,i_PHI+1, i_psi+1, i_lambda+1] = hconv
		#	}
			for(i_c in 1:3){
				if(is.nan(dconc)) hconc = NaN
				else hconc = hessian(gradientIntegralConcave,PHI_concave,psi_concave,lambda, modePHIConcave, modepsiConcave, modelambdaConcave,i_c)
				hessiansGradientsConcave[i_c,,,i_PHI+1, i_psi+1, i_lambda+1] = hconc
			}


			

		}
	}
}
print(c(100,"% completed"))
}


print("saving data")
#now, save to file

source("floatconversion.R")

#saveTable("dataConvex.csv",c(res_PHI,res_psi,res_lambda),headerPHIConvex, headerPsiConvex, headerLambdaConvex,dataConvex, gradientsConvex, hessiansConvex)
saveTable("dataConcave.csv",c(res_PHI,res_psi,res_lambda),headerPHIConcave, headerPsiConcave, headerLambdaConcave,dataConcave, gradientsConcave, hessiansConcave)

if(DODERIVATIVES){
	print("saving derivatives and hessians")
	for(i in 1:dimensions){
		#saveTable(paste0("dataConvex",(i-1),".csv"),c(res_PHI,res_psi,res_lambda),headerPHIConvex, headerPsiConvex, headerLambdaConvex,gradientsConvex[i,,,], gradientsGradientsConvex[i,,,,], hessiansGradientsConvex[i,,,,,])
		saveTable(paste0("dataConcave",(i-1),".csv"),c(res_PHI,res_psi,res_lambda),headerPHIConcave, headerPsiConcave, headerLambdaConcave,gradientsConcave[i,,,], gradientsGradientsConcave[i,,,,], hessiansGradientsConcave[i,,,,,])
	}
}


#saveTableRaw("dataConvex.raw",c(res_PHI,res_psi,res_lambda),headers,dataConvex)

