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





#the resolution of our grid
res_lambda = 10
res_psi = 10
res_PHI = 10

#limits for the parameters
max_lambda = pi/2
max_psi = pi
max_PHI = pi

dimensions=3
fd=0.00001


ex = matrix(c(1,0,0),3,1)
ey = matrix(c(0,1,0),3,1)
ez = matrix(c(0,0,1),3,1)

#normal vector to the plane connecting the origin, the integration origin and the interace center of the circular region
nOrigin = matrix(c(0,0,1),3,1)


saveTable <-function(filename, dimensions, header, tbl, tblgradients, tblhessians){

	headerdat = array(0,dim=c(8,length(header)))
	tbldat = array(0,dim=c(8,length(tbl)))
	tblgradientsdat = array(0,dim=c(8,length(tblgradients)))
	tblhessiansdat = array(0,dim=c(8,length(tblhessians)))
	write(file=filename,c(length(dimensions),nrow(header)),append=FALSE)
	write(file=filename,dimensions,append=TRUE)
	for(i in 1:length(header)){
		headerdat[,i] = double2CharArray(header[i])
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

	write(file=filename,headerdat,ncol=8,append=TRUE)
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
normalize <- function(a){
	a / sqrt(c(t(a) %*% a))
}

#gives angle between two normal vectors
normalAngle <-function(a,b){
	acos(c(t(a) %*% b))
}

#converts a PHI angle into a phi angle
PHI2phi <- function(PHI, psi, lambda){


	if(PHI==0){
		#stability
		res=0
	}
	else{
		border = mvMultiply(rotz(-lambda), ex)
		v = mvMultiply(rotz(-psi), ex)
		a = border-ex
		aRot = mvMultiply(rotx(PHI), a)
	
		ip = ex + aRot
	
		nIntersection = crossproduct(v,ip)
		nIntersection = normalize(nIntersection)
		phiIntersection = normalAngle(nIntersection,nOrigin)
		res = phiIntersection
	}



	#res = acos((-cos(PHI) * cos(psi) * sin(lambda)+cos(lambda) * sin(psi))/(sqrt(cos(psi)^2 * sin(lambda)^2-2 * cos(PHI) * cos(lambda) * cos(psi) * sin(lambda) * sin(psi)+(cos(lambda)^2+sin(PHI)^2 * sin(lambda)^2) * sin(psi)^2)))


	res


}

tmp=c(0,0)

maxPHI <- function(psi, lambda){
	tmp <<- c(psi,lambda)

	if(psi>=pi/2) res=pi/2
	else if(lambda==pi/2) res= pi/2
	else if(psi<lambda) res = pi/2
	else res = acos(cot(psi)*tan(lambda))

	res
}


halfSpherePHI <- function(psi,lambda)
{
	rho = acos(cos(lambda)*csc(psi))
	v = mvMultiply(rotx(rho),ey)
	n = mvMultiply(rotz(psi),ex) * cos(lambda)
	t0 = normalize(v-n)
	v2 = mvMultiply(rotz(psi-lambda),ex)
	t1 = normalize(v2-n)
	res = normalAngle(t0,t1)

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
	v =1+(-cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	#for(i in 1:length(v))
	#	v[i] = min(pi/2,v[i])
	v
}

#gives theta values dependent on phi, psi and lambda angles for a concave spherical arc
arcConcave <- function(phi, psi, lambda){
	v = 1-(cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	#for(i in 1:length(v))
	#	v[i] = min(pi/2,v[i])
	v
}


phiLimitAbs <-function(psi,lambda){
	#print(c(psi,lambda))
	#if(abs(cos(psi)) >= abs(cos(lambda))) limit=pi/2
	#else limit=acos(sqrt((cos(lambda)^2-cos(psi)^2) * csc(psi)^2))

	limit = PHI2phi(maxPHI(psi,lambda),psi,lambda)

	limit
}


integralConvex <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#calculate the maximal phi value to which we can integrate
	#phiLimit = PHI2phi(pi/2,psi,lambda)


	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<=lambda){

		phiLimit = phiLimitAbs(psi,lambda)
	
		if(PHI==0) phi=pi


		#calculations for convex arcs
		if(phi<=pi/2){
			A = A + abs(integrate(arcConvex,lower=0,upper=phi,psi=psi,lambda=lambda)$val)
		}
		else{
			A = A + abs(integrate(arcConvex,lower=0,upper=pi/2,psi=psi,lambda=lambda)$val)
			A = A + abs(integrate(arcConcave,lower=pi/2,upper=pi-phi,psi=psi,lambda=lambda)$val)
		}
		
	}
	else{
		#there are no convex regions above pi/2 (they are handled by the mirror-side)
		if(psi>pi/2){
			A=NaN
		}
		else{

			#we reject every PHI that is "on the other side of the half-sphere"
			if(psi-lambda>=pi/2 || (psi+lambda>=pi/2 && PHI > halfSpherePHI(psi,lambda))){
				A= NaN
			}
			else{
				
				phiLimit = phiLimitAbs(psi,lambda)
				PHILimit = maxPHI(psi,lambda)
				#for the convex case
				if(PHI < PHILimit){
					A = A - abs(integrate(arcConcave,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
					ip=phiLimit
				}
				else ip=phi
				
				A = A + abs(integrate(arcConvex,lower=0,upper=ip,psi=psi,lambda=lambda)$val)
			}
		}
	}
	A
}

integralConcave <-function(PHI, psi, lambda, modePHI, modepsi, modelambda, i){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#calculate the maximal phi value to which we can integrate
	
	

	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<=lambda){
		phiLimit = phiLimitAbs(psi,lambda)
		#this is empty on purpose, calculating these values makes only theoretically sense, we will not use them
	}
	else{
		if(psi-lambda>=pi/2 || (psi+lambda>=pi/2 && PHI > halfSpherePHI(psi,lambda))){
			A= NaN
		}
		else{


			phiLimit = phiLimitAbs(psi,lambda)
			PHILimit = maxPHI(psi,lambda)
			#for the concave case
			if(PHI > PHILimit){
				A = A + abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A + abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda)$val)
		}
	}
	A
		

}


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
	)
	switch(modepsi, 
		central={pfd_psi=fd; nfd_psi= -fd; cfd_psi=0},
		forward={pfd_psi=2*fd; nfd_psi=0; cfd_psi=fd},
		backward={pfd_psi=0; nfd_psi= -2*fd; cfd_psi= -fd},
	)
	switch(modelambda, 
		central={pfd_lambda=fd; nfd_lambda= -fd; cfd_lambda=0},
		forward={pfd_lambda=2*fd; nfd_lambda=0; cfd_lambda=fd},
		backward={pfd_lambda=0; nfd_lambda= -2*fd; cfd_lambda= -fd},
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
	)
	switch(modepsi, 
		central={pfd_psi=fd; nfd_psi= -fd; cfd_psi=0},
		forward={pfd_psi=2*fd; nfd_psi=0; cfd_psi=fd},
		backward={pfd_psi=0; nfd_psi= -2*fd; cfd_psi= -fd},
	)
	switch(modelambda, 
		central={pfd_lambda=fd; nfd_lambda= -fd; cfd_lambda=0},
		forward={pfd_lambda=2*fd; nfd_lambda=0; cfd_lambda=fd},
		backward={pfd_lambda=0; nfd_lambda= -2*fd; cfd_lambda= -fd},
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

headers=array(0,dim=c(dimensions,max(res_PHI,res_psi,res_lambda)))

#iterate over psi angles
for(i_psi in 0:(res_psi-1)){
	print(i_psi)
	psi = max_psi*i_psi/(res_psi-1)

	headers[2,i_psi+1]=psi

	#iterate over lambda angles
	for(i_lambda in 1:(res_lambda)){
		lambda = max_lambda*i_lambda/res_lambda

		headers[3,i_lambda]=lambda



		for(i_PHI in 0:(res_PHI-1)){
			PHI = max_PHI * i_PHI/(res_PHI-1)

			headers[1,i_PHI+1]=PHI



			dconv = integralConvex(PHI,psi,lambda)
			dconc = integralConcave(PHI,psi,lambda)

			dataConvex[i_PHI+1, i_psi+1, i_lambda] = dconv
			dataConcave[i_PHI+1, i_psi+1, i_lambda] = dconc
			modePHI="central"
			modepsi="central"
			modelambda="central"

			if(i_PHI==0) modePHI="forward"
			else if(i_PHI==(res_PHI-1)) modePHI="backward"
			else modePHI="central"
			if(i_psi==0) modepsi="forward"
			else if(i_psi==(res_psi-1)) modepsi="backward"
			else modepsi="central"
			if(i_lambda==1) modelambda="forward"
			else if(i_lambda==res_lambda) modelambda="backward"
			else modelambda="central"
				
			if(is.nan(dconv)) gconv = NaN
			else gconv = gradient(integralConvex,PHI,psi,lambda, modePHI, modepsi, modelambda)
			gradientsConvex[,i_PHI+1, i_psi+1, i_lambda] = gconv

			if(is.nan(dconc)) gconc = NaN
			else gconc = gradient(integralConcave,PHI,psi,lambda, modePHI, modepsi, modelambda)
			gradientsConcave[,i_PHI+1, i_psi+1, i_lambda] = gconc

			if(is.nan(dconv)) hconv = NaN
			else hconv = hessian(integralConvex, PHI, psi, lambda, modePHI, modepsi, modelambda)
			hessiansConvex[,,i_PHI+1, i_psi+1, i_lambda] = hconv

			if(is.nan(dconc)) hconc = NaN
			else hconc = hessian(integralConcave, PHI, psi, lambda, modePHI, modepsi, modelambda)
			hessiansConcave[,,i_PHI+1, i_psi+1, i_lambda] = hconc

			

		}
	}
}


#build tables for the components of the gradient
for(i_psi in 0:(res_psi-1)){
	print(i_psi)
	psi = max_psi*i_psi/(res_psi-1)

	#iterate over lambda angles
	for(i_lambda in 1:(res_lambda)){
		lambda = max_lambda*i_lambda/res_lambda

		for(i_PHI in 0:(res_PHI-1)){
			PHI = max_PHI * i_PHI/(res_PHI-1)

			dconv = dataConvex[i_PHI+1, i_psi+1, i_lambda]
			dconc = dataConcave[i_PHI+1, i_psi+1, i_lambda]


			modePHI="central"
			modepsi="central"
			modelambda="central"

			if(i_PHI==0) modePHI="forward"
			else if(i_PHI==(res_PHI-1)) modePHI="backward"
			else modePHI="central"
			if(i_psi==0) modepsi="forward"
			else if(i_psi==(res_psi-1)) modepsi="backward"
			else modepsi="central"
			if(i_lambda==1) modelambda="forward"
			else if(i_lambda==res_lambda) modelambda="backward"
			else modelambda="central"

			for(i_c in 1:3){
				if(is.nan(dconv)) gconv = NaN
				else gconv = gradient(gradientIntegralConvex,PHI,psi,lambda, modePHI, modepsi, modelambda,i_c)
				gradientsGradientsConvex[i_c,,i_PHI+1, i_psi+1, i_lambda] = gconv
			}
			for(i_c in 1:3){
				if(is.nan(dconv)) gconv = NaN
				else gconc = gradient(gradientIntegralConcave,PHI,psi,lambda, modePHI, modepsi, modelambda,i_c)
				gradientsGradientsConcave[i_c,,i_PHI+1, i_psi+1, i_lambda] = gconc
			}

			for(i_c in 1:3){
				if(is.nan(dconv)) hconv = NaN
				else gconv = hessian(gradientIntegralConvex,PHI,psi,lambda, modePHI, modepsi, modelambda,i_c)
				hessiansGradientsConvex[i_c,,,i_PHI+1, i_psi+1, i_lambda] = hconv
			}
			for(i_c in 1:3){
				if(is.nan(dconv)) hconv = NaN
				else hconc = hessian(gradientIntegralConcave,PHI,psi,lambda, modePHI, modepsi, modelambda,i_c)
				hessiansGradientsConcave[i_c,,,i_PHI+1, i_psi+1, i_lambda] = hconc
			}


			

		}
	}
}



#now, save to file

source("floatconversion.R")

saveTable("dataConvex.csv",c(res_PHI,res_psi,res_lambda),headers,dataConvex, gradientsConvex, hessiansConvex)
saveTable("dataConcave.csv",c(res_PHI,res_psi,res_lambda),headers,dataConcave, gradientsConcave, hessiansConcave)

for(i in 1:dimensions){
	saveTable(paste0("dataConvex",(i-1),".csv"),c(res_PHI,res_psi,res_lambda),headers,gradientsConvex[i,,,], gradientsGradientsConvex[i,,,,], hessiansGradientsConvex[i,,,,,])
	saveTable(paste0("dataConcave",(i-1),".csv"),c(res_PHI,res_psi,res_lambda),headers,gradientsConcave[i,,,], gradientsGradientsConcave[i,,,,], hessiansGradientsConcave[i,,,,,])
}


#saveTableRaw("dataConvex.raw",c(res_PHI,res_psi,res_lambda),headers,dataConvex)


#write(file="dataConvex.csv",c(res_PHI,res_psi,res_lambda),ncol=dimensions,append=FALSE)
#write(file="dataConvex.csv",headers,ncol=ncol(headers),append=TRUE)
#write(file="dataConcave.csv",c(res_PHI,res_psi,res_lambda),ncol=dimensions,append=FALSE)
#write(file="dataConcave.csv",headers,ncol=ncol(headers),append=TRUE)
#write(file="gradientsConvex.csv",c(dimensions,res_PHI,res_psi,res_lambda),ncol=dimensions+1,append=FALSE)
#write(file="gradientsConvex.csv",headers,ncol=ncol(headers),append=TRUE)
#write(file="gradientsConcave.csv",c(dimensions,res_PHI,res_psi,res_lambda),ncol=dimensions+1,append=FALSE)
#write(file="gradientsConcave.csv",headers,ncol=ncol(headers),append=TRUE)
#write(file="hessiansConvex.csv",c(dimensions,dimensions,res_PHI,res_psi,res_lambda),ncol=dimensions+2,append=FALSE)
#write(file="hessiansConvex.csv",headers,ncol=ncol(headers),append=TRUE)
#write(file="hessiansConcave.csv",c(dimensions,dimensions,res_PHI,res_psi,res_lambda),ncol=dimensions+2,append=FALSE)
#write(file="hessiansConcave.csv",headers,ncol=ncol(headers),append=TRUE)

#for(i_lambda in 1:res_lambda){
#	write(file="dataConvex.csv",dataConvex[,,i_lambda],ncol=res_psi,append=TRUE)
#	write(file="dataConcave.csv",dataConcave[,,i_lambda],ncol=res_psi,append=TRUE)
#
#	for(i_psi in 1:res_psi){
#		write(file="gradientsConvex.csv",gradientsConvex[,,i_psi,i_lambda],ncol=res_PHI,append=TRUE)
#		write(file="gradientsConcave.csv",gradientsConcave[,,i_psi,i_lambda],ncol=res_PHI,append=TRUE)
#
#		for(i_PHI in 1:res_PHI){
#			write(file="hessiansConvex.csv",hessiansConvex[,,i_PHI,i_psi,i_lambda],ncol=dimensions,append=TRUE)
#			write(file="hessiansConcave.csv",hessiansConcave[,,i_PHI,i_psi,i_lambda],ncol=dimensions,append=TRUE)
#		}
#	}
#}



