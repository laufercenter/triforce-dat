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
fd=0.0001

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
		border = mvMultiply(rotz(-lambda), n)
		v = mvMultiply(rotz(-psi), n)
		a = border-n
		aRot = mvMultiply(rotx(PHI), a)
	
		ip = n + aRot
	
		nIntersection = crossproduct(v,ip)
		nIntersection = normalize(nIntersection)
		phiIntersection = normalAngle(nIntersection,nOrigin)
		res = phiIntersection
	}
	res

}


#normal vector in x-direction
n = matrix(c(1,0,0),3,1)

#normal vector to the plane connecting the origin, the integration origin and the interace center of the circular region
nOrigin = matrix(c(0,0,1),3,1)


csc <-function(x){
	1/sin(x)
}

#gives theta values dependent on phi, psi and lambda angles for a convex spherical arc
arcConvex <- function(phi, psi, lambda){
	v = 1+(-cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	#for(i in 1:length(v))
	#	if(is.nan(v[i])) v[i]=0
	v
}

#gives theta values dependent on phi, psi and lambda angles for a concave spherical arc
arcConcave <- function(phi, psi, lambda){
	v = 1-(cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)	
	#for(i in 1:length(v))
	#	if(is.nan(v[i])) v[i]=0
	v
}


integralConvex <-function(PHI, psi, lambda){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#calculate the maximal phi value to which we can integrate
	phiLimit = PHI2phi(pi/2,psi,lambda)


	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<=lambda){
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
		#for the convex case
		if(PHI < pi/2){
			A = A + abs(integrate(arcConcave,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
			ip=phiLimit
		}
		else ip=phi
		
		A = A + abs(integrate(arcConvex,lower=0,upper=ip,psi=psi,lambda=lambda)$val)
	}
	A
}

integralConcave <-function(PHI, psi, lambda){
	#convert PHI to phi
	phi = PHI2phi(PHI,psi,lambda)
	#calculate the maximal phi value to which we can integrate
	phiLimit = PHI2phi(pi/2,psi,lambda)


	A=0

	#if the integration origin is inside the circular area, we have to be careful how the integration limits are set
	if(psi<=lambda){
		#this is empty on purpose, calculating these values makes only theoretically sense, we will not use them
	}
	else{
		#for the concave case
		if(PHI > pi/2){
			A = A + abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
			ip=phiLimit
		}
		else ip=phi
		
		A = A + abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda)$val)
	}
	A
		

}


gradient <- function(integrator, PHI, psi, lambda){
	res = array(0,dim=c(3))
	res[1] = (integrator(PHI+fd, psi, lambda) - integrator(PHI-fd, psi, lambda))/(2*fd)
	res[2] = (integrator(PHI, psi+fd, lambda) - integrator(PHI, psi-fd, lambda))/(2*fd)
	res[3] = (integrator(PHI, psi, lambda+fd) - integrator(PHI, psi, lambda-fd))/(2*fd)
	res
	
}



hessian <- function(integrator, PHI, psi, lambda, f){
	res = array(0,dim=c(3,3))
	#f PHI PHI
	res[1,1] = (integrator(PHI+fd, psi, lambda) - 2*f + integrator(PHI-fd, psi, lambda))/fd^2
	#f PHI psi
	res[1,2] = (integrator(PHI+fd, psi+fd, lambda) - integrator(PHI+fd, psi-fd, lambda) - integrator(PHI-fd, psi+fd, lambda) + integrator(PHI-fd, psi-fd, lambda))/(4*fd^2)
	#f PHI lambda
	res[1,3] = (integrator(PHI+fd, psi, lambda+fd) - integrator(PHI+fd, psi, lambda-fd) - integrator(PHI-fd, psi, lambda+fd) + integrator(PHI-fd, psi, lambda-fd))/(4*fd^2)

	#f psi PHI
	res[2,1] = res[1,2]
	#f psi psi
	res[2,2] = (integrator(PHI, psi+fd, lambda) - 2*f + integrator(PHI, psi-fd, lambda))/fd^2
	#f psi lambda
	res[2,3] = (integrator(PHI, psi+fd, lambda+fd) - integrator(PHI, psi+fd, lambda-fd) - integrator(PHI, psi-fd, lambda+fd) + integrator(PHI, psi-fd, lambda-fd))/(4*fd^2)


	#f lambda  PHI
	res[3,1] = res[1,3]
	#f lambda psi
	res[3,2] = res[2,3]
	#f psi psi
	res[3,3] = (integrator(PHI, psi, lambda+fd) - 2*f + integrator(PHI, psi, lambda-fd))/fd^2

	res
}

















#empty grid data
dataConvex = array(0,dim=c(res_PHI,res_psi,res_lambda))
dataConcave = array(0,dim=c(res_PHI,res_psi,res_lambda))

gradientsConvex = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda))
hessiansConvex = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))

gradientsConcave = array(0,dim=c(dimensions, res_PHI,res_psi,res_lambda))
hessiansConcave = array(0,dim=c(dimensions, dimensions, res_PHI,res_psi,res_lambda))

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
				
			gradientsConvex[,i_PHI+1, i_psi+1, i_lambda] = gradient(integralConvex,PHI,psi,lambda)
			gradientsConcave[,i_PHI+1, i_psi+1, i_lambda] = gradient(integralConcave,PHI,psi,lambda)

			hessiansConvex[,,i_PHI+1, i_psi+1, i_lambda] = hessian(integralConvex, PHI, psi, lambda, dconv)
			hessiansConcave[,,i_PHI+1, i_psi+1, i_lambda] = hessian(integralConcave, PHI, psi, lambda, dconc)

		}
	}
}