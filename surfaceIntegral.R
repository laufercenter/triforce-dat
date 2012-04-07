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



res_lambda = 10
res_psi = 10
res_PHI = 10

max_lambda = pi
max_psi = pi
max_PHI = pi


rotz <- function(theta){
	matrix(c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1),3,3)
}

rotx <- function(theta){
	matrix(c(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta)),3,3)
}


crossproduct <-function(a,b){
	c = matrix(0,3,1)
	c[1] = a[2]*b[3] - a[3]*b[2]
	c[2] = a[3]*b[1] - a[1]*b[3]
	c[3] = a[1]*b[2] - a[2]*b[1]
	c
}

mvMultiply <- function(A,b){
	matrix(c(A %*% b),3,1)
}

normalize <- function(a){
	a / sqrt(c(t(a) %*% a))
}

normalAngle <-function(a,b){
	acos(c(t(a) %*% b))
}

PHI2phi <- function(PHI, psi, lambda){

	#print("######################")

	border = mvMultiply(rotz(-lambda), n)
	#print(border)
	#print("------------------")
	a = border-n
	#print(a)
	#print("------------------")
	aRot = mvMultiply(rotx(PHI), a)
	#print(aRot)
	#print("------------------")


	nIntersection = crossproduct(aRot,n)
	#print(nIntersection)
	#print("------------------")
	nIntersection = normalize(nIntersection)
	#print(nIntersection)
	#print("------------------")
	phiIntersection = normalAngle(nIntersection,nOrigin)
	#print(phiIntersection)
	#print("------------------")
	
	#print("######################")
	phiIntersection
}


nl = matrix(c(0,0,1),3,1)

n = matrix(c(1,0,0),3,1)

#normal vector for the origin plane
nOrigin = matrix(c(0,0,1),3,1)


arcConvex <- function(phi, psi, lambda){
	1+(-cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
}

arcConcave <- function(phi, psi, lambda){
	1-(cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)	
}


dataConvex = array(0,dim=c(res_PHI+1,res_PHI+1,res_psi+1,res_lambda))
dataConcave = array(0,dim=c(res_PHI+1,res_PHI+1,res_psi+1,res_lambda))

for(i_psi in 0:res_psi){
	psi = max_psi*i_psi/res_psi
	for(i_lambda in 1:(res_lambda-1)){
		lambda = max_lambda*i_lambda/res_lambda

		phiLimit0 = PHI2phi(-pi/2,psi,lambda)
		phiLimit1 = PHI2phi(pi/2,psi,lambda)

		for(i_PHI0 in 0:res_PHI){
			PHI0 = max_PHI * (i_PHI0-(res_PHI/2))/res_PHI

			phi0 = PHI2phi(PHI0,psi,lambda)



			for(i_PHI1 in 0:res_PHI){
				PHI1 = max_PHI * (i_PHI1-(res_PHI/2))/res_PHI

				phi1 = PHI2phi(PHI1,psi,lambda)



				#for the concave case
				A=0
				if(PHI0 < pi/4){
					A = A + abs(integrate(arcConvex,lower=phi0,upper=phiLimit0,psi=psi,lambda=lambda)$val)
					ip0=phiLimit0
				}
				else ip0=phi0

				if(PHI1 > pi/4){
					A = A + abs(integrate(arcConvex,lower=phi1,upper=phiLimit1,psi=psi,lambda=lambda)$val)
					ip1=phiLimit1
				}
				else ip1=phi1
				
				A = A + abs(integrate(arcConcave,lower=ip0,upper=ip1,psi=psi,lambda=lambda)$val)

				dataConcave[i_PHI0+1, i_PHI1+1, i_psi+1, i_lambda+1] = A

				#for the convex case
				A=0
				if(PHI0 > pi/4){
					A = A + abs(integrate(arcConcave,lower=phi0,upper=phiLimit0,psi=psi,lambda=lambda)$val)
					ip0=phiLimit0
				}
				else ip0=phi0

				if(PHI1 < pi/4){
					A = A + abs(integrate(arcConcave,lower=phi1,upper=phiLimit1,psi=psi,lambda=lambda)$val)
					ip1=phiLimit1
				}
				else ip1=phi1
				
				A = A + abs(integrate(arcConvex,lower=ip0,upper=ip1,psi=psi,lambda=lambda)$val)

				dataConvex[i_PHI0+1, i_PHI1+1, i_psi+1, i_lambda+1] = A
				

				
			}
		}
	}
}