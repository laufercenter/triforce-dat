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



res_lambda = 35
res_psi = 35
res_PHI = 35

max_lambda = pi/2
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

	if(PHI==0){
		res=0
	}
	else{
	

	#print("######################")

	border = mvMultiply(rotz(-lambda), n)
	v = mvMultiply(rotz(-psi), n)
	#print(border)
	#print("------------------")
	a = border-n
	#print(a)
	#print("------------------")
	aRot = mvMultiply(rotx(PHI), a)
	#print(aRot)
	#print("------------------")

	ip = n + aRot



	nIntersection = crossproduct(v,ip)
	#print(nIntersection)
	#print("------------------")
	nIntersection = normalize(nIntersection)
	#print(nIntersection)
	#print("------------------")
	phiIntersection = normalAngle(nIntersection,nOrigin)
	#print(phiIntersection)
	#print("------------------")
	
	#print("######################")
	res = phiIntersection
	}
	res

}


nl = matrix(c(0,0,1),3,1)

n = matrix(c(1,0,0),3,1)

#normal vector for the origin plane
nOrigin = matrix(c(0,0,1),3,1)


csc <-function(x){
	1/sin(x)
}


arcConvex <- function(phi, psi, lambda){
	v = 1+(-cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
	#for(i in 1:length(v))
	#	if(is.nan(v[i])) v[i]=0
	v
}

arcConcave <- function(phi, psi, lambda){
	v = 1-(cos(lambda) * cos(psi)+sqrt(cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)	
	#for(i in 1:length(v))
	#	if(is.nan(v[i])) v[i]=0
	v
}

phiLimitNegative <-function(psi,lambda){
	-acos(sqrt((cos(lambda)^2-cos(psi)^2) * csc(psi)^2))
}

phiLimitPositive <-function(psi,lambda){
	acos(sqrt((cos(lambda)^2-cos(psi)^2) * csc(psi)^2))
}

findLimit <-function(psi,lambda){
	r = 0
	l = -pi

	for(i in 0:1000){
		m <- (l-r)/2
		t <- arcConvex(m,psi,lambda)
		if(is.nan(t)){
			l<-m
		}
		else{
			r<-m
		}
	}
	r
}

dataConvex = array(0,dim=c(res_PHI+1,res_psi+1,res_lambda))
dataConcave = array(0,dim=c(res_PHI+1,res_psi+1,res_lambda))


pconvex <- function(psi,lambda){
res=1000
at = array(0,dim=c(res))
h = array(0,dim=c(res))
for(i in 1:res){
	at[i]=arcConvex(max_PHI*2 * (i-(res/2))/res,psi,lambda)
	h[i]=max_PHI*2 * (i-(res/2))/res
}
plot(h,at)
}

for(i_psi in 0:res_psi){
	print(i_psi)
	psi = max_psi*i_psi/res_psi
	for(i_lambda in 1:(res_lambda-1)){
		lambda = max_lambda*i_lambda/res_lambda

		phiLimit = PHI2phi(pi/2,psi,lambda)

		for(i_PHI in 0:res_PHI){
			PHI = max_PHI * i_PHI/res_PHI

			phi = PHI2phi(PHI,psi,lambda)


			#for the concave case
			A=0
			if(PHI > pi/2){
				A = A + abs(integrate(arcConvex,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A + abs(integrate(arcConcave,lower=0,upper=ip,psi=psi,lambda=lambda)$val)

			dataConcave[i_PHI+1, i_psi+1, i_lambda+1] = A

			#for the convex case
			A=0
			if(PHI < pi/2){
				A = A + abs(integrate(arcConcave,lower=phi,upper=phiLimit,psi=psi,lambda=lambda)$val)
				ip=phiLimit
			}
			else ip=phi
			
			A = A + abs(integrate(arcConvex,lower=0,upper=ip,psi=psi,lambda=lambda)$val)

			dataConvex[i_PHI+1, i_psi+1, i_lambda+1] = A
				

				
			
		}
	}
}