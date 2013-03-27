
resolution= 10
grid=12

res_lambda = grid
res_psi = grid
res_phi = resolution

#limits for the parameters
max_lambda = pi/2
max_psi = pi
max_phi = pi

dimensions=3




saveTable <-function(filename, dimensions, headerPhi, headerPsi, headerLambda, tbl){

	DERIVATIVE_LEVEL = 0
	STORE_PHI = 0
	headerPhidat = array(0,dim=c(8,length(headerPhi)))
	headerPsidat = array(0,dim=c(8,length(headerPsi)))
	headerLambdadat = array(0,dim=c(8,length(headerLambda)))
	tbldat = array(0,dim=c(8,length(tbl)))
	for(i in 1:length(headerPhi)){
		headerPhidat[,i] = double2CharArray(headerPhi[i])
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

	write(file=filename,DERIVATIVE_LEVEL, append=FALSE)
	write(file=filename,STORE_PHI, append=TRUE)
	write(file=filename,3,append=TRUE)
	write(file=filename,c(1,res_phi,1,res_psi,1,res_lambda),ncol=3,append=TRUE)
	write(file=filename,headerPhidat,ncol=8,append=TRUE)
	write(file=filename,headerPsidat,ncol=8,append=TRUE)
	write(file=filename,headerLambdadat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
}




csc <-function(x){
	1/sin(x)
}

cot <- function(x){
	1/tan(x)
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



#gives theta values dependent on phi, psi and lambda angles for a convex spherical arc
thetaConvex <- function(phi, psi, lambda){
	phi=abs(phi)

	if(lambda==0){
		v = 0
	}
	else{
		s = cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		if(s<0){
			print(c("perturbation error: s=",s))
			s = 0
		}
		a = (cos(lambda) * cos(psi)-sqrt(s))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		if(abs(a)>1){
			print(c("boundary error: a=",a))
			a = sign(a)
		}
		v = acos(a)
	}

	v
}

#gives theta values dependent on phi, psi and lambda angles for a concave spherical arc
thetaConcave <- function(phi, psi, lambda){
	phi=abs(phi)
	if(lambda==0){
		v = 0
	}
	else{
		s = cos(phi)^2 * sin(psi)^2 * (-cos(lambda)^2+cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		if(s<0){
			print(c("perturbation error: s=",s))
			s = 0
		}
		a = (cos(lambda) * cos(psi)+sqrt(s))/(cos(psi)^2+cos(phi)^2 * sin(psi)^2)
		if(abs(a)>1){
			print(c("boundary error: a=",a))
			a=sign(a)
		}
		
		v = acos(a)
	}



	v
}





calculateDepth <- function(phi, psi, lambda){
	phi=abs(phi)

	if(psi<lambda){
		if(phi <= phiLimitPlain(psi,lambda)){
			forw = thetaConvex(phi, psi, lambda)
			backw = 2*pi - thetaConcave(phi+pi, psi, lambda)
		}
		else{
			forw = thetaConcave(phi, psi, lambda)
			backw = 2*pi - thetaConvex(phi+pi, psi, lambda)
		}
	}
	else if(psi+lambda >= pi){
		if(phi <= phiLimitPlain(psi,lambda)){
			forw = thetaConcave(phi, psi, lambda)
			backw = 2*pi - thetaConvex(phi+pi, psi, lambda)
		}
		else{
			forw = thetaConvex(phi, psi, lambda)
			backw = 2*pi - thetaConcave(phi+pi, psi, lambda)
		}
	}
	else{
		if(phi <= phiLimitPlain(psi,lambda)){
			forw = thetaConcave(phi, psi, lambda)
			backw = thetaConvex(phi, psi, lambda)
		}
		else{
			forw=NaN
			backw=NaN
		}
			
	}
	
	c(forw,backw)
}


draw <-function(psi,lambda){
	plot(dataForward[,psi,lambda],type="l",ylim=c(0,2*pi))
	lines(dataBackward[,psi,lambda],ylim=c(0,2*pi))
}

drawf <-function(psi,lambda){
	plot(dataForward[,psi,lambda],type="l",ylim=c(0,2*pi))
}

drawb <-function(psi,lambda){
	plot(dataBackward[,psi,lambda],type="l",ylim=c(0,2*pi))
}



dataForward = array(0,dim=c(res_phi,res_psi,res_lambda))
dataBackward = array(0,dim=c(res_phi,res_psi,res_lambda))

headerLambda=array(0,dim=c(res_lambda))
headerPsi=array(0,dim=c(res_psi))
headerPhi=array(0,dim=c(res_phi))



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

for(i_phi in 0:(res_phi-1)){

	phi = (i_phi * (max_phi)/(res_phi-1)) -(pi/2)
	headerPhi[i_phi+1]=phi
}


for(i_lambda in 0:(res_lambda-1)){
	lambda = max_lambda*i_lambda/(res_lambda-1)
	
	print(c(round(100*i_lambda/(res_lambda)),"% completed"))
	

	#iterate over psi angles
	for(i_psi in 0:(res_psi-1)){
		psi = (i_psi * (max_psi)/(res_psi-1))
		
		for(i_phi in 0:(res_phi-1)){
		phi = (i_phi * (max_phi)/(res_phi-1)) - (pi/2)
			d = calculateDepth(phi,psi,lambda)
			dataForward[i_phi+1, i_psi+1, i_lambda+1] = d[1]
			dataBackward[i_phi+1, i_psi+1, i_lambda+1] = d[2]
		}
		
		
	}
}
print(c(100,"% completed"))

source("floatconversion.R")
saveTable("depthBufferForward.csv",c(res_phi,res_psi,res_lambda),headerPhi, headerPsi, headerLambda, dataForward)
saveTable("depthBufferBackward.csv",c(res_phi,res_psi,res_lambda),headerPhi, headerPsi, headerLambda, dataBackward)






