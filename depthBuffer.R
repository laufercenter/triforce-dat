
resolution= 25
grid=12

res_lambda = grid
res_psi = grid
res_g = resolution

#limits for the parameters
max_lambda = pi/2
max_psi = pi
max_g = 1

dimensions=3




saveTable <-function(filename, dimensions, headerG, headerPsi, headerLambda, tbl){

	DERIVATIVE_LEVEL = 0
	STORE_PHI = 0
	headerGdat = array(0,dim=c(8,length(headerG)))
	headerPsidat = array(0,dim=c(8,length(headerPsi)))
	headerLambdadat = array(0,dim=c(8,length(headerLambda)))
	tbldat = array(0,dim=c(8,length(tbl)))
	for(i in 1:length(headerG)){
		headerGdat[,i] = double2CharArray(headerG[i])
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
	write(file=filename,c(1,res_g,1,res_psi,1,res_lambda),ncol=3,append=TRUE)
	write(file=filename,headerGdat,ncol=8,append=TRUE)
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


kappa <- function(g, psi, lambda){
	if(abs(g)==1) g=g*0.999
	a = ((cos(lambda) - g*cos(psi)) * csc(psi)) / (sqrt(1-g^2))
	if(abs(a)>1) res = NaN
	else res= acos(a)
	res
}


calculateDepth <-function(g,psi,lambda){
	k=kappa(g, psi, lambda)
	if(is.nan(k)){
		if(psi>pi/2) q = -1
		else q = 1
		m = q*cos(lambda)*cos(psi)
		if(psi<lambda && g>m) depth = pi
		else if(psi+lambda>pi && g<m) depth=pi
		else depth=0
			
	}
	else depth=k
	
	depth
}



draw <-function(psi,lambda){
	plot(data[,psi,lambda],type="l",ylim=c(0,pi))
}



data = array(0,dim=c(res_g,res_psi,res_lambda))

headerLambda=array(0,dim=c(res_lambda))
headerPsi=array(0,dim=c(res_psi))
headerG=array(0,dim=c(res_g))



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

for(i_g in 0:(res_g-1)){

	g = (i_g * (max_g*2)/(res_g-1)) -1
	headerG[i_g+1]=g
}


for(i_lambda in 0:(res_lambda-1)){
	lambda = max_lambda*i_lambda/(res_lambda-1)
	
	print(c(round(100*i_lambda/(res_lambda)),"% completed"))
	

	#iterate over psi angles
	for(i_psi in 0:(res_psi-1)){
		psi = (i_psi * (max_psi)/(res_psi-1))
		
		for(i_g in 0:(res_g-1)){
		g = (i_g * (max_g*2)/(res_g-1)) - 1
			d = calculateDepth(g,psi,lambda)
			data[i_g+1, i_psi+1, i_lambda+1] = d
		}
		
		
	}
}
print(c(100,"% completed"))

source("floatconversion.R")
saveTable("depthBuffer.csv",c(res_g,res_psi,res_lambda),headerG, headerPsi, headerLambda, data)






