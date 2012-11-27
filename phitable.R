res_gj = 12
res_gk = 12
res_cosrho = 12

min_gj = 0
max_gj = 1

min_gk = 0
max_gk = 1

min_cosrho = -1
max_cosrho = 1

fd = 0.001

dimensions = 3

data = array(0,dim=c(res_gj, res_gk, res_cosrho))
gradients = array(0,dim=c(dimensions,res_gj, res_gk, res_cosrho))
hessians = array(0,dim=c(dimensions, dimensions, res_gj, res_gk, res_cosrho))

headerGj = array(0,dim=c(res_gj))
headerGk = array(0,dim=c(res_gk))
headerCosrho = array(0,dim=c(res_cosrho))


ex = c(1,0,0)
ey = c(0,1,0)
ez = c(0,0,1)


THRESHOLD_NUMERICAL = 0.0001
MINISCULE = 0.01

isWithinNumericalLimits <- function(x,l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) r = TRUE
	else r = FALSE
	
	r
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

#gjves angle between two normal vectors
normalAngle <-function(a,b){
	acos(c(t(a) %*% b))
}

vlength <-function(v){
	sqrt(c(t(v) %*% v))

}

angle <-function(a,b){
	a1=normalise(a)
	b1=normalise(b)
	normalAngle(a1,b1)
}



isWithinLimits <-function(g_j, g_k, cosrho){
	r = TRUE

	if(g_j<=min_gj || g_j>=max_gj || g_k<=min_gk || g_k>=max_gk || cosrho < min_cosrho || cosrho > max_cosrho)
		r = FALSE
	else{
		lambda_j = acos(g_j)
		lambda_k = acos(g_k)

		rho = acos(cosrho)

		if(rho >= lambda_j+lambda_k || rho+lambda_j < lambda_k || rho+lambda_k < lambda_j){
			r = FALSE
		}
	}

	r
}

calculatePHI <- function(g_j, g_k, cosrho, modegj, modegk, modecosrho, i){
	#print(c("entering",gj,gk,cosrho))

	lambda_j = acos(g_j)
	lambda_k = acos(g_k)

	rho = acos(cosrho)
	if(isWithinNumericalLimits(rho,0)) rho=MINISCULE

	

	if(rho >= lambda_j+lambda_k || rho+lambda_j < lambda_k || rho+lambda_k < lambda_j){
		PHI = NaN
	}
	else{
		#print(c("starting",gj,gk,rho))
		#projection
		r_i = 1;

		mu_j = ey
		mu_k = mvMultiply(rotx(rho), mu_j)

		a_k = sqrt(r_i * r_i - g_k * g_k)
		a_j = sqrt(r_i * r_i - g_j * g_j)


		tau_kj = (g_k - g_j * cos(rho)) / ((sin(rho)^2))
		tau_jk = (g_j - g_k * cos(rho)) / ((sin(rho)^2))




		eta_jk = eta_kj = (mu_k * tau_kj) + (mu_j * tau_jk)



		omega_jk = crossproduct(mu_j, mu_k) / sin(rho)


		gamma_jk = sqrt(r_i * r_i - g_j*tau_jk - g_k*tau_kj)

		x0 = (omega_jk * gamma_jk) + eta_jk

		#print(x0)

		if(isWithinNumericalLimits(x0[1],1) && isWithinNumericalLimits(x0[2],0) && isWithinNumericalLimits(x0[3],0))
			PHI=0
		else{


			#measurement points
			o = ex
			v = ey * g_j
			s = normalise(crossproduct(v,o))
			s = a_j*s
			p0 = v+s
			
			s2 = normalise(-crossproduct(v,s))
			s2 = a_j*s2
			p1 = v+s2

			

			#PHI
			vx = x0-v
			vp0 = p0-v
			vp1 = p1-v

			eta = angle(vp1,vx)
			
			sn = sign(dot(vp0,vx))


			
			if(sn > 0) eta = -eta
			

			PHI = eta
		}
	}


	PHI



}






gradient <- function(integrator, gj, gk, cosrho, modegj, modegk, modecosrho, i){

	switch(modegj, 
		central={pfd_gj=fd; nfd_gj= -fd; cfd_gj=0},
		forward={pfd_gj=2*fd; nfd_gj=0; cfd_gj=fd},
		backward={pfd_gj=0; nfd_gj= -2*fd; cfd_gj= -fd},
		invalid={pfd_gj=0; nfd_gj= 0; cfd_gj= 0},
	)
	switch(modegk, 
		central={pfd_gk=fd; nfd_gk= -fd; cfd_gk=0},
		forward={pfd_gk=2*fd; nfd_gk=0; cfd_gk=fd},
		backward={pfd_gk=0; nfd_gk= -2*fd; cfd_gk= -fd},
		invalid={pfd_gk=0; nfd_gk= 0; cfd_gk= 0},
	)
	switch(modecosrho, 
		central={pfd_cosrho=fd; nfd_cosrho= -fd; cfd_cosrho=0},
		forward={pfd_cosrho=2*fd; nfd_cosrho=0; cfd_cosrho=fd},
		backward={pfd_cosrho=0; nfd_cosrho= -2*fd; cfd_cosrho= -fd},
		invalid={pfd_cosrho=0; nfd_cosrho= 0; cfd_cosrho= 0},
	)


	res = array(0,dim=c(3))
	res[1] = (integrator(gj+pfd_gj, gk, cosrho, modegj, modegk, modecosrho, i) - integrator(gj+nfd_gj, gk, cosrho, modegj, modegk, modecosrho, i))/(2*fd)
	res[2] = (integrator(gj, gk+pfd_gk, cosrho, modegj, modegk, modecosrho, i) - integrator(gj, gk+nfd_gk, cosrho, modegj, modegk, modecosrho, i))/(2*fd)
	res[3] = (integrator(gj, gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) - integrator(gj, gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i))/(2*fd)
	res
	
}




hessian <- function(integrator, gj, gk, cosrho, modegj, modegk, modecosrho, i){
	res = array(0,dim=c(3,3))

	switch(modegj, 
		central={pfd_gj=fd; nfd_gj= -fd; cfd_gj=0},
		forward={pfd_gj=2*fd; nfd_gj=0; cfd_gj=fd},
		backward={pfd_gj=0; nfd_gj= -2*fd; cfd_gj= -fd},
		invalid={pfd_gj=0; nfd_gj= 0; cfd_gj= 0},
	)
	switch(modegk, 
		central={pfd_gk=fd; nfd_gk= -fd; cfd_gk=0},
		forward={pfd_gk=2*fd; nfd_gk=0; cfd_gk=fd},
		backward={pfd_gk=0; nfd_gk= -2*fd; cfd_gk= -fd},
		invalid={pfd_gk=0; nfd_gk= 0; cfd_gk= 0},
	)
	switch(modecosrho, 
		central={pfd_cosrho=fd; nfd_cosrho= -fd; cfd_cosrho=0},
		forward={pfd_cosrho=2*fd; nfd_cosrho=0; cfd_cosrho=fd},
		backward={pfd_cosrho=0; nfd_cosrho= -2*fd; cfd_cosrho= -fd},
		invalid={pfd_cosrho=0; nfd_cosrho= 0; cfd_cosrho= 0},
	)


	#f gj gj
	res[1,1] = (integrator(gj+pfd_gj, gk, cosrho, modegj, modegk, modecosrho, i) - 2*integrator(gj+cfd_gj, gk, cosrho, modegj, modegk, modecosrho, i) + integrator(gj+nfd_gj, gk, cosrho, modegj, modegk, modecosrho, i))/fd^2
	#f gj gk
	res[1,2] = (integrator(gj+pfd_gj, gk+pfd_gk, cosrho, modegj, modegk, modecosrho, i) - integrator(gj+pfd_gj, gk+nfd_gk, cosrho, modegj, modegk, modecosrho, i) - integrator(gj+nfd_gj, gk+pfd_gk, cosrho, modegj, modegk, modecosrho, i) + integrator(gj+nfd_gj, gk+nfd_gk, cosrho, modegj, modegk, modecosrho, i))/(4*fd^2)
	#f gj cosrho
	res[1,3] = (integrator(gj+pfd_gj, gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) - integrator(gj+pfd_gj, gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i) - integrator(gj+nfd_gj, gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) + integrator(gj+nfd_gj, gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i))/(4*fd^2)

	#f gk gj
	res[2,1] = res[1,2]
	#f gk gk
	res[2,2] = (integrator(gj, gk+pfd_gk, cosrho, modegj, modegk, modecosrho, i) - 2*integrator(gj, gk+cfd_gk, cosrho, modegj, modegk, modecosrho, i) + integrator(gj, gk+nfd_gk, cosrho, modegj, modegk, modecosrho, i))/fd^2
	#f gk cosrho
	res[2,3] = (integrator(gj, gk+pfd_gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) - integrator(gj, gk+pfd_gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i) - integrator(gj, gk+nfd_gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) + integrator(gj, gk+nfd_gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i))/(4*fd^2)


	#f cosrho  gj
	res[3,1] = res[1,3]
	#f cosrho gk
	res[3,2] = res[2,3]
	#f gk gk
	res[3,3] = (integrator(gj, gk, cosrho+pfd_cosrho, modegj, modegk, modecosrho, i) - 2*integrator(gj, gk, cosrho+cfd_cosrho, modegj, modegk, modecosrho, i) + integrator(gj, gk, cosrho+nfd_cosrho, modegj, modegk, modecosrho, i))/fd^2

	res
}






saveTable <-function(filename, dimensions, headerGj, headerGk, headerCosrho, tbl, tblgradients, tblhessians){

	headerGjdat = array(0,dim=c(8,length(headerGj)))
	headerGkdat = array(0,dim=c(8,length(headerGk)))
	headerCosrhodat = array(0,dim=c(8,length(headerCosrho)))
	tbldat = array(0,dim=c(8,length(tbl)))
	tblgradientsdat = array(0,dim=c(8,length(tblgradients)))
	tblhessiansdat = array(0,dim=c(8,length(tblhessians)))
	for(i in 1:length(headerGj)){
		headerGjdat[,i] = double2CharArray(headerGj[i])
	}
	for(i in 1:length(headerGk)){
		headerGkdat[,i] = double2CharArray(headerGk[i])
	}
	for(i in 1:length(headerCosrho)){
		headerCosrhodat[,i] = double2CharArray(headerCosrho[i])
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
	write(file=filename,c(1,res_gj,1,res_gk,1,res_cosrho),ncol=3,append=TRUE)
	write(file=filename,headerGjdat,ncol=8,append=TRUE)
	write(file=filename,headerGkdat,ncol=8,append=TRUE)
	write(file=filename,headerCosrhodat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
	write(file=filename,tblgradientsdat,ncol=8,append=TRUE)
	write(file=filename,tblhessiansdat,ncol=8,append=TRUE)
}





for(i_gj in 0:(res_gj-1)){
	print(c(round(100*i_gj/(res_gj)),"% completed"))
	gj = max_gj*i_gj/(res_gj-1)

	headerGj[i_gj+1]=gj

	if(isWithinNumericalLimits(gj,0)) gj=MINISCULE

	#iterate over gk angles
	for(i_gk in 0:(res_gk-1)){

		gk = (i_gk * (max_gk)/(res_gk-1))

		headerGk[i_gk+1]=gk

		if(isWithinNumericalLimits(gk,0)) gk=MINISCULE

		for(i_cosrho in 0:(res_cosrho-1)){
			cosrho = (max_cosrho-min_cosrho) * i_cosrho/(res_cosrho-1) + min_cosrho

			headerCosrho[i_cosrho+1]=cosrho

			#if(isWithinNumericalLimits(cosrho,-1)) cosrho= -1+MINISCULE

			d = calculatePHI(gj,gk,cosrho)
			data[i_gj+1 ,i_gk+1 ,i_cosrho+1] = d

			modeGj="central"
			modeGk="central"
			modeCosrho="central"

			if(!isWithinLimits(gj+2*fd,gk,cosrho) && !isWithinLimits(gj-2*fd,gk,cosrho)) modeGj="invalid"
			else if(modeGj!="backward" && !isWithinLimits(gj+2*fd,gk,cosrho)) modeGj="backward"
			else if(modeGj!="forward" && !isWithinLimits(gj-2*fd,gk,cosrho)) modeGj="forward"

			if(!isWithinLimits(gj,gk+2*fd,cosrho) && !isWithinLimits(gj,gk-2*fd,cosrho)) modeGk="invalid"
			else if(modeGk!="backward" && !isWithinLimits(gj,gk+2*fd,cosrho)) modeGk="backward"
			else if(modeGk!="forward" && !isWithinLimits(gj,gk-2*fd,cosrho)) modeGk="forward"

			if(!isWithinLimits(gj,gk,cosrho+2*fd) && !isWithinLimits(gj,gk,cosrho-2*fd)) modeCosrho="invalid"
			else if(modeCosrho!="backward" && !isWithinLimits(gj,gk,cosrho+2*fd)) modeCosrho="backward"
			else if(modeCosrho!="forward" && !isWithinLimits(gj,gk,cosrho-2*fd)) modeCosrho="forward"



			if(is.nan(d)) g = NaN
			else g = gradient(calculatePHI,gj,gk,cosrho, modeGj, modeGk, modeCosrho)
			gradients[,i_gj+1, i_gk+1, i_cosrho+1] = g

			if(is.nan(d)) h = NaN
			else h = hessian(calculatePHI, gj, gk, cosrho, modeGj, modeGk, modeCosrho)
			hessians[,,i_gj+1, i_gk+1, i_cosrho+1] = h

		}
	}
}
print(c(100,"% completed"))


source("floatconversion.R")
saveTable("dataPHI.csv",c(res_gj,res_gk,res_cosrho),headerGj, headerGk, headerCosrho,data, gradients, hessians)

