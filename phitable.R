#
#	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
#	Email: nils.drechsel@gmail.com
#	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
#


REHASH = TRUE
DERIVATIVES = FALSE

res_gj = 12
res_gk = 12
res_cosrho = 12

min_gj = 0
max_gj = 1

min_gk = 0
max_gk = 1

min_cosrho = -1
max_cosrho = 1

fd = 0.01

parameter_gj_log = 1
parameter_gk_log = 1
parameter_cosrho_log = 50

dimensions = 3
if(REHASH){
data = array(0,dim=c(res_gj, res_gk, res_cosrho))
data_raw = array(0,dim=c(res_gj, res_gk, res_cosrho))
gradients = array(0,dim=c(dimensions,res_gj, res_gk, res_cosrho))
hessians = array(0,dim=c(dimensions, dimensions, res_gj, res_gk, res_cosrho))

headerGj = array(0,dim=c(res_gj))
headerGk = array(0,dim=c(res_gk))
headerCosrho = array(0,dim=c(res_cosrho))
}

ex = c(1,0,0)
ey = c(0,1,0)
ez = c(0,0,1)


THRESHOLD_NUMERICAL = 0.0001
MINISCULE = 0.00001


library("cubature")

isWithinNumericalLimits <- function(x,l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) r = TRUE
	else r = FALSE
	
	r
}


rotationalAngle <-function(nI, nJ){
	ez=c(0,0,1)
	ex=c(1,0,0)
	
	n0 = normalise(crossproduct(ex,nI));
	#n0=c(0,0,1)
	n2 = normalise(crossproduct(nI,n0));
	
	n1 = normalise(crossproduct(nI, nJ))

	sn = sign(dot(n1,n2))

	print(c("sn: ",sn))

	a = asin(dot(n0,n1))

	print(c("a: ",a))

	a = (1-sn)*pi + sn * asin(dot(n0,n1))

	
	a
	
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

	if(g_j<min_gj || g_j>max_gj || g_k<min_gk || g_k>max_gk || cosrho < min_cosrho || cosrho > max_cosrho){
		r = FALSE
	}
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

csc <- function(x){
	1/sin(x)
}

cot <- function(x){
	1/tan(x)
}



calculatePHI2 <- function(g_j, g_k, cosrho){
	lambda_j=acos(g_j)
	lambda_k=acos(g_k)
	rho = acos(cosrho)

	eta = pi/2 - acos(-cot(lambda_j)*cot(rho)+cos(lambda_k)*csc(lambda_j)*csc(rho))

	eta
}


calculatePHI <- function(g_j, g_k, cosrho, modegj, modegk, modecosrho, i,branch=0){
	if(isWithinNumericalLimits(g_j,0)) g_j=g_j+MINISCULE
	if(isWithinNumericalLimits(g_j,1)) g_j=g_j-MINISCULE
	if(isWithinNumericalLimits(g_k,0)) g_k=g_k+MINISCULE
	if(isWithinNumericalLimits(g_k,1)) g_k=g_k-MINISCULE

	lambda_j=acos(g_j)
	lambda_k=acos(g_k)
	
	min_value_cosrho = cos(lambda_j+lambda_k) + MINISCULE

	if(lambda_j > lambda_k) max_value_cosrho = cos(lambda_j-lambda_k) - MINISCULE
	else max_value_cosrho = cos(lambda_k-lambda_j) - MINISCULE

	

	#print(c(min_value_cosrho,max_value_cosrho, lambda_j,lambda_k))


	#gj = log(parameter_gj_log*max_gj*i_gj/(res_gj-1)+1) / log(parameter_gj_log*max_gj*(res_gj-1)/(res_gj-1)+1)


	#if(i_cosrho < res_cosrho/2){
	#	cosrho = (log(parameter_cosrho_log*(max_value_cosrho-min_value_cosrho) * i_cosrho/(res_cosrho-1)+1)+min_value_cosrho) / (log(parameter_cosrho_log*(max_value_cosrho-min_value_cosrho) * (res_cosrho-1)/(res_cosrho-1)+1)+min_value_cosrho)
	#}
	#else{
	#	cosrho = (log(parameter_cosrho_log*(max_value_cosrho-min_value_cosrho) * i_cosrho/(res_cosrho-1)+min_value_cosrho+1)) / (log(parameter_cosrho_log*(max_value_cosrho-min_value_cosrho) * (res_cosrho-1)/(res_cosrho-1)+min_value_cosrho+1))
	#}



	#print(c("entering",gj,gk,cosrho))


	if(cosrho>1){
		#print("perturbation error in cosrho",cosrho)
		cosrho=1
	}
	rho = acos(cosrho)

	if(isWithinNumericalLimits(rho,0)) rho = rho + MINISCULE
	if(isWithinNumericalLimits(rho,1)) rho = rho - MINISCULE


	#print(c("rho",cosrho,rho))

	#if(isWithinNumericalLimits(g_j,g_k) && isWithinNumericalLimits(rho,0)){
	#	PHI = pi/2
	#}
	#else
	{

			#print(c("starting",gj,gk,rho))
			#projection
			r_i = 1;

			mu_j = ey
			mu_k = mvMultiply(rotx(rho), mu_j)

			a_k = sqrt(r_i * r_i - g_k * g_k)
			a_j = sqrt(r_i * r_i - g_j * g_j)


			tau_kj = (g_k - g_j * cos(rho)) / ((sin(rho)^2))
			tau_jk = (g_j - g_k * cos(rho)) / ((sin(rho)^2))


			#print(c("tau:",tau_jk,tau_kj))


			eta_jk = eta_kj = (mu_k * tau_kj) + (mu_j * tau_jk)



			#omega_kj = crossproduct(mu_k, mu_j) / sin(rho)
			#omega_jk = crossproduct(mu_j, mu_k) / sin(rho)

			omega_jk = c(1,0,0)



			#gamma_kj = sqrt(r_i * r_i - g_k*tau_kj - g_j*tau_jk)
			s = r_i * r_i - g_j*tau_jk - g_k*tau_kj
			if(s<0){
				#print(c("perturbation error: ",s))
				s = 0
			}
			gamma_jk = sqrt(s)

			x0 = (omega_jk * gamma_jk) + eta_jk
			#else x0 = (omega_kj * gamma_kj) + eta_kj


			#print(x0)

			#if(isWithinNumericalLimits(x0[1],1) && isWithinNumericalLimits(x0[2],0) && isWithinNumericalLimits(x0[3],0))
			#	PHI=0
			#else{
				#measurement points

					o = ex
					v = ey * g_j
					s = normalise(crossproduct(v,o))
					s = a_j*s
					p0 = v+s
					
					s2 = normalise(-crossproduct(v,s))
					s2 = a_j*s2
					p1 = v+s2
					#print(c("-s-:",s,s2,o,v))

				

				#PHI
					vx = x0-v
					vp0 = p0-v
					vp1 = p1-v

				eta = angle(vp1,vx)

				#print(c("eta",eta))
				
				sn = sign(dot(vp0,vx))

				#print(c("sn",sn))

				
				if(sn > 0) eta = -eta
				#if(isWithinNumericalLimits(g_j,1)) eta = -eta
				

				PHI = eta

				#print(c("PHI",PHI))

			#}
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



interpolate <-function(sp,x){
	p = c(headerGj[sp[1]],headerGk[sp[2]],headerCosrho[sp[3]])
	f = data[sp[1],sp[2],sp[3]]
	g = gradients[,sp[1],sp[2],sp[3]]
	d = c(x-p)
	v =  f + dot(g,d)
	v
}

interpolate2 <-function(sp,f,g,x){
	p = c(headerGj[sp[1]],headerGk[sp[2]],headerCosrho[sp[3]])
	d = c(x-p)
	v =  f + dot(g,d);
	v
}


weights <-function(x,sps){
	lengths=c(abs(headerGj[2]-headerGj[1]),abs(headerGk[2]-headerGk[1]),abs(headerCosrho[2]-headerCosrho[1]))
	maxw = 0
	wghts=array(0,0)
	for(i in 1:nrow(sps)){
		sp = c(sps[i,])
		p_sp = c(headerGj[sp[1]],headerGk[sp[2]],headerCosrho[sp[3]])
		d=c(0,0,0)
		for(j in 1:3){
			d[j]=(abs(p_sp[j]-x[j]) / lengths[j])
			#print(abs(p_sp[j]-x[j]))
		}
		w = 1.0-max(d)
		maxw=maxw+w
		wghts = c(wghts,w)

	}
	wghts = wghts / maxw
	wghts
}

calculateInterpolationArea <- function(x){
	if(isWithinLimits(x[1],x[2],x[3])){
		a = 1
	}
	else a = 0
	a
}

calculateInterpolationError2 <- function(x, v, g, sps, base){
	#print(c(x[1],x[2],x[3]))
	#print("interpolation")
	if(isWithinLimits(x[1],x[2],x[3])){
		realv = calculatePHI(x[1],x[2],x[3])
		#print(realv)
		extsps = rbind(sps,matrix(base,1,3))

		wghts = weights(x,extsps)
		intv = 0
		for(i in 1:nrow(sps)){
			sp = c(sps[i,])
			v = interpolate(sp,x)
			intv = intv + wghts[i] * v
			
		}
		v = interpolate2(base,v, g, x)
		intv = intv + wghts[length(wghts)] * v

		#print(c(realv,intv))
		#print(wghts)
				

		
		e = (realv - intv)^2

		#e = (realv - intv)^2 + (realv - v)^2

		
	}
	else e = 0
	e
}


calculateInterpolationError <- function(x, v, g, sps, base){
	#print(c(x[1],x[2],x[3]))
	#print("interpolation")
	if(isWithinLimits(x[1],x[2],x[3])){
		realv = calculatePHI(x[1],x[2],x[3])
		intv = interpolate2(base,v, g, x)
		
		e = (realv - intv)^2

		
	}
	else e = 0
	e
}





integrateInterpolationError <-function(x, v, base, sps, octants, sps_numb){
	g = x[1:3]

	e = 0



	s = 1
	for(i in 1:length(sps_numb)){
		sps_octant = sps[s:sps_numb[i],]
		octant = c(octants[i,])
		s = s + sps_numb[i]

		lowerLimit = c(headerGj[base[1]], headerGk[base[2]], headerCosrho[base[3]])
		upperLimit = c(headerGj[base[1]+octant[1]], headerGk[base[2]+octant[2]], headerCosrho[base[3]+octant[3]])

		#print(c("lower: ",lowerLimit))
		#print(c("upper: ",upperLimit))

		#r = adaptIntegrate(calculateInterpolationArea, lower=lowerLimit, upper=upperLimit, maxEval=30)
		#a = abs(r$integral)
		#print(c("AREA: ",a))
		a = 1

		if(a > 0){
			r = adaptIntegrate(calculateInterpolationError, v=v, g=g, sps=sps_octant, base=base, lower=lowerLimit, upper=upperLimit, maxEval=30)
			e = e+abs(r$integral)/a
		}
	}
	#print(x)
	#print(e)
	e
}
	

optimalBrainSurgery <- function(base, sps, octants, sps_numb){

	print(sps)
	print(octants)
	print(sps_numb)

	v = getValueThroughLimitingPoint(base)

	if(!is.nan(v)){

		gr = c(0,0,0)
		for(i in 1:nrow(sps)){
			gr = gr + gradients[,sps[i,1],sps[i,2],sps[i,3]]
		}
		gr = gr / nrow(sps)

		start=gr
		print(c("START: ",start))
		print(c("BASE V: ",v))

		par = optim(start,fn=integrateInterpolationError, v = v, base=base, sps=sps, octants=octants, sps_numb=sps_numb)$par

		g=par[1:3]

		r = c(v,g)
	}
	else r = NaN

	r

	

}


getSupportPointsForOctant <- function(base, octant){
	sps = matrix(0,0,3)
	invalid = FALSE
	rawnan = 0
	for(i in 0:octant[1])
		for(j in 0:octant[2])
			for(k in 0:octant[3]){
				cgj = base[1] + i
				cgk = base[2] + j
				ccosrho = base[3] + k
				if(!(cgj >=1 && cgj <=res_gj && cgk >=1 && cgk <=res_gk && ccosrho >=1 && ccosrho <=res_cosrho)){
					invalid = TRUE
				}
				else
				if(!is.nan(data_raw[cgj,cgk,ccosrho])){
					sps = rbind(sps, matrix(c(cgj,cgk,ccosrho),1,3))
				}
			}

	if(invalid) sps = matrix(0,0,3)

	sps
}


getValueThroughLimitingPoint <- function(base){
	in_gj = base[1]
	in_gk = base[2]
	in_cosrho = base[3]

	lengths=c(abs(headerGj[2]-headerGj[1]),abs(headerGk[2]-headerGk[1]),abs(headerCosrho[2]-headerCosrho[1]))
	
	lambda_j = acos(headerGj[in_gj])
	lambda_k = acos(headerGk[in_gk])
	rho = acos(headerCosrho[in_cosrho])

	v = 0
	cv = 0
	w_sum = 0

	if(in_gj+1 <= res_gj && !is.nan(data_raw[in_gj+1,in_gk,in_cosrho])){
		limit = cos(rho+lambda_k)
		f = ((pi/2 - data_raw[in_gj+1,in_gk,in_cosrho]) / (limit - headerGj[in_gj+1])) * lengths[1] + data_raw[in_gj+1,in_gk,in_cosrho]
		w = abs(limit - headerGj[in_gj+1])
		v = v + f * w
		w_sum = w_sum + w
		cv = cv + 1
		print(c(f,w,w_sum))
	}

	if(in_gj-1 >= 1 && !is.nan(data_raw[in_gj-1,in_gk,in_cosrho]) && !(acos(gk)>acos(gj)+acos(cosrho)) && !(acos(gj)>acos(gk)+acos(cosrho))){
		limit = cos(rho-lambda_k)
		f = ((-pi/2 - data_raw[in_gj-1,in_gk,in_cosrho]) / (limit - headerGj[in_gj-1])) * lengths[1] + data_raw[in_gj-1,in_gk,in_cosrho]
		w = abs(limit - headerGj[in_gj-1])
		v = v + f * w
		w_sum = w_sum + w
		cv = cv + 1
		print(c(f,w,w_sum))
	}

	
	if(in_gk+1 <= res_gk && !is.nan(data_raw[in_gj,in_gk+1,in_cosrho])){
		limit = cos(rho+lambda_j)
		f = ((pi/2 - data_raw[in_gj,in_gk+1,in_cosrho]) / (limit - headerGk[in_gk+1])) * lengths[2] + data_raw[in_gj,in_gk+1,in_cosrho]
		w = abs(limit - headerGk[in_gk+1])
		v = v + f * w
		w_sum = w_sum + w
		cv = cv + 1
		print(c(f,w,w_sum))
	}

	if(in_gk-1 >= 1 && !is.nan(data_raw[in_gj,in_gk-1,in_cosrho]) && !(acos(gk)>acos(gj)+acos(cosrho)) && !(acos(gj)>acos(gk)+acos(cosrho))){
		limit = cos(rho-lambda_j)
		f = ((pi/2 - data_raw[in_gj,in_gk-1,in_cosrho]) / (limit - headerGk[in_gk-1])) * lengths[2] + data_raw[in_gj,in_gk-1,in_cosrho]
		w = abs(limit - headerGk[in_gk-1])
		v = v + f * w
		w_sum = w_sum + w
		cv = cv + 1
		print(c(f,w,w_sum))
	}

	if(in_cosrho+1 <= res_cosrho && !is.nan(data_raw[in_gj,in_gk,in_cosrho+1])){
		limit = cos(lambda_j+lambda_k)
		f = ((pi/2 - data_raw[in_gj,in_gk,in_cosrho+1]) / (limit - headerCosrho[in_cosrho+1])) * lengths[3] + data_raw[in_gj,in_gk,in_cosrho+1]
		w = abs(limit - headerCosrho[in_cosrho+1])
		v = v + f * w
		w_sum = w_sum + w
		cv = cv + 1
		print(c(f,w,w_sum))
	}

	if(in_cosrho-1 >= 1 && !is.nan(data_raw[in_gj,in_gk,in_cosrho-1])){
		if(lambda_j > lambda_k){
			limit = cos(lambda_j-lambda_k)
			f = ((pi/2 - data_raw[in_gj,in_gk,in_cosrho-1]) / (limit - headerCosrho[in_cosrho-1])) * lengths[3] + data_raw[in_gj,in_gk,in_cosrho-1]
			w = abs(limit - headerCosrho[in_cosrho-1])
			v = v + f * w
			w_sum = w_sum + w
			cv = cv + 1
			print(c(f,w,w_sum))
		}
		else{
			limit = cos(lambda_k-lambda_k)
			f = ((-pi/2 - data_raw[in_gj,in_gk,in_cosrho-1]) / (limit - headerCosrho[in_cosrho-1])) * lengths[3] + data_raw[in_gj,in_gk,in_cosrho-1]
			w = abs(limit - headerCosrho[in_cosrho-1])
			v = v + f * w
			w_sum = w_sum + w
			cv = cv + 1
			print(c(f,w,w_sum))
		}
	}


	if(cv >= 1){
		v = v / w_sum
	}
	else v = NaN
	v
}


insertLimitInterpolationPoints <-function(){

	gradients_fixed=gradients
	for(i_gj in 0:(res_gj-1)){
		in_gj = i_gj+1
		gj = max_gj*i_gj/(res_gj-1)

		for(i_gk in 0:(res_gk-1)){
			gk = (i_gk * (max_gk)/(res_gk-1))
			in_gk = i_gk+1

			for(i_cosrho in 0:(res_cosrho-1)){
				cosrho = (max_cosrho-min_cosrho) * i_cosrho/(res_cosrho-1) + min_cosrho
				in_cosrho = i_cosrho+1

				if(is.nan(data[in_gj, in_gk, in_cosrho])){
					print(c("FOUND NAN:",in_gj, in_gk, in_cosrho))

					base = c(in_gj, in_gk, in_cosrho)

					sps = matrix(0,0,3)
					octants = matrix(0,0,3)
					sps_numb = c()

					for(i in c(-1,1))
						for(j in c(-1,1))
							for(k in c(-1,1)){
								octant = c(i,j,k)
								tsps = getSupportPointsForOctant(base, octant)
								if(nrow(tsps)>0){
									sps_numb = c(sps_numb,nrow(tsps))
									octants = rbind(octants,matrix(octant,1,3))
									sps = rbind(sps,tsps)
								}
								
							}

					if(nrow(sps)>0){
						par = optimalBrainSurgery(base, sps, octants, sps_numb)
						print(par)
						if(!is.nan(par[1])){
							if(max(par)==1 || sqrt(dot(par[2:4],par[2:4]))>30){
								print("OPTIMISATION PROBLEM")
								#stop(2)
							}
							data[in_gj, in_gk, in_cosrho] <<- par[1]
							gradients[,in_gj, in_gk, in_cosrho] <<- par[2:4]
						}
					}




					
				}
			}
		}
	}

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




if(REHASH){
for(i_gj in 0:(res_gj-1)){
#i_gj=11
#{
	print(c(round(100*i_gj/(res_gj)),"% completed"))

	gj = max_gj*i_gj/(res_gj-1)
	#gj = log(parameter_gj_log*max_gj*i_gj/(res_gj-1)+1) / log(parameter_gj_log*max_gj*(res_gj-1)/(res_gj-1)+1)

	headerGj[i_gj+1]=gj




	#iterate over gk intrusions
	for(i_gk in 0:(res_gk-1)){
	#i_gk=3
	#{

		gk = max_gk*i_gk/(res_gk-1)
		#gk = log(parameter_gk_log*max_gk*i_gk/(res_gk-1)+1) / log(parameter_gk_log*max_gk*(res_gj-1)/(res_gk-1)+1)

		headerGk[i_gk+1]=gk



		for(i_cosrho in 0:(res_cosrho-1)){
		#i_cosrho=11
		#{

			#print(c("I:",i_gj,i_gk,i_cosrho))

			cosrho = max_cosrho-min_cosrho * i_cosrho/(res_cosrho-1) + min_cosrho



			headerCosrho[i_cosrho+1]=cosrho

			#if(isWithinNumericalLimits(cosrho,-1)) cosrho= -1+MINISCULE

			#d = calculatePHI(gj,gk,cosrho)
			d = calculatePHI(headerGj[i_gj+1],headerGk[i_gk+1],headerCosrho[i_cosrho+1])
			#print(c(i_gj,i_gk,i_cosrho,d))

			if(is.nan(d)){
				print(c(i_gj,i_gk,i_cosrho,gj,gk,cosrho))
				stop("ABORT")
			}
			data[i_gj+1 ,i_gk+1 ,i_cosrho+1] = d
			data_raw[i_gj+1 ,i_gk+1 ,i_cosrho+1] = d


			if(DERIVATIVES){
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
}
print(c(100,"% completed"))
}

#insertLimitInterpolationPoints()


#source("floatconversion.R")
#saveTable("dataPHI.csv",c(res_gj,res_gk,res_cosrho),headerGj, headerGk, headerCosrho,data, gradients, hessians)

