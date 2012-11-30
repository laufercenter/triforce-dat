res_gj = 12
res_gk = 12
res_cosrho = 12

min_gj = 0
max_gj = 1

min_gk = 0
max_gk = 1

min_cosrho = -1
max_cosrho = 1

fd = 0.0001

dimensions = 3
if(FALSE){
data = array(0,dim=c(res_gj, res_gk, res_cosrho))
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
MINISCULE = 0.01


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

	if(sn<0 && a>0) a = pi-a
	if(sn<0 && a<0) a = -pi-a

	
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

calculatePHI <- function(g_j, g_k, cosrho, modegj, modegk, modecosrho, i,branch=0){
	#print(c("entering",gj,gk,cosrho))

	lambda_j = acos(g_j)
	lambda_k = acos(g_k)

	rho = acos(cosrho)
	if(isWithinNumericalLimits(rho,0)) rho=MINISCULE

	

	if(rho >= lambda_j+lambda_k || rho+lambda_j < lambda_k || rho+lambda_k < lambda_j){
		PHI = NaN
	}
	else{
		if(g_j == g_k && (cosrho == -1 || cosrho == 1)){
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



			omega_kj = crossproduct(mu_k, mu_j) / sin(rho)
			omega_jk = crossproduct(mu_j, mu_k) / sin(rho)


			gamma_kj = sqrt(r_i * r_i - g_k*tau_kj - g_j*tau_jk)
			gamma_jk = sqrt(r_i * r_i - g_j*tau_jk - g_k*tau_kj)

			if(branch==0) x0 = (omega_jk * gamma_jk) + eta_jk
			else x0 = (omega_kj * gamma_kj) + eta_kj

			#print(x0)

			if(isWithinNumericalLimits(x0[1],1) && isWithinNumericalLimits(x0[2],0) && isWithinNumericalLimits(x0[3],0))
				PHI=0
			else{
				#measurement points

				if(isWithinNumericalLimits(g_j,0)){
					o = ex
					v = ey
					s = normalise(crossproduct(v,o))
					s = s
					p0 = s
					
					s2 = normalise(-crossproduct(v,s))
					s2 = s2
					p1 = s2
				}
				else{
					o = ex
					v = ey * g_j
					s = normalise(crossproduct(v,o))
					s = a_j*s
					p0 = v+s
					
					s2 = normalise(-crossproduct(v,s))
					s2 = a_j*s2
					p1 = v+s2
				}

				

				#PHI
				if(isWithinNumericalLimits(g_j,0)){
					vx = x0
					vp0 = p0
					vp1 = p1
				}
				else{
					vx = x0-v
					vp0 = p0-v
					vp1 = p1-v
				}

				eta = angle(vp1,vx)
				
				sn = sign(dot(vp0,vx))


				
				if(sn > 0) eta = -eta
				

				PHI = eta
			}
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



interpolate <-function(sp,x){
	p = c(headerGj[sp[1]],headerGk[sp[2]],headerCosrho[sp[3]])
	f = data[sp[1],sp[2],sp[3]]
	g = gradients[,sp[1],sp[2],sp[3]]
	d = x-p
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
		d=c(0,0,0)
		for(j in 1:3)
			d[j]=(abs(sp[j]-x[j]) / lengths[j])
		w = 1.0-max(d)
		maxw=maxw+w
		wghts = c(wghts,w)

	}
	wghts = wghts / maxw
	wghts
}


calculateInterpolationError <- function(xd, dimensions, v, g, sps, point, ul){
	p=point


	limit = c(headerGj[p[1]],headerGk[p[2]],headerCosrho[p[3]])
	x = c(headerGj[p[1]],headerGk[p[2]],headerCosrho[p[3]])
	j=1
	for(i in 1:3){
		if(dimensions[i] == 1){
			x[i] = xd[j]
			limit[i] = ul[j]
			j = j +1
		}
	}


	#print(c(x[1],x[2],x[3]))
	if(isWithinLimits(x[1],x[2],x[3])){
		realv = calculatePHI(x[1],x[2],x[3])
		extsps = rbind(sps,matrix(p,1,3))

		wghts = weights(x,extsps)
		intv = 0
		for(i in nrow(sps)){
			sp = c(sps[i,])
			v = interpolate(sp,x)
			intv = intv + wghts[i] * v
			
		}
		v = interpolate2(p,v, g, x)
		intv = intv + wghts[length(wghts)] * v
				

		
		e = (realv - intv)^2

		#sw = (1/dot(x-limit,x-limit)) ^2

		#e = e/sw



		
	}
	else e = 0
	#print(e)
	e
}

integrateInterpolationError <-function(x,dimensions, sps, point, lowerLimit, upperLimit){
	v = x[1]

	g = c(0,0,0)
	j=2
	for(i in 1:3){
		if(dimensions[i] == 1){
			g[i] = x[j]
			j = j +1
		}
	}

	e = 0

	l = array(dim=0)
	u = array(dim=0)

	for(i in 1:3){
		if(dimensions[i] == 1){
			l = c(l,lowerLimit[i])
			u = c(u,upperLimit[i])
		}
	}


	#print("INTEGRATION")

	r = adaptIntegrate(calculateInterpolationError, dimensions=dimensions, v=v, g=g, sps=sps, point=point, ul= u, lower = l, upper = u, maxEval=30)
	#print(r)
	e= r$integral
	
	e
}
	

optimalBrainSurgery <- function(sps,point,lowerLimit,upperLimit){

	start=c(1)
	dimensions=c(0,0,0)
	for(i in 1:3){
		if(!is.nan(upperLimit[i])){
			dimensions[i]=1
			start = c(start,1)
		}
	}

	print("parameters")
	print(sps)
	print(point)
	print(lowerLimit)
	print(upperLimit)
	print(start)


	par = optim(start,fn=integrateInterpolationError, dimensions=dimensions, sps=sps, point=point, lowerLimit=lowerLimit, upperLimit=upperLimit)$par

	v=par[1]
	g=c(0,0,0)
	j = 2
	for(i in 1:3){
		if(dimensions[i]==1){
			g[i] = par[j]
			j=j+1
		}
	}


	c(v,g)

}




insertLimitInterpolationPoints <-function(){
	data_fixed=data
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

				if(in_gj==1 && in_gk==9 && in_cosrho==11)
				if(is.nan(data[in_gj, in_gk, in_cosrho]) && !in_gj==in_gk && !(in_cosrho==1 || in_cosrho==res_cosrho)){
					print(c("FOUND NAN:",in_gj, in_gk, in_cosrho))

					upperLimit = c(NaN,NaN,NaN)
					lowerLimit = c(NaN,NaN,NaN)
					sps = array(dim=c(0,3))
					point = matrix(c(in_gj, in_gk, in_cosrho),1,3)
					case0=FALSE
					case1=FALSE
					if(in_gj+1 <= res_gj && !is.nan(data[in_gj+1,in_gk,in_cosrho])){
						sp = c(in_gj+1, in_gk, in_cosrho)
						sps = rbind(sps, matrix(sp,1,3))
						rho = acos(cosrho)
						if(isWithinNumericalLimits(rho,0)) rho=MINISCULE
						lambda_k = acos(gk)
						upperLimit[1] = cos(rho+lambda_k)
						lowerLimit[1]= headerGj[in_gj+1]
						case0 = TRUE
						
						
						print("CASE 0")
					}

					if(in_gj-1 >= 1 && !is.nan(data[in_gj-1,in_gk,in_cosrho]) && !(acos(gk)>acos(gj)+acos(cosrho)) && !(acos(gj)>acos(gk)+acos(cosrho))){
						sp = c(in_gj-1, in_gk, in_cosrho)
						sps = rbind(sps, matrix(sp,1,3))
						rho = acos(cosrho)
						if(isWithinNumericalLimits(rho,0)) rho=MINISCULE
						lambda_k = acos(gk)
						upperLimit[1] = cos(rho-lambda_k)
						lowerLimit[1]= headerGj[in_gj-1]
						case0 = TRUE
						print("CASE 1")
					}

					
					if(in_gk+1 <= res_gk && !is.nan(data[in_gj,in_gk+1,in_cosrho])){
						sp = c(in_gj, in_gk+1, in_cosrho)
						sps = rbind(sps, matrix(sp,1,3))
						rho = acos(cosrho)
						if(isWithinNumericalLimits(rho,0)) rho=MINISCULE
						lambda_j = acos(gj)
						upperLimit[2] = cos(rho+lambda_j)
						lowerLimit[2]= headerGk[in_gk+1]
						case1 = TRUE
						print("CASE 2")
					}

					if(in_gk-1 >= 1 && !is.nan(data[in_gj,in_gk-1,in_cosrho]) && !(acos(gk)>acos(gj)+acos(cosrho)) && !(acos(gj)>acos(gk)+acos(cosrho))){
						sp = c(in_gj, in_gk-1, in_cosrho)
						sps = rbind(sps, matrix(sp,1,3))
						rho = acos(cosrho)
						if(isWithinNumericalLimits(rho,0)) rho=MINISCULE
						lambda_j = acos(gj)
						upperLimit[2] = cos(rho-lambda_j)
						lowerLimit[2]= headerGk[in_gk-1]
						case0 = TRUE
						print("CASE 3")
					}

					if(in_cosrho+1 <= res_cosrho && !is.nan(data[in_gj,in_gk,in_cosrho+1])){
						sp = c(in_gj, in_gk, in_cosrho+1)
						sps = rbind(sps, matrix(sp,1,3))
						lambda_j = acos(gj)
						lambda_k = acos(gk)
						upperLimit[3] = cos(lambda_j+lambda_k)
						lowerLimit[3] = headerCosrho[in_cosrho+1]
						case0 = TRUE
						print("CASE 4")
					}

					if(in_cosrho-1 >= 1 && !is.nan(data[in_gj,in_gk,in_cosrho-1])){
						sp = c(in_gj, in_gk, in_cosrho-1)
						sps = rbind(sps, matrix(sp,1,3))
						lambda_j = acos(gj)
						lambda_k = acos(gk)
						if(lambda_j > lambda_k){
							upperLimit[3] = cos(lambda_j-lambda_k)
							lowerLimit[3] = headerCosrho[in_cosrho-1]
							case0 = TRUE
						print("CASE 5")
						}
						else{
							upperLimit[3] = cos(lambda_k-lambda_k)
							lowerLimit[3] = headerCosrho[in_cosrho-1]
							case1 = TRUE
						print("CASE 6")
						}
					}

					par = optimalBrainSurgery(sps,point, lowerLimit, upperLimit)
					print(par)
					if(max(par)==1){
						print("JAJA")
						stop(2)
					}
					data_fixed[in_gj, in_gk, in_cosrho] = par[1]
					gradients[,in_gj, in_gk, in_cosrho] = par[2:4]


					if(case0 && case1){
						print(c("ABORT", in_gj, in_gk, in_cosrho))
						stop("ABORT")
					}


					
				}
			}
		}
	}

	data <<- data_fixed

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




if(FALSE){
for(i_gj in 0:(res_gj-1)){
	print(c(round(100*i_gj/(res_gj)),"% completed"))
	gj = max_gj*i_gj/(res_gj-1)

	headerGj[i_gj+1]=gj

	#if(isWithinNumericalLimits(gj,0)) gj=MINISCULE

	#iterate over gk angles
	for(i_gk in 0:(res_gk-1)){

		gk = (i_gk * (max_gk)/(res_gk-1))

		headerGk[i_gk+1]=gk

		#if(isWithinNumericalLimits(gk,0)) gk=MINISCULE

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
}

insertLimitInterpolationPoints()



#source("floatconversion.R")
#saveTable("dataPHI.csv",c(res_gj,res_gk,res_cosrho),headerGj, headerGk, headerCosrho,data, gradients, hessians)

