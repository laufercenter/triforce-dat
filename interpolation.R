#
#	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
#	Email: nils.drechsel@gmail.com
#	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
#


#needs data from the surfaceIntegral.R script
#source("surfaceIntegral.R")


test = array(0,dim=c(3))

dist <- function(a,b){
	sqrt(c(t(a) %*% b))
}

TaylorConvex <- function(i_PHI, i_psi, i_lambda, x){
	#print(c(i_PHI,i_psi,i_lambda,x))
	p = c(headers[1,i_PHI], headers[2,i_psi], headers[3,i_lambda])
	#print(dataConvex[i_PHI, i_psi, i_lambda])
	#print(gradientsConvex[,i_PHI, i_psi, i_lambda])
	#print(hessiansConvex[,,i_PHI, i_psi, i_lambda])
	v = dataConvex[i_PHI, i_psi, i_lambda] + c(t(gradientsConvex[,i_PHI, i_psi, i_lambda]) %*% (x-p)) + 0.5 *(c(t((x-p)) %*% hessiansConvex[,,i_PHI, i_psi, i_lambda] %*% (x-p)))
	v
}

interpolate <- function(x){
	#find n points closest to x
	d_PHI=Inf
	p_PHI=0
	d_psi=Inf
	p_psi=0
	d_lambda=Inf
	p_lambda=0
	

#print("---------------------")

	
	for(i in 1:res_PHI){
		d = x[1]-headers[1,i]
		if(d >= 0 && d <= d_PHI){
			d_PHI=d
			p_PHI=i
		}
	}	
	for(i in 1:res_psi){
		d = x[2]-headers[2,i]
		if(d >= 0 && d <= d_psi){
			d_psi=d
			p_psi=i
		}
	}	
	for(i in 1:res_lambda){
		d = x[3]-headers[3,i]
		if(d >= 0 && d <= d_lambda){
			d_lambda=d
			p_lambda=i
		}
		#print(d)
	}


	#print("-----------------------_")

	closest=c(p_PHI,p_psi,p_lambda)
	#print(c("closest ",closest))
	all=array(0,dim=c(3,8))
	l=1
	for(i in 0:1)
		for(j in 0:1)
			for(k in 0:1){
				all[,l]=closest+c(i,j,k)
				l = l+1
			}

			
	v = array(0,dim=c(8))		
	w = array(0,dim=c(8))		
	dst=array(0,dim=c(3))
	for(i in 1:8){
		a=all[,i]
		if(a[1]>ncol(headers) || a[2]>ncol(headers) || a[3]>ncol(headers)){
			v[i]=0
			w[i]=0
		}
		else{
	
			v[i]=TaylorConvex(a[1],a[2],a[3],x)

			#calculate weights
			p = c(headers[1,a[1]], headers[2,a[2]], headers[3,a[3]])
			#test <<- p
			#print(p)
			for(j in 1:3)
				dst[j] = abs(p[j] - x[j])/abs(headers[j,2]-headers[j,1])

			#print(dst)
			w[i]=1-max(dst)
		}

	}
	#print(v)
	res=0
	s=0
	for(i in 1:8){
		s = s + w[i]
	}

	for(i in 1:8){
		w[i] = w[i] / s
		res = res + w[i] * v[i]
	}
	#res = res / 8

	#print("weights")
	#print(w)

	res

}

p=0.6
p2=0.5
#plb <- function(p,p2){
res=1000
a=array(0,dim=c(res))
b=array(0,dim=c(res))

for(i in 0:(res-1)){
	v = pi * i/(res)
	a[i+1]=interpolate(c(p,v,p2))
	b[i+1]=v
}
plot(b,a,type="l")

#}

#x=c(0.2*pi, 0.2*pi, 0.2*pi)


#s = interpolate(x)

#print(s)
