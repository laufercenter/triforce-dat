#
#	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
#	Email: nils.drechsel@gmail.com
#	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
#


library(MASS)
library(combinat)
source("floatconversion.R")


DERIVATIVE_LEVEL=0
STORE_PHI=1

#the resolution of our grid
grid=80
res_lambda = grid
res_psi = grid
res_PHI = grid

mode=0

dim=3
fd=0.00025


ex = matrix(c(1,0,0),3,1)
ey = matrix(c(0,1,0),3,1)
ez = matrix(c(0,0,1),3,1)



saveTable <-function(filename, dimensions, headerPHI, headerPsi, headerLambda, tbl,dataDimensions){

	
	headerPHIdat = array(0,dim=c(8,length(headerPHI)))
	headerPsidat = array(0,dim=c(8,length(headerPsi)))
	headerLambdadat = array(0,dim=c(8,length(headerLambda)))
	tbldat = array(0,dim=c(8,length(tbl)))
	for(i in 1:length(headerPHI)){
		headerPHIdat[,i] = double2CharArray(headerPHI[i])
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

	write(file=filename,0, append=FALSE)
	write(file=filename,0, append=TRUE)
	write(file=filename,3,append=TRUE)
	write(file=filename,c(1,res_PHI,1,res_psi,1,res_lambda),ncol=3,append=TRUE)
	write(file=filename,dataDimensions,append=TRUE)
	write(file=filename,headerPHIdat,ncol=8,append=TRUE)
	write(file=filename,headerPsidat,ncol=8,append=TRUE)
	write(file=filename,headerLambdadat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
}

dot <- function(a,b){
	c(t(a) %*% b)
}

vlen<-function(a){
	sqrt(dot(a,a))
}


supportPoints<-function(x){
	p0=c(0,0,0)
	x0 = x
	sp=array(0,dim=c(dim,dim))
	for(i in 1:dim){
		if(x[i]>0.5){
			p0[i]=1
			x0[i]=1-x0[i]
		}
	}
	sp[1,]=p0
	o=order(x0,decreasing=T)
	print(o)
	for(i in 2:dim){
		p=p0
		if(p[o[i-1]]==0) p[o[i-1]]=1
		else p[o[i-1]]=0
		sp[i,]=p
	}
	sp
}

bary <- function(x){
	
	d=array(0,dim=dim)
	for(l in 1:dim){
		d[l]=vlen(completeCrds[l,]-x)
	}
	o=order(d)
	
	crds[1,]=array(0.5,dim)
	for(l in 2:dim){
		crds[l,]=completeCrds[o[l-1],]
	}
	#print(crds)
	
	T=t(crds)
	Tinv = ginv(T)
	
	
	b = Tinv %*% x
	1-sum(b)
}


interp <- function(d){
	h0=c()
	h1=c()
	h2=c()
	for(i in 0:d){
		x=pmin(pmax(c(i/d+0.2, i/d-0.2, 0.3),0),1)
		h0=c(h0,bary2(x)[1])
		h1=c(h1,bary2(x)[2])
		h2=c(h2,bary2(x)[3])
	}
	plot(h0,col="Red",ylim=c(0,1))
	points(h1,col="Blue")
	points(h2,col="Green")
}

polyt <- function(x){
	d=array(0,dim=dim)
	for(l in 1:nrow(completeCrds)){
		d[l]=vlen(completeCrds[l,]-x)
	}
	o=order(d)
	dt=array(0,dim=dim)
	df=array(0,dim=dim)
	for(l in 1:dim){
		dt[l]=d[o[l]]
		df[l]=d[o[2*dim-l+1]]
	}
	w=pmax(df-dt,0)
	w
	
}


bary2 <- function(x){
	
	d=array(0,dim=dim)
	for(l in 1:nrow(completeCrds)){
		d[l]=vlen(completeCrds[l,]-x)
	}
	o=order(d)
	#o2=supportPoints(x)
	
	#crds[1,]=array(0,dim)
	for(l in 1:dim){
		crds[l,]=completeCrds[o[l],]
	}
	clast=as.vector(array(0.5,dim))
	
	
	#print("support:")
	#print(x)
	#print("original:")
	#print(crds)
	#print("new:")
	#print(o2)
	
	
	T=t(crds)
	T=T-clast
	Tinv = ginv(T)
	xlast=x-clast
	
	b = Tinv %*% xlast
	b
	pmin(pmax(b,0),1)^2
}

bary3 <- function(x){
	
	d=array(0,dim=dim)
	for(l in 1:dim){
		d[l]=vlen(completeCrds[l,]-x)
	}
	o=order(d)
	
	if(dim==3)
		h=which(HP[,1]==o[1] & HP[,2]==o[2] & HP[,3]==o[3])
	if(dim==6)
		h=which(HP[,1]==o[1] & HP[,2]==o[2] & HP[,3]==o[3] & HP[,4]==o[4] & HP[,5]==o[5] & HP[,6]==o[6])
	
	#b = Tinv[1,] %*% xlast
	#1-sum(b)
	#b[1]
	h[1]
}


buildHyperplanes<-function(){
	h=array(0,dim=c(0,dim+dim^2))
	a=c(1:dim)
	b=permn(a)
	for(i in 1:length(b)){
		o=b[[i]]
		
		crds[1,]=array(0,dim)
		for(l in 2:dim){
			crds[l,]=completeCrds[o[l-1],]
		}
		clast=as.vector(array(0.5,dim))
		
		T=t(crds)
		T=T-clast
		Tinv = ginv(T)
		for(j in 1:dim){
			g=c(o[1:dim], Tinv[j,])
			h[j,,]=rbind(h[j,,],g)
		}
		
	}
	h
}



#create corner nodes
zero=array(0,dim=dim)
crds=array(0,dim=c(dim,dim))
completeCrds=array(0,dim=c(dim,dim))
for(i in 1:(dim)){
	z=zero
	z[i]=1
	completeCrds[i,]=z
}

completeCrds=array(0,dim=c(0,dim))

for(i in 0:1){
	for(j in 0:1){
		for(k in 0:1){
			p=c(i,j,k)
			#if(!(i==0 && j==0 && k==0)){
				completeCrds=rbind(completeCrds,p)
			#}
		}
	}
}

print(completeCrds)




#T=array(0,dim=c(dim,dim))
#for(i in 1:dim){
#	for(j in 1:dim){
#		T[i,j] = crds[j,i]
#	}
#}



halfstep=(1/grid)/2

if(mode==0){
	griddim=array(grid,dim+1)
	griddim[1]=dim
	data=array(0,dim=griddim)
	header0=array((0:(grid-1))/grid+halfstep,dim=grid)
	header1=array((0:(grid-1))/grid+halfstep,dim=grid)
	header2=array((0:(grid-1))/grid+halfstep,dim=grid)
}else{
	griddim=array(grid,dim)
	HP=buildHyperplanes()
}



if(dim==3){
	for(i in 0:(grid-1)){
		print(c(round(100*i/(grid-1)),"% completed"))
		for(j in 0:(grid-1)){
			for(k in 0:(grid-1)){
				x=(c(i,j,k)/(grid))+halfstep
				
				if(mode==0){
					g=bary2(x)
					#g=pmax(0,g)
					data[,i+1,j+1,k+1]=g
					
				}
				else{
					g=bary3(x)
					data[,i+1,j+1,k+1]=g
				}
			}
		}
	}
}

if(dim==6){
	for(i in 0:(grid-1)){
		for(j in 0:(grid-1)){
			print(c(round(100*((j/(grid-1))/grid + i/(grid))),"% completed"))
			for(k in 0:(grid-1)){
				for(l in 0:(grid-1)){
					for(m in 0:(grid-1)){
						for(n in 0:(grid-1)){
							x=(c(i,j,k,l,m,n)/(grid-1))+halfstep
							
							if(mode==0){
								g=bary2(x)
								data[,i+1,j+1,k+1,l+1,m+1,n+1]=g
							}
						}
					}
				}
			}
		}
	}
}


if(T && mode==0){
	if(dim==3){
		saveTable("barycentricWeights.csv",c(grid,grid,grid),header0,header1,header2,data,dim)
	}
	if(dim==6){
		saveTable6D("barycentricWeights.csv",c(grid,grid,grid,grid,grid,grid),header0,header1,header2,header3,header4,header5,data,dim)
	}

}
