calculateVariance <-function(data,m){
	total = 0
	for(x in data){
		total = total + (x-m) * (x-m)
	}
	total / length(data)
}



saveTable <-function(filename, header, tbl, auxiliary){

	
	headerdat = array(0,dim=c(8,length(header)))
	tbldat = array(0,dim=c(8,length(tbl)))
	auxiliarydat = array(0,dim=c(8,length(auxiliary)))
	for(i in 1:length(header)){
		headerdat[,i] = double2CharArray(header[i])
	}
	for(i in 1:length(tbl)){
		tbldat[,i] = double2CharArray(tbl[i])
	}
	for(i in 1:length(auxiliary)){
		auxiliarydat[,i] = double2CharArray(auxiliary[i])
	}
	#derivatives=0, auxiliary=1, dimensions=1, dimensions header=1
	write(file=filename,c(0,1,1,1),append=FALSE)
	write(file=filename,length(header),append=TRUE)
	write(file=filename,headerdat,ncol=8,append=TRUE)
	write(file=filename,tbldat,ncol=8,append=TRUE)
	write(file=filename,auxiliarydat,ncol=8,append=TRUE)
}


bins=32
ex=array(0,bins)
ey=array(0,bins)

en0=array(0,bins)
en1=array(0,bins)
auxiliary0=array(0,bins)
auxiliary1=array(0,bins)
g0=array(0,bins)
g1=array(0,bins)



m=read.table("g.txt")
a=m[which(m[,1]==0),2]
b=m[which(m[,1]==1),2]

a=a[which(a<pi)]
b=b[which(b<pi)]

prior0 = length(a)/(length(a)+length(b))
prior1 = length(b)/(length(a)+length(b))

print(c("p(Occluded)= ",prior0))
print(c("p(Exposed)= ",prior1))

auxiliary0[1]=prior0
auxiliary1[1]=prior1


for(i in 1:bins) ex[i]=max(a)*(i-1)/(bins-1)

step = max(a)/(bins-1)

h=hist(a,breaks=ex,xlim=c(0,3))
t0=0
for(j in 1:(bins-1)) en0[j] = h$density[j]
en0[length(en0)] = en0[length(en0)-1]
for(j in 1:(bins-1)) t0 = t0 +en0[j]
for(j in 1:(bins-1)) en0[j] = en0[j]/(t0)
for(j in 1:(bins-1)) g0[j] = en0[j]/step



h=hist(b,breaks=ex,xlim=c(0,3))
t1=0
for(j in 1:(bins-1)) en1[j] = h$density[j]
en1[length(en1)] = en1[length(en1)-1]
for(j in 1:(bins-1)) t1 = t1 +en1[j]
for(j in 1:(bins-1)) en1[j] = en1[j]/(t1)
for(j in 1:(bins-1)) g1[j] = en1[j]/step


pdf("surfacePointProbabilityDistributions.pdf");
plot(ex,g1,type="l", ann=FALSE)
lines(ex,g0, pch=22, lty=2)
title(main="boundary segment probability distributions");
title(xlab="d (distance in radians)");
title(ylab="p(d|e)");
legend(2.1, 8.2, c("p(d|e=exposed)","p(d|e=buried)"), cex=0.8, lty=1:2);

dev.off()



source("floatconversion.R")

saveTable("occludedDistribution.csv",ex,en0, auxiliary0)
saveTable("exposedDistribution.csv",ex,en1, auxiliary1)





