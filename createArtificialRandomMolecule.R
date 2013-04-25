filename="artificialStructure.xyzr"
box=c(500,50,50)
radii=c(2.4,5.4)
wl=10

maxDensity=box[2]*box[3]/(mean(radii)^2*(wl*0.1))

step=maxDensity/box[1]
started=FALSE
cnt=0
for(i in 1:box[1]){
	for(j in 1:floor(step*i)){
		x = c(runif(1,min=(i-box[1]/2)-wl/2, max=(i-box[1]/2)+wl/2), runif(1,min=-box[2]/2, max=+box[2]/2), runif(1,min=-box[3]/2, max=+box[3]/2), runif(1,min=radii[1],max=radii[2]))
		if(started) write(file=filename,x, append=TRUE)
		else write(file=filename,x, append=FALSE)
		started=TRUE
		cnt=cnt+1
		
	}
}

print(cnt)
