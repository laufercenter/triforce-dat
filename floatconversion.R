#
#	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
#	Email: nils.drechsel@gmail.com
#	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
#


library("accuracy")
library("bitops")

double2CharArray <-function(x){
	if(is.nan(x)) res=c(255,255,255,255,255,255,255,255)
	else{
		a=frexp(x)
		significand = a[1,"Mantissa"]
		exponent = a[1,"Exponent"]
		
		significandInt32 = double2FixedSignedInt32(significand, 30)

		data0 = fixedSignedInt322CharArray(significandInt32)
		
		exponentInt32=int2SignedInt32(exponent)
		data1 = fixedSignedInt322CharArray(exponentInt32)
		
		res = c(data0,data1)
	}
	res
}
	
	
fixedSignedInt322CharArray <-function(x){
	INT32BYTEMASK=255

	if(abs(x)>2^31-1) print("ERROR, number is larger than representable range")

	data=c(0,0,0,0)
	
	data[1]=bitAnd(bitShiftR(x,24), INT32BYTEMASK)
	data[2]=bitAnd(bitShiftR(x,16), INT32BYTEMASK)
	data[3]=bitAnd(bitShiftR(x,8), INT32BYTEMASK)
	data[4]=bitAnd(x,INT32BYTEMASK)
	data
	
}
	
double2FixedSignedInt32 <-function(x, fraction){
	factor= bitShiftL(1,fraction)
	d = floor(x*factor)
	d
}

int2SignedInt32 <-function(x){
	x
}


