read.in.haul=function(x, nrows, name="Haul", nextra=0){
	y=read.table(paste0(name,x, ".txt"), header=TRUE, nrows=nrows)
	if(nextra>0)
	{
		for (i in 1:nextra)
		{
			temp=read.table(paste0(name,x, ".txt"),  skip=nrows+i, nrows=1)
			y[,as.character(temp[1])]=temp[2]
		}
	}
	#temp=read.table(paste0(name,x, ".txt"),  skip=nrows+1, nrows=nextra)
	#y2=data.frame(
	y$haul=x
	return(y)
}
