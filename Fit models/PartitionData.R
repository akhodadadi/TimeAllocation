#this script gets a vector x and partitioned it to
#numOfPart partitions. It returns a vector with the same size
#as x where each element indicates which partition the 
#corresponding element in x belongs to.


PartitionData = function(x,numOfPart){
  rg = range(x)
  limits=seq(rg[1]-.001,rg[2]+.001,length=(numOfPart+1))
  index=rep(1,length(x))
  
  for (i in 1:(length(limits)-1)){
    index[x>=limits[i] & x<limits[i+1]]=i
  }
  output=index
  
}