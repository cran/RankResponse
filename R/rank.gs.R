#    <rank.gs>
#    Copyright (C) <2014>  <Hsiuying Wang, Yu-Jun Lin>
#
#
#    This program is free software; you can redistribute it and/or modify
#
#    it under the terms of the GNU General Public License as published by
#
#    the Free Software Foundation; either version 2 of the License, or
#
#    (at your option) any later version.
#
# 
#
#    This program is distributed in the hope that it will be useful,
#
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#
#    GNU General Public License for more details.

rank.gs=function(data,alpha,type=2)
{
  data=as.matrix(data)
  data=data[!apply(apply(data,1,is.na),2,any),]
  n=dim(data)[1]
  k=dim(data)[2]

  m=apply(data,2,sum)
  names(m)=c(1:k)
  z=qnorm(1-alpha/2)
  pi=m/n
  pi_temp=sort(pi)

  t=as.numeric(names(pi_temp))
  data=data[,t]

  score=numeric(k) 
  p=k
  score[k]=k

  for(i in 1:(k-1))
  {
    q=k-i
    a=sum(data[,p])/n;b=sum(data[,q])/n

    if(type==1)
    {
       x=abs(a-b)/sqrt((a+b)/n)
    }else{

       c=sum(data[data[,p]==1&data[,q]==1,1])/n
       x=abs(a-b)/sqrt((a+b-2*c)/n)   
    }  

    if(x>=z)
    {
       p=k-i
       score[k-i]=q
    }else{

       p=k-i 
       score[k-i]=p     
    }
  }

  rank_temp=numeric(k)
  rank_temp[k]=1
  for(i in (k-1):1)
  {
      if(score[i+1]>score[i])
      {
         rank_temp[i]=rank_temp[i+1]+1
      }else{

         rank_temp[i]=rank_temp[i+1]
      }
   }

   rank=numeric(k)
   for(i in 1:k)
   {
     rank[t[i]]=rank_temp[i]
   }
   probability=pi
   result=rbind(probability,rank)
   return(result)

}
