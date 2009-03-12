#########################################################################
############ Some basic code for outputting the results #################
############ from Illuminus 
############ Taane Clark: last updated 1/8/2007
#########################################################################

# if you are using files with pert scores, you will need to uncomment 
# some of the code below to remove the extra column when assessing the 
# calls

########################
### read in the data ###
########################

example <- as.matrix(read.table("example.txt", header=T)) ## raw data 
calls1 <- as.matrix(read.table("out_calls", header=F))  ## calls 

########################
### ploting function ###
########################

plotthe <- function(i)
{

ncol3<- dim(example)[2]
nrow3 <-dim(calls1)[1]

par(ask = TRUE)

while(i<=nrow3)
{
#print(table(calls1[i,-c(1,2)]))
print(table(calls1[i,-c(1,2,3)]))
xx <- as.numeric(unlist(example[i,seq(3,ncol3-1,by=2)]))
yy <- as.numeric(unlist(example[i,seq(4,ncol3,  by=2)]))
#print(table(xx<0.0)); print(table(yy<0.0))
mx <- max(xx,na.rm=T)
my <- max(yy,na.rm=T)
cc <- (yy-xx)/(xx+yy)
ss <- log(as.numeric(xx)+as.numeric(yy))
# print(table(xx<=0 &  yy<=0))
missing <- as.numeric(xx)<=0 & as.numeric(yy) <= 0
cc[ missing ] <- 0.0
ss[ missing ] <- 0.0
#print(summary(cc))
#print(summary(ss))

############ plot 1: raw scale ###################

#plot(xx[!missing],yy[!missing],main=paste(as.character(example[i,1]),as.character(example[i,2]),
#"miss = ",sum(1*missing)),xlab="A",ylab="B", xlim=c(0,mx),ylim=c(0,my),type="n")
#points(xx[!missing],yy[!missing],pch=20, col = as.numeric(calls1[i,-c(1,2)])[!missing],cex=0.5)
plot(xx[!missing],yy[!missing],main=paste(as.character(example[i,1]),as.character(example[i,2]),
as.character(calls1[i,3]),"miss = ",sum(1*missing)),xlab="A",ylab="B",
xlim=c(0,mx),ylim=c(0,my),type="n")
points(xx[!missing],yy[!missing],pch=20, col = as.numeric(calls1[i,-c(1,2,3)])[!missing],cex=0.5) 

######## plot 2: contrast / strength scale ########

#plot(cc,ss,main=paste(as.character(example[i,1]),as.character(example[i,2])),
#xlab="Contrast",ylab="Strength",type="n")
#points(cc,ss,pch=20, col = as.numeric(calls1[i,-c(1,2)]),cex=0.5)
plot(cc,ss,main=paste(as.character(example[i,1]),as.character(example[i,2]),
as.character(calls1[i,3])),xlab="Contrast",ylab="Strength",type="n") 
points(cc,ss,pch=20, col = as.numeric(calls1[i,-c(1,2,3)]),cex=0.5)

i <- i+1
}
}

par(mfrow=c(1,2))
plotthe(1)


