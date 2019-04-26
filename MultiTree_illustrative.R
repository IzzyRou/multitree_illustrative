library(parallel)
library(lubridate)
#illustrating methods in Routledge, Isobel, et al. "Estimating spatiotemporally varying malaria reproduction numbers in a near elimination setting." Nature communications 9.1 (2018): 2476.


##################
#load fake line list
###################

#you will need to change the directory here to wherever dummy data is stored
d <- read.csv("[insertdirectoryname]/dummy_data.csv",stringsAsFactors=FALSE)
head(d)

#######################
#formatting data
######################

d<-d[complete.cases(d),]

d$Date.symptoms<-decimal_date(as.Date(d$Date, format="%d/%m/%Y"))
d$Date.symptoms[duplicated(d$Date.symptoms)]=d$Date.symptoms[duplicated(d$Date.symptoms)]+0.1

library(data.table)
d$Month<-as.numeric(t(as.data.table(strsplit(d$Date,'/')))[,2])

d<-d[order(d$Date.symptoms),]
#d$Date.symptoms[duplicated(d$Date.symptoms)]=d$Date.symptoms[duplicated(d$Date.symptoms)]+1e-1
d<-d[order(d$Date.symptoms),]
head(d)

d$I=0
d$I[d$Imported=='Y']=1


NI<-which(d$I==0)
d$Date.symptoms.year<-d$Date.symptoms
d$Date.symptoms=d$Date.symptoms-d$Date.symptoms[1]
d$Date.symptoms<-d$Date.symptoms*365

data<-list(n=nrow(d),q=length(NI),NI=as.numeric(NI),t=d$Date.symptoms)

#######################
# define wc function
######################

wc<-function(t,shift,a,epsilon) {
  t=t-shift
  likelihood = a*t*exp(-0.5*a*t*t)/epsilon
  likelihood[likelihood<0]=0
  return(likelihood)
}

n_prior=500

# choose priors
a=runif(n_prior,0.001,0.01)
e=1e-10
shift<-runif(n_prior,10,15)
ftol<-function(f,f2) (f - f2)/f


A_list<-list()

###########################
##run multitree algorithm
##########################


registerDoParallel(40)
Nedges<-5000
A_list<-foreach(xx = 1:n_prior) %dopar% {
  A<-matrix(0,nrow=length(data$NI),ncol=data$n)
  delta_cji<-0
  R<-matrix(1,nrow=length(data$NI),ncol=data$n)
  score<-rep(NA,Nedges)
  score<-rep(NA,Nedges)
  tol<-1
  for(K in 1:Nedges){
    delta_ji<-matrix(0,nrow=length(data$NI),ncol=data$n)
    k=1
    for(i in NI){
      scores<-wc(data$t[i]-data$t[1:(i-1)],shift[xx],a[xx],e)
      scores<-scores*R[k,1:(i-1)]			
      for(j in 1:(i-1)){
        delta_cji = delta_cji + sum(scores[-j],na.rm=TRUE)
        delta_ji[k,j] = log(delta_cji + scores[j]) - log(delta_cji + 1)
      }
      k=k+1
    }
    greedy_index<-which(delta_ji == max(delta_ji,na.rm=TRUE), arr.ind = TRUE)
    A[greedy_index]<-delta_ji[greedy_index]
    R[greedy_index]<-NA
    score[K]=delta_cji
    if(K!=1) tol<-ftol(score[K],score[K-1]); print(tol)
    if(tol<1e-5) break
  }
  return(list(A,score))
}

########################
#posterior mean of matrix
#########################

score<-list()
A<-matrix(0,nrow=length(data$NI),ncol=data$n)
for(i in 1:length(A_list)){
  A=A+A_list[[i]][[1]]
  score[[i]]<-A_list[[i]][[2]]
}
A=A/length(A_list)

A_avg=A
#score analysis
index<-450
plot(1:index,score[[1]][1:index],pch=16,cex=0.5,type='l',ylim=c(min(unlist(score),na.rm=TRUE),max(unlist(score),na.rm=TRUE)))
for(i in 1:length(A_list)){
  lines(score[[i]][1:index],pch=16,col=i,cex=0.5)
}

score

#################################
##results and visualisation
################################

A=A_avg


#mean alpha values for locally acquried cases
mean(colMeans(A[,NI]))]

#mean alpha values for imported cases
mean(colMeans(A[,-NI]))


#
par(mfrow=c(1,2))

As<-A/rowSums(A)
#As<-As[complete.cases(As),]
As[is.nan(As)]=0




#define reproductive numbers for each case
d$Rt<-colSums(As)
mean(d$Rt)

######################################
#plot reproduction numbers over time
#######################################

f=data.frame(x=d$Date.symptoms.year,y=colSums(As,na.rm=TRUE))
library(mgcv)
g=gam(y~s(x),data=f)  

plot(f,type='b',pch=16,col='red',xlab='Date',ylab='R')
lines(f$x,g$fitted.values,col='blue')
abline(h=1)

mean(colSums(As),na.rm=T)



match.cols<-function(val,n){
  colfunc <- colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick","darkred"), space = "rgb",bias=5)#colors
  #	colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"), space = "rgb",bias=1)
  
  col<-data.frame(val=seq(min(val),max(val),length.out=n),col=colfunc(n))
  out<-rep(NA,length(col))
  for(i in 1:length(val)){
    out[i]<-as.character(col[which.min(abs(col$val-val[i])),'col'])
  }	
  return(out)
}
cols<-match.cols(As,1000)
cols[As==0]='white'
mat<-matrix(cols,nrow=nrow(As),ncol=ncol(As))

#par(mfrow=c(1,2))
x=1:ncol(As)
y=1:nrow(As)
plot(x[1],y[1],col=mat[1,1],pch=15,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),bty="n" , xaxt="n", yaxt="n",xlab='Infector',ylab='Infectee')

for(i in 1:nrow(As)){
  for(j in 1:ncol(As)){
    points(x[j],y[i],col=mat[i,j],pch=15)
  }
}
for(i in 1:nrow(As)){
  if(sum(As[i,])==0){
    for(j in 1:ncol(As)){
      points(x[j],y[i],col='yellow',pch=22)
    }
  }
}


text(x=x,y=-1 , cex=0.5,
     labels =d$Date , srt = 90, pos = 1, xpd = TRUE)
text(y=y,x=-2, cex=0.5,
     labels = d$Date[d$I==0],  pos = 2, xpd = TRUE,offset=-0.5)  


# perform distribution analysis
Rt<-colSums(As)
Rt[Rt==0]<-1e-4

library(fitdistrplus)
fitg <- fitdist(Rt, "gamma")
fitg$aic
fitg <- fitdist(Rt, "exp")
fitg$aic
fitg <- fitdist(Rt, "lnorm")
fitg$aic

rbinom(1000,10,0.21)

fitg <- fitdist(Rt, "gamma")
mean(rgamma(10000,shape=fitg$estimate[1],rate=fitg$estimate[2]))
pgamma(1,fitg$estimate[1],fitg$estimate[2])

rate=fitg$estimate[2]
prob = pgamma(1,shape=fitg$estimate[1],rate=rate)

while(prob<=0.95){
  rate=rate+0.01
  prob = pgamma(1,shape=fitg$estimate[1],rate=rate)
}
mean(rgamma(10000,shape=fitg$estimate[1],rate=rate))

cbind(2016:2030,predict.gam(g,newdata=data.frame(x=2016:2030)))

#### Izzy month analysis
Rt_mat<-matrix(nrow=12,ncol=length(A_list))
for(i in 1:length(A_list)){
  A2=A_list[[i]][[1]]
  A2s<-A2/rowSums(A2)
  A2s[is.nan(A2s)]=0
  Rts<-colSums(A2s)
  Rt_mat[,i]<-as.matrix(aggregate(list(Rts),by=list(d$Month),FUN=mean))[,2]
}
avg<-apply(Rt_mat,1,mean)
q<-apply(Rt_mat,1,quantile,probs=c(0.05,0.95))

plot(1:12,avg,pch=16,ylim=c(0,2))
arrows(x, q[1,], x, q[2,], length=0.05, angle=90, code=3)
abline(h=mean(Rt_mat),col='blue')

