###this is just me playing around trying to use a simple markov model 
# to represent the hypothesized dynamic data generating process
# this is used as a reality check to ensure that our understanding of the 
#results of the analysis are supported
make.J<-function(x1){
  J<-matrix(0,3,3)
  colnames(J)<-rownames(J)<-c('active','dormant','dead')
  J[1,]<-c(1.1-.1*x1,.01+.1*x1,.19) #prop of active that go to each category each time step
  J[2,]<-c(.1,.9,0) #prop of dorm the more to each time step (10% become active, none die)
  J[3,]<-c(0,0,.85) #some dead are lost from system
  return(J)}
#simulate this for 100 subjects with cue of no cue
x=rep(c(0,1),each=50)#cues for each in
#simulate the model for each individual and capture the last time step
subs<-matrix(0,100,3)
colnames(subs)<-c('active','dormant','dead')
nsim<-50
for(j in 1:100){
  #set up matrix to capture time series
  J<-make.J(x[j])
  t=matrix(0,3,nsim)
  t[,1]<-c(rpois(1,1000),0,0)
  for(i in 2:nsim)t[,i]<-t(J)%*%t[,(i-1)]
  subs[j,]<-round(t[,nsim])
}
cor(subs/apply(subs,1,sum))

#look at stable distribtion among categories in each case
##1) no cue
J<-make.J(0)
eigen(t(J))
#stable distrbution for structure of cell proportions
print(s0<-eigen(t(J))$vectors[,1]/sum(eigen(t(J))$vectors[,1]))

##2) with dormance cue
J<-make.J(1)
eigen(t(J))
#stable distrbution for structure of cell proportions
print(s1<-eigen(t(J))$vectors[,1]/sum(eigen(t(J))$vectors[,1]))

#try our analysis on this should get the stable age dist results
sim.dat<-data.frame(subs,'total'=apply(subs,1,sum))
head(sim.dat)
f.sim.act<-glm(active~x,offset = log(total),data=sim.dat,family=poisson)
f.sim.dor<-glm(dormant~x,offset = log(total),data=sim.dat,family=poisson)
f.sim.ded<-glm(dead~x,offset = log(total),data=sim.dat,family=poisson)

coef.mat<-rbind(coef(f.sim.act),coef(f.sim.dor),coef(f.sim.ded))
coef.mat #x has predicted form of having neg coef for active and dead and 
#positive for dormant
#no cue case i.e. x=0
exp(coef.mat[,1]) #predicted stable distribution based on our analysis
s0 #actual stable distrubtion based on data generating process 
#with cue case x=1
exp(apply(coef.mat,1,sum)) #predicted stable distribution based on our analysis
s1 #predicted stable distribution based on our analysis
