###########################################################################
#This page shows the the logic behind the analysis, that it works and the
#inferences and interpretations that can be drawn from it
############################################################################

#####simuate data####
# Imagine that x1 is a process that cues cells to shift 
#from active to dormant state and the active cells are more likely
#to die than dormant cells
tot<-rpois(1200,1000)
#create a factor to act as the cue for dormancy
x1<-rep(c(0,1,0),each=400)
#put cell into categories
dor<-rpois(1200,tot*(.1+.3*x1))
act<-rpois(1200,tot*(.5-.1*x1))
ded<-rpois(1200,tot*(.4-.2*x1))

####look at raw correlation sturcture of data####
#the total number of dor cells should be negative correlated
#with the numbers of active and dead, while number of active and dead should
#be positively correlated
#count cor:
cor(cbind(dor,act,ded))
#prop cor: should be the same story
cor(cbind(dor/tot,act/tot,ded/tot))

####fit the model derived on analysis details.pdf
#fit the model we derived 
f.dor<-glm(dor~x1,offset = log(tot),family=poisson)
f.act<-glm(act~x1,offset = log(tot),family=poisson)
f.ded<-glm(ded~x1,offset = log(tot),family=poisson)
#prop in each category with x1 absent (should recover above values (.1,.5,.4))
print(z1<-c(exp(coef(f.dor)[1]),exp(coef(f.act)[1]),exp(coef(f.ded)[1])))
#should sum to about 1
sum(z1)
#prop in each category with x1 present (should be .4,.4,.2)
print(z2<-c(c(exp(sum(coef(f.dor))),exp(sum(coef(f.act))),exp(sum(coef(f.ded))))))
#still should sum to about 1
sum(z2)

#since we have completely captured the mechanism for cell counts, the 
#correlations of residuals should be near zero
cor(cbind(resid(f.dor),resid(f.act),resid(f.ded)))
pairs(cbind('dor'=resid(f.dor),'act'=resid(f.act),'ded'=resid(f.ded)))
#if we fail to capture the cue mechanism then we will observed residual correlation among 
#the categories

#don't include cue
f.dor1<-glm(dor~1,offset = log(tot),family=poisson)
f.act1<-glm(act~1,offset = log(tot),family=poisson)
f.ded1<-glm(ded~1,offset = log(tot),family=poisson)
#look at props
print(z3<-c(exp(coef(f.dor1)),exp(coef(f.act1)),exp(coef(f.ded1))))
#should sum to about 1
sum(z3)
#numbers in each catgory still correlated
cor(cbind(resid(f.dor1),resid(f.act1),resid(f.ded1)))
pairs(cbind('dor'=resid(f.dor1),'act'=resid(f.act1),'ded'=resid(f.ded1)))

