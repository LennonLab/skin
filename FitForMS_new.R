###############################################################
####This is the new anaysis based on the analysis details.pdf
###############################################################
setwd("~/Documents/Work/LennonSkin/skin/")
setwd("~/GitHub/skin")
#read in data 
ages.raw <- read.table("data/skin.age.txt", sep = "\t", header = TRUE)
#remove subject "G1"
ages <- ages.raw[ ! ages.raw$subject == "G1", ]
#create dummy variables for multinomial variables
#shower frequency
show.f2<-ifelse(ages$show.freq==2,1,0)
show.f3<-ifelse(ages$show.freq==3,1,0)
#shower length
show.l2<-ifelse(ages$show.leng==2,1,0)
show.l3<-ifelse(ages$show.leng==3,1,0)
#last shower
show.last2<-ifelse(ages$last.show==2,1,0)
show.last3<-ifelse(ages$last.show==3,1,0)
#create new data.frame with all goodies.
ages.2<-data.frame(ages[,c(1:6,11:16)],show.f2,show.f3,show.l2,show.l3,show.last2,show.last3)
#set up a design matrix that we might use later
X<-data.frame(ages.2[,-c(1,3:5)])
#center age for later interp
mean(ages$age)
X$age<-X$age-mean(X$age)

####look at raw correlations####
cor(ages.2[,c(3,4,5)]/ages.2$total)
#this uphold our idea that all cells are usually in an active or dormant state
#and that dormant cells are more likely to die

####Do the analysis####
#Active
summary(f.act<-glm.nb(ages$active~.+offset(log(X$total)),data=X[,-2]))
#Dormant
summary(f.dor<-glm.nb(ages$dormant~.+offset(log(X$total)),data=X[,-2]))
#Dead
summary(f.ded<-glm.nb(ages$dead~.+offset(log(X$total)),data=X[,-2]))
#Total
summary(f.tot<-glm.nb(ages$total~.,data=X[,-2],maxit=1000))
#look at the predicted prop in category for an individual of mean age and reference 
#state for all categorical variables
#this should be very close to 1
(exp(c(coef(f.act)[1],coef(f.dor)[1],coef(f.ded)[1])))
sum((exp(c(coef(f.act)[1],coef(f.dor)[1],coef(f.ded)[1]))))

####look at which variables are significant####
#Active
rownames(summary(f.act)$coef)[which(summary(f.act)$coef[,4]<0.05)]
#Dormant
rownames(summary(f.dor)$coef)[which(summary(f.dor)$coef[,4]<0.05)]
#Dead
rownames(summary(f.ded)$coef)[which(summary(f.ded)$coef[,4]<0.05)]
#Total
rownames(summary(f.tot)$coef)[which(summary(f.tot)$coef[,4]<0.05)]
summary(f.tot)$coef[which(summary(f.tot)$coef[,4]<0.05)]
#age and shower length both significantly reduce cells and affect distribution of cells
#eczema increases total cell number but it does not affect the distribution across metabolic classes

##look at coef: factors that are of one sign for active and dead and the opposite for
##dormancy are the potential dormance cues
coef.matrix<-rbind(coef(f.act),coef(f.dor),coef(f.ded))
rownames(coef.matrix)<-c('active','dormant','dead')
coef.matrix
coef.matrix[,"age"]#age is a dormancy inducing mechanism (and reduces total cell number)
coef.matrix[,"show.f3"]#show frequncy is a dormancy reducing mechanism
coef.matrix[,"show.l3"]#show length is not a dormancy mechanism, but 
#works to preferentially remove dead cells--this interp is supported by the 
#significant effects on total cell count see detail below

#Notice that shower freq and last show should be on average inversly related and the 
#coefs show that they tell the same story shower freq tends to reduce dormancy
#so, the longer its been since showering the higher proportion is in the dormant category

#look at the distibution of cells into each metabolic category of a subject of average age
#and all other factors set to lowest (reference) level
exp(coef(f.tot)[1])#total cell numbers
exp(coef.matrix[,1]) #almost all cells are in dead category
#remember there is some error which is why the proportions don't sum to exactly 1

#look at the effect of long showers, it reduces the total number of cells A LOT, disproportionately 
#at the expence of dead ones
exp(coef(f.tot)[1])*exp(coef(f.tot)["show.l3"])
exp(coef.matrix[,1])*exp(coef.matrix[,"show.l3"])
#thus result of long showers is fewer total cells and a smaller proportion of them
#in the dead category

#eczema has a large effect on the total number of cells, but each metabolic category seems 
#to be affected proportionaltely
exp(coef(f.tot)[1])#total cell numbers
exp(coef.matrix[,1]) #almost all cells are in dead category

exp(coef(f.tot)[1])*exp(coef(f.tot)["ecze"])
exp(coef.matrix[,1])*exp(coef.matrix[,"ecze"])

#Residual correlation
pairs(data.frame(active=resid(f.act),dormant=resid(f.dor),dead=resid(f.ded)))
cor(data.frame(active=resid(f.act),dormant=resid(f.dor),dead=resid(f.ded)))
##The strong residual correlations suggests that we are missing some
##major cues the cells are using to determine activity

####fit the rest of the SEM
#set up list for results
pred.fits<-list()
#fit models
dim(X)
for(i in 1:12)pred.fits[[i]]<-glm(X[,i+2]~age,data=X,family=binomial)
#print fits
for(i in 1:12)print(summary(pred.fits[[i]]))
#which ones have sig age effect?
which(sapply(pred.fits,function(x)summary(x)$coef[2,"Pr(>|z|)"])<0.052)
colnames(X)[which(sapply(pred.fits,function(x)summary(x)$coef[2,"Pr(>|z|)"])<0.052)+2]
##look at what we found
summary(pred.fits[[6]])#older people are more likely to drink alcohol
summary(pred.fits[[10]])#older people are less likely to take long showers
#for shower length 15<x<30 calculate Tjur's R^2
mean(predict(pred.fits[[10]],type='response')[which(show.l3==1)])-
  mean(predict(pred.fits[[10]],type='response')[which(show.l3==0)])

#test structure
fits<-c(list(f.act,f.ded,f.dor,f.tot),pred.fits)
sapply(fits,function(x)summary(x)$coefficients)
names(fits)<-names(ages.2)[3:18]
#set up adjacency matrix. This will allow us to find all the pairs of nodes 
#that are no connected directly by an arrow and therefore predicted by the model structure 
#to be conditionally indepdentant (i.e. indep. after accounting for all the arrows pointing at each)
A<-matrix(0,17,16)
rownames(A)<-names(ages.2)[c(2,7:18,6,3:5)]
colnames(A)<-names(ages.2)[c(7:18,6,3:5)]

#use fits list to fill in the adjacency matrix
for(i in 1:length(fits))A[names(coef(fits[[i]])[-1]),names(fits)[i]]<-1
#add effect of total
A["total",c("active","dead","dormant")]<-1
print(A)
#find missing links between our target variables 
miss.links<-NULL
for(i in 13:16)for(j in 1:i)if(A[j,i]==0)miss.links<-cbind(miss.links,
                                                             c(rownames(A)[j],colnames(A)[i]))
print(miss.links)
#test for significant correlations among them. Remember, our initial model said these should equal zero
pvals<-rep(0,ncol(miss.links))
for(i in 1:3)pvals[i]<-cor.test(resid(fits[miss.links[2,i]][[1]]),resid(fits[miss.links[1,i]][[1]]))$p.value
data.frame(t(miss.links),pvals)
#all three are significant. let's look at them
resid.cors<-list()
for(i in 1:3)resid.cors[[i]]<-cor.test(resid(fits[miss.links[2,i]][[1]]),resid(fits[miss.links[1,i]][[1]]))
for(i in 1:3)names(resid.cors)[i]<-paste(miss.links[1,i],miss.links[2,i],sep='-')
 #use fisher's combined test as a test of the structure of our model
1-pchisq(-2*sum(log(pvals[1:3])),2*3) #it fails...we know why.

#if we add those correlations to the model, we have no more degrees of freedom to test the 
#structure of the model. Our conclusion is that although we have found some
#mechanisms, we know now that we have we have missed some mechanisms responible
#for determining how cells are partitioned across the metabolic classes. 
