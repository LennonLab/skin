#create function to erase whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#redo fits with better formulas
for(i in 1:12){
   fmla<-as.formula(paste(names(X)[i+2],'~age',sep=''))
   pred.fits[[i]]<-glm(fmla,data=X,family=binomial)
   }
#grab all the fits
fits<-list(f.act,f.dor,f.ded,f.tot)
#create graphviz file

sink(file = "~/Path/To/WorkingDir/network.dot")
cat("digraph{")
cat("\n")
#go thru each metabolic class fit, grab "sig" predictors and add to the graph
for(j in 1:4){
  fit<-fits[[j]]
  f<-formula(fit)
  f1<-as.character(f)
  rhs<-strsplit(f1[2],"\\$")[[1]][2];rhs
  #lhs<-trim(strsplit(f1[3],"\\+")[[1]])
  who<-names(coef(fit))[summary(fit)$coef[,4]<0.055];
 if(who[1]=="(Intercept)")who<-who[-1];
  for(i in 1:length(who)){
    if(coef(fit)[who[i]]>0)cat(paste(who[i],"->",rhs,"[style=dashed]")) else cat(paste(who[i],"->",rhs))
    cat("\n")
    }
#add total for each of the metabolic classes
    if(j<4)cat(paste("total","->",rhs))
  cat("\n")
}
#add the fit for show.l2~age
f1<-as.character(formula(pred.fits[[10]]))
cat(paste(f1[3],"->",f1[2]))
cat("\n")
cat("}")
sink()
#run dot(graphviz) to create the graphs from the input files we just created
system(paste("dot -Tpdf ~/Path/To/WorkingDir/network.dot -o ~/Path/To/WorkingDir/network.pdf"))

