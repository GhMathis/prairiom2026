
###graphical representation of communities
########################
library(igraph)

func.plot.matrix<-function(x,s.col=ls.col,x.title=NULL,y.lab=NULL,opposed.y=FALSE,no.ylab=FALSE){
  min<-min(x)
  max<-max(x)
  if(is.null(y.lab)){
	yLabels<-rownames(x)
  }
  else{
	yLabels<-y.lab
  }
  xLabels<-colnames(x)
  title<-x.title
  if(is.null(xLabels)){
    xLabels<-c(1:ncol(x))
  }
  if(is.null(yLabels)){
    yLabels<-c(1:nrow(x))
  }
  #layout(matrix(data=c(1,2),nrow=1,ncol=2),widths=c(4,1),heights=c(1,1))
  
  reverse<-nrow(x):1
  yLabels<-yLabels[reverse]
  x<-x[reverse,]
  par(mar=c(10,10,5,5))
  image(1:length(xLabels),1:length(yLabels),t(x),col=s.col, xlab="", ylab="",axes=FALSE,zlim=c(min,max))
  if(!no.ylab){
	title(ylab="Daphnia pop replicates", line=7, cex.lab=1.2)
  }
  title(xlab="", line=7, cex.lab=1.2)
  if(!is.null(title)){
    title(x.title)
  }
  if(opposed.y){
  axis(RIGHT<-4,at=1:length(yLabels), labels=as.factor(as.character(yLabels)),las= 2,cex.axis=0.5)
  axis(RIGHT<-4,at=1:length(yLabels),labels=rep("",length(yLabels)),las=2,cex.axis<-0.5)
  }
  else{
  #axis(BELOW<-1,at=1:length(xLabels),labels=rep("",length(xLabels)),las =2,cex.axis=0.5)
  axis(LEFT<-2,at=1:length(yLabels), labels=as.factor(as.character(yLabels)),las= 2,cex.axis=0.5)
  #axis(BELOW<-1,at=1:length(xLabels),labels=rep("",length(xLabels)),las =2,cex.axis=0.5)
  axis(LEFT<-2,at=1:length(yLabels),labels=rep("",length(yLabels)),las=2,cex.axis<-0.5)
  }
}

func.plot.matrix.byGroup<-function(x,group,s.col=ls.col,x.title=NULL,y.lab=NULL,opposed.y=FALSE,no.ylab=FALSE){
  x<-x+matrix(2*(as.numeric(group)-1),nrow=dim(x)[1],ncol=dim(x)[2],byrow=F)
  min<-min(x)
  max<-max(x)
  if(is.null(y.lab)){
	yLabels<-rownames(x)
  }
  else{
	yLabels<-y.lab
  }
  xLabels<-colnames(x)
  title<-x.title
  if(is.null(xLabels)){
    xLabels<-c(1:ncol(x))
  }
  if(is.null(yLabels)){
    yLabels<-c(1:nrow(x))
  }
  #layout(matrix(data=c(1,2),nrow=1,ncol=2),widths=c(4,1),heights=c(1,1))
  
  reverse<-nrow(x):1
  yLabels<-yLabels[reverse]
  x<-x[reverse,]
  par(mar=c(10,10,5,5))
  image(1:length(xLabels),1:length(yLabels),t(x),col=s.col, xlab="", ylab="",axes=FALSE,zlim=c(min,max))
  if(!no.ylab){
	title(ylab="Daphnia pop replicates", line=7, cex.lab=1.2)
  }
  title(xlab="", line=7, cex.lab=1.2)
  if(!is.null(title)){
    title(x.title)
  }
  if(opposed.y){
  #axis(RIGHT<-4,at=1:length(yLabels), labels=as.factor(as.character(yLabels)),las= 2,cex.axis=0.5)
  #axis(RIGHT<-4,at=1:length(yLabels),labels=rep("",length(yLabels)),las=2,cex.axis<-0.5)
  }
  else{
  #axis(BELOW<-1,at=1:length(xLabels),labels=rep("",length(xLabels)),las =2,cex.axis=0.5)
  #axis(LEFT<-2,at=1:length(yLabels), labels=as.factor(as.character(yLabels)),las= 2,cex.axis=0.5)
  #axis(BELOW<-1,at=1:length(xLabels),labels=rep("",length(xLabels)),las =2,cex.axis=0.5)
  axis(LEFT<-2,at=1:length(yLabels),labels=rep("",length(yLabels)),las=2,cex.axis<-0.5)
  }
}

###CCA: extracting chi2, F, etc. and speeding up comparisons #deprecated functions to remove?
########################

nulltozero<-function(x){
	if(is.null(x)){
		x<-0
	}
	x
}

F.computation<-function(ccaobj,size){
	if(is.null(ccaobj)){
		0
	}
	else{
		ccaobj$CCA$tot.chi*(size-ccaobj$CA$rank-nulltozero(ccaobj$pCCA$rank))/(ccaobj$CA$tot.chi*ccaobj$CCA$qrank)
	}
}

chi2.computation<-function(ccaobj){
ccaobj$CCA$tot.chi
}

extractTerms<-function(model){
    paste0(formula(model$terms))[3]

}

only.F<-function(incid,var1,var2,var3){
  n<-dim(incid)[1]
  modules<-tab.disjonctif(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
  one.1<-cca(modules ~ var1, na.action=na.omit)
  one.2<-cca(modules ~ var2, na.action=na.omit)
  one.3<-cca(modules ~ var3, na.action=na.omit)
  two.1<-cca(modules ~ var1 + var2, na.action=na.omit)
  two.2<-cca(modules ~ var1 + var3, na.action=na.omit)
  two.3<-cca(modules ~ var2 + var3, na.action=na.omit)
  three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
  one.1.alone<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
  one.2.alone<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
  one.3.alone<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
  one.1.vs.3<-cca(modules ~ var1 +  Condition(var3), na.action=na.omit)
  one.1.vs.2<-cca(modules ~ var1 +  Condition(var2), na.action=na.omit)
  one.2.vs.3<-cca(modules ~ var2 +  Condition(var3), na.action=na.omit)
  one.2.vs.1<-cca(modules ~ var2 +  Condition(var1), na.action=na.omit)
  one.3.vs.1<-cca(modules ~ var3 +  Condition(var1), na.action=na.omit)
  one.3.vs.2<-cca(modules ~ var3 +  Condition(var2), na.action=na.omit)
  xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
  
  sapply(1:16,function(x) F.computation(xxx[[x]],n))
}

only.F.alt <- function(memberships, var1, var2, var3) {
  n <- length(memberships)
  
  modules <-
    tab.disjonctif(as.factor(memberships))
  
  one.1 <- cca(modules ~ var1, na.action = na.omit)
  one.2 <- cca(modules ~ var2, na.action = na.omit)
  one.3 <- cca(modules ~ var3, na.action = na.omit)
  two.1 <- cca(modules ~ var1 + var2, na.action = na.omit)
  two.2 <- cca(modules ~ var1 + var3, na.action = na.omit)
  two.3 <- cca(modules ~ var2 + var3, na.action = na.omit)
  three <- cca(modules ~ var1 + var2 + var3, na.action = na.omit)
  one.1.alone <-
    cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action = na.omit)
  one.2.alone <-
    cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action = na.omit)
  one.3.alone <-
    cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action = na.omit)
  one.1.vs.3 <-
    cca(modules ~ var1 +  Condition(var3), na.action = na.omit)
  one.1.vs.2 <-
    cca(modules ~ var1 +  Condition(var2), na.action = na.omit)
  one.2.vs.3 <-
    cca(modules ~ var2 +  Condition(var3), na.action = na.omit)
  one.2.vs.1 <-
    cca(modules ~ var2 +  Condition(var1), na.action = na.omit)
  one.3.vs.1 <-
    cca(modules ~ var3 +  Condition(var1), na.action = na.omit)
  one.3.vs.2 <-
    cca(modules ~ var3 +  Condition(var2), na.action = na.omit)
  
  xxx <-
    list(
      one.1,
      one.2,
      one.3,
      two.1,
      two.2,
      two.3,
      three,
      one.1.alone,
      one.2.alone,
      one.3.alone,
      one.1.vs.3,
      one.1.vs.2,
      one.2.vs.3,
      one.2.vs.1,
      one.3.vs.1,
      one.3.vs.2
    )
  
  sapply(1:16, function(x)
    F.computation(xxx[[x]], n))
}

only.F.alt2 <- function(memberships, var1, var2) {
  n <- length(memberships)
  
  modules <-
    dummy(as.factor(memberships), levelsToKeep = levels(as.factor(memberships)))
  
  one.1 <- cca(modules ~ var1, na.action = na.omit)
  one.2 <- cca(modules ~ var2, na.action = na.omit)
  two.1 <- cca(modules ~ var1 + var2, na.action = na.omit)
  one.1.alone <-
    cca(modules ~ var1 + Condition(var2), na.action = na.omit)
  one.2.alone <-
    cca(modules ~ Condition(var1) + var2, na.action = na.omit)
  
  xxx <-
    list(
      one.1,
      one.2,
      two.1,
      one.1.alone,
      one.2.alone)
  
  sapply(1:5, function(x)
    F.computation(xxx[[x]], n))
}

only.F2<-function(incid,incid1,var2,var3){
n<-dim(incid)[1]
modules<-dummy(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
var1<-as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid1)))$mem[1:n])
n<-n-1
three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
two.1<-cca(modules ~ Condition(var1) + var2 + var3, na.action=na.omit)
two.2<-cca(modules ~ var1 + Condition(var2) + var3, na.action=na.omit)
two.3<-cca(modules ~ var1 + var2 + Condition(var3), na.action=na.omit)
one.1<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
one.2<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
one.3<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)

F.three<-three$CCA$tot.chi*(n-three$CCA$qrank)/(three$CA$tot.chi*three$CCA$qrank)
F.two.1<-two.1$CCA$tot.chi*(n-two.1$CCA$qrank-two.1$pCCA$rank)/(two.1$CA$tot.chi*two.1$CCA$qrank)
F.two.2<-two.2$CCA$tot.chi*(n-two.2$CCA$qrank-two.2$pCCA$rank)/(two.2$CA$tot.chi*two.2$CCA$qrank)
F.two.3<-two.3$CCA$tot.chi*(n-two.3$CCA$qrank-two.3$pCCA$rank)/(two.3$CA$tot.chi*two.3$CCA$qrank)
F.one.1<-one.1$CCA$tot.chi*(n-one.1$CCA$qrank-one.1$pCCA$rank)/(one.1$CA$tot.chi*one.1$CCA$qrank)
F.one.2<-one.2$CCA$tot.chi*(n-one.2$CCA$qrank-one.2$pCCA$rank)/(one.2$CA$tot.chi*one.2$CCA$qrank)
F.one.3<-one.3$CCA$tot.chi*(n-one.3$CCA$qrank-one.3$pCCA$rank)/(one.3$CA$tot.chi*one.3$CCA$qrank)

c(F.three,F.two.1,F.two.2,F.two.3,F.one.1,F.one.2,F.one.3)
}

only.F3<-function(cle,cle1,var2,var3){
n<-length(cle)
nl<-nlevels(as.factor(cle))
if(nl==1){
	rep(0,16)
}
else{
	modules<-dummy(as.factor(cle))
	var1<-as.factor(cle1)
	nl1<-nlevels(var1)
	one.1<-cca(modules ~ var1, na.action=na.omit)
	one.2<-cca(modules ~ var2, na.action=na.omit)
	one.3<-cca(modules ~ var3, na.action=na.omit)
	two.1<-cca(modules ~ var1 + var2, na.action=na.omit)
	two.2<-cca(modules ~ var1 + var3, na.action=na.omit)
	two.3<-cca(modules ~ var2 + var3, na.action=na.omit)
	three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
	one.1.alone<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)	
	one.1.vs.3<-cca(modules ~ var1 +  Condition(var3), na.action=na.omit)
	one.1.vs.2<-cca(modules ~ var1 +  Condition(var2), na.action=na.omit)
	one.2.vs.3<-cca(modules ~ var2 +  Condition(var3), na.action=na.omit)
	one.3.vs.2<-cca(modules ~ var3 +  Condition(var2), na.action=na.omit)
	if(nl1>1){
		one.2.alone<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)		
		one.3.alone<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
		one.2.vs.1<-cca(modules ~ var2 +  Condition(var1), na.action=na.omit)
		one.3.vs.1<-cca(modules ~ var3 +  Condition(var1), na.action=na.omit)
	}
	else{
		one.2.alone<-one.2.vs.3
		one.3.alone<-one.3.vs.2
		one.2.vs.1<-one.2
		one.3.vs.1<-one.3
		one.1<-NULL
		one.1.alone<-NULL
		one.1.vs.3<-NULL
		one.1.vs.2<-NULL
		
	}	
	xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
	sapply(1:16,function(x) F.computation(xxx[[x]],n))
}
}

analysis.function<-function(incid,var1,var2,var3,configs=NULL){#industrial processing of CCA, two permutations
n<-dim(incid)[1]
modules<-dummy(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
one.1<-cca(modules ~ var1, na.action=na.omit)
one.2<-cca(modules ~ var2, na.action=na.omit)
one.3<-cca(modules ~ var3, na.action=na.omit)
two.1<-cca(modules ~ var1 + var2, na.action=na.omit)
two.2<-cca(modules ~ var1 + var3, na.action=na.omit)
two.3<-cca(modules ~ var2 + var3, na.action=na.omit)
three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
one.1.alone<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
one.2.alone<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
one.3.alone<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
one.1.vs.3<-cca(modules ~ var1 +  Condition(var3), na.action=na.omit)
one.1.vs.2<-cca(modules ~ var1 +  Condition(var2), na.action=na.omit)
one.2.vs.3<-cca(modules ~ var2 +  Condition(var3), na.action=na.omit)
one.2.vs.1<-cca(modules ~ var2 +  Condition(var1), na.action=na.omit)
one.3.vs.1<-cca(modules ~ var3 +  Condition(var1), na.action=na.omit)
one.3.vs.2<-cca(modules ~ var3 +  Condition(var2), na.action=na.omit)
xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)

Forms<-sapply(1:16,function(x) extractTerms(xxx[[x]]))
F<-sapply(1:16,function(x) F.computation(xxx[[x]],n))
P<-sapply(1:16,function(x) anova(xxx[[x]],permutations = how(nperm=NPERM-1))$P[1])

if(is.null(configs)) {
	results<-data.frame("Formulas"=Forms,"F"=F,"P"=P)
}
else {
	depth<-dim(configs)[3]
	all.F<-sapply(1:depth,function(x) only.F(configs[,,x],var1,var2,var3))
	all.ecdf<-apply(all.F,1,ecdf)
	P2<-sapply(1:7,function(x) 1-((all.ecdf[[x]])(F[x])))
	results<-data.frame("Formulas"=Forms,"F"=F,"P"=P,"P2"=P2)
}
results
}

DSFP.cca<-function(anovaobjectcca){
  c(anovaobjectcca$Df[1],anovaobjectcca$ChiSquare[1],anovaobjectcca[[3]][1],anovaobjectcca$P[1])
}

analysis.function.alt<-function(incid,memb_obs,var1,var2,var3,NPERM =1000, configs=NULL, plt = NULL){#industrial processing of CCA, two permutations
  n <- dim(incid)[1]
  modules <- tab.disjonctif(as.factor(memb_obs))
  
  one.1 <- cca(modules ~ var1, na.action=na.omit)
  one.2 <- cca(modules ~ var2, na.action=na.omit)
  one.3 <- cca(modules ~ var3, na.action=na.omit)
  two.1 <- cca(modules ~ var1 + var2, na.action=na.omit)
  two.2 <- cca(modules ~ var1 + var3, na.action=na.omit)
  two.3 <- cca(modules ~ var2 + var3, na.action=na.omit)
  three <- cca(modules ~ var1 + var2 + var3, na.action=na.omit)
  one.1.alone <- cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
  one.2.alone <- cca(modules ~  var2 + Condition(var1) +Condition(var3), na.action=na.omit)
  one.3.alone <- cca(modules ~ var3 + Condition(var1) + Condition(var2) , na.action=na.omit)
  one.1.vs.3 <- cca(modules ~ var1 +  Condition(var3), na.action=na.omit)
  one.1.vs.2 <- cca(modules ~ var1 +  Condition(var2), na.action=na.omit)
  one.2.vs.3 <- cca(modules ~ var2 +  Condition(var3), na.action=na.omit)
  one.2.vs.1 <- cca(modules ~ var2 +  Condition(var1), na.action=na.omit)
  one.3.vs.1 <- cca(modules ~ var3 +  Condition(var1), na.action=na.omit)
  one.3.vs.2 <- cca(modules ~ var3 +  Condition(var2), na.action=na.omit)
  
  xxx <- list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
  
  chi2 <- sapply(1:16, function(x) chi2.computation(xxx[[x]]))
  Forms <- sapply(1:16, function(x) extractTerms(xxx[[x]]))
  F_stat <- sapply(1:16, function(x) F.computation(xxx[[x]],n))
  if(is.null(plt)){
    anova_perm <- t(sapply(1:16, function(x) DSFP.cca(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))
  }else{
    anova_perm <- t(sapply(1:16, function(x) DSFP.cca(anova(xxx[[x]],
                                                          permutations = how(plots = Plots(strata = plt, type = "free") , nperm = NPERM-1)))))
  }
  
  if(is.null(configs)) {
    results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F_stat,"P"=anova_perm[,4])
  } else {
    depth <- dim(configs)[2]
    all.F <- sapply(1:depth,function(x) only.F.alt(memberships = configs[,x],var1,var2,var3))
    all.ecdf <- apply(all.F,1,ecdf)
    P2 <- sapply(1:16,function(x) 1-((all.ecdf[[x]])(F_stat[x])))
    results <- data.frame("Formulas"=Forms,df =anova_perm[,1] ,"chi2"=chi2,"F"=F_stat,"P"=anova_perm[,4],"P2"=P2)
  }
  results
}




analysis.function.alt.4var<-function(incid,memb_obs,X,NPERM =1000, configs=NULL){#industrial processing of CCA, two permutations
  n <- dim(incid)[1]
  modules <- tab.disjonctif(as.factor(memb_obs))
  
  # Define the variables
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  X4 <- X[, 4]
  
  # Perform Canonical Correspondence Analysis
  aeghklno <- cca(modules ~ X1, na.action = na.omit)
  befiklmo <- cca(modules ~ X2, na.action = na.omit)
  cfgjlmno <- cca(modules ~ X3, na.action = na.omit)
  dhijkmno <- cca(modules ~ X4, na.action = na.omit)
  abefghiklmno <- cca(modules ~ X1 + X2, na.action = na.omit)
  acefghjklmno <- cca(modules ~ X1 + X3, na.action = na.omit)
  adeghijklmno <- cca(modules ~ X1 + X4, na.action = na.omit)
  bcefgijklmno <- cca(modules ~ X2 + X3, na.action = na.omit)
  bdefhijklmno <- cca(modules ~ X2 + X4, na.action = na.omit)
  cdfghijklmno <- cca(modules ~ X3 + X4, na.action = na.omit)
  abcefghijklmno <- cca(modules ~ X1 + X2 + X3, na.action = na.omit)
  abdefghijklmno <- cca(modules ~ X1 + X2 + X4, na.action = na.omit)
  acdefghijklmno <- cca(modules ~ X1 + X3 + X4, na.action = na.omit)
  bcdefghijklmno <- cca(modules ~ X2 + X3 + X4, na.action = na.omit)
  abcdefghijklmno <- cca(modules ~ X1 + X2 + X3 + X4, na.action = na.omit)
  
  # Define individual fractions
  a <- cca(modules ~ X1 + Condition(X2) + Condition(X3) + Condition(X4), na.action = na.omit)
  b <- cca(modules ~ X2 + Condition(X1) + Condition(X3) + Condition(X4), na.action = na.omit)
  c <- cca(modules ~ X3 + Condition(X1) + Condition(X2) + Condition(X4), na.action = na.omit)
  d <- cca(modules ~ X4 + Condition(X1) + Condition(X2) + Condition(X3), na.action = na.omit)
  # Define other fractions if needed...
  
  # Controlling 2 tables X
  ae <- cca(modules ~ X1 + Condition(X3) + Condition(X4), na.action = na.omit)
  ag <- cca(modules ~ X1 + Condition(X2) + Condition(X4), na.action = na.omit)
  ah <- cca(modules ~ X1 + Condition(X2) + Condition(X3), na.action = na.omit)
  be <- cca(modules ~ X2 + Condition(X3) + Condition(X4), na.action = na.omit)
  bf <- cca(modules ~ X2 + Condition(X1) + Condition(X4), na.action = na.omit)
  bi <- cca(modules ~ X2 + Condition(X1) + Condition(X3), na.action = na.omit)
  cf <- cca(modules ~ X3 + Condition(X1) + Condition(X4), na.action = na.omit)
  cg <- cca(modules ~ X3 + Condition(X2) + Condition(X4), na.action = na.omit)
  cj <- cca(modules ~ X3 + Condition(X1) + Condition(X2), na.action = na.omit)
  dh <- cca(modules ~ X4 + Condition(X2) + Condition(X3), na.action = na.omit)
  di <- cca(modules ~ X4 + Condition(X1) + Condition(X3), na.action = na.omit)
  dj <- cca(modules ~ X4 + Condition(X1) + Condition(X2), na.action = na.omit)
  
  # Controlling 1 table X
  aghn <- cca(modules ~ X1 + Condition(X2), na.action = na.omit)
  aehk <- cca(modules ~ X1 + Condition(X3), na.action = na.omit)
  aegl <- cca(modules ~ X1 + Condition(X4), na.action = na.omit)
  bfim <- cca(modules ~ X2 + Condition(X1), na.action = na.omit)
  beik <- cca(modules ~ X2 + Condition(X3), na.action = na.omit)
  befl <- cca(modules ~ X2 + Condition(X4), na.action = na.omit)
  cfjm <- cca(modules ~ X3 + Condition(X1), na.action = na.omit)
  cgjn <- cca(modules ~ X3 + Condition(X2), na.action = na.omit)
  cfgl <- cca(modules ~ X3 + Condition(X4), na.action = na.omit)
  dijm <- cca(modules ~ X4 + Condition(X1), na.action = na.omit)
  dhjn <- cca(modules ~ X4 + Condition(X2), na.action = na.omit)
  dhik <- cca(modules ~ X4 + Condition(X3), na.action = na.omit)
  
  xxx <- list( aeghklno = aeghklno, befiklmo = befiklmo, cfgjlmno = cfgjlmno,
                              dhijkmno = dhijkmno, abefghiklmno = abefghiklmno, acefghjklmno = acefghjklmno,
                              adeghijklmno = adeghijklmno, bcefgijklmno = bcefgijklmno, bdefhijklmno = bdefhijklmno,
                              cdfghijklmno = cdfghijklmno, abcefghijklmno = abcefghijklmno, abdefghijklmno = abdefghijklmno, 
                              acdefghijklmno = acdefghijklmno, bcdefghijklmno = bcdefghijklmno, abcdefghijklmno = abcdefghijklmno,
                              a = a, b = b, c = c, d = d, ae = ae, ag = ag, ah = ah, be = be, bf = bf, bi = bi, cf = cf, cg = cg, cj = cj,
                              dh = dh, di = di, dj = dj, aghn = aghn, aehk = aehk, aegl = aegl, bfim = bfim, beik = beik, befl = befl, cfjm = cfjm,
                              cgjn = cgjn, cfgl = cfgl, dijm = dijm, dhjn = dhjn, dhik = dhik )
  
  chi2 <- sapply(1:length(xxx), function(x) chi2.computation(xxx[[x]]))
  Forms <- sapply(1:length(xxx), function(x) extractTerms(xxx[[x]]))
  F_stat <- sapply(1:length(xxx), function(x) F.computation(xxx[[x]],n))
  anova_perm <- t(sapply(1:length(xxx), function(x) DSFP.cca(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))

  if(is.null(configs)) {
    results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F_stat,"P"=anova_perm[,4])
  } else {
    depth <- dim(configs)[2]
    all.F <- sapply(1:depth,function(x) only.F.alt(memberships = configs[,x],var1,var2,var3))
    all.ecdf <- apply(all.F,1,ecdf)
    P2 <- sapply(1:length(xxx),function(x) 1-((all.ecdf[[x]])(F_stat[x])))
    results <- data.frame("Formulas"=Forms,df =anova_perm[,1] ,"chi2"=chi2,"F"=F_stat,"P"=anova_perm[,4],"P2"=P2)
  }
  results
}
#ajout Virginie pour cas où on regarde que 2 variables
analysis.function.alt2<-function(incid,memb_obs,var1,var2,configs=NULL){#industrial processing of CCA, two permutations
  n <- dim(incid)[1]
  modules <- dummy(as.factor(memb_obs), levelsToKeep = levels(as.factor(memb_obs)))
  
  one.1 <- cca(modules ~ var1, na.action=na.omit)
  one.2 <- cca(modules ~ var2, na.action=na.omit)
  two.1 <- cca(modules ~ var1 + var2, na.action=na.omit)
  one.1.alone <- cca(modules ~ var1 + Condition(var2), na.action=na.omit)
  one.2.alone <- cca(modules ~ Condition(var1) + var2, na.action=na.omit)
  
  xxx <- list(one.1,one.2,two.1,one.1.alone,one.2.alone)
  
  chi2 <- sapply(1:5, function(x) chi2.computation(xxx[[x]]))
  Forms <- sapply(1:5, function(x) extractTerms(xxx[[x]]))
  F <- sapply(1:5, function(x) F.computation(xxx[[x]],n))
  P <- sapply(1:5, function(x) anova(xxx[[x]],permutations = how(nperm=NPERM-1))$P[1])
  
  if(is.null(configs)) {
    results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P)
    } else {
    depth <- dim(configs)[2]
    all.F <- sapply(1:depth,function(x) only.F.alt2(configs[,x],var1,var2))
    all.ecdf <- apply(all.F,1,ecdf)
    P2 <- sapply(1:5,function(x) 1-((all.ecdf[[x]])(F[x])))
    results <- data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P,"P2"=P2)
  }
  results
}



analysis.function2<-function(incid,incid1,var2,var3,configs=NULL,configs1=NULL){#industrial processing of CCA, two permutations
n<-dim(incid)[1]
modules<-dummy(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
var1<-as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid1)))$mem[1:n])
n<-n-1
three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
two.1<-cca(modules ~ Condition(var1) + var2 + var3, na.action=na.omit)
two.2<-cca(modules ~ var1 + Condition(var2) + var3, na.action=na.omit)
two.3<-cca(modules ~ var1 + var2 + Condition(var3), na.action=na.omit)
one.1<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
one.2<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
one.3<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)

Forms<- c(extractTerms(three),extractTerms(two.1),extractTerms(two.2),extractTerms(two.3),extractTerms(one.1),extractTerms(one.2),extractTerms(one.3))

F.three<-three$CCA$tot.chi*(n-three$CCA$qrank)/(three$CA$tot.chi*three$CCA$qrank)
F.two.1<-two.1$CCA$tot.chi*(n-two.1$CCA$qrank-two.1$pCCA$rank)/(two.1$CA$tot.chi*two.1$CCA$qrank)
F.two.2<-two.2$CCA$tot.chi*(n-two.2$CCA$qrank-two.2$pCCA$rank)/(two.2$CA$tot.chi*two.2$CCA$qrank)
F.two.3<-two.3$CCA$tot.chi*(n-two.3$CCA$qrank-two.3$pCCA$rank)/(two.3$CA$tot.chi*two.3$CCA$qrank)
F.one.1<-one.1$CCA$tot.chi*(n-one.1$CCA$qrank-one.1$pCCA$rank)/(one.1$CA$tot.chi*one.1$CCA$qrank)
F.one.2<-one.2$CCA$tot.chi*(n-one.2$CCA$qrank-one.2$pCCA$rank)/(one.2$CA$tot.chi*one.2$CCA$qrank)
F.one.3<-one.3$CCA$tot.chi*(n-one.3$CCA$qrank-one.3$pCCA$rank)/(one.3$CA$tot.chi*one.3$CCA$qrank)
F<-c(F.three,F.two.1,F.two.2,F.two.3,F.one.1,F.one.2,F.one.3)

P.three<-anova(three,permutations = how(nperm=NPERM-1))$P[1]
P.two.1<-anova(two.1,permutations = how(nperm=NPERM-1))$P[1]
P.two.2<-anova(two.2,permutations = how(nperm=NPERM-1))$P[1]
P.two.3<-anova(two.3,permutations = how(nperm=NPERM-1))$P[1]
P.one.1<-anova(one.1,permutations = how(nperm=NPERM-1))$P[1]
P.one.2<-anova(one.2,permutations = how(nperm=NPERM-1))$P[1]
P.one.3<-anova(one.3,permutations = how(nperm=NPERM-1))$P[1]
P<-c(P.three,P.two.1,P.two.2,P.two.3,P.one.1,P.one.2,P.one.3)

if(is.null(configs)) {
	results<-data.frame("Formulas"=Forms,"F"=F,"P"=P)
}
else {
	depth<-dim(configs)[3]
	all.F<-sapply(1:depth,function(x) only.F2(configs[,,x],configs1[,,x],var2,var3))
	all.ecdf<-apply(all.F,1,ecdf)
	P2<-sapply(1:7,function(x) 1-((all.ecdf[[x]])(F[x])))
	results<-data.frame("Formulas"=Forms,"F"=F,"P"=P,"P2"=P2)
}
results
}

analysis.function3<-function(incid,cle1,var2,var3,cleconfigs=NULL,cleconfigs1=NULL){#industrial processing of CCA, two permutations
n<-dim(incid)[1]
modules<-dummy(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
var1<-as.factor(cle1)
one.1<-cca(modules ~ var1, na.action=na.omit)
one.2<-cca(modules ~ var2, na.action=na.omit)
one.3<-cca(modules ~ var3, na.action=na.omit)
two.1<-cca(modules ~ var1 + var2, na.action=na.omit)
two.2<-cca(modules ~ var1 + var3, na.action=na.omit)
two.3<-cca(modules ~ var2 + var3, na.action=na.omit)
three<-cca(modules ~ var1 + var2 + var3, na.action=na.omit)
one.1.alone<-cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
one.2.alone<-cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
one.3.alone<-cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
one.1.vs.3<-cca(modules ~ var1 +  Condition(var3), na.action=na.omit)
one.1.vs.2<-cca(modules ~ var1 +  Condition(var2), na.action=na.omit)
one.2.vs.3<-cca(modules ~ var2 +  Condition(var3), na.action=na.omit)
one.2.vs.1<-cca(modules ~ var2 +  Condition(var1), na.action=na.omit)
one.3.vs.1<-cca(modules ~ var3 +  Condition(var1), na.action=na.omit)
one.3.vs.2<-cca(modules ~ var3 +  Condition(var2), na.action=na.omit)
xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
chi2<-sapply(1:16,function(x) chi2.computation(xxx[[x]]))
Forms<-sapply(1:16,function(x) extractTerms(xxx[[x]]))
F<-sapply(1:16,function(x) F.computation(xxx[[x]],n))
P<-sapply(1:16,function(x) anova(xxx[[x]],permutations = how(nperm=NPERM-1))$P[1])

if(is.null(cleconfigs)) {
	results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P)
}
else {
	depth<-dim(cleconfigs)[2]
	all.F<-sapply(1:depth,function(x) only.F3(cleconfigs[,x],cleconfigs1[,x],var2,var3))
	all.ecdf<-apply(all.F,1,ecdf)
	P2<-sapply(1:16,function(x) 1-((all.ecdf[[x]])(F[x])))
	results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P,"P2"=P2)
}
results
}

DSFP<-function(anovaobject){
c(anovaobject$D[1],anovaobject$Variance[1],anovaobject$F[1],anovaobject$P[1])
}

only.F.rda<-function(memberships,var1,var2,var3,dista="jaccard"){
n<-length(memberships)
modules<-dummy(as.factor(memberships))
one.1<-dbrda(modules ~ var1, distance = dista, na.action=na.omit)
one.2<-dbrda(modules ~ var2, distance = dista, na.action=na.omit)
one.3<-dbrda(modules ~ var3, distance = dista, na.action=na.omit)
two.1<-dbrda(modules ~ var1 + var2, distance = dista, na.action=na.omit)
two.2<-dbrda(modules ~ var1 + var3, distance = dista, na.action=na.omit)
two.3<-dbrda(modules ~ var2 + var3, distance = dista, na.action=na.omit)
three<-dbrda(modules ~ var1 + var2 + var3, distance = dista, na.action=na.omit)
one.1.alone<-dbrda(modules ~ var1 + Condition(var2) + Condition(var3), distance = dista, na.action=na.omit)
one.2.alone<-dbrda(modules ~ Condition(var1) + var2 + Condition(var3), distance = dista, na.action=na.omit)
one.3.alone<-dbrda(modules ~ Condition(var1) + Condition(var2) + var3, distance = dista, na.action=na.omit)
one.1.vs.3<-dbrda(modules ~ var1 +  Condition(var3), distance = dista, na.action=na.omit)
one.1.vs.2<-dbrda(modules ~ var1 +  Condition(var2), distance = dista, na.action=na.omit)
one.2.vs.3<-dbrda(modules ~ var2 +  Condition(var3), distance = dista, na.action=na.omit)
one.2.vs.1<-dbrda(modules ~ var2 +  Condition(var1), distance = dista, na.action=na.omit)
one.3.vs.1<-dbrda(modules ~ var3 +  Condition(var1), distance = dista, na.action=na.omit)
one.3.vs.2<-dbrda(modules ~ var3 +  Condition(var2), distance = dista, na.action=na.omit)
xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
res<-t(sapply(1:16,function(x) DSFP(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))
res[,3]
}

analysis.function.dbrda<-function(incid,var1,var2,var3,dista="jaccard",configs=NULL){#industrial processing of dbrda
n<-dim(incid)[1]
modules<-dummy(as.factor(cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(incid)))$mem[1:n]))
one.1<-dbrda(modules ~ var1, distance = dista, na.action=na.omit)
one.2<-dbrda(modules ~ var2, distance = dista, na.action=na.omit)
one.3<-dbrda(modules ~ var3, distance = dista, na.action=na.omit)
two.1<-dbrda(modules ~ var1 + var2, distance = dista, na.action=na.omit)
two.2<-dbrda(modules ~ var1 + var3, distance = dista, na.action=na.omit)
two.3<-dbrda(modules ~ var2 + var3, distance = dista, na.action=na.omit)
three<-dbrda(modules ~ var1 + var2 + var3, distance = dista, na.action=na.omit)
one.1.alone<-dbrda(modules ~ var1 + Condition(var2) + Condition(var3), distance = dista, na.action=na.omit)
one.2.alone<-dbrda(modules ~ Condition(var1) + var2 + Condition(var3), distance = dista, na.action=na.omit)
one.3.alone<-dbrda(modules ~ Condition(var1) + Condition(var2) + var3, distance = dista, na.action=na.omit)
one.1.vs.3<-dbrda(modules ~ var1 +  Condition(var3), distance = dista, na.action=na.omit)
one.1.vs.2<-dbrda(modules ~ var1 +  Condition(var2), distance = dista, na.action=na.omit)
one.2.vs.3<-dbrda(modules ~ var2 +  Condition(var3), distance = dista, na.action=na.omit)
one.2.vs.1<-dbrda(modules ~ var2 +  Condition(var1), distance = dista, na.action=na.omit)
one.3.vs.1<-dbrda(modules ~ var3 +  Condition(var1), distance = dista, na.action=na.omit)
one.3.vs.2<-dbrda(modules ~ var3 +  Condition(var2), distance = dista, na.action=na.omit)
xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
res<-t(sapply(1:16,function(x) DSFP(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))
F<-res[,3]
Forms<-sapply(1:16,function(x) extractTerms(xxx[[x]]))

if(is.null(configs)) {
	results<-data.frame("Formulas"=Forms,"df"=res[,1],"SSq"=res[,2],"F"=res[,3],"P"=res[,4])
}
else {
	depth<-dim(configs)[2]
	all.F<-sapply(1:depth,function(x) only.F.dbrda(configs[,x],var1,var2,var3,dista))
	all.ecdf<-apply(all.F,1,ecdf)
	P2<-sapply(1:16,function(x) 1-((all.ecdf[[x]])(F[x])))
	results<-data.frame("Formulas"=Forms,"df"=res[,1],"SSq"=res[,2],"F"=res[,3],"P"=res[,4],"P2"=P2)
}
results

}


###SVD
########################

incidence_to_traits<-function(incid,nc=NULL){
gr.svd<-svd(incid)
svd.L <- gr.svd$u
svd.S <- diag(gr.svd$d)
svd.Ssqrt <- structure(vapply(svd.S, sqrt, numeric(1)),dim=dim(svd.S))
svd.F<-svd.L %*% svd.Ssqrt
if(is.null(nc)){
	svd.F
}
else{
	svd.F[,1:nc]
}
}

probabilize<-function(predictions){
	mip<-min(predictions)
	map<-max(predictions)
	(predictions-mip)/(map-mip)
}

approx.from.svd<-function(incid, svd.L, svd.R, nvect){
	d<-dim(incid)
	appr<-svd.L[,1:nvect] %*% t(svd.R[,1:nvect])
	obs<-c(unlist(incid))
	probs<-probabilize(c(unlist(appr)))
	#thresh<-optim.thresh(obs,probs)$'max.sensitivity+specificity'
	roc_obj <- ROCR::prediction(probs, obs)
	roc_perf <- ROCR::performance(roc_obj, "tpr", "fpr")

	AUC <- ROCR::performance(roc_obj, measure = "auc")@y.values
	max_sum <- which.max(roc_perf@y.values[[1]] + (1 - roc_perf@x.values[[1]]))
	thresh <- roc_perf@alpha.values[[1]][max_sum]
	pred<-probs
	pred[probs>=thresh]<-1
	pred[probs<thresh]<-0
	array(pred,dim=d)
}

###RDA
########################

analysis.function.rdpg <-
  function(incid, var1, var2, var3, nc = NULL) {
    #industrial processing of varpart results
    gr.svd <- svd(incid)
    svd.L <- gr.svd$u
    svd.S <- diag(gr.svd$d)
    svd.Ssqrt <-
      structure(vapply(svd.S, sqrt, numeric(1)), dim = dim(svd.S))
    svd.F <- svd.L %*% svd.Ssqrt
    if (is.null(nc)) {
      c(
        varpart(svd.F, var1, var2, var3)$part$fract[["Adj.R.squared"]],
        varpart(svd.F, var1, var2, var3)$part$indfract[["Adj.R.squared"]],
        varpart(svd.F, var1, var2, var3)$part$contr1[["Adj.R.squared"]]
      )
    }
    else{
      c(
        varpart(svd.F[, 1:nc], var1, var2, var3)$part$fract[["Adj.R.squared"]],
        varpart(svd.F[, 1:nc], var1, var2, var3)$part$indfract[["Adj.R.squared"]],
        varpart(svd.F[, 1:nc], var1, var2, var3)$part$contr1[["Adj.R.squared"]]
      )
    }
  }

analysis.function.rdpg2 <-
  function(incid, var1, var2, nc = NULL) {
    #industrial processing of varpart results
    gr.svd <- svd(incid)
    svd.L <- gr.svd$u
    svd.S <- diag(gr.svd$d)
    svd.Ssqrt <-
      structure(vapply(svd.S, sqrt, numeric(1)), dim = dim(svd.S))
    svd.F <- svd.L %*% svd.Ssqrt
    if (is.null(nc)) {
      c(
        varpart(svd.F, var1, var2)$part$fract[["Adj.R.squared"]],
        varpart(svd.F, var1, var2)$part$indfract[["Adj.R.squared"]],
        varpart(svd.F, var1, var2)$part$contr1[["Adj.R.squared"]]
      )
    }
    else{
      c(
        varpart(svd.F[, 1:nc], var1, var2)$part$fract[["Adj.R.squared"]],
        varpart(svd.F[, 1:nc], var1, var2)$part$indfract[["Adj.R.squared"]],
        varpart(svd.F[, 1:nc], var1, var2)$part$contr1[["Adj.R.squared"]]
      )
    }
  }

analysis.function.rdpg.fromtraits<-function(traits,var1,var2,var3,nc=NULL){
svd.F<-traits
if(is.null(nc)){
	c(varpart(svd.F,var1,var2,var3)$part$fract[["Adj.R.square"]],varpart(svd.F,var1,var2,var3)$part$indfract[["Adj.R.square"]],varpart(svd.F,var1,var2,var3)$part$contr1[["Adj.R.square"]])
}
else{
	c(varpart(svd.F[,1:nc],var1,var2,var3)$part$fract[["Adj.R.square"]],varpart(svd.F[,1:nc],var1,var2,var3)$part$indfract[["Adj.R.square"]],varpart(svd.F[,1:nc],var1,var2,var3)$part$contr1[["Adj.R.square"]])
}
}

rsquared.function.rdpg<-function(incid,var1,var2,var3,nc=NULL){
gr.svd<-svd(incid)
svd.L <- gr.svd$u
svd.S <- diag(gr.svd$d)
svd.Ssqrt <- structure(vapply(svd.S, sqrt, numeric(1)),dim=dim(svd.S))
svd.F <- svd.L %*% svd.Ssqrt
if(is.null(nc)){
	c(varpart(svd.F,var1,var2,var3)$part$fract[["R.square"]])
}
else{
	c(varpart(svd.F[,1:nc],var1,var2,var3)$part$fract[["R.square"]])
}
}

rsquared.function.rdpg.fromtraits<-function(traits,var1,var2,var3,nc=NULL){
svd.F<-traits
if(is.null(nc)){
	c(varpart(svd.F,var1,var2,var3)$part$fract[["R.square"]])
}
else{
	c(varpart(svd.F[,1:nc],var1,var2,var3)$part$fract[["R.square"]])
}
}

