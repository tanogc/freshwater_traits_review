
# Setting working directory
setwd("your_folder")
res_fold<-paste(getwd(),"/results/",sep="")

library(reshape)
library(sqldf)
library(stringr)
library(MuMIn)
library(nnet)
library(plyr)
library(vegan)
library(ineq)
library(viridis)
library(wesanderson)

# loading data
wos<-read.table("wos.txt",h=T,sep="\t", encoding = "UFT-8") # literature
dat<-read.table("survey.txt",h=T,sep="\t", encoding = "UFT-8") # Survey

# Renaming responses

for (i in 1:ncol(dat)) dat[which(dat[,i]=="dont know"),i] <- "don't know"
for (i in 1:ncol(dat)) dat[which(dat[,i]=="ecosystem function"),i] <- "function"
for (i in 1:ncol(dat)) dat[which(dat[,i]=="structural measure"),i] <- "structural"

sapply(dat, class)

for (i in 1:ncol(dat)) as.factor(dat[,i])->dat[,i]

sapply(dat, class)

library(RColorBrewer)
coul <- viridis(5)

# Median number of trait terms used by researchers

median(apply(dat[,9:14],1,function(x) length(table(x[which(x!="don't know")]))))

# Mean+SD "I don't know"

length(which(dat[,9:14]=="don't know"))/length(unlist(dat[,9:14]))

mean(apply(dat[,9:14],2,function(x) length(which(x=="don't know")))/nrow(dat))
sd(apply(dat[,9:14],2,function(x) length(which(x=="don't know")))/nrow(dat))
range(apply(dat[,9:14],2,function(x) length(which(x=="don't know")))/nrow(dat))

# Are dominant responses?

even_res<-rep(NA, 7)

a=0

for (i in 8:14) { 
  a=a+1
  #gini_res2[a]<-Gini(dat[,i])
  # Pielou evenness index, 0=uneven distribution, 1=perfectly even distribution
  # https://en.wikipedia.org/wiki/Species_evenness
  even_res[a]<-diversity(table(dat[,i]))/log(specnumber(table(dat[,i])))
  }
#

names(even_res)<-c("Q7 Trait definition", "Q8 Body size", "Q9 Enzymatic activity", "Q10 Chlorophyll-a", "Q11 Plant growth form", 
  "Q12 Nutrient uptake", "Q13 OM decomposition")

write.table(even_res, paste(res_fold,"even_res.txt",sep=""), sep="\t")

res.list<-apply(dat[,9:13],1,table)

dat$num_terms<-unlist(lapply(res.list, function(x) specnumber(x)))

unlist(lapply(res.list, function(x) diversity(table(x))/log(specnumber(x))))

dont_know<-rep(NA, length(c(8:14)))

a=0

for (i in c(8:14)) {
  a=a+1
dont_know[a]<-table(dat[,i])[2]/length(dat[,i])
}

mean(dont_know)

sel.var<-names(dat)[c(8:14)] # selecting biodiversity metrics

r2_order<-res_av<-res_d<-res_d_full<-res<-list() # to store data

a=0

for (i in c(8:14)) {
  
  a=a+1
  
  y=dat[,i]
  
  # Models to test
  mod0<-multinom(y ~ 1, dat)
  mod1<-multinom(y ~ q1_experience, dat)
  mod2<-multinom(y ~ q2_region, dat)
  mod3<-multinom(y ~ q3_ecosystem, dat)
  mod4<-multinom(y ~ q4_field, dat)
  mod5<-multinom(y ~ q5_organisms, dat)
  mod6<-multinom(y ~ q6_functions, dat)

  mod.list<-list(mod0, mod1, mod2, mod3, mod4, mod5, mod6)
  
  res_d_full[[a]]<-mod_d<-model.sel(mod.list, rank = "AICc")

}

res_d_full.df<-data.frame(var=rep(sel.var, unlist(lapply(res_d_full,nrow))), do.call(rbind.data.frame, res_d_full))

# Respondent characteristics explaining variability in responses Q7-Q13

# Response tables aggregated by researcher features
q7_org <- ddply(dat,.(q5_organisms),function(x) round(100*table(x$q7_trait)/nrow(x),1))
q7_field <- ddply(dat,.(q4_field),function(x) round(100*table(x$q7_trait)/nrow(x),1))
q8_field <- ddply(dat,.(q4_field),function(x) round(100*table(x$q8_)/nrow(x),1))
q9_region <- ddply(dat,.(q2_region),function(x) round(100*table(x$q9_)/nrow(x),1))
q9_function <- ddply(dat,.(q6_functions),function(x) round(100*table(x$q9_)/nrow(x),1))
q10_field <- ddply(dat,.(q4_field),function(x) round(100*table(x$q10_)/nrow(x),1))
q11_region <- ddply(dat,.(q2_region),function(x) round(100*table(x$q11_)/nrow(x),1))
q12_field <- ddply(dat,.(q4_field),function(x) round(100*table(x$q12_)/nrow(x),1))
q13_function <- ddply(dat,.(q6_functions),function(x) round(100*table(x$q13_)/nrow(x),1))

# Exporting results
write.table(res_d_full.df, paste(res_fold,"mod_res.txt",sep=""), sep="\t",row.names = F)
write.table(q7_org, paste(res_fold,"q7_org.txt",sep=""), sep="\t",row.names = F)
write.table(q7_field, paste(res_fold,"q7_field.txt",sep=""), sep="\t",row.names = F)
write.table(q8_field, paste(res_fold,"q8_field.txt",sep=""), sep="\t",row.names = F)
write.table(q9_region, paste(res_fold,"q9_region.txt",sep=""), sep="\t",row.names = F)
write.table(q9_function, paste(res_fold,"q9_function.txt",sep=""), sep="\t",row.names = F)
write.table(q10_field, paste(res_fold,"q10_field.txt",sep=""), sep="\t",row.names = F)
write.table(q11_region, paste(res_fold,"q11_region.txt",sep=""), sep="\t",row.names = F)
write.table(q12_field, paste(res_fold,"q12_field.txt",sep=""), sep="\t",row.names = F)
write.table(q13_function, paste(res_fold,"q13_function.txt",sep=""), sep="\t",row.names = F)


  imp<-c(res_d_full.df$weight[which(res_d_full.df[,3]=="+")],
        res_d_full.df$weight[which(res_d_full.df[,4]=="+")],
        res_d_full.df$weight[which(res_d_full.df[,5]=="+")],
        res_d_full.df$weight[which(res_d_full.df[,6]=="+")],
        res_d_full.df$weight[which(res_d_full.df[,7]=="+")],
        res_d_full.df$weight[which(res_d_full.df[,8]=="+")])

  imp.df<-data.frame(imp, var=rep(c("Experience", "Region", "Ecosystem",
                          "Discipline", "Biotic group", "Function"),each=7))
  pie.col=viridis(5)

  # Fig 1. Use of trait-related terms in literature (old)
  
  pdf(file=paste(res_fold,"Fig 1 (old).pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5,height=7)  
  
  par(mfrow=c(2,1),cex.lab=1.3, cex.axis=1.1, mar=c(3,3,3,3))
  
  barplot(as.matrix(wos[,c(2,4,7,8)]),beside=T, col=viridis(3),ylab="Papers",ylim=c(0,12000))
  legend("topright", c("Biological trait", "Functional trait", "Species trait"), fill=viridis(3),inset=.02,box.lty=0)
  barplot(as.matrix(wos[,c(3,5,6,9)]), col=c(pie.col[c(1,2)],"brown"),ylab="Papers",ylim=c(0,3000))
  
  
  dev.off()
  
  # New Fig 1.
  
  pdf(file=paste(res_fold,"Fig 1.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5,height=5)  
  
  par(mfrow=c(1,1),cex.lab=1.1, cex.axis=1, mar=c(5,5,5,5))
  
  ord.field<-order(colSums(wos[,-1]),decreasing = F)+1
  bp<-barplot(as.matrix(wos[,ord.field]), horiz=T, las=1, col=c(pie.col[c(1,5)],"#3CBC75FF"),xlab="Papers", xlim=c(0,17000))
  
  legend("bottomright", c("Biological trait", "Functional trait", "Species trait"), fill=c(pie.col[c(1,5)],"#3CBC75FF"),inset=.02,box.lty=0)
  
  per<-as.numeric(round(100*as.matrix(wos[,ord.field])[2,]/colSums(as.matrix(wos[,ord.field])),0))
  
  text(as.numeric(colSums(as.matrix(wos[,ord.field])))+1000, bp,labels = paste(per,"%",sep=""))
  
  dev.off()
  
  round(100*as.matrix(wos[,ord.field])[2,]/colSums(as.matrix(wos[,ord.field])),0)
  
  # Fig. 2 Respondent characteristics
  
  sort(round(100*table(dat$q1_)/sum(table(dat$q1_)),1),decreasing=T)
  sort(round(100*table(dat$q2_region)/sum(table(dat$q2_region)),1),decreasing=T)
  sort(round(100*table(dat$q3)/sum(table(dat$q3)),1),decreasing=T)
  sort(round(100*table(dat$q4)/sum(table(dat$q4)),1),decreasing=T)
  sort(round(100*table(dat$q5)/sum(table(dat$q5)),1),decreasing=T)
  sort(round(100*table(dat$q6)/sum(table(dat$q6)),1),decreasing=T)
  
  pdf(file=paste(res_fold,"Fig 2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=10.5,height=9)  
  
  par(mfrow=c(3,2),mar=c(2,2,2,2),cex=1.2)
  
  mytable=round(100*table(dat$q1_)/sum(table(dat$q1_)),0)
  lbls <- paste(names(table(dat$q1_)), " (", mytable,"%)", sep="")
  pie(table(dat$q1_experience), labels=lbls, col = rainbow(5), border="black", main="Q1 Experience")
  
  mytable=round(100*table(dat$q2_region)/sum(table(dat$q2_region)),0)
  lbls <- paste(names(table(dat$q2_region)), " (", mytable,"%)", sep="")
  pie(table(dat$q2_region), labels=lbls, col = rainbow(5), border="black", main="Q2 Region")
  
  mytable=round(100*table(dat$q3)/sum(table(dat$q3)),0)
  lbls <- paste(names(table(dat$q3)), " (", mytable,"%)", sep="")
  pie(table(dat$q3), labels=lbls, col = rainbow(5), border="black", main="Q3 Ecosystem")
  
  mytable=round(100*table(dat$q4)/sum(table(dat$q4)),0)
  lbls <- paste(names(table(dat$q4)), " (", mytable,"%)", sep="")
  pie(table(dat$q4), labels=lbls, col = rainbow(5), border="black", main="Q4 Field")
  
  mytable=round(100*table(dat$q5)/sum(table(dat$q5)),0)
  lbls <- paste(names(table(dat$q5)), " (", mytable,"%)", sep="")
  pie(table(dat$q5), labels=lbls, col = rainbow(7), border="black", main="Q5 Organism")
  
  mytable=round(100*table(dat$q6)/sum(table(dat$q6)),0)
  lbls <- paste(names(table(dat$q6)), " (", mytable,"%)", sep="")
  pie(table(dat$q6), labels=lbls, col = rainbow(6), border="black", main="Q6 Functions")
  
  dev.off()
  
  # Fig 3. Responses to Q7 and evenness of responses for each Q
  
  pdf(file=paste(res_fold,"fig 3a.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5.25,height=4)  
  
  par(mfrow=c(1,1),mar=c(4,4,4,4),cex=1.2)
  
  mytable=round(100*table(dat$q7)/sum(table(dat$q7)),0)
  lbls <- paste(names(table(dat$q7)), " (", mytable,"%)", sep="")
  pie(table(dat$q7), labels=lbls, col = pie.col[c(1,2,5,3)], border="black", main="Q7 Trait definition")
  
  dev.off()
  
  pdf(file=paste(res_fold,"fig 3b.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5.25,height=4)  
  
  par(mfrow=c(1,1),mar=c(4,9,4,4),cex=1.2)
  
  barplot(even_res[7:1], horiz = T, las=1, xlab="Evenness", col= magma(7),xlim = c(0,1))
  
  dev.off()
  
  # Fig 4. Percentage of respondents selecting each term to 
  # best represent trait measures 
  
  pdf(file=paste(res_fold,"fig 4.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=10.5,height=9)  
  
  par(mfrow=c(3,2),mar=c(2,2,2,2),cex=1.2)
  
  col.f2<-c(pie.col[c(1,4,2,5)],"#1F9E89FF")
  
  mytable=round(100*table(dat$q8)/sum(table(dat$q8)),1)
  lbls <- paste(names(table(dat$q8)), " (", mytable,"%)", sep="")
  pie(table(dat$q8), labels=lbls, col = col.f2, border="black", main="Q8 Body size")
  
  mytable=round(100*table(dat$q9)/sum(table(dat$q9)),1)
  lbls <- paste(names(table(dat$q9)), " (", mytable,"%)", sep="")
  pie(table(dat$q9), labels=lbls, col=col.f2, border="black", main="Q9 Enzymatic activity")
  
  mytable=round(100*table(dat$q10)/sum(table(dat$q10)),1)
  lbls <- paste(names(table(dat$q10)), " (", mytable,"%)", sep="")
  pie(table(dat$q10), labels=lbls, col=col.f2, border="black", main="Q10 Chlorophyll-a")
  
  mytable=round(100*table(dat$q11)/sum(table(dat$q11)),1)
  lbls <- paste(names(table(dat$q11)), " (", mytable,"%)", sep="")
  pie(table(dat$q11), labels=lbls, col=col.f2, border="black", main="Q11 Plant growth form")
  
  mytable=round(100*table(dat$q12)/sum(table(dat$q12)),1)
  lbls <- paste(names(table(dat$q12)), " (", mytable,"%)", sep="")
  pie(table(dat$q12), labels=lbls, col=col.f2, border="black", main="Q12 Nutrient uptake")
  
  mytable=round(100*table(dat$q13)/sum(table(dat$q13)),1)
  lbls <- paste(names(table(dat$q13)), " (", mytable,"%)", sep="")
  pie(table(dat$q13), labels=lbls, col=col.f2, border="black", main="Q13 OM decomposition")
  
  dev.off()
  
  # Fig5. Importance of respondent characteristics
  
  pdf(file=paste(res_fold,"Fig 5.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=5.5,height=4)  
  
  par(mfrow=c(1,1),mar=c(4,8,4,4),cex.lab=1.3, cex.axis=1.25)
  
  boxplot(imp~var, data=imp.df, horizontal = T, las=1, col=magma(6)[6:1],ylab="",xlab="Variable importance",ylim=c(0,1))
  
  dev.off()
