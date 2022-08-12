library(nortest)

library(AID)

setwd('/Users/LAPT0087/Desktop/ME_CFS/Epigenetics/')

mDNA.data<-read.csv('Hannun_etal_2013/methDNA_data_Hannun_etal_2013.csv')

dim(mDNA.data)

pheno.data<-read.table('Hannun_etal_2013/pheno_data_Hannun_et_al_2013.txt',header=T,sep='\t')

dim(pheno.data)

colnames(pheno.data)

aux<-which(mDNA.data$ProbeID=='cg27193080')

mdata[aux,]

colnames(mDNA.data)[-1]==pheno.data$GEO.ID

output<-vector('numeric',dim(mDNA.data)[1])

x<-pheno.data$Age

for(i in 1:dim(mDNA.data)[1]){
		
	y<-t(mDNA.data[i,-1])
	est<-cor(x,y,method='spearman')
	output[i]<-est
	
}

summary(output)

aux<-which.max(output)

age<-pheno.data$Age

x<-100*t(mDNA.data[aux,-1])

log.odds<-log(x/(100-x))

plot(x,age,las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='methylation level (%)',xlim=c(10,60),main=mDNA.data$ProbeID[aux])

lm.fit1<-lm(age~x)

summary(lm.fit1)

qqnorm(lm.fit1$residuals)

qqline(lm.fit1$residuals)

lillie.test(lm.fit1$residuals)

shapiro.test(lm.fit1$residuals)

plot(log.odds,age,las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='log odds of methylation level (%)',main=mDNA.data$ProbeID[aux])

lm.fit2<-lm(age~log.odds)

summary(lm.fit2)

qqnorm(lm.fit2$residuals)

qqline(lm.fit2$residuals)

lillie.test(lm.fit2$residuals)

shapiro.test(lm.fit2$residuals)


plot(x,log(age),las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='log odds of methylation level (%)',main=mDNA.data$ProbeID[aux])

lm.fit3<-lm(log(age)~x)

summary(lm.fit3)

qqnorm(lm.fit3$residuals)

qqline(lm.fit3$residuals)

lillie.test(lm.fit3$residuals)

shapiro.test(lm.fit3$residuals)

plot(log.odds,log(age),las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='log odds of methylation level (%)',main=mDNA.data$ProbeID[aux])

lm.fit4<-lm(log(age)~log.odds)

summary(lm.fit4)

qqnorm(lm.fit4$residuals)

qqline(lm.fit4$residuals)

lillie.test(lm.fit4$residuals)

shapiro.test(lm.fit4$residuals)

boxcoxlm(x,log(age))

boxcoxlm(log.odds,age)

out<-boxcoxlm(x,age,method='sw')

fit<-lm(out$tf.data~x)

summary(fit)

lillie.test(fit$residuals)

shapiro.test(fit$residuals)

plot(x,out$tf.data,las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='log odds of methylation level (%)',main=mDNA.data$ProbeID[aux])

out2<-boxcoxlm(log.odds,age,method='sw')

plot(log.odds,out2$tf.data,las=1,pch=21,bg='pink',col='black',cex=0.75,xlab='log odds of methylation level (%)',main=mDNA.data$ProbeID[aux])

fit<-lm(out2$tf.data~x)

summary(fit)

lillie.test(fit$residuals)

shapiro.test(fit$residuals)
