# =========================================================
# Function: Multivariable linear regression
# Output: Beta (95% CI) and P value
# =========================================================

generate_multi_linear <- function(tdata,tcov,tcov.class,tcov.name,tphe,tadd=NULL){

	tfom<-paste0(paste0(tphe,"~",paste(c(tcov,tadd),collapse="+")))
    tmod<-lm(formula(tfom),data=tdata)
    tmod_coef<-matrix(coef(summary(tmod))[-1,c(1,2,4)],ncol=3)
    rownames(tmod_coef)<-rownames(coef(summary(tmod)))[-1]
    
	tmod_coef<-data.frame(Beta=tmod_coef[,1],SE=tmod_coef[,2],RR=paste0(sprintf("%1.2f",tmod_coef[,1])," (",sprintf("%1.2f",tmod_coef[,1]-1.96*tmod_coef[,2]),", ",sprintf("%1.2f",tmod_coef[,1]+1.96*tmod_coef[,2]),")"),P=formatC(tmod_coef[,3],format="e",digits=2))
    
    output<-data.frame(rbindlist(lapply(1:length(tcov),function(tcovi){ # "c" = continuous, "f" = categorical
        if (tcov.class[tcovi]=="c"){
            data.frame(Factors=tcov.name[tcovi],Levels="-",Mean_sd="-",n=length(tdata[,tcov[tcovi]][!is.na(tdata[,tcov[tcovi]])]),matrix(unlist(tmod_coef[tcov[tcovi],]),ncol=4))
        } else if (tcov.class[tcovi]=="f"){
            tlev<-levels(factor(tdata[,tcov[tcovi]]))
            data.frame(Factors=c(tcov.name[tcovi],rep("",length(tlev)-1)),Levels=tlev,Mean_sd=tapply(tdata[,tphe],tdata[,tcov[tcovi]],function(x) paste0(sprintf("%1.2f",mean(x)),"±",sprintf("%1.2f",sd(x))))[tlev],n=unclass(table(tdata[,tcov[tcovi]])[tlev]),rbind(c("ref","ref","ref","-"),matrix(unlist(tmod_coef[paste0(tcov[tcovi],tlev[-1]),]),ncol=4)))
        }
    })))
    colnames(output)<-c("Factors","Levels","Mean±sd","n","Beta","SE","BetaCI","P")
	output_final <- subset(output,Levels!="Unknown")
    output_final[,c(1:4,7:8)]
}

generate_multi_linear_digit3 <- function(tdata,tcov,tcov.class,tcov.name,tphe,tadd=NULL){

	tfom<-paste0(paste0(tphe,"~",paste(c(tcov,tadd),collapse="+")))
    tmod<-lm(formula(tfom),data=tdata)
    tmod_coef<-matrix(coef(summary(tmod))[-1,c(1,2,4)],ncol=3)
    rownames(tmod_coef)<-rownames(coef(summary(tmod)))[-1]
    
	tmod_coef<-data.frame(Beta=tmod_coef[,1],SE=tmod_coef[,2],RR=paste0(sprintf("%1.3f",tmod_coef[,1])," (",sprintf("%1.3f",tmod_coef[,1]-1.96*tmod_coef[,2]),", ",sprintf("%1.3f",tmod_coef[,1]+1.96*tmod_coef[,2]),")"),P=formatC(tmod_coef[,3],format="e",digits=2))
    
    output<-data.frame(rbindlist(lapply(1:length(tcov),function(tcovi){
        if (tcov.class[tcovi]=="c"){
            data.frame(Factors=tcov.name[tcovi],Levels="-",Mean_sd="-",n=length(tdata[,tcov[tcovi]][!is.na(tdata[,tcov[tcovi]])]),matrix(unlist(tmod_coef[tcov[tcovi],]),ncol=4))
        } else if (tcov.class[tcovi]=="f"){
            tlev<-levels(factor(tdata[,tcov[tcovi]]))
            data.frame(Factors=c(tcov.name[tcovi],rep("",length(tlev)-1)),Levels=tlev,Mean_sd=tapply(tdata[,tphe],tdata[,tcov[tcovi]],function(x) paste0(sprintf("%1.3f",mean(x)),"±",sprintf("%1.3f",sd(x))))[tlev],n=unclass(table(tdata[,tcov[tcovi]])[tlev]),rbind(c("ref","ref","ref","-"),matrix(unlist(tmod_coef[paste0(tcov[tcovi],tlev[-1]),]),ncol=4)))
        }
    })))
    colnames(output)<-c("Factors","Levels","Mean±sd","n","Beta","SE","BetaCI","P")
	output_final <- subset(output,Levels!="Unknown")
    output_final[,c(1:4,7:8)]
}

# =========================================================
# Function: Multivariable logistic regression
# Output: Odds Ratio (95% CI)
# =========================================================

generate_multi_logit <- function(tdata,tcov,tcov.class,tcov.name,tphe,tadd=NULL){

	tfom<-paste0(paste0(tphe,"~",paste(c(tcov,tadd),collapse="+")))
    tmod<-glm(formula(tfom),data=tdata,fam=binomial)
    tmod_coef<-matrix(coef(summary(tmod))[-1,c(1,2,4)],ncol=3)
    rownames(tmod_coef)<-rownames(coef(summary(tmod)))[-1]
    tmod_coef<-data.frame(Beta=tmod_coef[,1],SE=tmod_coef[,2],RR=paste0(sprintf("%1.2f",exp(tmod_coef[,1]))," (",sprintf("%1.2f",exp(tmod_coef[,1]-1.96*tmod_coef[,2])),", ",sprintf("%1.2f",exp(tmod_coef[,1]+1.96*tmod_coef[,2])),")"),P=formatC(tmod_coef[,3],format="e",digits=2))
	
	output<-data.frame(rbindlist(lapply(1:length(tcov),function(tcovi){
        if (tcov.class[tcovi]=="c"){
            data.frame(Factors=tcov.name[tcovi],Levels="-",Mean_sd="-",n=length(tdata[,tcov[tcovi]][!is.na(tdata[,tcov[tcovi]])]),matrix(unlist(tmod_coef[tcov[tcovi],]),ncol=4))
        } else if (tcov.class[tcovi]=="f"){
            tlev<-levels(factor(tdata[,tcov[tcovi]]))
            data.frame(Factors=c(tcov.name[tcovi],rep("",length(tlev)-1)),Levels=tlev,Mean_sd=tapply(tdata[,tphe],tdata[,tcov[tcovi]],function(x) paste0(sprintf("%1.2f",mean(x)),"±",sprintf("%1.2f",sd(x))))[tlev],n=unclass(table(tdata[,tcov[tcovi]])[tlev]),rbind(c("ref","ref","ref","-"),matrix(unlist(tmod_coef[paste0(tcov[tcovi],tlev[-1]),]),ncol=4)))
        }
    })))
    colnames(output)<-c("Factors","Levels","Mean±sd","n","Beta","SE","RRCI","P")
    output_final <- subset(output,Levels!="Unknown")
    output_final[,c(1:4,7:8)]
}
