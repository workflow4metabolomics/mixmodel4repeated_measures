#######  R functions to perform linear mixed model for repeated measures 
#######  on a multi var dataset using 3 files as used in W4M
##############################################################################################################
lmRepeated2FF <- function(ids, ifixfact, itime, isubject, ivd, ndim, nameVar=colnames(ids)[[ivd]],dffOption, 
                          pvalCutof=0.05, visu , tit = "", least.confounded = FALSE, outlier.limit =3) {
   ### function to perform linear mixed model with 1 Fixed factor + Time + random factor subject
   ### based on lmerTest package providing functions giving the same results as SAS proc mixed
   
   if (!is.numeric(ids[[ivd]]))     {stop("Dependant variable is not numeric")}
   if (!is.factor(ids[[ifixfact]])) {stop("fixed factor is not a factor")}
   if (!is.factor(ids[[itime]]))    {stop("Repeated factor is not a factor")}
   if (!is.factor(ids[[isubject]])) {stop("Random factor is not a factor")}
   # a ce stade, il faudrait pr?voir des tests sur la validit? du plan d'exp?rience
   
   time <- ids[[itime]]
   fixfact <- ids[[ifixfact]]
   subject <- ids[[isubject]]
   vd <- ids[[ivd]]
   
   # argument of the function instead of re re-running ndim <- defColRes(ids,ifixfact,itime)
   # nfp : number of main factors + model infos (REML, varSubject) + normality test 
   nfp <- ndim[1];
   # ncff number of comparison of the fixed factor
   nlff <- ndim[2];  ncff <- ndim[3]
   # nct number of comparison of the time factor
   nlt <- ndim[4] ; nct <- ndim[5]
   # nci number of comparison of the interaction
   nli <- ndim[6];  nci <- ndim[7]
   # number of all lmer results
   nresf <- ncff+nct+nci
   ## initialization of the result vector (1 line)
   res <- data.frame(array(rep(NA,(nfp+2*nresf))))
   colnames(res)[1] <- "resultLM"
   
   ### if at least one subject have data for only 1 time, mixed model is not possible and variable must be skip
   ### after excluding NA, table function is used to seek subjects with only 1 data
   ids <- ids[!is.na(ids[[ivd]]),]
   skip <- length(which(table(ids[[isubject]])==1))
   
   if (skip==0) {
      
      mfl <- lmer( vd ~ time + fixfact + time:fixfact + (1| subject), ids) # lmer remix
      
      # ## NL add   
      # ###  DEPLACE APRES CALCUL PVALUES AJUSTEES ET NE FAIRE QUE SI AU MOINS 1 FACTEUR SIGNIFICATIF
      #  if(visu) diagmflF(mfl, title = tit, least.confounded = least.confounded, outlier.limit = outlier.limit)
      # ## end of NL add
      
      
      rsum <- summary(mfl,ddf = dffOption)
      ## test Shapiro Wilks on the residus of the model 
      rShapiro <- shapiro.test(rsum$residuals)
      raov <- anova(mfl,ddf = dffOption)
      dlsm1  <- difflsmeans(mfl,test.effs=NULL)
      ddlsm1 <- dlsm1$diffs.lsmeans.table
      
      ## writing the results on a single line
      namesFactEstim <- paste("estimate ",rownames(ddlsm1)[c(1:(nct+ncff))],sep="")
      namesFactPval <- paste("pvalue ",rownames(ddlsm1)[c(1:(nct+ncff))],sep="")
      namesInter <- rownames(ddlsm1)[-c(1:(nct+ncff))]
      ncI <- nchar(namesInter)
      namesEstimate <- paste("estimate ",substr(namesInter,15,ncI),sep="")
      namespvalues <- paste("pvalue ",substr(namesInter,15,ncI),sep="")
      namesFactprinc <- c("pval_time","pval_trt","pval_inter")
      
      ### lmer results on 1 vector row
      # pvalue of shapiro Wilks test of the residuals
      res[1,] <- rShapiro$p.value; rownames(res)[1] <- "Shapiro.pvalue.residuals" 
      res[2,] <- rsum$varcor$subject[1] ;rownames(res)[2] <- "Subject.Variance"
      res[3,] <- rsum$devcomp$cmp[7] ; rownames(res)[3] <- "REML"
      ### 3 principal factors pvalues results + shapiro test =>  nfp <- 4
      res[c((nfp-2):nfp),] <- raov[,6]; rownames(res)[c((nfp-2):nfp)] <- namesFactprinc
      
      
      ####################  Residuals diagnostics for significants variables #########################
      ### Il at least 1 factor is significant and visu=TRUE NL graphics add to pdf
       if (length(which(raov[,6]<pvalCutof))>0 & visu == 'yes')  {
          diagmflF(mfl, title = tit, least.confounded = least.confounded, outlier.limit = outlier.limit)
          cat("Signif")
       }
      
      # pvalue of fixed factor comparisons
      res[(nfp+1):(nfp+nct),] <- ddlsm1[c(1:nct),7]
      res[(nfp+nct+1):(nfp+nct+ncff),] <- ddlsm1[(nct+1):(nct+ncff),7]
      rownames(res)[(nfp+1):(nfp+nct+ncff)] <- namesFactPval
      res[(nfp+nct+ncff+1):(nfp+nresf),] <- ddlsm1[(nct+ncff+1):(nresf),7]
      rownames(res)[(nfp+nct+ncff+1):(nfp+nresf)] <- namespvalues
      # Estimate of the difference between levels of factors
      res[(nfp+nresf+1):(nfp+nresf+nct),] <- ddlsm1[c(1:nct),1]
      res[(nfp+nresf+nct+1):(nfp+nresf+nct+ncff),] <- ddlsm1[(nct+1):(nct+ncff),1]
      rownames(res)[(nfp+nresf+1):(nfp+nresf+nct+ncff)] <- namesFactEstim
      res[(nfp+nresf+nct+ncff+1):(nfp+2*nresf),] <- ddlsm1[(nct+ncff+1):(nresf),1]
      rownames(res)[(nfp+nresf+nct+ncff+1):(nfp+2*nresf)] <- namesEstimate
   }
   else
      ## one of the subject has only one time, subject can't be a random variable
      ## A repeated measure could be run instead function lme of package nlme, next version
   {   res[1,] <- NA
      #cat("impossible computing\n")
   
   #  # ## NL add (useless)
   #    if(visu){
   #       grid.arrange(ggplot(data.frame()) + geom_point() + xlim(-1, 1) + ylim(-1, 1)+
   #                    annotate("text", x = 0, y = 0, label = "impossible computing")+
   #                    xlab(NULL) +  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
   #                    ylab(NULL) +  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
   #                    theme(panel.grid.minor = element_blank() ,
   #                          panel.grid.major = element_blank() ,
   #                          panel.background = element_rect(fill = "white"))
   #                 , top = textGrob(tit,gp=gpar(fontsize=40,font=4)))
   #    
   #    }
   # # ## end of NL add
   
   }
   tres <- data.frame(t(res)); rownames(tres)[1] <- nameVar
   return(tres)
}

##############################################################################################################
lmRepeated1FF <- function(ids, ifixfact=0, itime, isubject, ivd, ndim, nameVar=colnames(ids)[[ivd]], 
                          dffOption,pvalCutof=0.05) 
   {
   ### function to perform linear mixed model with factor Time + random factor subject
   ### based on lmerTest package providing functions giving the same results as SAS proc mixed
   
   if (!is.numeric(ids[[ivd]]))     {stop("Dependant variable is not numeric")}
   if (!is.factor(ids[[itime]]))    {stop("Repeated factor is not a factor")}
   if (!is.factor(ids[[isubject]])) {stop("Random factor is not a factor")}
   # a ce stade, il faudrait pr?voir des tests sur la validit? du plan d'exp?rience
   
   time <- ids[[itime]]
   subject <- ids[[isubject]]
   vd <- ids[[ivd]] ## dependant variables (quatitative)
   
   # ndim <- defColRes(ids,0,itime)
   # nfp : nombre de facteurs principaux + model infos + normality test 
   nfp <- ndim[1]
   # nct number of comparison of the time factor
   nlt <- ndim[4] ; nct <- ndim[5]
   # number of all lmer results
   nresf <- nct
   ## initialization of the result vector (1 line)
   res <- data.frame(array(rep(NA,(nfp+2*nresf))))
   colnames(res)[1] <- "resultLM"
   
   ### if at least one subject have data for only 1 time, mixed model is not possible and variable must be skip
   ### after excluding NA, table function is used to seek subjects with only 1 data
   ids <- ids[!is.na(ids[[ivd]]),]
   skip <- length(which(table(ids[[isubject]])==1))
   
   if (skip==0) {
      
      mfl <- lmer( vd ~ time + (1| subject), ids) # lmer remix
      rsum <- summary(mfl,ddf = dffOption)
      ## test Shapiro Wilks on the residus of the model 
      rShapiro <- shapiro.test(rsum$residuals)
      raov <- anova(mfl,ddf = dffOption)
      ## Sum of square : aov$'Sum Sq', Mean square : aov$`Mean Sq`, proba : aov$`Pr(>F)`
      
      ## Test of all differences estimates between levels as SAS proc mixed. 
      ## results are in diffs.lsmeans.table dataframe
      ## test.effs=NULL perform all comparisons 2 a 2 including interaction effect
      dlsm1  <- difflsmeans(mfl,test.effs=NULL)
      ddlsm1 <- dlsm1$diffs.lsmeans.table
      
      ## writing the results on a single line
      namesFactEstim <- paste("estimate ",rownames(ddlsm1)[c(1:(nct))],sep="")
      namesFactPval <- paste("pvalue ",rownames(ddlsm1)[c(1:(nct))],sep="")
      namesFactprinc <- "pval_time"
      
      ### lmer results on 1 vector
      # pvalue of shapiro Wilks test of the residuals
      res[1,] <- rShapiro$p.value; rownames(res)[1] <- "Shapiro.pvalue.residuals" 
      res[2,] <- rsum$varcor$subject[1] ;rownames(res)[2] <- "Subject.Variance"
      res[3,] <- rsum$devcomp$cmp[7] ; rownames(res)[3] <- "REML"
      
      ### principal factor time pvalue results + shapiro test 
      res[nfp,] <- raov[,6]; rownames(res)[nfp] <- namesFactprinc
      # pvalue of fixed factor comparisons
      res[(nfp+1):(nfp+nct),] <- ddlsm1[c(1:nct),7]
      rownames(res)[(nfp+1):(nfp+nct)] <- namesFactPval
      
      # Estimate of the difference between levels of factors
      res[(nfp+nresf+1):(nfp+nresf+nct),] <- ddlsm1[c(1:nct),1]
      rownames(res)[(nfp+nresf+1):(nfp+nresf+nct)] <- namesFactEstim
   }
   else
      ## one of the subject has only one time, subject can't be a random variable
      ## A repeated measure could be run instead function lme of package nlme, next version
   {   res[1,] <- NA
   #cat("traitement impossible\n")
   }
   tres <- data.frame(t(res)); rownames(tres)[1] <- nameVar
   return(tres)
}

##############################################################################################################
defColRes <- function(ids, ifixfact, itime) {
   ## define the size of the result file depending on the numbers of levels of the fixed and time factor.
   ## Numbers of levels define the numbers of comparisons with pvalue and estimate of the difference.
   ## The result file also contains the pvalue of the fixed factor, time factor and interaction
   ## plus Shapiro normality test. This is define by nfp 
   ## subscript of fixed factor=0 means no other fixed factor than "time"
   if (ifixfact>0){
      nfp <- 6 # shapiro+time+fixfact+interaction+ others....
      time <- ids[[itime]]
      fixfact <- ids[[ifixfact]]
      
      cat("\n levels fixfact",levels(fixfact))
      cat("\n levels time",levels(time))
      
      # ncff number of comparisons of the fixed factor (nlff number of levels of fixed factor)
      nlff <- length(levels(fixfact)); ncff <- (nlff*(nlff-1))/2
      # nct number of comparison of the time factor (nlt number of levels of time factor)
      nlt <- length(levels(time)); nct <- (nlt*(nlt-1))/2
      # nci number of comparison of the interaction
      nli <- nlff*nlt;  nci <- (nli*(nli-1))/2
      ndim <- c(NA,NA,NA,NA,NA,NA,NA)
      
      ndim[1] <- nfp   # pvalues of fixed factor, time factor and interaction (3columns) and shapiro test pvalue
      ndim[2] <- nlff  # number of levels of fixed factor
      ndim[3] <- ncff  # number of comparisons (2by2) of the fixed factor
      ndim[4] <- nlt   # number of levels of time factor
      ndim[5] <- nct   # number of comparisons (2by2) of the time factor
      ndim[6] <- nli   # number of levels of interaction
      ndim[7] <- nci   # number of comparisons (2by2) of the interaction
      
   } 
   else {
      nfp <- 4 # shapiro+time
      time <- ids[[itime]]
      # nct number of comparison of the time factor
      nlt <- length(levels(time)); nct <- (nlt*(nlt-1))/2
      ndim <- c(NA,NA,NA,NA,NA,NA,NA)
      
      ndim[1] <- nfp   # pvalues of time factor and shapiro test pvalue
      ndim[4] <- nlt   # number of levels of time factor
      ndim[5] <- nct   # number of comparisons (2by2) of the time factor
   }
   return(ndim)
}

##############################################################################################################
lmixedm <- function(datMN, 
                    samDF,
                    varDF,
                    fixfact, time, subject,
                    logtr = "none", 
                    pvalCutof = 0.05,
                    pvalcorMeth = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[7],
                    dffOption,
                    visu = "no",
                    least.confounded = FALSE,
                    outlier.limit = 3,
					     pdfC
                    )
   {
   sampids <- samDF
   dataMatrix <- datMN
   varids <- varDF
   cat("dff computation method=",dffOption,"\n")
   ### Function running lmer function on a set of variables described in 
   ### 3 different dataframes as used by W4M
   ### results are merge with the metadata variables varids
   ### ifixfact, itime, isubject are subscripts of the dependant variables
   if (fixfact=="none")  ifixfact <-0 else ifixfact <- which(colnames(sampids)==fixfact)
   itime    <- which(colnames(sampids)==time)
   isubject <- which(colnames(sampids)==subject)
   
    #lmmds <- dataMatrix[,-1]
   
   lmmds <- dataMatrix
   if (logtr!="log10" & logtr!="log2") logtr <- "none"
   if (logtr=="log10") lmmds <- log10(lmmds+1)
   if (logtr== "log2") lmmds <- log2(lmmds+1)
   
   #idsamp <- dataMatrix[,1]
   #lmmds <- t(lmmds)
   dslm <- cbind(sampids,lmmds)

   nvar <- ncol(lmmds); firstvar <- ncol(sampids)+1; lastvar <- firstvar+ncol(lmmds)-1

   dslm[[ifixfact]] <- factor(dslm[[ifixfact]])
   dslm[[itime]]    <- factor(dslm[[itime]])  
   dslm[[isubject]] <- factor(dslm[[isubject]])
   ## call defColres to define the numbers of test and so the number of columns of results
   ## depends on whether or not there is a fixed factor with time. If only time factor ifixfact=0
   if (ifixfact>0) {
      ndim <- defColRes(dslm[,c(ifixfact,itime)],ifixfact=1,itime=2)
      nColRes <- ndim[1]+(2*(ndim[3]+ndim[5]+ndim[7]))
      firstpval <- ndim[1]-2
      lastpval <- ndim[1]+ndim[3]+ndim[5]+ndim[7]
  } else 
   {
      ndim <- defColRes(dslm[,itime],ifixfact=0,itime=1) 
      nColRes <- ndim[1]+(2*(ndim[5]))
      firstpval <- ndim[1]
      lastpval <- ndim[1]+ndim[5]
   }
   ## initialisation of the  result file 
   resLM <- data.frame(array(rep(NA,nvar*nColRes),dim=c(nvar,nColRes))); colnames(resLM)[1] <- "res"

   ###############  test ecriture dans pdf
   if(visu == "yes") {
      pdf(pdfC, onefile=TRUE, height = 15, width = 30)
      par(mfrow=c(1,3))
   }
   ###############  fin test ecriture dans pdf
   ## pour test : lastvar <- 15
   for (i in firstvar:lastvar) {
    
     ## NL modif
     cat("\n[",colnames(dslm)[i],"] ")
     ## end of NL modif
     
     subds <- dslm[,c(ifixfact,itime,isubject,i)]

      ## NL modif
     tryCatch({
      if (ifixfact>0)
        reslmer <- lmRepeated2FF(subds,ifixfact=1,2,3, ivd=4, ndim=ndim, visu = visu,  tit = varids[i-firstvar+1,1], pvalCutof,dffOption,
                                 least.confounded = least.confounded, outlier.limit = outlier.limit) else 
        reslmer <- lmRepeated1FF(subds,ifixfact=0,1,2, ivd=3, ndim=ndim, pvalCutof,dffOption)
      ## end of NL modif
      resLM[i-firstvar+1,] <- reslmer
     }, error=function(e){cat("ERROR : ",conditionMessage(e), "\n")})
      if (i==firstvar) {colnames(resLM) <- colnames(reslmer)}
   }
   
   
   ## pvalue correction with p.adjust library multtest
   ## Possible methods of pvalue correction
   AdjustMeth <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr","none")
   if (length(which(pvalcorMeth == AdjustMeth))==0) pvalcorMeth <- "none"

   if (pvalcorMeth !="none") {
      for (k in firstpval:lastpval){
         resLM[[k]]=p.adjust(resLM[[k]], method=pvalcorMeth, n=dim(resLM[k])[[1]])
      }
   }

   ## for each variables, set pvalues and estimates to NA when pvalue of factor > Cut-Off value define by user (pvalCutof) 
   if (ifixfact>0) {
      ## time effect
      resLM[which(resLM[,firstpval]> pvalCutof),c((lastpval+1):(lastpval+ndim[5]))] <- NA
      resLM[which(resLM[,firstpval]> pvalCutof),c((ndim[1]+1):(ndim[1]+ndim[5]))] <- NA
   ## treatment effect
      resLM[which(resLM[,firstpval+1]> pvalCutof),c((lastpval+ndim[5]+1):(lastpval+ndim[5]+ndim[3]))] <- NA
      resLM[which(resLM[,firstpval+1]> pvalCutof),c((ndim[1]+ndim[5]+1):(ndim[1]+ndim[5]+ndim[3]))] <- NA
   ## interaction effect
      resLM[which(resLM[,firstpval+2]> pvalCutof),c((lastpval+ndim[5]+ndim[3]+1):(lastpval+ndim[5]+ndim[3]+ndim[7]))] <- NA
      resLM[which(resLM[,firstpval+2]> pvalCutof),c((ndim[1]+ndim[5]+ndim[3]+1):(ndim[1]+ndim[5]+ndim[3]+ndim[7]))] <- NA
   } else {
      ## time effect only
      resLM[which(resLM[,firstpval]> pvalCutof),c((lastpval+1):(lastpval+ndim[5]))] <- NA
      resLM[which(resLM[,firstpval]> pvalCutof),c((firstpval+1):(firstpval+ndim[5]))] <- NA
   }
   
   
   ## NL add
   if(visu == "yes") dev.off()
   ## end of NL add
   
   
   
   ## return result file with pvalues ans estimates
   resLM <- cbind(varids,resLM)
   return(resLM)
}

