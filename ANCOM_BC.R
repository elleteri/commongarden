# ANCOM with Chem data ####
## Create ANCOM Function ####

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL,check.names=FALSE)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T & repeated==F ){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  if(adjusted==T & repeated==T ){
    adjusted.formula = paste0(base.formula," + ", adj.formula," | ", repeat.var)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  ######################## Function 
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  ######################## Function 
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}


#####End of functions

#FULL GC ANCOM ####
##Remove entries with insufficient areas
OCG_GC_subset[OCG_GC_subset < 10] <- 0
OCG_GC_subset <- OCG_GC_subset[rowSums(OCG_GC_subset) > 0,] # each observation needs at least 10

summary(rowSums(OCG_GC_subset)) #179
summary(colSums(OCG_GC_subset)) #2500

OCG_GC_subset <- OCG_GC_subset[,colSums(OCG_GC_subset) > 10] # each sample needs at least 10 

md.OCG.GC.a <- subset(md.OCG.GC, row.names(md.OCG.GC) %in% row.names(OCG_GC_subset)) 

OCG_GC_subset.t <- t(OCG_GC_subset) 
OCG_GC_subset_t <- OCG_GC_subset.t[, colnames(OCG_GC_subset.t) %in% row.names(md.OCG.GC.a), drop = FALSE]
OCG_GC_subset <- t(OCG_GC_subset_t)
OCG_GC_subset <- as.data.frame(OCG_GC_subset) 

OCG_GC_subset[OCG_GC_subset < 10] <- 0  # repeat cleaning after trimming
OCG_GC_subset <- OCG_GC_subset[rowSums(OCG_GC_subset) > 0,] # each observation needs at least 10

summary(rowSums(OCG_GC_subset)) #179
summary(colSums(OCG_GC_subset)) #2500

OCG_GC_subset <- OCG_GC_subset[,colSums(OCG_GC_subset) > 10]

#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
md.OCG.GC_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC), md.OCG.GC)
OCG_GC_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_subset), OCG_GC_subset, check.names = F)
row.names(OCG_GC_sbst) == row.names(md.OCG.GC_sbst) #TRUE

### ANCOM SUBSPECIES GC ####
ANCOM_subspecies_GC <- ANCOM.main(OCG_GC_sbst,md.OCG.GC_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant ASVs
sigGC_subspecies <- subset(ANCOM_subspecies_GC$W.taxa, ANCOM_subspecies_GC$W.taxa$W_stat > 0)[,1]
sigGC_subspecies <- as.data.frame(sigGC_subspecies) 
row.names(sigGC_subspecies) <- sigGC_subspecies[54:1,1] 

sigGC_subspecies[,1] <- c(54:1) 
sigGC_subspecies_t <- t(sigGC_subspecies) 
sigGC_subspecies_t <- as.data.frame(sigGC_subspecies_t)
colnames(sigGC_subspecies_t) <- as.character(colnames(sigGC_subspecies_t))
print(colnames(sigGC_subspecies_t))   
sigGC_subspecies_t <- sigGC_subspecies_t[,order(colnames(sigGC_subspecies_t))]
rownames(sigGC_subspecies_t) <- c("sig_rank")

GC.OCG_sig_subspecies <-  t(subset(t(OCG_GC_sbst), colnames(OCG_GC_sbst) %in% row.names(sigGC_subspecies)))
GC.OCG_sig_subspecies <- GC.OCG_sig_subspecies[,order(colnames(GC.OCG_sig_subspecies))]
colnames(sigGC_subspecies_t) == colnames(GC.OCG_sig_subspecies) #TRUE

GC.OCG_sig_subspecies <- rbind(GC.OCG_sig_subspecies, sigGC_subspecies_t)
GC.OCG_sig_subspecies_t <- as.data.frame(t(GC.OCG_sig_subspecies))
GC.OCG_sig_subspecies_t$sig_rank <- as.numeric(GC.OCG_sig_subspecies_t$sig_rank) 
GC.OCG_sig_subspecies_t <- GC.OCG_sig_subspecies_t[order(GC.OCG_sig_subspecies_t$sig_rank),] 
GC.OCG_sig_subspecies_t <- subset(GC.OCG_sig_subspecies_t, select=-c(sig_rank))
GC.OCG_sig_subspecies_t <- as.data.frame(t(GC.OCG_sig_subspecies_t))

#Build objects for plotting
sbsplotGC <- data.frame(GC.OCG_sig_subspecies_t[,1:10], "subspecies" = md.OCG.GC_sbst$Subspecies, check.names = FALSE)
sbsplotGC[,1:10] <- lapply(sbsplotGC[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC[,1:10] <- log(sbsplotGC[,1:10]+1)
sbsplotGC_sub <- data.frame(sample=rownames(sbsplotGC),sbsplotGC, check.names = F)
sbsplotGC_sublong <- melt(sbsplotGC_sub)

ANCOM_subspecies_GC$W.taxa

sbsplotGC_sublong$subspecies <-  factor(sbsplotGC_sublong$subspecies, levels = c("T", "V", "W"))

##### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS ALL GC ####
sp1 <- ggplot(sbsplotGC_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("Full GC ANCOM for Subspecies")+ 
  theme_classic()

# 2012 GC ANCOM ####
OCG_GC_2012_subset[OCG_GC_2012_subset < 10] <- 0
OCG_GC_2012_subset <- OCG_GC_2012_subset[rowSums(OCG_GC_2012_subset) > 0,] 

summary(rowSums(OCG_GC_2012_subset)) #7962
summary(colSums(OCG_GC_2012_subset)) #2358

OCG_GC_2012_subset <- OCG_GC_2012_subset[,colSums(OCG_GC_2012_subset) > 10] 
md.OCG.GC.2012.a <- subset(md.OCG.GC.2012, row.names(md.OCG.GC.2012) %in% row.names(OCG_GC_2012_subset)) 

OCG_GC_2012_subset.t <- t(OCG_GC_2012_subset) 
OCG_GC_2012_subset_t <- OCG_GC_2012_subset.t[, colnames(OCG_GC_2012_subset.t) %in% row.names(md.OCG.GC.2012.a), drop = FALSE]
OCG_GC_2012_subset <- t(OCG_GC_2012_subset_t)
OCG_GC_2012_subset <- as.data.frame(OCG_GC_2012_subset) 

OCG_GC_2012_subset[OCG_GC_2012_subset < 10] <- 0  
OCG_GC_2012_subset <- OCG_GC_2012_subset[rowSums(OCG_GC_2012_subset) > 0,] 

summary(rowSums(OCG_GC_2012_subset)) #7962
summary(colSums(OCG_GC_2012_subset)) #2358

OCG_GC_2012_subset <- OCG_GC_2012_subset[,colSums(OCG_GC_2012_subset) > 10]

#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
md.OCG.GC.2012_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.2012.a), md.OCG.GC.2012.a)
OCG_GC_2012_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_2012_subset), OCG_GC_2012_subset, check.names = F)
row.names(OCG_GC_2012_sbst) == row.names(md.OCG.GC.2012_sbst) #TRUE

## SUBSPECIES 2012 GC ####
ANCOM_subspecies_GC.2012 <- ANCOM.main(OCG_GC_2012_sbst,md.OCG.GC.2012_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant ASVs
sigGC.2012_subspecies <- subset(ANCOM_subspecies_GC.2012$W.taxa, ANCOM_subspecies_GC.2012$W.taxa$W_stat > 0)[,1]
sigGC.2012_subspecies <- as.data.frame(sigGC.2012_subspecies) 
row.names(sigGC.2012_subspecies) <- sigGC.2012_subspecies[47:1,1] 

sigGC.2012_subspecies[,1] <- c(47:1) 
sigGC.2012_subspecies_t <- t(sigGC.2012_subspecies) 
sigGC.2012_subspecies_t <- as.data.frame(sigGC.2012_subspecies_t)
colnames(sigGC.2012_subspecies_t) <- as.character(colnames(sigGC.2012_subspecies_t))
print(colnames(sigGC.2012_subspecies_t))

sigGC.2012_subspecies_t <- sigGC.2012_subspecies_t[,order(colnames(sigGC.2012_subspecies_t))]
rownames(sigGC.2012_subspecies_t) <- c("sig_rank")

GC.2012.OCG_sigsbst_subspecies <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_subspecies)))

GC.2012.OCG_sigsbst_subspecies <- GC.2012.OCG_sigsbst_subspecies[,order(colnames(GC.2012.OCG_sigsbst_subspecies))]
colnames(sigGC.2012_subspecies_t) == colnames(GC.2012.OCG_sigsbst_subspecies) #sanity check:TRUE

GC.2012.OCG_sigsbst_subspecies <- rbind(GC.2012.OCG_sigsbst_subspecies, sigGC.2012_subspecies_t)
GC.2012.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspecies))
GC.2012.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_subspecies_t$sig_rank) 
GC.2012.OCG_sigsbst_subspecies_t <- GC.2012.OCG_sigsbst_subspecies_t[order(GC.2012.OCG_sigsbst_subspecies_t$sig_rank),] 
GC.2012.OCG_sigsbst_subspecies_t <- subset(GC.2012.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.GC.2012_sbst$Subspecies, check.names = FALSE)
sbsplotGC2012[,1:10] <- lapply(sbsplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC2012[,1:10] <- log(sbsplotGC2012[,1:10]+1)
sbsplotGC2012_sub <- data.frame(sample=rownames(sbsplotGC2012),sbsplotGC2012, check.names = F)
sbsplotGC2012_sublong <- melt(sbsplotGC2012_sub)

ANCOM_subspecies_GC.2012$W.taxa

sbsplotGC2012_sublong$subspecies <-  factor(sbsplotGC2012_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS 2012 GC ####
sp2 <- ggplot(sbsplotGC2012_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Subspecies")+
  labs(color = "Compounds") +
  theme_classic()

## PLOIDY 2012 GC ####
#Run ANCOM, specify variable
ANCOM_ploidy_GC <- ANCOM.main(OCG_GC_2012_sbst,md.OCG.GC.2012_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigGC.2012_ploidy <- subset(ANCOM_ploidy_GC$W.taxa, ANCOM_ploidy_GC$W.taxa$W_stat > 0)[,1]
sigGC.2012_ploidy <- as.data.frame(sigGC_ploidy)
row.names(sigGC.2012_ploidy) <- sigGC.2012_ploidy[47:1,1] 

sigGC.2012_ploidy[,1] <- c(47:1) 
sigGC.2012_ploidy_t <- t(sigGC.2012_ploidy) 
sigGC.2012_ploidy_t <- as.data.frame(sigGC.2012_ploidy_t)
colnames(sigGC.2012_ploidy_t) <- as.character(colnames(sigGC.2012_ploidy_t))
print(colnames(sigGC.2012_ploidy_t))

sigGC.2012_ploidy_t <- sigGC.2012_ploidy_t[,order(colnames(sigGC.2012_ploidy_t))]
rownames(sigGC.2012_ploidy_t) <- c("sig_rank")

GC.2012.OCG_sigsbst_ploidy <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_ploidy)))

GC.2012.OCG_sigsbst_ploidy <- GC.2012.OCG_sigsbst_ploidy[,order(colnames(GC.2012.OCG_sigsbst_ploidy))]
colnames(sigGC.2012_subspecies_t) == colnames(GC.2012.OCG_sigsbst_ploidy) #sanity check:TRUE

GC.2012.OCG_sigsbst_ploidy <- rbind(GC.2012.OCG_sigsbst_ploidy, sigGC.2012_ploidy_t)
GC.2012.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_ploidy))
GC.2012.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_ploidy_t$sig_rank) 
GC.2012.OCG_sigsbst_ploidy_t <- GC.2012.OCG_sigsbst_ploidy_t[order(GC.2012.OCG_sigsbst_ploidy_t$sig_rank),] 
GC.2012.OCG_sigsbst_ploidy_t <- subset(GC.2012.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.GC.2012_sbst$Ploidy, check.names = FALSE)
plplotGC2012[,1:10] <- lapply(plplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
plplotGC2012[,1:10] <- log(plplotGC2012[,1:10]+1)
plplotGC2012_sub <- data.frame(sample=rownames(plplotGC2012),plplotGC2012, check.names = F)
plplotGC2012_sublong <- melt(plplotGC2012_sub)

ANCOM_subspecies_GC.2012$W.taxa

plplotGC2012_sublong$ploidy<-  factor(plplotGC2012_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS 2012 GC ####
pl1 <- ggplot(sbsplotGC2012.ploidy_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


## SUBSPECIES PLOIDY 2012 GC ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.GC12 <- ANCOM.main(OCG_GC_2012_sbst, md.OCG.GC.2012_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2012_subsp_ploidy <- subset(ANCOM_subsp_ploidy.GC12$W.taxa, ANCOM_subsp_ploidy.GC12$W.taxa$W_stat > 0)[,1]
sigGC.2012_subsp_ploidy <- as.data.frame(sigGC.2012_subsp_ploidy)
row.names(sigGC.2012_subsp_ploidy) <- sigGC.2012_subsp_ploidy[47:1,1]
sigGC.2012_subsp_ploidy[,1] <- c(47:1)

sigGC.2012_subsp_ploidy_t <- t(sigGC.2012_subsp_ploidy)
sigGC.2012_subsp_ploidy_t <- as.data.frame(sigGC.2012_subsp_ploidy_t)
colnames(sigGC.2012_subsp_ploidy_t) <- as.character(colnames(sigGC.2012_subsp_ploidy_t))
print(colnames(sigGC.2012_subsp_ploidy_t))
sigGC.2012_subsp_ploidy_t <- sigGC.2012_subsp_ploidy_t[,order(colnames(sigGC.2012_subsp_ploidy_t))]
rownames(sigGC.2012_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

GC.2012.OCG_sigsbst_subspploidy <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_subsp_ploidy)))

GC.2012.OCG_sigsbst_subspploidy <- GC.2012.OCG_sigsbst_subspploidy[,order(colnames(GC.2012.OCG_sigsbst_subspploidy))]
colnames(sigGC.2012_subsp_ploidy_t) == colnames(GC.2012.OCG_sigsbst_subspploidy) #sanity check:TRUE

GC.2012.OCG_sigsbst_subspploidy <- rbind(GC.2012.OCG_sigsbst_subspploidy, sigGC.2012_subsp_ploidy_t)
GC.2012.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspploidy))
GC.2012.OCG_sigsbst_subspploidy_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_subspploidy_t$sig_rank) 
GC.2012.OCG_sigsbst_subspploidy_t <- GC.2012.OCG_sigsbst_subspploidy_t[order(GC.2012.OCG_sigsbst_subspploidy_t$sig_rank),] 
GC.2012.OCG_sigsbst_subspploidy_t <- subset(GC.2012.OCG_sigsbst_subspploidy_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspploidy_t))

#Build objects for plotting
sbsp.plplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.GC.2012_sbst$Subsp_ploidy, check.names = FALSE)
sbsp.plplotGC2012[,1:10] <- lapply(sbsp.plplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
sbsp.plplotGC2012[,1:10] <- log(sbsp.plplotGC2012[,1:10]+1)
sbsp.plplotGC2012_sub <- data.frame(sample=rownames(sbsp.plplotGC2012),sbsp.plplotGC2012, check.names = F)
sbsp.plplotGC2012_sublong <- melt(sbsp.plplotGC2012_sub)

ANCOM_subsp_ploidy.GC12$W.taxa

sbsp.plplotGC2012$subsp_ploidy <-  factor(sbsp.plplotGC2012$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

#### FIGURE FOR SUBSPECIES PLOIDY ####
spl1 <- ggplot(sbsp.plplotGC2012_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


# 2021 GC ANCOM ####
OCG_GC_2021_subset[OCG_GC_2021_subset < 10] <- 0
OCG_GC_2021_subset <- OCG_GC_2021_subset[rowSums(OCG_GC_2021_subset) > 0,] 

summary(rowSums(OCG_GC_2021_subset)) #179
summary(colSums(OCG_GC_2021_subset)) #842.5

OCG_GC_2021_subset <- OCG_GC_2021_subset[,colSums(OCG_GC_2021_subset) > 10] 
md.OCG.GC.2021.a <- subset(md.OCG.GC.2021, row.names(md.OCG.GC.2021) %in% row.names(OCG_GC_2021_subset)) 

OCG_GC_2021_subset.t <- t(OCG_GC_2021_subset) 
OCG_GC_2021_subset_t <- OCG_GC_2021_subset.t[, colnames(OCG_GC_2021_subset.t) %in% row.names(md.OCG.GC.2021.a), drop = FALSE]
OCG_GC_2021_subset <- t(OCG_GC_2021_subset_t)
OCG_GC_2021_subset <- as.data.frame(OCG_GC_2021_subset) 

OCG_GC_2021_subset[OCG_GC_2021_subset < 10] <- 0  
OCG_GC_2021_subset <- OCG_GC_2021_subset[rowSums(OCG_GC_2021_subset) > 0,] 

summary(rowSums(OCG_GC_2021_subset)) #179
summary(colSums(OCG_GC_2021_subset)) #842.5

OCG_GC_2021_subset <- OCG_GC_2021_subset[,colSums(OCG_GC_2021_subset) > 10]

md.OCG.GC.2021_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.2021.a), md.OCG.GC.2021.a)
OCG_GC_2021_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_2021_subset), OCG_GC_2021_subset, check.names = F)
row.names(OCG_GC_2021_sbst) == row.names(md.OCG.GC.2021_sbst) #TRUE

## SUBSPECIES 2021 GC ####
ANCOM_subspecies_GC.2021 <- ANCOM.main(OCG_GC_2021_sbst,md.OCG.GC.2021_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant compounds
sigGC.2021_subspecies <- subset(ANCOM_subspecies_GC.2021$W.taxa, ANCOM_subspecies_GC.2021$W.taxa$W_stat > 0)[,1]
sigGC.2021_subspecies <- as.data.frame(sigGC.2021_subspecies) 
row.names(sigGC.2021_subspecies) <- sigGC.2021_subspecies[29:1,1] 

sigGC.2021_subspecies[,1] <- c(29:1) 
sigGC.2021_subspecies_t <- t(sigGC.2021_subspecies) 
sigGC.2021_subspecies_t <- as.data.frame(sigGC.2021_subspecies_t)
colnames(sigGC.2021_subspecies_t) <- as.character(colnames(sigGC.2021_subspecies_t))
print(colnames(sigGC.2021_subspecies_t))

sigGC.2021_subspecies_t <- sigGC.2021_subspecies_t[,order(colnames(sigGC.2021_subspecies_t))]
rownames(sigGC.2021_subspecies_t) <- c("sig_rank")

GC.2021.OCG_sigsbst_subspecies <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_subspecies)))

GC.2021.OCG_sigsbst_subspecies <- GC.2021.OCG_sigsbst_subspecies[,order(colnames(GC.2021.OCG_sigsbst_subspecies))]
colnames(sigGC.2021_subspecies_t) == colnames(GC.2021.OCG_sigsbst_subspecies) #sanity check:TRUE

GC.2021.OCG_sigsbst_subspecies <- rbind(GC.2021.OCG_sigsbst_subspecies, sigGC.2021_subspecies_t)
GC.2021.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspecies))
GC.2021.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_subspecies_t$sig_rank) 
GC.2021.OCG_sigsbst_subspecies_t <- GC.2021.OCG_sigsbst_subspecies_t[order(GC.2021.OCG_sigsbst_subspecies_t$sig_rank),] 
GC.2021.OCG_sigsbst_subspecies_t <- subset(GC.2021.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotGC2021 <- data.frame(GC.2021.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.GC.2021_sbst$Subspecies, check.names = FALSE)
sbsplotGC2021[,1:10] <- lapply(sbsplotGC2021[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC2021[,1:10] <- log(sbsplotGC2021[,1:10]+1)
sbsplotGC2021_sub <- data.frame(sample=rownames(sbsplotGC2021),sbsplotGC2021, check.names = F)
sbsplotGC2021_sublong <- melt(sbsplotGC2021_sub)

ANCOM_subspecies_GC.2021$W.taxa

sbsplotGC2021_sublong$subspecies <-  factor(sbsplotGC2021_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS 2021 GC ####
sp3 <- ggplot(sbsplotGC2021_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2021 GC ANCOM for Subspecies")+theme_classic()

## PLOIDY 2021 GC ####
#Run ANCOM, specify variable
ANCOM_ploidy_GC.21 <- ANCOM.main(OCG_GC_2021_sbst,md.OCG.GC.2021_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2021_ploidy <- subset(ANCOM_ploidy_GC.21$W.taxa, ANCOM_ploidy_GC.21$W.taxa$W_stat > 0)[,1]
sigGC.2021_ploidy <- as.data.frame(sigGC.2021_ploidy)
row.names(sigGC.2021_ploidy) <- sigGC.2021_ploidy[37:1,1] 

sigGC.2021_ploidy[,1] <- c(37:1) 
sigGC.2021_ploidy_t <- t(sigGC.2021_ploidy) 
sigGC.2021_ploidy_t <- as.data.frame(sigGC.2021_ploidy_t)
colnames(sigGC.2021_ploidy_t) <- as.character(colnames(sigGC.2021_ploidy_t))
print(colnames(sigGC.2021_ploidy_t))

sigGC.2021_ploidy_t <- sigGC.2021_ploidy_t[,order(colnames(sigGC.2021_ploidy_t))]
rownames(sigGC.2021_ploidy_t) <- c("sig_rank")

GC.2021.OCG_sigsbst_ploidy <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_ploidy)))

GC.2021.OCG_sigsbst_ploidy <- GC.2021.OCG_sigsbst_ploidy[,order(colnames(GC.2021.OCG_sigsbst_ploidy))]
colnames(sigGC.2021_ploidy_t) == colnames(GC.2021.OCG_sigsbst_ploidy) #sanity check:TRUE

GC.2021.OCG_sigsbst_ploidy <- rbind(GC.2021.OCG_sigsbst_ploidy, sigGC.2021_ploidy_t)
GC.2021.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_ploidy))
GC.2021.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_ploidy_t$sig_rank) 
GC.2021.OCG_sigsbst_ploidy_t <- GC.2021.OCG_sigsbst_ploidy_t[order(GC.2021.OCG_sigsbst_ploidy_t$sig_rank),] 
GC.2021.OCG_sigsbst_ploidy_t <- subset(GC.2021.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotGC2021.ploidy <- data.frame(GC.2021.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.GC.2021_sbst$Ploidy, check.names = FALSE)
plplotGC2021.ploidy[,1:10] <- lapply(plplotGC2021.ploidy[,1:10], function(x) as.numeric(as.character(x)))
plplotGC2021.ploidy[,1:10] <- log(plplotGC2021.ploidy[,1:10]+1)
plplotGC2021.ploidy_sub <- data.frame(sample=rownames(plplotGC2021.ploidy),plplotGC2021.ploidy, check.names = F)
plplotGC2021.ploidy_sublong <- melt(plplotGC2021.ploidy_sub)

ANCOM_subspecies_GC.2021$W.taxa

plplotGC2021.ploidy_sublong$ploidy<-  factor(plplotGC2021.ploidy_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS 2021 GC ####
pl2 <- ggplot(sbsplotGC2021.ploidy_sublong, aes(y = value, x = ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Ploidy") + 
  ggtitle("2021 GC ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

## SIGNIFICANCE BY SUBSPECIES PLOIDY 2021 GC ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.GC21 <- ANCOM.main(OCG_GC_2021_sbst, md.OCG.GC.2021_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2021_subsp_ploidy <- subset(ANCOM_subsp_ploidy.GC21$W.taxa, ANCOM_subsp_ploidy.GC21$W.taxa$W_stat > 0)[,1]
sigGC.2021_subsp_ploidy <- as.data.frame(sigGC.2021_subsp_ploidy)
row.names(sigGC.2021_subsp_ploidy) <- sigGC.2021_subsp_ploidy[36:1,1]
sigGC.2021_subsp_ploidy[,1] <- c(36:1)

sigGC.2021_subsp_ploidy_t <- t(sigGC.2021_subsp_ploidy)
sigGC.2021_subsp_ploidy_t <- as.data.frame(sigGC.2021_subsp_ploidy_t)
colnames(sigGC.2021_subsp_ploidy_t) <- as.character(colnames(sigGC.2021_subsp_ploidy_t))
print(colnames(sigGC.2021_subsp_ploidy_t))
sigGC.2021_subsp_ploidy_t <- sigGC.2021_subsp_ploidy_t[,order(colnames(sigGC.2021_subsp_ploidy_t))]
rownames(sigGC.2021_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

GC.2021.OCG_sigsbst_subspploidy <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_subsp_ploidy)))

GC.2021.OCG_sigsbst_subspploidy <- GC.2021.OCG_sigsbst_subspploidy[,order(colnames(GC.2021.OCG_sigsbst_subspploidy))]
colnames(sigGC.2021_subsp_ploidy_t) == colnames(GC.2021.OCG_sigsbst_subspploidy) #sanity check:TRUE

GC.2021.OCG_sigsbst_subspploidy <- rbind(GC.2021.OCG_sigsbst_subspploidy, sigGC.2021_subsp_ploidy_t)
GC.2021.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspploidy))
GC.2021.OCG_sigsbst_subspploidy_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_subspploidy_t$sig_rank) 
GC.2021.OCG_sigsbst_subspploidy_t <- GC.2021.OCG_sigsbst_subspploidy_t[order(GC.2021.OCG_sigsbst_subspploidy_t$sig_rank),] 
GC.2021.OCG_sigsbst_subspploidy_t <- subset(GC.2021.OCG_sigsbst_subspploidy_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspploidy_t))

#Build objects for plotting
sbsp.plplotGC2021 <- data.frame(GC.2021.OCG_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.GC.2021_sbst$Subsp_ploidy, check.names = FALSE)
sbsp.plplotGC2021[,1:10] <- lapply(sbsp.plplotGC2021[,1:10], function(x) as.numeric(as.character(x)))
sbsp.plplotGC2021[,1:10] <- log(sbsp.plplotGC2021[,1:10]+1)
sbsp.plplotGC2021_sub <- data.frame(sample=rownames(sbsp.plplotGC2021),sbsp.plplotGC2021, check.names = F)
sbsp.plplotGC2021_sublong <- melt(sbsp.plplotGC2021_sub)

ANCOM_subsp_ploidy.GC21$W.taxa

sbsp.plplotGC2021$subsp_ploidy <-  factor(sbsp.plplotGC2021$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

#### FIGURE FOR SUBSPECIES PLOIDY 2021 GC ####
spl2 <- ggplot(sbsp.plplotGC2021_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2021 GC ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

# Full GC ANCOM ####
OCG_GC_subset[OCG_GC_subset < 10] <- 0
OCG_GC_subset <- OCG_GC_subset[rowSums(OCG_GC_subset) > 0,] 

summary(rowSums(OCG_GC_subset)) #179
summary(colSums(OCG_GC_subset)) #2500

OCG_GC_subset <- OCG_GC_subset[,colSums(OCG_GC_subset) > 10] 
md.OCG.GC.a <- subset(md.OCG.GC, row.names(md.OCG.GC) %in% row.names(OCG_GC_subset)) 

OCG_GC_subset.t <- t(OCG_GC_subset) 
OCG_GC_subset_t <- OCG_GC_subset.t[, colnames(OCG_GC_subset.t) %in% row.names(md.OCG.GC.a), drop = FALSE]
OCG_GC_subset <- t(OCG_GC_subset_t)
OCG_GC_subset <- as.data.frame(OCG_GC_subset) 

OCG_GC_subset[OCG_GC_subset < 10] <- 0  
OCG_GC_subset <- OCG_GC_subset[rowSums(OCG_GC_subset) > 0,] 

summary(rowSums(OCG_GC_subset)) #179
summary(colSums(OCG_GC_subset)) #2500

OCG_GC_subset <- OCG_GC_subset[,colSums(OCG_GC_subset) > 10]

md.OCG.GC_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.a), md.OCG.GC.a)
OCG_GC_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_subset), OCG_GC_subset, check.names = F)
row.names(OCG_GC_sbst) == row.names(md.OCG.GC_sbst) #TRUE

## YEAR Full GC ####
#Run ANCOM, specify variable
ANCOM_yr.GC <- ANCOM.main(OCG_GC_sbst,md.OCG.GC_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)

sigGC_yr <- subset(ANCOM_yr.GC$W.taxa, ANCOM_yr.GC$W.taxa$W_stat > 0)[,1]
sigGC_yr <- as.data.frame(sigGC_yr)
row.names(sigGC_yr) <- sigGC_yr[54:1,1]
sigGC_yr[,1] <- c(54:1)

sigGC_yr_t <- t(sigGC_yr)
sigGC_yr_t <- as.data.frame(sigGC_yr_t)
colnames(sigGC_yr_t) <- as.character(colnames(sigGC_yr_t))
print(colnames(sigGC_yr_t))
sigGC_yr_t <- sigGC_yr_t[,order(colnames(sigGC_yr_t))]
rownames(sigGC_yr_t) <- c("sig_rank")

GC_sig_yr <-  t(subset(t(OCG_GC_sbst), colnames(OCG_GC_sbst) %in% row.names(sigGC_yr)))

GC_sig_yr <- GC_sig_yr[,order(colnames(GC_sig_yr))]
colnames(sigGC_yr_t) == colnames(GC_sig_yr) #sanity check:TRUE

GC_sig_yr <- rbind(GC_sig_yr, sigGC_yr_t)
GC_sig_yr_t <- as.data.frame(t(GC_sig_yr))
GC_sig_yr_t$sig_rank <- as.numeric(GC_sig_yr_t$sig_rank) 
GC_sig_yr_t <- GC_sig_yr_t[order(GC_sig_yr_t$sig_rank),] 
GC_sig_yr_t <- subset(GC_sig_yr_t, select=-c(sig_rank))
GC_sig_yr_t <- as.data.frame(t(GC_sig_yr_t))

#Build objects for plotting
plotGC.yr <- data.frame(GC_sig_yr_t[,1:10], "year" = md.OCG.GC_sbst$Year, check.names = FALSE)
plotGC.yr[,1:10] <- lapply(plotGC.yr[,1:10], function(x) as.numeric(as.character(x)))
plotGC.yr[,1:10] <- log(plotGC.yr[,1:10]+1)
plotGC.yr_sub <- data.frame(sample=rownames(plotGC.yr),plotGC.yr, check.names = F)
plotGC.yr_sublong <- melt(plotGC.yr_sub)

plotGC.yr$year <-  factor(plotGC.yr$year, levels = c("2012","2021"))

### FIGURE FOR SUBSPECIES PLOIDY ####
ggplot(plotGC.yr_sublong, aes(y = value, x = year, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("GC ANCOM for year")+
  labs(color = "Compounds") +
  theme_classic()

# LCMS ANCOM ####
OCG_LCMS_3uL_subset[OCG_LCMS_3uL_subset < 10] <- 0
OCG_LCMS_3uL_subset <- OCG_LCMS_3uL_subset[rowSums(OCG_LCMS_3uL_subset) > 0,] 

summary(rowSums(OCG_LCMS_3uL_subset)) #31016795
summary(colSums(OCG_LCMS_3uL_subset)) #384380

OCG_LCMS_3uL_subset <- OCG_LCMS_3uL_subset[,colSums(OCG_LCMS_3uL_subset) > 10] 
md.OCG.LCMS.3.a <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_3uL_subset)) 

OCG_LCMS_3uL_subset.t <- t(OCG_LCMS_3uL_subset) 
OCG_LCMS_3uL_subset_t <- OCG_LCMS_3uL_subset.t[, colnames(OCG_LCMS_3uL_subset.t) %in% row.names(md.OCG.LCMS.3.a), drop = FALSE]
OCG_LCMS_3uL_subset <- t(OCG_LCMS_3uL_subset_t)
OCG_LCMS_3uL_subset <- as.data.frame(OCG_LCMS_3uL_subset) 

OCG_LCMS_3uL_subset[OCG_LCMS_3uL_subset < 10] <- 0  
OCG_LCMS_3uL_subset <- OCG_LCMS_3uL_subset[rowSums(OCG_LCMS_3uL_subset) > 0,] 

summary(rowSums(OCG_LCMS_3uL_subset)) #31016795
summary(colSums(OCG_LCMS_3uL_subset)) #384380

OCG_LCMS_3uL_subset <- OCG_LCMS_3uL_subset[,colSums(OCG_LCMS_3uL_subset) > 10]

md.OCG.LCMS.3_sbst <- data.frame("Sample.ID" = row.names(md.OCG.LCMS.3.a), md.OCG.LCMS.3.a)
OCG_LCMS_sbst <- data.frame("Sample.ID" = row.names(OCG_LCMS_3uL_subset), OCG_LCMS_3uL_subset, check.names = F)
row.names(OCG_LCMS_sbst) == row.names(md.OCG.LCMS.3_sbst) #TRUE

## SUBSPECIES LCMS ####
ANCOM_subspecies_LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant compounds
sigLCMS_subspecies <- subset(ANCOM_subspecies_LCMS$W.taxa, ANCOM_subspecies_LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_subspecies <- as.data.frame(sigLCMS_subspecies) 
row.names(sigLCMS_subspecies) <- sigLCMS_subspecies[302:1,1] 

sigLCMS_subspecies[,1] <- c(302:1) 
sigLCMS_subspecies_t <- t(sigLCMS_subspecies) 
sigLCMS_subspecies_t <- as.data.frame(sigLCMS_subspecies_t)
colnames(sigLCMS_subspecies_t) <- as.character(colnames(sigLCMS_subspecies_t))
print(colnames(sigLCMS_subspecies_t))

sigLCMS_subspecies_t <- sigLCMS_subspecies_t[,order(colnames(sigLCMS_subspecies_t))]
rownames(sigLCMS_subspecies_t) <- c("sig_rank")

LCMS.OCG_sigsbst_subspecies <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_subspecies)))

LCMS.OCG_sigsbst_subspecies <- LCMS.OCG_sigsbst_subspecies[,order(colnames(LCMS.OCG_sigsbst_subspecies))]
colnames(sigLCMS_subspecies_t) == colnames(LCMS.OCG_sigsbst_subspecies) #sanity check:TRUE

LCMS.OCG_sigsbst_subspecies <- rbind(LCMS.OCG_sigsbst_subspecies, sigLCMS_subspecies_t)
LCMS.OCG_sigsbst_subspecies_t <- as.data.frame(t(LCMS.OCG_sigsbst_subspecies))
LCMS.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(LCMS.OCG_sigsbst_subspecies_t$sig_rank) 
LCMS.OCG_sigsbst_subspecies_t <- LCMS.OCG_sigsbst_subspecies_t[order(LCMS.OCG_sigsbst_subspecies_t$sig_rank),] 
LCMS.OCG_sigsbst_subspecies_t <- subset(LCMS.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
LCMS.OCG_sigsbst_subspecies_t <- as.data.frame(t(LCMS.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotLCMS <- data.frame(LCMS.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.LCMS.3_sbst$Subspecies, check.names = FALSE)
sbsplotLCMS[,1:10] <- lapply(sbsplotLCMS[,1:10], function(x) as.numeric(as.character(x)))
sbsplotLCMS[,1:10] <- log(sbsplotLCMS[,1:10]+1)
sbsplotLCMS_sub <- data.frame(sample=rownames(sbsplotLCMS),sbsplotLCMS, check.names = F)
sbsplotLCMS_sublong <- melt(sbsplotLCMS_sub)

ANCOM_subspecies_LCMS$W.taxa

sbsplotLCMS_sublong$subspecies <-  factor(sbsplotLCMS_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS ALL LCMS ####
ggplot(sbsplotLCMS_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Subspecies")+ 
  theme_classic()

## PLOIDY LCMS ####
#Run ANCOM, specify variable
ANCOM_ploidy_LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_ploidy <- subset(ANCOM_ploidy_LCMS$W.taxa, ANCOM_ploidy_LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_ploidy <- as.data.frame(sigLCMS_ploidy)
row.names(sigLCMS_ploidy) <- sigLCMS_ploidy[302:1,1] 

sigLCMS_ploidy[,1] <- c(302:1) 
sigLCMS_ploidy_t <- t(sigLCMS_ploidy) 
sigLCMS_ploidy_t <- as.data.frame(sigLCMS_ploidy_t)
colnames(sigLCMS_ploidy_t) <- as.character(colnames(sigLCMS_ploidy_t))
print(colnames(sigLCMS_ploidy_t))

sigLCMS_ploidy_t <- sigLCMS_ploidy_t[,order(colnames(sigLCMS_ploidy_t))]
rownames(sigLCMS_ploidy_t) <- c("sig_rank")

LCMS.OCG_sigsbst_ploidy <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_ploidy)))

LCMS.OCG_sigsbst_ploidy <- LCMS.OCG_sigsbst_ploidy[,order(colnames(LCMS.OCG_sigsbst_ploidy))]
colnames(sigLCMS_ploidy_t) == colnames(LCMS.OCG_sigsbst_ploidy) #sanity check:TRUE

LCMS.OCG_sigsbst_ploidy <- rbind(LCMS.OCG_sigsbst_ploidy, sigLCMS_ploidy_t)
LCMS.OCG_sigsbst_ploidy_t <- as.data.frame(t(LCMS.OCG_sigsbst_ploidy))
LCMS.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(LCMS.OCG_sigsbst_ploidy_t$sig_rank) 
LCMS.OCG_sigsbst_ploidy_t <- LCMS.OCG_sigsbst_ploidy_t[order(LCMS.OCG_sigsbst_ploidy_t$sig_rank),] 
LCMS.OCG_sigsbst_ploidy_t <- subset(LCMS.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
LCMS.OCG_sigsbst_ploidy_t <- as.data.frame(t(LCMS.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotLCMS.ploidy <- data.frame(LCMS.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.LCMS.3_sbst$Ploidy, check.names = FALSE)
plplotLCMS.ploidy[,1:10] <- lapply(plplotLCMS.ploidy[,1:10], function(x) as.numeric(as.character(x)))
plplotLCMS.ploidy[,1:10] <- log(plplotLCMS.ploidy[,1:10]+1)
plplotLCMS.ploidy_sub <- data.frame(sample=rownames(plplotLCMS.ploidy),plplotLCMS.ploidy, check.names = F)
plplotLCMS.ploidy_sublong <- melt(plplotLCMS.ploidy_sub)

ANCOM_ploidy_LCMS$W.taxa

plplotLCMS.ploidy_sublong$ploidy<-  factor(plplotLCMS.ploidy_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS LCMS ####
ggplot(plplotLCMS.ploidy_sublong, aes(y = value, x = ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


## SUBSPECIES PLOIDY LCMS ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_subsp_ploidy <- subset(ANCOM_subsp_ploidy.LCMS$W.taxa, ANCOM_subsp_ploidy.LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_subsp_ploidy <- as.data.frame(sigLCMS_subsp_ploidy)
row.names(sigLCMS_subsp_ploidy) <- sigLCMS_subsp_ploidy[302:1,1]
sigLCMS_subsp_ploidy[,1] <- c(302:1)

sigLCMS_subsp_ploidy_t <- t(sigLCMS_subsp_ploidy)
sigLCMS_subsp_ploidy_t <- as.data.frame(sigLCMS_subsp_ploidy_t)
colnames(sigLCMS_subsp_ploidy_t) <- as.character(colnames(sigLCMS_subsp_ploidy_t))
print(colnames(sigLCMS_subsp_ploidy_t))
sigLCMS_subsp_ploidy_t <- sigLCMS_subsp_ploidy_t[,order(colnames(sigLCMS_subsp_ploidy_t))]
rownames(sigLCMS_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

LCMS_sigsbst_subspploidy <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_subsp_ploidy)))

LCMS_sigsbst_subspploidy <- LCMS_sigsbst_subspploidy[,order(colnames(LCMS_sigsbst_subspploidy))]
colnames(sigLCMS_subsp_ploidy_t) == colnames(LCMS_sigsbst_subspploidy) #sanity check:TRUE

LCMS_sigsbst_subspploidy <- rbind(LCMS_sigsbst_subspploidy, sigLCMS_subsp_ploidy_t)
LCMS_sigsbst_subspploidy_t <- as.data.frame(t(LCMS_sigsbst_subspploidy))
LCMS_sigsbst_subspploidy_t$sig_rank <- as.numeric(LCMS_sigsbst_subspploidy_t$sig_rank) 
LCMS_sigsbst_subspploidy_t <- LCMS_sigsbst_subspploidy_t[order(LCMS_sigsbst_subspploidy_t$sig_rank),] 
LCMS_sigsbst_subspploidy_t <- subset(LCMS_sigsbst_subspploidy_t, select=-c(sig_rank))
LCMS_sigsbst_subspploidy_t <- as.data.frame(t(LCMS_sigsbst_subspploidy_t))

#Build objects for plotting
plotLCMS.subspploidy <- data.frame(LCMS_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.LCMS.3_sbst$Subsp_ploidy, check.names = FALSE)
plotLCMS.subspploidy[,1:10] <- lapply(plotLCMS.subspploidy[,1:10], function(x) as.numeric(as.character(x)))
plotLCMS.subspploidy[,1:10] <- log(plotLCMS.subspploidy[,1:10]+1)
plotLCMS.subspploidy_sub <- data.frame(sample=rownames(plotLCMS.subspploidy),plotLCMS.subspploidy, check.names = F)
plotLCMS.subspploidy_sublong <- melt(plotLCMS.subspploidy_sub)

ANCOM_subsp_ploidy$W.taxa

plotLCMS.subspploidy$subsp_ploidy <-  factor(plotLCMS.subspploidy$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

### FIGURE FOR SUBSPECIES PLOIDY ####
ggplot(plotLCMS.subspploidy_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

## YEAR LCMS ####
#Run ANCOM, specify variable
ANCOM_yr.LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_yr <- subset(ANCOM_yr.LCMS$W.taxa, ANCOM_yr.LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_yr <- as.data.frame(sigLCMS_yr)
row.names(sigLCMS_yr) <- sigLCMS_yr[302:1,1]
sigLCMS_yr[,1] <- c(302:1)

sigLCMS_yr_t <- t(sigLCMS_yr)
sigLCMS_yr_t <- as.data.frame(sigLCMS_yr_t)
colnames(sigLCMS_yr_t) <- as.character(colnames(sigLCMS_yr_t))
print(colnames(sigLCMS_yr_t))
sigLCMS_yr_t <- sigLCMS_yr_t[,order(colnames(sigLCMS_yr_t))]
rownames(sigLCMS_yr_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

LCMS_sig_yr <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_yr)))

LCMS_sig_yr <- LCMS_sig_yr[,order(colnames(LCMS_sig_yr))]
colnames(sigLCMS_yr_t) == colnames(LCMS_sig_yr) #sanity check:TRUE

LCMS_sig_yr <- rbind(LCMS_sig_yr, sigLCMS_yr_t)
LCMS_sig_yr_t <- as.data.frame(t(LCMS_sig_yr))
LCMS_sig_yr_t$sig_rank <- as.numeric(LCMS_sig_yr_t$sig_rank) 
LCMS_sig_yr_t <- LCMS_sig_yr_t[order(LCMS_sig_yr_t$sig_rank),] 
LCMS_sig_yr_t <- subset(LCMS_sig_yr_t, select=-c(sig_rank))
LCMS_sig_yr_t <- as.data.frame(t(LCMS_sig_yr_t))

#Build objects for plotting
plotLCMS.yr <- data.frame(LCMS_sig_yr_t[,1:10], "year" = md.OCG.LCMS.3_sbst$Year, check.names = FALSE)
plotLCMS.yr[,1:10] <- lapply(plotLCMS.yr[,1:10], function(x) as.numeric(as.character(x)))
plotLCMS.yr[,1:10] <- log(plotLCMS.yr[,1:10]+1)
plotLCMS.yr_sub <- data.frame(sample=rownames(plotLCMS.yr),plotLCMS.yr, check.names = F)
plotLCMS.yr_sublong <- melt(plotLCMS.yr_sub)

plotLCMS.yr$year <-  factor(plotLCMS.yr$year, levels = c("2012","2021"))

### FIGURE FOR SUBSPECIES PLOIDY ####
ggplot(plotLCMS.yr_sublong, aes(y = value, x = year, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for year")+
  labs(color = "Compounds") +
  theme_classic()

# ANCOM BC #### 
## Make phyloseq objects ####
row.names(OCG_GC_2012_subset) == row.names(md.OCG.GC.2012) #TRUE
physeq <- phyloseq(OCG_GC_2012_subset, md.OCG.GC.2012)
ancombc2(OCG_GC_2012_subset)
