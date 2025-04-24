#sov_EOLO data
library(ggplot2)
library(EpiDISH)
library(parallel)
library(dplyr)
library(carData)
library(car)
library(lumi)

#load results from EOLO_analysis.R
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/pd_EOLO.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/beta_EOLO.RData")

# https://bioconductor.org/packages/release/data/experiment/html/FlowSorted.CordBloodCombined.450k.html
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=beta_EOLO, ref.m=ref, method="CP")

pd_EOLO$eope = ifelse(pd_EOLO$`group:ch1` == "Early-onset preeclampsia", 1, 0)
pd = as.data.frame(cbind(out$estF, pd_EOLO))
pd = pd%>%dplyr::select("group","CD4T",
                        "CD8T","Bcell","Gran","NK","Mono")

beta_EOLO[beta_EOLO == 0] <- 1e-6  # Replace 0 with a small value
beta_EOLO[beta_EOLO == 1] <- 1 - 1e-6  # Replace 1 with a value slightly below 1
m = beta2m(beta_EOLO)
sovdat = data.frame(t(m))
Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)

for(i in 2:ncol(pd)){
  varname <- names(pd)[i]
  Ftab[varname] <- numeric()
}

calF <- function(probe = probecol){
  newdata <- pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  contrasts_list <- lapply(newdata, function(col) {
    if (is.factor(col)) return(contr.sum(length(levels(col))))
    return(NULL)
  })
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  fit <- lm(formstr, data = newdata)
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  F <- aovfit$`F value`
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  Ftab <- rbind(Ftab, F)
  return(Ftab)
}

Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 16)
Ftab <- do.call(rbind, Ftab)
Fmean <- colMeans(Ftab)
Fmean <- Fmean[order(-Fmean)]
Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)
finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))

save(Ftab, Fmean, finalvars, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/sov_res_m.RData")

#load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/sov_res_m.RData")
sovplot <- function(restab = MSSmean, clustername = 'Case', plottype = 'MSS',
                    textsize = 20){
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))+
        geom_hline(yintercept=1, linetype="dashed", color = "red")
    )
  }
}

pdf("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/SOV_EOLO_Gervin_m.pdf")
sovplot(restab = Fmean, plottype = 'F', textsize = 14)
dev.off()

Fmean
