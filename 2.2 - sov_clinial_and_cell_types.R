library(ggplot2)
library(parallel)
library(dplyr)
library(car)
library(lumi)

#load myCombat from 1 - Raw data pre-processing.Rmd 
#load pd_all from 2.0 - cell_type_deconvolution_Gervin.R
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/pd_all.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/myCombat.RData")

#convert beta to m values
myCombat = myCombat[, rownames(pd)]
myCombat[myCombat == 0] <- 1e-6  # Replace 0 with a small value
myCombat[myCombat == 1] <- 1 - 1e-6  # Replace 1 with a value slightly below 1
m = beta2m(myCombat)

pd = pd%>%dplyr::select("Sample_Group",  "Age", "BMI", "Eth2", "GA", "Parity","Smoker","CD4T",
                        "CD8T","Bcell","Gran","NK","Mono")#excluded nrbc to avoid perfect collinearity
pd$Eth2 = as.numeric(as.factor(pd$Eth2))

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
  
  #use sum contrast for type III ANOVA
  contrasts_list <- lapply(newdata, function(col) {
    if (is.factor(col)) return(contr.sum(length(levels(col))))
    return(NULL)
  })
  options(contrasts = c("contr.sum", "contr.poly"))
  
  #design matrix
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

save(Ftab, Fmean, finalvars, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/sov_res_0212.RData")


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

pdf("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/SOV_Gervin_0212.pdf")
sovplot(restab = Fmean, plottype = 'F', textsize = 7)
dev.off()

Fmean
