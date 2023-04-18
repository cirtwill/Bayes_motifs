#### TEST MOTIFS AND PERSISTENCE 
require(ggpmisc)
require(tidyverse)

# Libraries to fit the non-plotted lms
library(lme4)
library(lmerTest)
# For permanovas and betadisper
library(vegan)
# For R2 GLMM
library(MuMIn)


data<-read.table(file = '../../data/3sp_roles_participation.tsv', sep = '\t', header = TRUE)
data$Disturbance<-0.1+0.4*data$Basal_p
# ONLY FIT MODELS ON CONSUMERS!!!
consumers<-data[which(data$in_Degree>0),]
# Check that all predictors are numeric in case of R shenanigans
# All should be numeric except for species ID
for(c in c((1:4),(6:25))){
    if(!is.numeric(consumers[,c])){
        consumers[,c]=as.numeric(as.character(consumers[,c]))
        print(paste0(colnames(consumers)[c],' made numeric'))
    }
}

consumers$total_motifs=consumers$m12+consumers$m6+consumers$m36+consumers$m38
consumers$prop_chain=consumers$m12/consumers$total_motifs
consumers$prop_apparent=consumers$m6/consumers$total_motifs
consumers$prop_direct=consumers$m36/consumers$total_motifs
consumers$prop_omni=consumers$m38/consumers$total_motifs
consumers$netID=paste(consumers$Size,consumers$Connectance,consumers$Network,sep=':')

per_network_results=matrix(nrow=0,ncol=6)
colnames(per_network_results)=c("Network","S","C","Disturbance","Motif","Motif_slope")

for(n in levels(as.factor(consumers$Network))){
    print(n)
    netdat=consumers[which(consumers$Network==n),]
    S=netdat$Size[1]
    C=netdat$Connectance[1]
    for (d in c(0,0.2,0.4,0.6,0.8,1)){
        Ddat=netdat[which(netdat$Basal_p==d),]
        omnimod=with(Ddat,glm(Persistence~prop_omni,family='binomial'))
        chainmod=with(Ddat,glm(Persistence~prop_chain,family='binomial'))
        ACmod=with(Ddat,glm(Persistence~prop_apparent,family='binomial'))
        DCmod=with(Ddat,glm(Persistence~prop_direct,family='binomial'))
        per_network_results=rbind(per_network_results,
            c(n,S,C,d,'Omnivory',summary(omnimod)$coefficients[2,1]))
        per_network_results=rbind(per_network_results,
            c(n,S,C,d,'Chain',summary(chainmod)$coefficients[2,1]))
        per_network_results=rbind(per_network_results,
            c(n,S,C,d,'Apparent',summary(ACmod)$coefficients[2,1]))
        per_network_results=rbind(per_network_results,
            c(n,S,C,d,'Direct',summary(DCmod)$coefficients[2,1]))
    }
}

write.table(per_network_results,file='../../data/per_network_regressions.tsv',row.names=FALSE)
save.image('per_network_regressions.Rdata')

