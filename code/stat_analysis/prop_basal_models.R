load.image('all_tests.Rdata')

# Libraries to fit the non-plotted lms
library(lme4)
library(lmerTest)
# For permanovas and betadisper
library(vegan)
# For R2 GLMM
library(MuMIn)

# data is full data, including basals
data$network_ID=paste(data$Size,data$Connectance,data$Network,sep='_')
for(r in 1:nrow(data)){
    netdat=data[which(data$network_ID==data$network_ID[r]),]
    nBasal=length(which(netdat$in_Degree==0))
    n=data$Size[r]
    data$prop_Basal[r]=nBasal/n
}

# consumers is heterotrophs only


# A: Does persistence vary with proportion of basal resources?
SCper=with(consumers,lm(Persistence~scale(prop_Basal)*scale(Size)*scale(Connectance)*scale(Disturbance)))


# B: How do motif profiles vary with prop_Basal and C?
    motdata=read.table('../../data/3sp_motif_profiles.tsv',header=TRUE)
    motdata$mot_total=rowSums(motdata[,4:7])
    motdata$prop_chain=motdata$m12/motdata$mot_total
    motdata$prop_apparent=motdata$m6/motdata$mot_total
    motdata$prop_direct=motdata$m36/motdata$mot_total
    motdata$prop_omni=motdata$m38/motdata$mot_total

    # Expanding motif profiles for each level of disturbance
    motdata$netID=paste(motdata$Size,motdata$Connectance,motdata$Network,sep=':')
    consumers$netID=paste(consumers$Size,consumers$Connectance,consumers$Network,sep=':')
    consumers$net_chain=0
    consumers$net_apparent=0
    consumers$net_direct=0
    consumers$net_omni=0
    for(r in 1:nrow(consumers)){
        net=consumers$netID[r]
        consumers$net_chain[r]=motdata$prop_chain[which(motdata$netID==net)]
        consumers$net_apparent[r]=motdata$prop_apparent[which(motdata$netID==net)]
        consumers$net_direct[r]=motdata$prop_direct[which(motdata$netID==net)]
        consumers$net_omni[r]=motdata$prop_omni[which(motdata$netID==net)]
    }

    save.image('analysis.Rdata')

    total_lm=with(motdata,lm(mot_total~prop_Basal*Connectance))
    summary(total_lm)

    # PERMANOVA and variability of motif profiles
    propdist=vegdist(motdata[,10:13],method="bray")
    propperm=adonis(propdist~motdata$prop_Basal*motdata$Connectance,perm=9999)
    motdata$interact=paste0(motdata$prop_Basal,motdata$Connectance,sep=':')
    propinterdisp=betadisper(propdist,motdata$interact,type="centroid")
    Spropdisp=betadisper(propdist,motdata$prop_Basal,type="centroid")
    Cpropdisp=betadisper(propdist,motdata$Connectance,type="centroid")

    propintanova=anova(propinterdisp)
    Spropanova=anova(Spropdisp)
    Cpropanova=anova(Cpropdisp)

    # Printing variabilities of motif profiles for plotting
    interdata2=cbind(propinterdisp$distances,as.character(propinterdisp$group))
    colnames(interdata2)=c("Distance","BC")
    interdata2=as.data.frame(interdata2)    

    interdata2$SC=as.character(interdata2$BC)
    for(r in 1:nrow(interdata2)){
        levelcode=interdata2$"BC"[r]
        interdata2$S[r]=strsplit(levelcode,':')[[1]][1]
        interdata2$C[r]=strsplit(levelcode,':')[[1]][2]
    }
    write.table(interdata2,file='proportion_variability_BC.tsv',sep='\t')

    # LMs
    pchainlm=with(motdata,lm(prop_chain~prop_Basal*Connectance))
    pomnilm=with(motdata,lm(prop_omni~prop_Basal*Connectance))
    papplm=with(motdata,lm(prop_apparent~prop_Basal*Connectance))
    pdirlm=with(motdata,lm(prop_direct~prop_Basal*Connectance))


# C: How does motif participation vary with other properties?
    # Using subdata, just one level of basal_p, to not sextuple-count webs
    subdata=consumers[which(consumers$Disturbance==0.1),]

    # prop_Basal and C
    omni_prop=with(subdata,lm(prop_omni~prop_Basal*Connectance))
    chain_prop=with(subdata,lm(prop_chain~prop_Basal*Connectance))
    apparent_prop=with(subdata,lm(prop_apparent~prop_Basal*Connectance))
    direct_prop=with(subdata,lm(prop_direct~prop_Basal*Connectance))
    results=matrix(nrow=4,ncol=5)
    results[1,]=c("Omnivory",summary(omni_prop)$coefficients[,1])
    results[2,]=c("Chain",summary(chain_prop)$coefficients[,1])
    results[3,]=c("Apparent",summary(apparent_prop)$coefficients[,1])
    results[4,]=c("Direct",summary(direct_prop)$coefficients[,1])
    write.table(results,file='roles_vs_BC.tsv',sep='\t')


# D: How do frequencies of different positions vary with prop_Basal, C?

    for(c in 1:ncol(posprops)){
        posname=colnames(posprops)[c]
        SClm=lm(posprops[,c]~consumers$prop_Basal*consumers$Connectance)
        write.table(summary(SClm)$coefficients,file=paste0('positions_vs_other/BC_',posname,'.tsv',sep=''),sep='\t')
    }

save.image('basal_tests.Rdata')

