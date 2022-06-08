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

# # Proportions
# Creating proportions. Note that because proportions always sum to 1, doesn't make sense to do the models with all motifs together. 
# Have to remove one motif so that no predictor is a linear combination of the others.
consumers$total_motifs=consumers$m12+consumers$m6+consumers$m36+consumers$m38
consumers$prop_chain=consumers$m12/consumers$total_motifs
consumers$prop_apparent=consumers$m6/consumers$total_motifs
consumers$prop_direct=consumers$m36/consumers$total_motifs
consumers$prop_omni=consumers$m38/consumers$total_motifs
consumers$netID=paste(consumers$Size,consumers$Connectance,consumers$Network,sep=':')


# Section 2: How do motif profiles vary with S and C?
    motdata=read.table('../../data/3sp_motif_profiles.tsv',header=TRUE)
    motdata$mot_total=rowSums(motdata[,4:7])
    motdata$prop_chain=motdata$m12/motdata$mot_total
    motdata$prop_apparent=motdata$m6/motdata$mot_total
    motdata$prop_direct=motdata$m36/motdata$mot_total
    motdata$prop_omni=motdata$m38/motdata$mot_total

    # Expanding motif profiles for each level of disturbance
    motdata$netID=paste(motdata$Size,motdata$Connectance,motdata$Network,sep=':')
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

    total_lm=with(motdata,lm(mot_total~Size*Connectance))
    summary(total_lm)

    # PERMANOVA and variability of motif profiles
    propdist=vegdist(motdata[,10:13],method="bray")
    propperm=adonis(propdist~motdata$Size*motdata$Connectance,perm=9999) # Should I stratify within a network? Probably doesn't matter so much since dispersion not normal
    motdata$interact=paste0(motdata$Size,':',motdata$Connectance)
    propinterdisp=betadisper(propdist,motdata$interact,type="centroid")
    Spropdisp=betadisper(propdist,motdata$Size,type="centroid")
    Cpropdisp=betadisper(propdist,motdata$Connectance,type="centroid")

    propintanova=anova(propinterdisp)
    Spropanova=anova(Spropdisp)
    Cpropanova=anova(Cpropdisp)

    # Printing variabilities of motif profiles for plotting
    interdata2=cbind(propinterdisp$distances,as.character(propinterdisp$group))
    colnames(interdata2)=c("Distance","SC")
    interdata2=as.data.frame(interdata2)    

    interdata2$SC=as.character(interdata2$SC)
    for(r in 1:nrow(interdata2)){
        levelcode=interdata2$"SC"[r]
        interdata2$S[r]=strsplit(levelcode,':')[[1]][1]
        interdata2$C[r]=strsplit(levelcode,':')[[1]][2]
    }
    write.table(interdata2,file='proportion_variability_SC.tsv',sep='\t')

    # LMs - no need for network-level random effect since it's within-network motif profiles
    pchainlm=with(motdata,lm(prop_chain~Size*Connectance))
    pomnilm=with(motdata,lm(prop_omni~Size*Connectance))
    papplm=with(motdata,lm(prop_apparent~Size*Connectance))
    pdirlm=with(motdata,lm(prop_direct~Size*Connectance))


# Section 3: How does persistence vary with network motif profile?
    # Needs a separate data frame for profiles repeated per disturbance, mean persistence per network/disturbance
    simdata=matrix(nrow=0,ncol=ncol(motdata)+2)
    for(net in levels(as.factor(consumers$netID))){
        netdat=consumers[which(consumers$netID==as.character(net)),]
        for(dist in levels(as.factor(consumers$Disturbance))){
            minidat=motdata[which(motdata$netID==net),]
            per=mean(netdat$Persistence[which(round(netdat$Disturbance,2)==as.numeric(as.character(dist)))])
            rowly=cbind(minidat,per,dist)
            simdata=rbind(simdata,rowly)
        }
    }

    simdata=as.data.frame(simdata)
    simdata$per=as.numeric(simdata$per)
    simdata$dist=as.numeric(as.character(simdata$dist))

    # Mean persistence ~ chain, omni, apparent, direct with disturbance
    netchain_mean=with(simdata,lm(per~scale(prop_chain)*scale(dist)))
    netomni_mean=with(simdata,lm(per~scale(prop_omni)*scale(dist)))
    netapparent_mean=with(simdata,lm(per~scale(prop_apparent)*scale(dist)))
    netdirect_mean=with(simdata,lm(per~scale(prop_direct)*scale(dist)))


    # PERMANOVA needs to be on mean persistence across levels of disturbance, for size constraints.
    for(r in 1:nrow(motdata)){
        motdata$persistence[r]=mean(simdata$per[which(simdata$netID==motdata$netID[r])])
    }

    # Check cols
    propdist=vegdist(motdata[,10:13], method='bray')
    persist_perm=adonis(propdist~motdata$persistence,perm=9999)
    # Betadisper of motif profiles at each level of mean persistence
    # Binning mean_persistence since levels are generally unique.
    persist_disp=betadisper(propdist,round(motdata$persistence,3),type="centroid")
    # ANOVA of dispersion~mean persistence
    persist_anova=anova(persist_disp)
    # LM dispersion ~ mean persistence
    persist_lm=lm(persist_disp$distances~round(motdata$persistence,3))

# Section 1: How does persistence vary with global network properties?

    # SCper=with(consumers,lm(Persistence~scale(Size)*scale(Connectance)*scale(Disturbance)))
    SCper_mean=with(simdata,lm(per~scale(Size)*scale(Connectance)*scale(dist)))

# Section 4: How does persistence vary with proportions?
    # # lmer with random intercepts for each level of S:C
    # consumers$Global=paste0(consumers$Size,':',consumers$Connectance)
    # propchain_lmer1<-with(consumers,lmer(Persistence~scale(prop_chain)*scale(Disturbance)+(1|Global)))
    # propapparent_lmer1<-with(consumers,lmer(Persistence~scale(prop_apparent)*scale(Disturbance)+(1|Global)))
    # propdirect_lmer1<-with(consumers,lmer(Persistence~scale(prop_direct)*scale(Disturbance)+(1|Global)))
    # propomni_lmer1<-with(consumers,lmer(Persistence~scale(prop_omni)*scale(Disturbance)+(1|Global)))

    # R2chain=r.squaredGLMM(propchain_lmer1,null=lmer(consumers$Persistence~(1|Global)))
    # R2app=r.squaredGLMM(propapparent_lmer1,null=lmer(consumers$Persistence~(1|Global)))
    # R2dir=r.squaredGLMM(propdirect_lmer1,null=lmer(consumers$Persistence~(1|Global)))
    # R2omni=r.squaredGLMM(propomni_lmer1,null=lmer(consumers$Persistence~(1|Global)))
    # Random intercepts for S:C and network ID. Conclusions same as above.
    propchain_lmer1_randnet<-with(consumers,lmer(Persistence~scale(prop_chain)*scale(Disturbance)+(1|Global)+(1|netID)))
    write.table(summary(propchain_lmer1_randnet)$coefficients,file='chain_lm.tsv',sep='\t')
    propapparent_lmer1_randnet<-with(consumers,lmer(Persistence~scale(prop_apparent)*scale(Disturbance)+(1|Global)+(1|netID)))
    write.table(summary(propapparent_lmer1_randnet)$coefficients,file='apparent_lm.tsv',sep='\t')
    propdirect_lmer1_randnet<-with(consumers,lmer(Persistence~scale(prop_direct)*scale(Disturbance)+(1|Global)+(1|netID)))
    write.table(summary(propdirect_lmer1_randnet)$coefficients,file='direct_lm.tsv',sep='\t')
    propomni_lmer1_randnet<-with(consumers,lmer(Persistence~scale(prop_omni)*scale(Disturbance)+(1|Global)+(1|netID)))
    write.table(summary(propomni_lmer1_randnet)$coefficients,file='omnivory_lm.tsv',sep='\t')

    R2chain_randnet=r.squaredGLMM(propchain_lmer1_randnet,null=lmer(consumers$Persistence~(1|Global)+(1|netID)))
    R2app_randnet=r.squaredGLMM(propapparent_lmer1_randnet,null=lmer(consumers$Persistence~(1|Global)+(1|netID)))
    R2dir_randnet=r.squaredGLMM(propdirect_lmer1_randnet,null=lmer(consumers$Persistence~(1|Global)+(1|netID)))
    R2omni_randnet=r.squaredGLMM(propomni_lmer1_randnet,null=lmer(consumers$Persistence~(1|Global)+(1|netID)))


# Section 5. How does persistence vary with degree and TL?
    # Persistence vs. STL and Degree
    # TLper=with(consumers,lmer(Persistence~scale(STL)*scale(Disturbance)+(1|Global)))
    # Degper=with(consumers,lmer(Persistence~scale(in_Degree)*scale(Disturbance)+(1|Global)))

    TLper_randnet=with(consumers,lmer(Persistence~scale(STL)*scale(Disturbance)+(1|Global)+(1|netID)))
    Degper_randnet=with(consumers,lmer(Persistence~scale(in_Degree)*scale(Disturbance)+(1|Global)+(1|netID)))

    write.table(summary(TLper_randnet)$coefficients,file='persistence_vs_TL.tsv',sep='\t')
    write.table(summary(Degper_randnet)$coefficients,file='persistence_vs_Deg.tsv',sep='\t')


# Section 6. How does motif participation vary with other properties?
    # Using subdata, just one level of basal_p, to not sextuple-count webs
    subdata=consumers[which(consumers$Disturbance==0.1),]

    # # S and C
    # omni_prop=with(subdata,lm(prop_omni~Size*Connectance))
    # chain_prop=with(subdata,lm(prop_chain~Size*Connectance))
    # apparent_prop=with(subdata,lm(prop_apparent~Size*Connectance))
    # direct_prop=with(subdata,lm(prop_direct~Size*Connectance))
    # results=matrix(nrow=4,ncol=5)
    # results[1,]=c("Omnivory",summary(omni_prop)$coefficients[,1])
    # results[2,]=c("Chain",summary(chain_prop)$coefficients[,1])
    # results[3,]=c("Apparent",summary(apparent_prop)$coefficients[,1])
    # results[4,]=c("Direct",summary(direct_prop)$coefficients[,1])
    # write.table(results,file='roles_vs_SC.tsv',sep='\t')
    # # Results qualitatively similar with network-level ranef
    omni_prop_randnet=with(subdata,lmer(prop_omni~Size*Connectance+(1|netID)))
    chain_prop_randnet=with(subdata,lmer(prop_chain~Size*Connectance+(1|netID)))
    apparent_prop_randnet=with(subdata,lmer(prop_apparent~Size*Connectance+(1|netID)))
    direct_prop_randnet=with(subdata,lmer(prop_direct~Size*Connectance+(1|netID)))
    results=matrix(nrow=4,ncol=5)
    results[1,]=c("Omnivory",summary(omni_prop_randnet)$coefficients[,1])
    results[2,]=c("Chain",summary(chain_prop_randnet)$coefficients[,1])
    results[3,]=c("Apparent",summary(apparent_prop_randnet)$coefficients[,1])
    results[4,]=c("Direct",summary(direct_prop_randnet)$coefficients[,1])
    write.table(results,file='roles_vs_SC.tsv',sep='\t')

    # # Deg and TL
    # # relating counts and proportions to degree, TL (repeats analysis in Cirwill & Wootton, in prep.)
    # # Degree, prop
    # chain_deg_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_chain~in_Degree+(1|Global)))
    # omni_deg_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_omni~in_Degree+(1|Global)))
    # apparent_deg_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_apparent~in_Degree+(1|Global)))
    # direct_deg_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_direct~in_Degree+(1|Global)))
    # # TL, prop
    # chain_TL_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_chain~STL+(1|Global)))
    # omni_TL_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_omni~STL+(1|Global)))
    # apparent_TL_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_apparent~STL+(1|Global)))
    # direct_TL_prop=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_direct~STL+(1|Global)))

    chain_deg_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_chain~in_Degree+(1|Global)+(1|netID)))
    omni_deg_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_omni~in_Degree+(1|Global)+(1|netID)))
    apparent_deg_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_apparent~in_Degree+(1|Global)+(1|netID)))
    direct_deg_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_direct~in_Degree+(1|Global)+(1|netID)))
    # TL, prop
    chain_TL_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_chain~STL+(1|Global)+(1|netID)))
    omni_TL_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_omni~STL+(1|Global)+(1|netID)))
    apparent_TL_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_apparent~STL+(1|Global)+(1|netID)))
    direct_TL_prop_randnet=with(consumers[which(consumers$Disturbance==0.1),],lmer(prop_direct~STL+(1|Global)+(1|netID)))

    results=matrix(nrow=8,ncol=7)
    colnames(results)=c("Predctor","Role","Motif","Intercept","Intercept_SD","Pred_slope","Pred_SD")
    results[1,]=c("Deg","Prop","Chain",summary(chain_deg_prop_randnet)$coefficients[1,1:2],summary(chain_deg_prop_randnet)$coefficients[2,1:2])
    results[2,]=c("Deg","Prop","Omnivory",summary(omni_deg_prop_randnet)$coefficients[1,1:2],summary(omni_deg_prop_randnet)$coefficients[2,1:2])
    results[3,]=c("Deg","Prop","Apparent",summary(apparent_deg_prop_randnet)$coefficients[1,1:2],summary(apparent_deg_prop_randnet)$coefficients[2,1:2])
    results[4,]=c("Deg","Prop","Direct",summary(direct_deg_prop_randnet)$coefficients[1,1:2],summary(direct_deg_prop_randnet)$coefficients[2,1:2])
    results[5,]=c("TL","Prop","Chain",summary(chain_TL_prop_randnet)$coefficients[1,1:2],summary(chain_TL_prop_randnet)$coefficients[2,1:2])
    results[6,]=c("TL","Prop","Omnivory",summary(omni_TL_prop_randnet)$coefficients[1,1:2],summary(omni_TL_prop_randnet)$coefficients[2,1:2])
    results[7,]=c("TL","Prop","Apparent",summary(apparent_TL_prop_randnet)$coefficients[1,1:2],summary(apparent_TL_prop_randnet)$coefficients[2,1:2])
    results[8,]=c("TL","Prop","Direct",summary(direct_TL_prop_randnet)$coefficients[1,1:2],summary(direct_TL_prop_randnet)$coefficients[2,1:2])
    write.table(results,file='roles_vs_TL_Deg.tsv',sep='\t')


# Not repeating position materials with network-level random effect since we're not using them.

# # Section 7. How does persistence vary with positions? ## Highlight number of prey, predators in figure

#    positions=consumers[,15:24] 
#    posprops=positions/rowSums(positions)

#    poslms=matrix(nrow=10,ncol=9)
#    for(c in 1:ncol(posprops)){
#     lm=lmer(consumers$Persistence~consumers$Disturbance*posprops[,c]+(1|consumers$Global))
#     poslms[c,1]=colnames(posprops)[c]
#     poslms[c,2:5]=summary(lm)$coefficients[,1]
#     poslms[c,6:9]=summary(lm)$coefficients[,5]
#    }
#    colnames(poslms)=c("Position","Intercept","Disturbance","Position","Interaction","Intercept_p","Disturbance_p","Position_p","Interaction_p")
#    write.table(poslms,file='persistence_vs_positions.tsv',sep='\t',row.names=FALSE)


# # Section 8. How do frequencies of different positions vary with S, C, in-degree, TL?

#     for(c in 1:ncol(posprops)){
#         posname=colnames(posprops)[c]
#         SClm=lm(posprops[,c]~consumers$Size*consumers$Connectance)
#         DTlm=lm(posprops[,c]~consumers$in_Degree*consumers$STL)
#         write.table(summary(SClm)$coefficients,file=paste0('positions_vs_other/SC_',posname,'.tsv',sep=''),sep='\t')
#         write.table(summary(DTlm)$coefficients,file=paste0('positions_vs_other/DT_',posname,'.tsv',sep=''),sep='\t')
#     }

# # Other, supplemental tests:\
save.image('all_tests.Rdata')


# # Do proportions of motifs vary due to network processing?

# original_data=read.table('../../data/3sp_motif_profiles_cyclicwebs.tsv',header=TRUE)
# original_data$prop_chain=original_data$m12/original_data$Total
# original_data$prop_apparent=original_data$m6/original_data$Total
# original_data$prop_direct=original_data$m36/original_data$Total
# original_data$prop_omni=original_data$m38/original_data$Total
# original_data$prop_other=(original_data$Total-rowSums(original_data[,4:7]))/original_data$Total

# data$label=paste(data$Size,data$Connectance,data$Network,sep=':')
# original_data$label=paste(original_data$Size,original_data$Connectance,original_data$Network,sep=':')

# used=original_data[which(original_data$label%in%data$label),]



# # # Z-scores (comparing a species' role to the average role within its web)
# # Making a unique label for each web, to calculate means and SDs
# # Means, SDs based on just one level of Basal_p (roles don't change with basal_p, makes sample size reflect n(webs))
# consumers$webID=as.factor(paste(as.character(consumers$Size),as.character(consumers$Connectance),as.character(consumers$Network),sep='_'))
# # Calculate mean and SD for each motif count, in each web
# subdata=consumers[which(consumers$Disturbance==0.1),]
# # Chain
# webmean_chain=tapply(subdata$m12,subdata$webID,mean)
# webSD_chain=tapply(subdata$m12,subdata$webID,sd)
# consumers$Zchain=(consumers$m12-webmean_chain[consumers$webID])/webSD_chain[consumers$webID]
# # Apparent competition
# webmean_app=tapply(subdata$m6,subdata$webID,mean)
# webSD_app=tapply(subdata$m6,subdata$webID,sd)
# consumers$Zapparent=(consumers$m6-webmean_app[consumers$webID])/webSD_app[consumers$webID]
# # Direct competition
# webmean_dir=tapply(subdata$m36,subdata$webID,mean)
# webSD_dir=tapply(subdata$m36,subdata$webID,sd)
# consumers$Zdirect=(consumers$m36-webmean_dir[consumers$webID])/webSD_dir[consumers$webID]
# # Omnivory
# webmean_omni=tapply(subdata$m38,subdata$webID,mean)
# webSD_omni=tapply(subdata$m38,subdata$webID,sd)
# consumers$Zomni=(consumers$m38-webmean_omni[consumers$webID])/webSD_omni[consumers$webID]


