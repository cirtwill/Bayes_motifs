require(tidyverse) ## for efficient data manipulation and plotting

web_nos=as.character(list.files(path="../data/empirical/global_verts/matrix/", full.names=FALSE, pattern="*.tsv"))

# for (k in 1:length(web_nos)){
#     web_no <- web_nos[k]
#     load(file=paste0("Results_restricted/Linear/",web_no, ".RData"))
#     final <- mutate(collect, web = web_no)
#    # combines the result from all webs (of a particular S and C) in a single df (tab)
#      if (k==1){
#         tab <- final
#     } else {
#         tab <- bind_rows(tab, final)
#     }                   
# }
# tab$per=unlist(tab$per)
# write.table(tab,file='../data/all_empirical_linear.tsv')
# save(tab, file="Results_restricted/Linear/empirical.RData")


for (k in 1:length(web_nos)){
    web_no <- strsplit(web_nos[k],'.tsv')[[1]][1]
    load(file=paste0("Results_restricted/Nonlinear/global_verts/",web_no, ".RData"))
    final <- mutate(collect, web = web_no)
   # combines the result from all webs (of a particular S and C) in a single df (tab)
     if (k==1){
        tab <- final
    } else {
        tab <- bind_rows(tab, final)
    }                   
}
tab$per=unlist(tab$per)
write.table(tab,file='../data/global_verts_nonlinear.tsv')
save(tab, file="Results_restricted/Nonlinear/global_verts.RData")

for (k in 1:length(web_nos)){
    web_no <- strsplit(web_nos[k],'.tsv')[[1]][1]
    load(file=paste0("Results_restricted/Linear/global_verts/",web_no, ".RData"))
    final <- mutate(collect, web = web_no)
   # combines the result from all webs (of a particular S and C) in a single df (tab)
     if (k==1){
        tab <- final
    } else {
        tab <- bind_rows(tab, final)
    }                   
}
tab$per=unlist(tab$per)
write.table(tab,file='../data/global_verts_linear.tsv')
save(tab, file="Results_restricted/Linear/global_verts.RData")

