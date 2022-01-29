require(tidyverse) ## for efficient data manipulation and plotting

sp_nos <- seq(50, 100, by=10)    # simulated sizes of networks (S)
c <- seq(0.02, 0.18, by=0.04)   # simulated connectances of networks (C)
web_nos <- c(1:100)              # web identity 


for (i in 1:length(sp_nos)){
    web_s <- sp_nos[i]
    for (j in 1:length(c)){
        connectance <- c[j]
        for (k in 1:length(web_nos)){
            web_no <- web_nos[k]
            load(file=paste0("Results_restricted/Linear/", web_s, "/", connectance, "/TESTweb_",web_no, ".RData"))
            final <- mutate(collect, web_size = web_s, c = connectance, web = web_no)
           # combines the result from all webs (of a particular S and C) in a single df (tab)
             if (k==1){
                tab <- final
            } else {
                tab <- bind_rows(tab, final)
            }                   
        }
        # combines the result from all C (and all webs 1-100) in a single df (data2)
        if (j==1){
            data2 <- tab
        } else {
            data2 <- bind_rows(data2, tab)
        }
    }
    # combines the result from all S (and all C and webs) in a single df (data) which is saved
    if (i==1){
        data <- data2
    } else {
        data <- bind_rows(data, data2)
    }
    save(data, file=paste0("Results_restricted/Linear/", "all.RData"))  
}

