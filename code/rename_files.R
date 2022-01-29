sp_nos <- seq(50, 100, by=10)
c <- seq(0.02, 0.18, by=0.04)
for (i in 1:length(sp_nos)){
    i <- sp_nos[i]
    for (j in 1:length(c)){
        j <- c[j]
        old_files <- list.files(paste0("/Users/annak52/Documents/motifs_and_stability/webs_restricted_tl/",i,"/",j), pattern="*.csv", full.names = TRUE)
        new_files <- paste0("/Users/annak52/Documents/motifs_and_stability/websR/", i, "/", j, "/initial_matrix_",1:length(old_files),".csv")
        file.copy(from = old_files, to = new_files)
    }
}
