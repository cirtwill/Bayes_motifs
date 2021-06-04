require(tidyverse) ## for efficient data manipulation and plotting
require(igraph)
source("getTL.R") ## to get prey average TL
dir.create('../data/TL',showWarnings=FALSE) # Storing the degree/TL files in the same way as the motif files
sp_nos <- seq(50,100, by=10)
c <- seq(0.02,0.18,by=0.04)
basal_proportion <- seq(0,1,by=0.1)

for (i in 1:length(sp_nos)){
    sp_no <- sp_nos[i]
    print(sp_no)
    dir.create(paste0('../data/TL/',sp_no,sep=''),showWarnings=FALSE)
    for (j in 1:length(c)){
        connectance <- c[j]
        print(connectance)
        dir.create(paste0('../data/TL/',sp_no,'/',connectance,sep=''),showWarnings=FALSE)
        matrixfiles<-as.character(list.files(path=paste0("../data/networks/acyclic_webs/",as.character(sp_no),"/",as.character(connectance)), full.names=TRUE, pattern="*.csv"))
        if(!length(matrixfiles)==100){
            print(c(sp_no,connectance,length(matrixfiles)))
            print('Error: Too many or too few webs')
        }
        for (k in 1:length(matrixfiles)){
            infile=matrixfiles[k]
            webnum=strsplit(strsplit(infile,'_')[[1]][4],'.csv')[[1]][1]
            A <- read.table(infile, sep=",", header=TRUE)
            S <- dim(A)[1]
            A <- as.matrix(A)
            dd <- distances(graph_from_adjacency_matrix(A), mode="out")
            basals <- rownames(A[rowSums(A)==0,])
            short_tl <- 1 + as.numeric(apply(dd[,as.numeric(basals)], 1, min))
            # Getting a subscript out of bounds error from the dd[,basals] withtout the as.numeric()
            FW <- list()
            FW$M <- t(A)
            FW <- getTL(FW)
            
            df <- data_frame(
                sp_no = rownames(A), 
                in_degree = rowSums(A),
                out_degree = colSums(A),
                shortest_tl = short_tl,
                preyav_tl = round(FW$TL,2)
            )
            write.table(df,file=paste0("../data/TL/", sp_no, "/", connectance, "/web_",webnum, ".tsv"),sep='\t')
            # save(df, file=paste0("../data/TL/", sp_no, "/", connectance, "/web_",web_no, ".RData"))
        }
    }
}



