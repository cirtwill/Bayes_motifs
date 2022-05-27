require(XRJulia) ## package allowing for Julia interfacing
require(igraph)
require(tidyverse) ## for efficient data manipulation and plotting
source("BN-DFS.R") ## to make the adjacency matrix acyclic
source("getTL.R") ## to get prey average TL, filter good webs

## sort topological and remove cycles
topo_sort_acyclic <- function(A) {
    A<-as.matrix(A)
    ## Make sure producers are the first rows; this is needed for the DFS
    ## algorithm to work (which removes all cycles from the digraph)
    A <- A[order(rowSums(A)),order(rowSums(A))]
    ## Now make the digraph acyclic
    A <- DFS(A)
    
    ## Find a sorting of A's rows and columns to make A lower triangular
    o <- graph_from_adjacency_matrix(A) %>% topo_sort(mode="in")
    ## Return sorted, lower triangular adjacency matrix A
    return(A[o,o])
}

goodfiles=as.character(list.files(path="../data/empirical/global_verts/matrix/", full.names=FALSE, pattern="*.tsv"))
for (g in goodfiles){
    print(g)
    A <- read.table(paste("../data/empirical/global_verts/matrix/",g,sep=''), sep="\t", row.names=1, header=TRUE)
    rownames(A) <- colnames(A) <- 1:nrow(A) 
    # image(t(apply(A,2,rev)))
    
    A <- topo_sort_acyclic(A)
    # image(t(apply(A,2,rev)))
    print(dim(A))
    print(sum(A)/(nrow(A)*(nrow(A)-1)))
    # Calculate the PATL to sort good webs
    FW <- list()
    FW$M <- t(A)
    FW <- getTL(FW)
    preyav_tl=round(FW$TL,2)
    # Let's only keep the ones with max PATL<6
    # Also indexing by "good webs" so I know how many to re-run

    write.table(A, file=paste("../data/empirical/global_verts/acyclic_matrix/",g,sep=''), sep=",")
    edges=get.edgelist(graph.adjacency(as.matrix(A)))
    write.table(edges,file=paste0('../data/empirical/global_verts/edgelists/',g,sep=''),sep='\t',row.names=FALSE,col.names=FALSE)
}



