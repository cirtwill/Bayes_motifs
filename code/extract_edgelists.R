library(igraph)
# Convert graphs to edgelists

for(s in seq(50,100,10)){
	print(paste0(as.character(s),' species'))
	dir.create(paste0('../data/edgelists/acyclic/',as.character(s)))
	for(c in seq(0.02,0.18,0.04)){
		print(paste0('connectance=',as.character(c)))
		dir.create(paste0('../data/edgelists/acyclic/',as.character(s),'/',as.character(c)))
        matrixfiles<-as.character(list.files(path=paste0("../data/networks/acyclic_webs/",as.character(s),"/",as.character(c)), full.names=TRUE, pattern="*.csv"))
        if(!length(matrixfiles)==100){
        	print(c(s,c,length(matrixfiles)))
        	print('Error: Too many or too few webs')
        }
		for(i in 1:length(matrixfiles)){
			# Select the json file to work with
			infile=matrixfiles[i]
			# Keeping the same number as the matrices
			webnum=strsplit(strsplit(infile,'_')[[1]][4],'.csv')[[1]][1]
			A=read.csv(infile,header=TRUE,sep=',')
			# convert to edgelist and save
			edges=get.edgelist(graph.adjacency(as.matrix(A)))
			write.table(edges,file=paste0('../data/edgelists/acyclic/',as.character(s),'/',as.character(c),'/initial_edges_',as.character(webnum),'.tsv',sep=''),sep='\t',row.names=FALSE,col.names=FALSE)
		}
	}
}

