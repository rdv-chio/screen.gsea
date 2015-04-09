## to reconstruct geneset database (gmt file) from correlation matrix

make_gmt = function(input,output='') {
	matrix = read.table(input,header=T,sep='\t')
	max.size = max(rowSums(matrix))
	gs.db = matrix(rep("null",nrow(matrix)*(max.size+2)), nrow=nrow(matrix), ncol=max(rowSums(matrix))+2)
	for(i in 1:nrow(matrix)) {
		gs.db[i,] = c(rownames(matrix)[i],'',colnames(matrix)[which(matrix[i,] == 1)],rep('null',max.size - sum(matrix[i,])))
	}
	write.table(gs.db,output,row.names=F,col.names=F,sep='\t')
}