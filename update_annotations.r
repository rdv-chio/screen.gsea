update_annotations = function(input,annotations,output) {
	data = read.table(input,header=T,sep='\t',colClasses=c('character','numeric','numeric'))
	annotations = read.table(annotations,header=T,sep='\t',colClasses=c('character','character'))
	genes = vector('character',length=length(data[[1]]))
	for(i in 1:length(data[[1]])) {
		genes[i] = annotations[annotations[[1]] == data[[1]][i],2]
	}
	data = data.frame(cbind(data,genes),stringsAsFactors=F)
	data = data[,c(1,4,2,3)]
	data[[3]] = as.numeric(data[[3]])
	data[[4]] = as.numeric(data[[4]])
	colnames(data) = c('shRNA','Gene','Effect','P Value')
	write.table(data,output,row.names=F,sep='\t')
	data
}