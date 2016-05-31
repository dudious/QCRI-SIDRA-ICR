Sample.Meta.data<-Sample.Meta.data[-nrow(Sample.Meta.data),]
Expression.Data<-Expression.Data[,-ncol(Expression.Data)]
Expression.Data<-as.matrix(Expression.Data)
mode (Expression.Data) <- "numeric"
colnames(Sample.Meta.data)<-c("DMFS_10y_time","DMFS_10y_event","PAM50","IBE","IBD")
Sample.Meta.data$DMFS_10y_time <-  as.numeric(as.character(Sample.Meta.data$DMFS_10y_time))
save(Expression.Data,Gene.Meta.data,P.metagene,Sample.Meta.data,file="./2 DATA/LM.BRCA/LM.Dataset.split.fixed.Rdata")
