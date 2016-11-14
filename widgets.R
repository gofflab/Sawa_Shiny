###################
# Utility Functions
###################
hclust2 <- function(x, method="ward", ...)
  hclust(x, method="ward.D", ...)

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  res <- unique(res)
  res
}

meltCDS<-function(cds,geneset,logMode=F){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-exprs(sub)
  if(logMode){
    sub.expr<-log10(sub.expr+1)
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","Cell_ID","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="Cell_ID",by.y="Cell_ID")
  res<-merge(res,fData(sub),by.x="gene_id",by.y="gene_id")
  res
}

myBarMap<-function(cds,geneset,facet_by="Treatment",color_by="Treatment",cluster="both",showSummary=T,scale="free",space="free_x",...){
  sub.melt<-meltCDS(cds,geneset,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  sub.melt.summary<-sub.melt %>%
    dplyr::group_by_(.dots=c("gene_short_name",facet_by_melt)) %>%
    dplyr::summarise(mean=mean(value),median=median(value),sd=sd(value),upper_bound=mean+sd,lower_bound=max(mean-sd,0))
  
  if(cluster %in% c("row","both",T)){
    sub.sum.mat<-sub.melt.summary %>%
      recast(as.formula(paste("gene_short_name ~",facet_by)),measure.var="mean",fun.aggregate=mean)
    sub.sum.hclust<-hclust2(dist(sub.sum.mat[,-1]))
    gene.order.idx<-order.dendrogram(as.dendrogram(sub.sum.hclust))
    gene.order<-sub.sum.mat$gene_short_name[gene.order.idx]
    sub.melt$gene_short_name<-factor(sub.melt$gene_short_name, levels=gene.order)
  }
  
  if(cluster %in% c("column","both",T)){
    sub.mat<-sub.melt %>%
      recast(as.formula("gene_short_name ~ Cell_ID"),measure.var="value",fun.aggregate=mean)
    sub.hclust<-hclust2(dist(t(sub.mat[,-1])))
    cell.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
    cell.order<-colnames(sub.mat[,-1])[cell.order.idx]
    #print(cell.order)
    sub.melt$Cell_ID<-factor(sub.melt$Cell_ID,levels=cell.order)
  }
  
  p<-ggplot(sub.melt)
  p<-p + geom_bar(aes_string(x="Cell_ID",y="value",fill=color_by,color=color_by),stat="identity")
  
  if(showSummary){
    p<-p + geom_hline(aes(yintercept=mean),data=sub.melt.summary,size=1.0)
    p<-p + geom_hline(aes(yintercept=upper_bound),data=sub.melt.summary,linetype="dashed")
    p<-p + geom_hline(aes(yintercept=lower_bound),data=sub.melt.summary,linetype="dashed")
  }
  p<-p +
    facet_grid(as.formula(paste("gene_short_name ~", facet_by)),scale=scale,space=space,labeller=labeller(.default=label_both,gene_short_name=label_value)) +
    theme_bw() + guides(color=FALSE) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.margin = unit(0, "lines"),
          panel.grid = element_blank()
    )
  p
}

mySummaryBarPlot<-function(cds,geneset,metric="mean",facet_by="Treatment",color_by="Treatment",cluster="both",showSummary=T,scale="free",space="free_x",...){
  sub.melt<-meltCDS(cds,geneset,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  sub.melt.summary<-sub.melt %>%
    dplyr::group_by_(.dots=c("gene_short_name",facet_by_melt)) %>%
    dplyr::summarise(mean=mean(value),median=median(value),sd=sd(value),upper_bound=mean+sd,lower_bound=max(mean-sd,0))
  
  p<-ggplot(sub.melt.summary)
  p<-p + geom_bar(aes_string(x="gene_short_name",group=color_by,y=metric,fill=color_by),position="dodge",stat="identity") +
    geom_errorbar(aes_string(x="gene_short_name",group=color_by,ymax="upper_bound",ymin="lower_bound"),position="dodge")
  
  p + monocle:::monocle_theme_opts()
}


myTSNEPlot<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("Cell_ID","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="Cell_ID")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,size="value")) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,size="value",shape=shape_by)) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}

myTSNEPlotAlpha<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=2){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes)
    genes<-melt(genes)
    colnames(genes)<-c("Cell_ID","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="Cell_ID")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value"),stroke=0,size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value",stroke=0,shape=shape_by),size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}

############
# Color Palette
############
celltype_colors<-c("chocolate1","blue4")
label_colors<-c("darkgreen","red")