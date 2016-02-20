##module is a function to get the count value information , the group information and the module information.##
#'count' is a data frame or matrix for case and control, with rows as variables (genes) and columns as samples#
#'module' is A matrix of modules which consist of function-related gene lists.##
##'case' is the sample labels for one phenotype.
#'control' is the sample labels for another phenotype.#
DERM<- function(count, module, case, control) 
{
  # get the counts value from the data frame
  group.data  <- read.table(file=count,header=T,sep="\t")
  group.data1 <-  group.data[,c(case,control)] 
  counts  <- as.matrix(group.data1)
  Count <- as.numeric(counts)
  
  #get the group information from the data
  Group <- matrix(rep(1:0,c(length(case)*nrow(group.data),length(control)*nrow(group.data))),nrow=nrow(group.data),ncol=ncol(group.data1))
  Group <- as.factor(Group)
 
  #get the Module information and construct module
  module.data <- read.table(file=module,header=F,quote="")
  modulematrix <- matrix(nrow=nrow(group.data), ncol=nrow(module.data)+1)
  
  for (i in 1:nrow(group.data)) {
    gene <- as.character(group.data[i,1])
    
    for (j in 1:nrow(module.data)) {
      genelist <- as.character(module.data[j,2])
      genes <- strsplit(genelist, ";")
      sum <- 0
      
      for (k in 1:length(genes[[1]])) 
      {
        if (genes[[1]][k] == gene)
        {
          sum <- sum + 1 
        }
        if (sum == 1) 
        {
          next
        }
      }
      
      modulematrix[i,j] <- sum
      
    }
    
  }
  
  for (i in 1:nrow(group.data)) {
    if (sum(modulematrix[i,1:nrow(module.data)]) == 0) {
      modulematrix[i,nrow(module.data)+1] <- 1
    } else {
      modulematrix[i,nrow(module.data)+1] <- 0
    }  
  }
 
#################################################################################################
  # the GLM model.
  #################################################################################################
  #calculate p-values
  library(MASS)
  pvalue.CRC <- matrix(nrow=ncol(modulematrix),ncol=1)
  for (i in 1:(ncol(modulematrix)-1))  {
    Module <- as.matrix(modulematrix[,i])
    Module <- Module[,rep(1,ncol(count))]
    Module <- as.factor(Module)  
    fm1 <- glm.nb(Count ~ Group * Module )
    pvalue.CRC[i] <- summary(fm1)$coefficients[4,4]
  }
  
  other <- matrix(nrow=1, ncol=2)
  other[1] <- "other"
  allgenes <- data.frame(count[,1], modulematrix[,ncol(modulematrix)])
  colnames(allgenes) <- c("gene", "other")
  genes.other <- subset(allgenes, other == 1)
  genelist.other <- paste(sort(genes.other$gene), collapse=";")
  other[2] <- genelist.other
  
  data.module2 <- rbind(module, other)
  
  #adjust p-values:
  pvalue.fdr <- p.adjust(pvalue.CRC, method='fdr')
  
  output.lml <- data.frame(data.module2[,1], pvalue.CRC, pvalue.fdr, data.module2[,2])
  colnames(output.lml) <- c("Module", "pvalue.group", "pvalue.fdr", "Gene.Symbol")
  result <- output.lml[order(output.lml[,3]),]
  write.table(result, "output", quote=FALSE, sep="\t", row.names=FALSE)
}
