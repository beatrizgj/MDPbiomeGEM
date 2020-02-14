library(caret)
library(tidyverse)
library(dplyr)
library(randomForest)
library(fpc) # pamk
library(cluster) # pam
library(ggplot2,quietly=TRUE)
suppressMessages(library(ellipse))
suppressMessages(library(RColorBrewer))
library(reshape2,quietly=TRUE) # Load the reshape2 package (for the melt() function)

source('robust.clustering.metagenomics.functions.r')


#####################################
### Function process.one.fluxes.file
#####################################
# a) To generate subjectID (suffix from file, between 'biomass_' and '.tsv')
# b) To generate sample ID: subjectID_[0*]seqNum
process.one.fluxes.file <- function(file){
  suffix = gsub('.tsv$','',gsub('^fluxes_','',file))
  table=read.table(file,sep="\t",header=TRUE)
  table[,grep('time',colnames(table))] <- NULL # Remove repeated time columns (from different strains)
  id.subject=rep(suffix,nrow(table))
  id.num=seq(1,nrow(table),1)
  id.num=str_pad(id.num,5, pad = '0') # Add leading zeros, to easy sort
  id.vector=paste(id.subject, id.num,sep = '_')
  rownames(table)=id.vector
  return(table)
} # end-function process.one.fluxes.file

#####################################
### Function rbind.with.rownames
#####################################
rbind.with.rownames <- function(datalist) {
  require(plyr)
  temp <- rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  return(temp)
} # end-function rbind.with.rownames

#####################################
### Function cluster names
#####################################
cluster_names <- function(df) {
  require(dplyr)
  means <- df %>% 
    select(starts_with("Biomass"), cluster) %>%
    group_by(cluster) %>%
    summarize_all(mean)
  signed <- c(as.character(means[which.max(rowSums(means[,2:ncol(means)])),]$cluster),
              as.character(means[which.min(rowSums(means[,2:ncol(means)])),]$cluster))
  df$cluster <- mapvalues(df$cluster, from = signed, to = c("growth", "stationary"))
  return(df)
} # end-function cluster names

#####################################
### Function rf.10CV
#####################################
# Build and evaluate (with 10 Cross-Validation) a model with Random Forest, using caret.
# Args:
#   df: data frame with all data (including predictors and class variable)
#   formula: object of class formula (i.e. class~ predictors (separated by '+'))
#   classVar: string/number of class variable in 'df'. Default: the last one.
#   ntree: number of trees in the forest
# Output:
#   finalModel: object with the model in randomForest format, for further evaluations or plots.
rf.10CV <- function(df,formula,classVar=length(df),ntree=1000){
  require(doMC) # parallel processing; in Windows, maybe doParallel would work
  registerDoMC(cores = 5)
  set.seed(123)
  train_control <- trainControl(method="cv", number=10, savePredictions = TRUE)
  model <- train(form=formula,data=df, trControl=train_control, method="rf", metric='Accuracy',
                 ntree=ntree, importance=TRUE, localImp=TRUE, na.action=na.omit,
                 allowParallel = T, ncores = 5, nodesize = 42)
  # Define output
  output <- list("model"=model$finalModel, "pred"=model$pred)
  # summarize results
  print(model)
  # plot model error
  print(plot(model))
  # Variable Importance
  pdf('varImportance_rf10CV.pdf')
  print(varImpPlot(model$finalModel))
  dev.off()
  print(plot(margin(model$finalModel,df[,classNum])))
  
  return(output)
}

####################################################
### plot.rf.var.importance.by.class.andMean.dotplot
####################################################
# Plot dotplot with variable importance mean over all classes
# Args:
#   model: random forest model already build
#   predVar: string of column ID with predictor/variables names values
#   classVar: string of class variable in 'df'
#   title: header of the plot
#   colorVector: vector of colors
#   nBestFeatures: number of top relevant features to show in the plot.
plot.rf.var.importance.by.class.andMean.dotplot <- function(model,predVar,classVar,title='',colorVector=NULL,nBestFeatures=NULL){
  library(reshape2)
  imp.df <- melt(importance(model)[,1:(length(model$classes)+1)])
  colnames(imp.df)=c(predVar,classVar,'value')
  # a.-Order features
  pred.order=names(sort(importance(model)[,'MeanDecreaseAccuracy'])) # My order according to global MeandDecreaseAccuracy
  imp.df[,predVar] <- factor(imp.df[,predVar], levels = pred.order)
  class.names=levels(imp.df[,classVar])
  levels(imp.df[,classVar]) <- c(class.names[1:(length(class.names)-1)],"MEAN")
  imp.df[,classVar] <- factor(imp.df[,classVar])
  # b.- Subset features to show
  if(!is.null(nBestFeatures)){
    imp.df=subset(imp.df,subset=(imp.df[,predVar] %in% tail(pred.order,n=nBestFeatures)))
  }
  p <- ggplot(imp.df, aes_string(x = 'value', y = predVar, group = predVar, colour = classVar)) +
    geom_segment(aes_string(yend=predVar), xend=0, colour="grey50") +
    geom_point( size = 1) +
    theme_bw() +
    facet_grid(reformulate(classVar)) +
    theme(panel.grid.major.y = element_blank()) +
    #theme(text = element_text(size=16)) +
    xlab(paste(predVar," importance",sep='')) +
    theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
    theme(legend.position="none") +
    ggtitle(title)
  if(!is.null(colorVector)){
    p +  scale_color_manual(values=colorVector)
  }else{
    p
  }
  return(p)
}

##########################################################################
# END FUNCTIONS
##########################################################################

# 1. READ DATA
# Read fluxes tables from MMODES
file <- list.files(pattern="^fluxes\\_.*\\.tsv$")
# Add subject and sample ID per table
tables.list=lapply(file, process.one.fluxes.file)
# Concatenate all fluxes data.frames (one per subject time series) in a unique one
table.all <- rbind.with.rownames(tables.list)
# Save
fluxes_raw  <- table.all
saveRDS(fluxes, "fluxes.rds")
write.table(fluxes,file='allTimeSeries_fluxes.tsv',sep='\t',quote=FALSE)
fluxes_raw <- readRDS("fluxes.rds")

# 1.2. FEATURE SELECTION
# Check columns with NA values
# fluxes.isna=colnames(fluxes)[colSums(is.na(fluxes)) > 0]
# NA values are generated if some Perturbations are missing in some simulation.
# We can afford to remove the rows of these simulations (generally 1 or 0 in the dataset).
fluxes <- fluxes_raw[complete.cases(fluxes_raw),]
# Check fluxes with constant values (sd=0)
output_sd=apply(fluxes, MARGIN=2, sd)
fluxes.sd0=colnames(fluxes)[output_sd==0]
# Remove those fluxes
if(length(fluxes.sd0)>0){
  fluxes = fluxes[, -which(names(fluxes) %in% fluxes.sd0)]
}
# Remove first all of the transporters and exchanges
# It does not remove knowledge, but interpretations of the inner metabolism is possible
fluxes <- fluxes %>%
  select(-starts_with("EX"), -starts_with("Ex"), -contains("tex_"), -contains("tpp_"), -contains("abcpp"))
# Some columns shares identical values across rows. That is because reactions are coupled,
# belonging to the same pathway
fluxes <- fluxes[,-which(duplicated(t(fluxes)) == TRUE)]
cat("Predictos variables reduced sized from", ncol(fluxes_raw), "to", ncol(fluxes), 
    "by feature selection.")

# normalization
library(scales)
fluxes.rescale <- apply(fluxes[,1:ncol(fluxes)-1], MARGIN = 2, rescale, to=c(-1,1))

# 2. RANDOM SUFFLE
iniSeed <- 1234
set.seed(iniSeed)
selected <- sample(nrow(fluxes.rescale))
fluxes.rescale.random <- fluxes.rescale[selected,] # for clustering
fluxes.random <- fluxes[selected,] # for random forest, without normalization

# 3. CLUSTERING
maxClus=10
eval.array2d <- array(0,dim=c(maxClus-1,3),dimnames=list(as.character(seq(2,maxClus,1)), list('SI','PS','Jaccard')))
# Silhouette
fitPamBest <- pamk(fluxes.rescale.random,krange=2:maxClus)
save(fitPamBest,file='fitPamBest.Rdata')
eval.array2d[,'SI']=fitPamBest$crit[2:maxClus]
# Prediction Strength
out.pred.str <- prediction.strength(fluxes.rescale.random, Gmin=2, Gmax=maxClus, M=50, clustermethod=claraCBI, classification="centroid")
eval.array2d[,'PS']=out.pred.str$mean.pred[2:maxClus]
# Jaccard
for(k in 2:maxClus){
  cf <- clusterboot(fluxes.rescale.random,B=100,bootmethod="boot",clustermethod=claraCBI,k=k,seed=iniSeed,count=FALSE)
  #print(mean(cf$bootmean))
  eval.array2d[as.character(k),'Jaccard'] <- mean(cf$bootmean)
}
save(eval.array2d,file='eval.arrays.RData')
plot.robust.clustering.PAM(eval.array2d)

# Select the best k (with the higher SI and (PS or Jac) > their thresholds)
kBest <- robust.clustering.decision(eval.array2d)
fit <- pam(fluxes.rescale.random,kBest)
print.clustering.results(eval.array2d,fit,kBest)

# Getting a list <sampleID,clusterID>
labels <-  as.data.frame(as.factor(fit$cluster))
colnames(labels) <- c('cluster')
df.fluxes <- as.data.frame(fluxes.random)
df <- merge(df.fluxes,labels,by='row.names')
row.names(df) <- df$Row.names
df$Row.names <- NULL
df.out <- subset(df,select=cluster)
write.table(df.out,'sampleId-cluster_pairs_fluxes.txt',quote=FALSE,sep=',',row.names=TRUE)
rm(df.fluxes,df.out)

# 4. RANDOM FOREST
# Print medoids with their values in the most relevant features
# Very important to intepretate the clusters!!
as_tibble(df[rownames(fit$medoids),])

# Supervised learning after feature selection
# Random forest with the clusters defined by PAM
formula <- formula("cluster~.")
model <- rf.10CV(df,formula)
saveRDS(model, "final_RF_bf_model.rsd")
sink('model_RF_bf.txt')
cat('Growth in medoids')
as_tibble(df[rownames(fit$medoids),]) %>% select(starts_with("Biomass"))
print(model)
sink()
mod <- model$model
pred <- model$pred
# fluxes importance independent by predicted class
plot.rf.var.importance.by.class.andMean.dotplot(mod,'flux','cluster',title='Most relevant fluxes',nBestFeatures=25)
# see predictions (and in which CV fold it was tested)

