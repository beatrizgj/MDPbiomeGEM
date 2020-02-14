# Bifidobacterium-Faecalibacterium data from MMODES Config
library(phyloseq)

# read otu and mapping table
phyloObject_datafile <- paste('data.normAndDist_definitiveClustering_BifFae.RData', sep = "")
load(phyloObject_datafile)

# Rename clusters for clarity
levels(sample_data(data.norm)$cluster) <- c("dysbiosis","risky","healthy")

# the set of perturbations
# rename perturbations field, for clarity
sample_data(data.norm)$diet <- as.character(sample_data(data.norm)$Perturbations)
Perturbations <- c('diet') # Only 1 perturbation, with different values

# sample attribute for time points
stepVar <- "time"

# sample attribute for subject 
subjectVar <- "subject"

# Associate perturbation to the sample it results in, rather than to the state where it is applied. To be the same than in experimental sampling, where the mapping value of the perturbations is associated to the sample after the perturbation was applied.
for(subject in unique(get_variable(data.norm,subjectVar))){
  subject.data <- phyloSubset(data.norm, subjectVar, subject)
  vectPert <- NULL
  vectPert <- get_variable(subject.data,'diet')
  newVecPert <- c('NA',vectPert[1:(length(vectPert)-1)])
  sample_data(data.norm)[sample_names(subject.data), "diet"] <- newVecPert
} # end-for move perturbation

concMetabolite <- function(phyloObj,subjectVar, subject, met){
  subject.data <- phyloSubset(phyloObj, subjectVar, subject)
  # concentrations <- get_variable(subject.data,met) # vector
  concentrations <- sample_data(subject.data)[,c(met)] # sample_data structure
  return(concentrations)
}

########################################################
### ANALYSIS OF METABOLOME CHANGES DUE TO PERTURBATIONS
########################################################
# Function get_iniMedium
  # What need to know the initial concentration of other metabolites different from butyrate (we do not know if their initial concentration is 0), to substract the initial concentration in the medium. Time=0 is not available (removed in MMODEStoMDPbiome procedure). So, to take from a text file, such as initial_medium.tsv, with column 1: name of metabolite and column 2: initial concentration in mM. 
get_iniMediumTable <- function(phylo, fileMedium='initial_concentration_metabolites.tsv'){
  firstMet=1
  lastMet=(grep('Perturbations',sample_variables(phylo)))-1
  metNameList=sample_variables(data.norm)[firstMet:lastMet]
  iniMediumTable_file=read.table(fileMedium,sep='\t',header=TRUE,stringsAsFactors=FALSE,row.names = 1)
  iniMediumTableTemp=data.frame(row.names=metNameList, check.names = FALSE, stringsAsFactors = FALSE)
  iniMediumTable=merge(iniMediumTableTemp,iniMediumTable_file,by="row.names",all.x=TRUE)
  iniMediumTable[is.na(iniMediumTable)] <- 0 # Fill-in with initial concentration=0 for metabolites not present in medium list
  row.names(iniMediumTable)=iniMediumTable[,1]
  iniMediumTable=subset(iniMediumTable,select=names(iniMediumTable)[2])
  return(iniMediumTable)
}

# Compute increase in fermentable products
initialMedium <- get_iniMediumTable(data.norm,'initial_concentration_metabolites_BifFae.tsv')
for (met in c('Butyrate','Lactate','Acetate')){
  metName <- sample_variables(data.norm)[grep(paste(met,".*_",sep=''),sample_variables(data.norm))]
  subjects <- unique(get_variable(data.norm,'subject'))
  for (isubject in subjects){
    vecConc <- concMetabolite(data.norm,subjectVar,isubject,metName)
    newVec <- as.numeric(rep(0,nsamples(vecConc)))
    for(pos in 2:(nsamples(vecConc))){
      newVec[pos] <- get_variable(vecConc[pos],metName) - get_variable(vecConc[pos-1],metName)
    } # end-for pos
    #print(vecConc)
    #print(newVec)
    sample_data(data.norm)[sample_names(vecConc), paste(met,'Increase',sep='')] <- newVec
  } # end-for isubject
  # Print computed values about change in metabolite concentration
  for (state in levels(sample_data(data.norm)$cluster)){
    #ss=subset_samples(data.norm,(cluster==state))
    ss=subset_samples(subset_samples(data.norm,(cluster==state)),diet!='NA')
    print(ss)
    #avg=mean(get_variable(ss,paste(met,'Increase',sep='')))
    all=get_variable(ss,paste(met,'Increase',sep=''))
    allPos=all[all>0]
    allNeg=all[all<0]
    avg=(mean(allPos)*length(allPos)/length(all))+(mean(allNeg)*length(allNeg)/length(all))
    print(paste(met,state,':',avg))
  } # end-for state
} # end-for met


# Compute ratio increase/decrease in ALL metabolites
initialMedium <- get_iniMediumTable(data.norm,'initial_concentration_metabolites_BifFae.tsv')
firstMet=1
lastMet=(grep('Perturbations',sample_variables(data.norm)))-1 # Compute position last metabolite in sample_variables
lastMet=lastMet-1 # The last one is Unnamed..129, giving error
for (posMet in seq(firstMet,lastMet)){
  metName <- sample_variables(data.norm)[posMet]
  print(metName)
  subjects <- unique(get_variable(data.norm,'subject'))
  for (isubject in subjects){
    vecConc <- concMetabolite(data.norm,subjectVar,isubject,metName)
    vecAbsInc <- as.numeric(rep(0,nsamples(vecConc)))
    vecRatio <- as.numeric(rep(0,nsamples(vecConc)))
    for(pos in 2:(nsamples(vecConc))){
      oldConc <- get_variable(vecConc[pos-1],metName)
      vecAbsInc[pos] <- get_variable(vecConc[pos],metName) - oldConc
      # When the concentration is the same than before: ratio=1 (conc_t1=conc_t2)
      # When previous oldConc=0: ratio=2.5 (it really increases) 
      if(oldConc==0){
        vecRatio[pos] <- 2.5 
      }else{
        vecRatio[pos] <- get_variable(vecConc[pos],metName)/oldConc
        if(vecRatio[pos]>2.5){ # To define a maximum
          vecRatio[pos] <- 2.5
        } 
      } # if-oldConc
    } # end-for pos
    sample_data(data.norm)[sample_names(vecConc), paste('AbsIncrease_',metName,sep='')] <- vecAbsInc
    sample_data(data.norm)[sample_names(vecConc), paste('RatioIncrease_',metName,sep='')] <- vecRatio
  } # end-for isubject
  save(data.norm,file='data.norm_withAbsAndRatioIncreaseMetabolites.Rdata')
} # end-for met
# Write metabolite increase (abs and ratio) per sample
namesVar=sample_variables(data.norm)[grep('^AbsIncrease_*',sample_variables(data.norm))]
tabMetAbsInc=sample_data(data.norm)[,namesVar]
write.table(tabMetAbsInc,file='table_metaboliteIncrease_perSample_Abs.tsv',sep='\t',col.names=NA)
namesVar=sample_variables(data.norm)[grep('^RatioIncrease_*',sample_variables(data.norm))]
tabMetRatioInc=sample_data(data.norm)[,namesVar]
write.table(tabMetRatioInc,file='table_metaboliteIncrease_perSample_Ratio.tsv',sep='\t',col.names=NA)


# Compute averages (per perturbation and/or state) and write tables
dfAbsIncreasePert=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfRatioIncreasePert=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfAbsIncreaseState=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfRatioIncreaseState=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfAbsIncreaseStatePert=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfRatioIncreaseStatePert=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
dfRatioIncreaseList=list()
for(countList in 1:length(levels(sample_data(data.norm)$cluster))){
  dfRatioIncreaseList[[countList]]=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])}
for (posMet in seq(firstMet,lastMet)){
  metName <- sample_variables(data.norm)[posMet]
  # Print computed values about ratio change in metabolite concentration, per PERTURBATION
  for (pert in unique(sample_data(data.norm)$diet)){
    if(pert!='NA'){
      ss=subset_samples(subset_samples(data.norm,(diet==pert)),diet!='NA')
      avg=mean(get_variable(ss,paste('AbsIncrease_',metName,sep='')))
      dfAbsIncreasePert[metName,pert]=avg
      avg=mean(get_variable(ss,paste('RatioIncrease_',metName,sep='')))
      dfRatioIncreasePert[metName,pert]=avg
      #print(paste(metName,pert,':',avg))
    } # end-if no-NA
  } # end-for pert
  # Print computed values about ratio change in metabolite concentration, per STATE
  for (state in levels(sample_data(data.norm)$cluster)){
    #ss=subset_samples(data.norm,cluster==state)
    ss=subset_samples(subset_samples(data.norm,(cluster==state)),diet!='NA')
    avg=mean(get_variable(ss,paste('AbsIncrease_',metName,sep='')))
    dfAbsIncreaseState[metName,state]=avg
    avg=mean(get_variable(ss,paste('RatioIncrease_',metName,sep='')))
    dfRatioIncreaseState[metName,state]=avg
  } # end-for state
  # Print computed values about ratio change in metabolite concentration, per STATE AND PERTURBATION
  for (pert in unique(sample_data(data.norm)$diet)){
    if(pert!='NA'){
      #ss=subset_samples(data.norm,diet==pert)
      ssp=subset_samples(subset_samples(data.norm,(diet==pert)),diet!='NA')
      for (state in levels(sample_data(ssp)$cluster)){
        label=paste(pert,state,sep='_')
        ss=subset_samples(ssp,cluster==state)
        print(label)
        avg=mean(get_variable(ss,paste('AbsIncrease_',metName,sep='')))
        dfAbsIncreaseStatePert[metName,label]=avg
        avg=mean(get_variable(ss,paste('RatioIncrease_',metName,sep='')))
        dfRatioIncreaseStatePert[metName,label]=avg
      } # end-for state
    } # end-if no-NA
  } # end-for pert
  # Print computed values about ratio change in metabolite concentration, per PERTURBATION, BY INDEPENDENT STATE
  countList=1
  for (state in levels(sample_data(data.norm)$cluster)){
    ssp=subset_samples(subset_samples(data.norm,(cluster==state)),diet!='NA')
    for (pert in unique(sample_data(data.norm)$diet)){
      if(pert!='NA'){
        ss=subset_samples(subset_samples(ssp,(diet==pert)),diet!='NA')
        avg=mean(get_variable(ss,paste('RatioIncrease_',metName,sep='')))
        dfRatioIncreaseList[[countList]][metName,pert] = avg
      } # end-for pert
    } # end-if no-NA
    countList=countList+1
  } # end-for state
} # end-for met
write.table(dfAbsIncreasePert,file='dfMetabolitesAbsIncreasePert.tsv',sep='\t',col.names=NA)
write.table(dfRatioIncreasePert,file='dfMetabolitesRatioIncreasePert.tsv',sep='\t',col.names=NA)
write.table(dfAbsIncreaseState,file='dfMetabolitesAbsIncreaseState.tsv',sep='\t',col.names=NA)
write.table(dfRatioIncreaseState,file='dfMetabolitesRatioIncreaseState.tsv',sep='\t',col.names=NA)
write.table(dfAbsIncreaseStatePert,file='dfMetabolitesAbsIncreaseStatePert.tsv',sep='\t',col.names=NA)
write.table(dfRatioIncreaseStatePert,file='dfMetabolitesRatioIncreaseStatePert.tsv',sep='\t',col.names=NA)
countList=1
for (state in levels(sample_data(data.norm)$cluster)){
  write.table(dfRatioIncreaseList[[countList]],file=paste('dfMetabolitesRatioIncreasePert_',state,'.tsv',sep=''),sep='\t',col.names=NA)
  countList=countList+1  
} # end-for
save(dfAbsIncreasePert,dfAbsIncreaseState,dfAbsIncreaseStatePert,dfRatioIncreasePert,dfRatioIncreaseState,dfRatioIncreaseStatePert,dfRatioIncreaseList,file='dfPerPertAndOrState.Rdata')



# PLOT HEATMAPS
############################################################################
### Function plot.heatmap.metabolite.x.perturbation
############################################################################
# Plot a heatmap with metabolite in rows and perturbations in columns, with
# the ratio average of increase of concentration of metabolite.
# Args:
#   dfIn: data frame with average ratios
#   metDelete: vector of metabolite names to not shown in plot
#   title: plot title
plot.heatmap.metabolite.x.perturbation  <- function(dfIn, metDelete,title){
  # Colors: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
  library(reshape2)
  library(viridis)
  # Remove metabolites
  dfHeatmap=data.frame()
  dfHeatmap=dfIn[!(row.names(dfIn) %in% metDelete),]
  # Rename metabolites
  dfHeatmap$met=unlist(lapply(rownames(dfHeatmap), function(x) gsub('_.*','',x)))
  df=melt(dfHeatmap,id=6,measure.vars=1:5)
  colnames(df)=c('metabolite','perturbation','IncConcRatio')
  pdf(paste('heatmap_',title,'.pdf',sep=''),width=3.5,height=7)
  p <- ggplot(df, aes(x = perturbation, y = metabolite, fill = IncConcRatio)) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE,option='plasma',direction=-1) +
    #scale_fill_gradient(low="white", high="blue") +
    theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
    theme(legend.position="top") +
    ggtitle(title)
  print(p)
  dev.off()
} # end-function plot.heatmap.metabolite.x.perturbation

############################################################################
### Function plot.heatmap.metabolite.x.perturbation.horizontal
############################################################################
# The same as plot.heatmap.metabolite.x.perturbation but in horizontal
plot.heatmap.metabolite.x.perturbation.horizontal  <- function(dfIn, metDelete,title,legendPosition='top'){
  # Colors: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
  library(reshape2)
  library(viridis)
  # Remove metabolites
  dfHeatmap=data.frame()
  dfHeatmap=dfIn[!(row.names(dfIn) %in% metDelete),]
  # Rename metabolites
  dfHeatmap$met=unlist(lapply(rownames(dfHeatmap), function(x) gsub('_.*','',x)))
  df=melt(dfHeatmap,id=6,measure.vars=1:5)
  colnames(df)=c('metabolite','perturbation','IncConcRatio')
  pdf(paste('heatmap_',title,'_horizontal.pdf',sep=''),width=7,height=3.5)
  p <- ggplot(df, aes(x = perturbation, y = metabolite, fill = IncConcRatio)) +
    geom_tile() +
    coord_flip() +
    scale_fill_viridis(discrete=FALSE,option='plasma',direction=-1) +
    #scale_fill_gradient(low="white", high="blue") +
    theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
    theme(legend.position=legendPosition) +
    ggtitle(title)
  print(p)
  dev.off()
} # end-function plot.heatmap.metabolite.x.perturbation.horizontal


# 1.- Filtering:
# 1.1.- Metabolites with 0 concentration in all samples --> 2.5 ratio
met0Conc=c()
for (posMet in seq(firstMet,lastMet)){
  metName <- sample_variables(data.norm)[posMet]
  if(sum(get_variable(data.norm, metName))==0){
    met0Conc=c(met0Conc,metName)
  }
} # 72 metabolites
# 1.2.- Metabolites with the same concentration in all samples --> 1 ratio
metConstant=c()
for (posMet in seq(firstMet,lastMet)){
  metName <- sample_variables(data.norm)[posMet]
  if(mean(get_variable(data.norm, metName))==get_variable(data.norm,metName)[1]&&(sd(get_variable(data.norm, metName))==0)){
    metConstant=c(metConstant,metName)
  }
} # 89, including also met0Conc
metDeleteHeatmap=union(met0Conc,metConstant)
metName <- sample_variables(data.norm)[grep(paste(met,".*_",sep=''),sample_variables(data.norm))]

# All samples
plot.heatmap.metabolite.x.perturbation(dfRatioIncreasePert, metDeleteHeatmap,'AllSamples')
plot.heatmap.metabolite.x.perturbation(dfRatioIncreaseList[[1]], metDeleteHeatmap,'Dysbiosis')
plot.heatmap.metabolite.x.perturbation(dfRatioIncreaseList[[2]], metDeleteHeatmap,'Risky')
plot.heatmap.metabolite.x.perturbation(dfRatioIncreaseList[[3]], metDeleteHeatmap,'Healthy')

colnames(dfRatioIncreasePert)=c('Inulin','FOS','Starch','Pro-high','Pro-low','met')
plot.heatmap.metabolite.x.perturbation.horizontal(dfRatioIncreasePert, metDeleteHeatmap,'',legendPosition='bottom')
colnames(dfRatioIncreaseList[[1]])=c('Inulin','FOS','Starch','Pro-high','Pro-low')
colnames(dfRatioIncreaseList[[2]])=c('Inulin','FOS','Starch','Pro-high','Pro-low')
colnames(dfRatioIncreaseList[[3]])=c('Inulin','FOS','Starch','Pro-high','Pro-low')
plot.heatmap.metabolite.x.perturbation.horizontal(dfRatioIncreaseList[[1]], metDeleteHeatmap,'A   Dysbiosis',legendPosition='bottom')
plot.heatmap.metabolite.x.perturbation.horizontal(dfRatioIncreaseList[[2]], metDeleteHeatmap,'C   Risky',legendPosition='bottom')
plot.heatmap.metabolite.x.perturbation.horizontal(dfRatioIncreaseList[[3]], metDeleteHeatmap,'B   Healthy',legendPosition='bottom')




#################################################
### FIGURES OF COMPARISON WITH Bauer et al. 2018
#################################################
# Plot to compare experimental vs simulated data
# Fig. 3b from Bauer2018: barplot experimental vs simulated

############################################################################
### Function get.mean.met.conc.state
############################################################################
# Get average concentration of a given metabolite in all samples of a 
# given state (not ratio, but concentration in mM)
# Args:
#   physeq: phyloseq object with all samples
#   state: state/cluster to select the samples (e.g. dysbiosis, healthy, etc.)
#   met: metabolite simple name, not complete name with formula
get.mean.met.conc.state <- function(physeq,state,met){
  avg=0
  #ss=subset_samples(physeq,cluster==state)  # subset_samples does not work in a loop (https://github.com/joey711/phyloseq/issues/752)
  ss=prune_samples(sample_data(physeq)$cluster==state,physeq)
  metName <- sample_variables(ss)[grep(paste('^',met,sep=''),sample_variables(ss))]
  if(length(metName)>0){
    avg=mean(get_variable(ss,metName))
  }
  return(avg)
} # end-funcion get.mean.met.conc.state

array3d <- array(0,dim=c(5,3,3),dimnames=list(list('Butyrate','Propionate','Isobutyrate','L.Lactate','Acetate'), list('exp','Bauer','MDPbiomeGEM'), list('control','case','ratio')))

# 1) Experimental values (metabolite concentrations) from: Table2 (control) and Table6/column 2 (disease) from Hove, H. & Mortensen, P. B. Influence of intestinal inflammation (IBD) and small and large bowel length on fecal short-chain fatty acids and lactate. Dig. Dis. Sci. 40, 1372–1380 (1995).
array3d[,'exp','control']=c(11,15,2,2,60)
array3d[,'exp','case']=c(7,11,1,6,38)
array3d[,'exp','ratio']=array3d[,'exp','control']/array3d[,'exp','case']

# 2) Bauer2018 simulated data from Fig.3b, grey bar, by-hand
# Only ratio available.
array3d[,'Bauer','ratio']=c(3.45,3.51,1.45,0.44,1.53)

# 3) MDPbiomeGEM simulated data
for(met in c('Butyrate','Propionate','Isobutyrate','L.Lactate','Acetate')){
  array3d[met,'MDPbiomeGEM','control']=get.mean.met.conc.state(data.norm,'healthy',met)
  array3d[met,'MDPbiomeGEM','case']=get.mean.met.conc.state(data.norm,'dysbiosis',met)
} # end-for
array3d[,'MDPbiomeGEM','ratio']=array3d[,'MDPbiomeGEM','control']/array3d[,'MDPbiomeGEM','case']

df=melt(array3d[,,'ratio'])
colnames(df)=c('metabolite','source','value')
pdf('barplot_expVsSimulated_MetRatioCtrlCase.pdf',height=3.5)
p <- ggplot(df, aes(x = metabolite, y = value, fill = source)) +
  geom_bar(stat='identity', position='dodge', color='black') +
  scale_fill_manual(values = c("#ff9900", "grey", "blue3")) +
  geom_hline(yintercept = 1, linetype=2, size=1) + 
  ylab('Controls to CD patients ratio') +
  theme_bw() +
  theme(legend.position="top", legend.direction="horizontal",  legend.title = element_blank())
print(p)
dev.off()

# Only butyrate and Acetate
df2=subset(df,metabolite=='Butyrate' | metabolite=='Acetate')
pdf('barplot_expVsSimulated_MetRatioCtrlCase_onlyButAct.pdf',height=3.5,width=3)
p2 <- ggplot(df2, aes(x = metabolite, y = value, fill = source)) +
  geom_bar(stat='identity', position='dodge', color='black') +
  scale_fill_manual(values = c("#ff9900", "grey", "blue3")) +
  geom_hline(yintercept = 1, linetype=2, size=1) + 
  ylab('Controls to CD patients ratio') +
  #ylab('Healthy to Dysbiosis patients ratio') +
  theme_bw() +
  theme(legend.position="top", legend.direction="horizontal",  legend.title = element_blank())
print(p2)
dev.off()


#############################
# Plot to compare metabolite concentration distribution, Bauer vs MDPbiomeGEM, per state and treatment
# boxplot
# To compare with Fig.2e-5e Bauer2018 (concentration (not ratio) in dysbiosis(not inulin), healthy and treatment (dysbiosis-inulin))

# To split which ones are following the treatment (diet='Inulin', cluster='dysbiosis') from those ones do not following it (diet!='Inulin', cluster='dysbiosis')
data.norm.temp=data.norm
ss=subset_samples(subset_samples(data.norm.temp,cluster=='dysbiosis'),diet=='Inulin')
ids=get_variable(ss,'id')
levels(sample_data(data.norm.temp)$cluster)=c(levels(sample_data(data.norm.temp)$cluster),"treatment")
sample_data(data.norm.temp)[sample_names(ss), "cluster"] <- rep("treatment",nsamples(ss))
# Filter out metabolites to not show in the distribution: Constant and Conc=0. Saved above in 'metDeleteHeatmap'
metNames <- sample_variables(data.norm.temp)[seq(firstMet,lastMet)]
metNames <- metNames[!(metNames %in% metDeleteHeatmap)]
data.norm.temp=subset_samples(data.norm.temp,cluster!='risky')
df=as.data.frame(sample_data(data.norm.temp))
df=subset(df,select=c(metNames,'cluster'))
names(df)=unlist(lapply(names(df), function(x) gsub('_.*','',x))) # Rename metabolites
# Convert to ggplot
dfBoxplot=melt(df,id='cluster')
rm(data.norm.temp,df)
colnames(dfBoxplot)=c('state','metabolite','concentration')
stateOrder=c("dysbiosis","treatment","healthy")
dfBoxplot[,'state'] <- factor(dfBoxplot[,'state'], levels = stateOrder)


pdf('boxplot_metabolitePerStateAndTreatment_MDPbiomeGEM.pdf',width=14)
p <- ggplot(dfBoxplot, aes(x=metabolite, y=concentration, color=state)) + 
  geom_boxplot(lwd=0.5, outlier.size = 0.3) +
  ylim(0,5) + # limit to metabolite concentration < 5 mM
  scale_color_manual(values = c("red", "purple", "dodgerblue2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=16), legend.position=c(0.95,0.91), legend.title = element_blank())
print(p)
dev.off()

dfBoxplot3=dfBoxplot
dfBoxplot3[,'concentration'] = log(1+dfBoxplot[,'concentration'])
pdf('boxplot_metabolitePerStateAndTreatment_MDPbiomeGEM_logPlus1.pdf',width=14)
p <- ggplot(dfBoxplot3, aes(x=metabolite, y=concentration, color=state)) + 
  geom_boxplot(lwd=0.8, outlier.size = 1.2) +
  ylim(0,1) + # limit to metabolite concentration < 1 mM
  scale_color_manual(values = c("red", "purple", "dodgerblue2")) +
  labs(y="log(concentration)+1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=24), legend.position=c(0.93,0.89), legend.title = element_blank())
print(p)
dev.off()


# Select a subset of metabolites
dfBoxplot2=subset(dfBoxplot,metabolite %in% c('Butyrate','Acetate'))
pdf('boxplot_metabolitePerStateAndTreatment_MDPbiomeGEM_only2met.pdf',width=5)
p <- ggplot(dfBoxplot2, aes(x=metabolite, y=concentration, color=state)) + 
  geom_boxplot(lwd=1) +
  #ylim(0,5) +
  scale_color_manual(values = c("red", "purple", "dodgerblue2")) +
  theme_bw() +
  theme(text = element_text(size=24), legend.position=c(0.8,0.91), legend.title = element_blank())
print(p)
dev.off()


# To compute the ratio of increase per independent state
# Print computed values about ratio change in metabolite concentration, per perturbation
for (state in levels(sample_data(data.norm)$cluster)){
  ss=subset_samples(data.norm,cluster==state)
  dfRatioIncreaseState=data.frame(row.names=sample_variables(data.norm)[firstMet:lastMet])
  for (posMet in seq(firstMet,lastMet-1)){
    # (lastMet-1): Because there isn't RatioIncrease_Unnamed..129 (the last attribute)
    metName <- sample_variables(data.norm)[posMet]
    avg=mean(get_variable(ss,paste('RatioIncrease_',metName,sep='')))
    dfRatioIncreaseState[metName,'RatioIncrease']=avg
    #print(paste(metName,pert,':',avg))
  } # end-for met
  write.table(dfRatioIncreaseState,file=paste('dfRatioIncreaseState_',state,'.tsv',sep=''),sep='\t',col.names=NA)
} # end-for state

