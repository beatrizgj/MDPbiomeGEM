library(plyr) # ldply, to join several data.frames
library(tidyverse)
library(stringr) # str_pad
library(phyloseq)


#####################################
### Function complete.tax.table
#####################################
# To complete tax_table with the unknown upper categories name. Using taxize library.
# Input parameters:
#   taxmat: taxonomy matriz, with names of known taxa
#   uniqueKnownTaxaLevel: the level of the taxonomy of the values given.
complete.tax.table <- function(taxmat,uniqueKnownTaxaLevel="Species")
{
  colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # To complete the remainder parent ranks of the taxonomy classification
  # To parser some of the hand-written class rank taxonomy, to make possible to find their parent taxonomy automatically
  taxmat[,uniqueKnownTaxaLevel] <- unlist(lapply(rownames(taxmat),function(y) gsub("*_no_class","\\",y)))
  require(taxize)
  for (i in 1:nrow(taxmat)){
    # With taxmat()
    # rows=1: To select the first result if there are several ones. ask=FALSE, not interactive for selecting in results.
    if(!is.na(get_uid(taxmat[i,uniqueKnownTaxaLevel],verbose=FALSE,rows=1,ask=FALSE)[1])){
      ## Phylum
      out.phylum <- tax_name(query=taxmat[i,uniqueKnownTaxaLevel],get='Phylum',verbose=FALSE,ask=FALSE)$phylum
      # If null, second possibility to compute phylum, with rows=1:  First, it is better without rows=1, to avoid errors, for example, in bacilli, returning 'Nematoda' instead of 'Firmicutes'. If 'NA', it means several rows, and someone should be chosen, for example, the first one. The rows=1 is neccesary when the class level is really a phylum level, for example, with 'Actinobacteria', when several valid rows are returned.
      if(is.null(out.phylum)){
        out.phylum <- tax_name(query=taxmat[i,uniqueKnownTaxaLevel],get='Phylum',verbose=FALSE,ask=FALSE,rows=1)$phylum  
      }  
      if(!(is.null(out.phylum))){
        taxmat[i,"Phylum"] <- out.phylum
      }
      ## Kingdom
      out.kingdom <- tax_name(query=taxmat[i,uniqueKnownTaxaLevel],get='Kingdom',verbose=FALSE,ask=FALSE)$kingdom
      # If null, second possibility to compute phylum, with rows=1; as in Phylum.
      if(is.null(out.kingdom)){
        out.kingdom <- tax_name(query=taxmat[i,uniqueKnownTaxaLevel],get='Kingdom',verbose=FALSE,ask=FALSE,rows=1)$kingdom
      }
      if(!(is.null(out.kingdom))){
        taxmat[i,"Kingdom"] <- out.kingdom
      }
    } #if it exists UID
    # TO-IMPROVE: to generalize for all parent levels!!!
    # Also see classification()
    # Additional info: tax_rank(query="Mollicutes") --> rank: It returns the rank in the taxonomy
  }
  return(taxmat)
} # end-function complete.tax.table


#####################################
### Function process.one.biomass.file
#####################################
# a) To generate subjectID (suffix from file, between '.')
# b) To generate sample ID: subjectID.[0*]seqNum
process.one.biomass.file <- function(file){
  suffix = gsub('.tsv$','',gsub('^biomass_','',file))
  table=read.table(file,sep="\t",header=TRUE,row.names=1)
  id.subject=rep(suffix,nrow(table))
  id.num=seq(1,nrow(table),1)
  id.num=str_pad(id.num,5,pad='0') # Add leading zeros, to easy sort
  id.vector=paste(id.subject,id.num,sep='_')
  table$time=rownames(table)
  rownames(table)=id.vector
  table$subject=id.subject
  table$id=id.vector
  return(table)
} # end-function process.one.biomass.file

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
### Function rename.biomass.file
#####################################
# prefixExp: complete, including '.' or '_'. Ej: 'atr2.'
rename.biomass.file <- function(folder,prefixExp){
  suffix=gsub(prefixExp,'',folder)
  fin=paste(folder,'filtered_plot.tsv',sep='/')
  prefix=unlist(strsplit(prefixExp,'[.]'))
  fout=paste('biomass_',prefix,'_',suffix,'.tsv',sep='')
  system(paste("cp -p",fin,fout))
} # end-function rename.biomass.file

#################################################################
# END FUNCTIONS
#################################################################
