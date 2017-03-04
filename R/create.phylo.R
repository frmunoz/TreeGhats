create.phylo <- function(names)
# The function creates a phylogeny of taxa based on species, genus and APGIII
# family information
# names is a table with columns "Binome", "Genus" and "Family_APGIII" 
{
  names <- tolower(names)
  if(!any(sapply(names,function(x) strsplit(x,split=" "))>1)) stop("Genus and species names to be separated with space")
  
  require(ape)
  require(Hmisc)
  require(phytools)
  
  #phylo.all.qian <-read.tree("C:/PhytoPhylo.tre") 
  data(phylo.all)
  phylo.all$tip.label <- unlist(sapply(phylo.all$tip.label,function(x) gsub(x,pattern="-",replacement="")))
  
  # Utilisation de l'algorithme de Qian and Jin (2015)
  #source("D:/François/Documents de travail/Thèses/Ruksan Bose/S.PhyloMaker.R")
  
  #nodes <- read.csv("C:/nodes.csv",header=T,stringsAsFactors =F)     # read in the nodes information of the megaphylogeny.
  data(nodes)
  
  # For diagnostic only
  #sum(tolower(unique(reftaxo$Family_APGIII))%in%tolower(nodes$family))/length(unique(reftaxo$Family_APGIII))  # 99.2%
  #unique(reftaxo$Family_APGIII)[!tolower(unique(reftaxo$Family_APGIII))%in%tolower(nodes$family)] # Centroplacaceae
  #sum(unique(reftaxo$Genus)%in%nodes$genus)/length(unique(reftaxo$Genus))  # 79.9%
  #sum(unique(rownames(reftaxo))%in%tolower(phylo.all.qian$tip.label))/length(unique(rownames(reftaxo)))  # 18.1%
  
  #phylomaker <- cbind(gsub(capitalize(names$Binome),pattern="_",replacement=" "),sapply(rownames(reftaxo),function(x) capitalize(reftaxo[x,"Genus"])),sapply(rownames(reftaxo),function(x) capitalize(reftaxo[x,"Family_APGIII"])));
  phylomaker <- cbind(gsub(capitalize(names$Binome),pattern="_",replacement=" "),sapply(rownames(reftaxo),function(x) capitalize(reftaxo[x,"Genus"])),sapply(rownames(reftaxo),function(x) capitalize(reftaxo[x,"Family_APGIII"])));
  colnames(phylomaker) <- c("species","genus","family"); phylomaker <- as.data.frame(phylomaker)
  
  result<-S.PhyloMaker(splist=phylomaker, tree=phylo.all.qian , nodes=nodes, scenarios = "S3")      # run the function S.PhyloMaker
  
  return(result$Scenario.3)
}