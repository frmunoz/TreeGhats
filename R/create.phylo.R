create.phylo <- function(names = NULL, scenarios = "S3")
# The function creates a phylogenic tree of taxa based on species, genus and APGIII
# family information
# names is a table with columns "Binome", "Genus" and "Family_APGIII" 
{
  if(!is.null(names))
  {
    if(any(!sapply(c("Binome","Genus","Family_APGIII"),function(x) x%in%colnames(names)))) stop("Input names table must include Binome, Genus and Family_APGIII columns")
    if(!any(sapply(names$Binome,function(x) length(strsplit(x,split=" ")[[1]]))>1)) stop("Genus and species names to be separated with space")
  } else names <- reftaxo;
  
  require(phytools)
  
  data(qian)
  #qian$phylo.all$tip.label <- unlist(sapply(qian$phylo.all$tip.label,function(x) gsub(x,pattern="-",replacement="")))
  data(reftaxo.phylo)
  
  phylomaker <- names[,c("Binome", "Genus", "Family_APGIII")]
  colnames(phylomaker) <- c("species","genus","family");
  rownames(phylomaker) <- capitalize(tolower(gsub(phylomaker$species,pattern=" ",replacement="_")))

  
  # For diagnostic only
  #sum(tolower(unique(reftaxo$Family_APGIII))%in%tolower(nodes$family))/length(unique(reftaxo$Family_APGIII))  # 99.2%
  #unique(reftaxo$Family_APGIII)[!tolower(unique(reftaxo$Family_APGIII))%in%tolower(nodes$family)] # Centroplacaceae
  #sum(unique(reftaxo$Genus)%in%nodes$genus)/length(unique(reftaxo$Genus))  # 79.9%
  #sum(unique(rownames(reftaxo))%in%tolower(phylo.all.qian$tip.label))/length(unique(rownames(reftaxo)))  # 18.1%
  
  # With whole Qian phylogeny
  #result<-S.PhyloMaker(splist=phylomaker, tree=qian$phylo.all , nodes=qian$nodes, scenarios = scenarios)      # run the function S.PhyloMaker
  # With phylogeny generated for WG trees
  result<-S.PhyloMaker(splist=phylomaker[rownames(phylomaker)%in%reftaxo.phylo$tip.label,], tree=reftaxo.phylo, nodes=qian$nodes, scenarios = scenarios)      # run the function S.PhyloMaker
  
  return(result)
}