create.phylo <- function(names = NULL, scenarios = "S3")
  # The function creates a phylogenic tree of taxa based on species, genus and APGIII
  # family information
  # names is a table with columns "Binome", "Genus" and "Family_APGIII"
{
  data('TreeGhatsData', package='TreeGhats', envir=environment())
  TreeGhatsData <- get("TreeGhatsData", envir=environment())
  data('qian', package='TreeGhats', envir=environment())
  qian <- get("qian", envir=environment())
  
  if(!is.null(names))
  {
    if(!"Binome"%in%colnames(names) & "Genus"%in%colnames(names) & "Species"%in%colnames(names))
      names$Binome <- paste(names$Genus,names$Species,sep=" ")
    names$Binome[is.na(names$Genus) | is.na(names$Species)] <- NA
    if(any(!sapply(c("Binome","Genus","Family_APGIII"),function(x) x%in%colnames(names)))) stop("Input names table must include Binome, Genus and Family_APGIII columns")
    if(!any(sapply(names$Binome,function(x) length(strsplit(x,split=" ")[[1]]))>1)) stop("Genus and species names to be separated with space")
  } else names <- TreeGhatsData;
  
  #qian$phylo.all$tip.label <- unlist(sapply(qian$phylo.all$tip.label,function(x) gsub(x,pattern="-",replacement="")))
  
  phylomaker <- names[,c("Binome", "Genus", "Family_APGIII")]
  phylomaker <- na.omit(phylomaker)
  phylomaker <- unique(phylomaker)
  colnames(phylomaker) <- c("species","genus","family");
  rownames(phylomaker) <- capitalize(tolower(gsub(phylomaker$species,pattern=" ",replacement="_")))
  
  # For diagnostic only
  #sum(tolower(unique(TreeGhatsData$Family_APGIII))%in%tolower(nodes$family))/length(unique(TreeGhatsData$Family_APGIII)) # 99.2%
  #unique(TreeGhatsData$Family_APGIII)[!tolower(unique(TreeGhatsData$Family_APGIII))%in%tolower(nodes$family)] # Centroplacaceae
  #sum(unique(TreeGhatsData$Genus)%in%nodes$genus)/length(unique(TreeGhatsData$Genus)) # 79.9%
  #sum(unique(rownames(TreeGhatsData))%in%tolower(phylo.all.qian$tip.label))/length(unique(rownames(TreeGhatsData))) # 18.1%
  
  # With whole Qian phylogeny
  #result<-S.PhyloMaker(splist=phylomaker, tree=qian$phylo.all , nodes=qian$nodes, scenarios = scenarios)      # run the function S.PhyloMaker
  # With phylogeny generated for WG trees
  result<-S.PhyloMaker(splist=na.omit(phylomaker), tree=qian$phylo.all, nodes=qian$nodes, scenarios = scenarios) # run the function S.PhyloMaker
  
  return(result)
}