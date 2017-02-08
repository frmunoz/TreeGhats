taxocheck <- function(names, apg = F)
# names = vector of taxa names (genus species, with space separation)
# apg = should we provide apg3 family names
# IMPORTANT reftaxo = base de données à utiliser, doit être chargée préalablement et dénommée ainsi
{
  names <- na.omit(names)
  names <- tolower(unique(names))
  orig.names <- names;
  # FoundName is the name found in the database, which can differ from the original name if there are typos
  tab <- data.frame(FoundName=rep(NA,length(names)),Typo=rep(NA,length(names)),StatusTBGRI=rep(NA,length(names)), StatusPlantList=rep(NA,length(names)), AcceptedTBGRI=rep(NA,length(names)), AcceptedPlantList=rep(NA,length(names)));

  # Research in Western Ghats database, reftaxo
  rownames(tab) <- names;
  sel <- intersect(reftaxo$Full.name,rownames(tab));
  tab[sel,]$StatusTBGRI <- unlist(sapply(sel,function(x) reftaxo[which(reftaxo$Full.name==x),"Status_TBGRI"])); 
  tab[sel,]$StatusPlantList <- unlist(sapply(sel,function(x) reftaxo[which(reftaxo$Full.name==x),"Status_thePlantlist"]))
  tab[sel,]$FoundName <- sel
  tab[sel,]$Typo <- F
  tab[!is.na(tab$StatusTBGRI) & tab$StatusTBGRI=="Synonym",]$AcceptedTBGRI <- unlist(sapply(rownames(tab[!is.na(tab$StatusTBGRI) & tab$StatusTBGRI=="Synonym",]),function(x) reftaxo[which(reftaxo$Full.name==x),]$AcceptedNames_TBGRI) )
  tab[!is.na(tab$StatusPlantList) & tab$StatusPlantList=="Synonym",]$AcceptedPlantList <- unlist(sapply(rownames(tab[!is.na(tab$StatusPlantList) & tab$StatusPlantList=="Synonym",]),function(x) reftaxo[which(reftaxo$Full.name==x),]$ReferenceNames_ThePlantList) )
  tab$habit <- NA;
  tab[!is.na(tab$StatusTBGRI),]$habit <- unlist(sapply(rownames(tab[!is.na(tab$StatusTBGRI),]),function(x) reftaxo[which(reftaxo$Full.name==x),]$Habit_TBGRI) )
  rownames(tab) <- orig.names;
  
  # For taxa absent from reftaxo, check in PlantList
  require(Taxonstand)
  tab.plantlist <- c();
  for(i in rownames(tab)[is.na(tab$StatusPlantList)]) {res=TPL(i); if(nrow(res)==1) tab.plantlist <- rbind(tab.plantlist,res) else tab.plantlist <- c(tab.plantlist,NA)}
  rownames(tab.plantlist) <- rownames(tab)[is.na(tab$StatusPlantList)]
  tab$StatusPlantList[is.na(tab$StatusPlantList)] <- unlist(tab.plantlist[rownames(tab)[is.na(tab$StatusPlantList)],"Taxonomic.status"])
  tab$StatusPlantList[is.na(tab$StatusPlantList)] <- unlist(tab.plantlist[rownames(tab)[is.na(tab$StatusPlantList)],"Typo"])
  tab$StatusPlantList[tab$StatusPlantList==""] <- NA
  tab[!is.na(tab$StatusPlantList) & tab$StatusPlantList=="Synonym" & is.na(tab$AcceptedPlantList),]$AcceptedPlantList <- unlist(sapply(rownames(tab[!is.na(tab$StatusPlantList) & tab$StatusPlantList=="Synonym" & is.na(tab$AcceptedPlantList),]),function(x) paste(tab.plantlist[x,"New.Genus"],tab.plantlist[x,"New.Species"],sep=" ")))
  tab[rownames(tab.plantlist),]$Typo <- tab.plantlist$Typo
  # Reallocate names based on accepted name in PlantList (which correct possible typos)
  sel <- is.na(tab$StatusTBGRI) & !is.na(tab$StatusPlantList) & tab$StatusPlantList=="Accepted" & tab$Typo==T
  tab[sel,]$FoundName <- tolower(paste(tab.plantlist[rownames(tab)[sel],"New.Genus"],tab.plantlist[rownames(tab)[sel],"New.Species"],sep=" "))
  sel <- is.na(tab$StatusTBGRI) & !is.na(tab$StatusPlantList) & tab$Typo==F
  tab[sel,]$FoundName <- rownames(tab)[sel]
  sel <- is.na(tab$StatusTBGRI) & !is.na(tab$StatusPlantList) & tab$StatusPlantList!="Accepted" & tab$Typo==T
  tab[sel,]$FoundName <- unlist(sapply(tab.plantlist[rownames(tab)[sel],"ID"],function(x) tolower(paste0(read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=",x,"&csv=true",sep=""), header = TRUE, sep = ",", fill = TRUE, colClasses = "character", as.is = TRUE)[1,c("Genus","Species")],collapse=" "))))
  
  # For species found in PlantList but absent from reftaxo, recheck in reftaxo
  sel <- is.na(tab$StatusTBGRI) & !is.na(tab$StatusPlantList) 
  tab[sel,]$StatusTBGRI <- unlist(sapply(tab[sel,]$FoundName,function(x) ifelse(sum(reftaxo$Full.name==x)!=0,reftaxo[which(reftaxo$Full.name==x),"Status_TBGRI"],NA))); 
  sel <- !is.na(tab$StatusTBGRI) & tab$StatusTBGRI=="Synonym" & is.na(tab$AcceptedTBGRI)
  tab[sel,]$AcceptedTBGRI <- unlist(sapply(tab[sel,]$FoundName,function(x) reftaxo[which(reftaxo$Full.name==x),]$AcceptedNames_TBGRI))
  
  # APG3 family
  if(apg)
  {
    tab$Genus <- sapply(rownames(tab),function(x) strsplit(x,split=" ")[[1]][1])
    require(taxize)
    gen.list <- unique(tab$Genus)
    fam.list <- c();
    for(x in gen.list)
    {
      res <- classification(get_uid(x,ask=F),db="ncbi");
      res <- res[[1]]
      fam.list <- c(fam.list,ifelse(is.null(dim(res)),NA,res[res[,2]=="family",1]));
    }
    tab$Family.apg <- unlist(sapply(tab$Genus,function(x) fam.list[which(gen.list==x)]))
  }
  
  # Return a table with original names in first columns, and information on these taxa in other columns
  if(spelling) return(cbind(orig.names,tab)) else return(tab)
}