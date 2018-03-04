
taxocheck <- function(names, otherinfo = T, max.distance = 2, phylo = F)

{
  # names = vector of taxa names (genus species, with space separation)
  if(!is.vector(names))
    if(!"NAMES"%in%toupper(colnames(names)) & !"BINOME"%in%toupper(colnames(names))) 
    {
      stop("input should be a vector of names")
    } else 
      if("NAMES"%in%toupper(colnames(names))) 
      {
        names <- names[,which(toupper(colnames(names))=="NAMES")[1]]
      } else if("BINOME"%in%toupper(colnames(names))) names <- names[,which(toupper(colnames(names))=="BINOME")[1]]
        
  # TreeGhatsData  must be use as the database
  data('TreeGhatsData', package='TreeGhats', envir=environment())
  TreeGhatsData <- get("TreeGhatsData", envir=environment())

  # Pb with definition of sp to be checked (see below)
  #sp <- NULL
  
  names <- na.omit(names)
  names <- tolower(unique(str_trim(names)))
  names<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)", "", names, perl=T)
  orig.names <- names;
  tab<-data.frame(FoundName=rep(NA, length(names)), Typo=rep(NA, length(names)), 
                  Genus=rep(NA, length(names)), 
                  Species=rep(NA,length(names)),
                  InfrataxonRank=rep(NA,length(names)), 
                  InfrataxonName=rep(NA,length(names)))
  rownames(tab) <- orig.names;
  # Detect incomplete names or names with number:
  num<-c()
  for(i in 0:9) {num<-rbind(num, str_detect(names, as.character(i)))}
  num<-apply(num, 2, function(x) any(x))
  num<-cbind(num, sapply(rownames(tab), function(x) length(unlist(strsplit(x, split=" ")))==1),
             unlist(lapply(rownames(tab), function(x) strsplit(x, split=" ")[[1]][2]%in%c("sp.", "sp"," species"))))
  num<-apply(num, 1, function(x) any(x))
  tab$FoundName<-ifelse(num==T, "IncompleteName", NA)
  sel<-tab$FoundName!= "IncompleteName" | is.na(tab$FoundName)
  if(sum(sel)!=0)
  {
    tab[sel,]$Genus <- capitalize(do.call(rbind, strsplit(as.vector(names[sel]), " "))[,1])
    tab[sel,]$Species <-  sapply(names[sel], function(x) ifelse(length(unlist(strsplit(x, " "))) > 1, strsplit(x, " ")[[1]][2], ""))
  }
  
  # Detect infrataxon
  vec0 <- c( "nothossp.", " nothossp ", "nothosubsp.", " nothosubsp ", "cultivar.", 
             " cultivar ",  " subfo ",  "subf."," subf ", " subproles ",  "cf.", " cf ", "aff.", " aff ",  "s.l.", "s.l ",  
             "s.str.", "s.str ", "x.", " x ", "X.", " X ",  "f.", " f ",  "fo.", " fo ", 
             " forma ", "subvar.", " subvar ",  "var.", " var ",  "subsp.", " subsp ",  
             "ssp.", " ssp ", " gama ", " grex ", "lus.", " lus ", " lusus ", " monstr ",  
             "nm.", " nm ", "prol.", " prol ", " proles ", " race ", "subvar.",  "cv.", " cv ")
  InfrataxonRank<-apply(sapply(names, function(names) 
    sapply(vec0, function(x) 
      ifelse(length(grep(x, names, fixed = TRUE)) > 0, T, NA))), 2, function(x) 
        ifelse(all(is.na(x)), NA, names(x[!is.na(x)])))
  InfrataxonRank<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)", "", InfrataxonRank, perl=T)
  
  if(length(unique(InfrataxonRank))>1)
  {
    for(j in 1:length(unique(InfrataxonRank[!is.na(InfrataxonRank)]))){
      names<-as.vector(sapply(names, function(x) gsub(unique(InfrataxonRank[!is.na(InfrataxonRank)])[j]," ", x, fixed = TRUE)))}
    names<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)", "", names, perl=T)
    ## Problem here because sp is undefined
    #names <- ifelse(substr(names, 1, 1) == " ", substr(sp, 2, nchar(names)), names)
    InfrataxonName <- sapply(names, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][3], ""))
    InfrataxonRank<-replace(InfrataxonRank, InfrataxonRank%in%c("subsp", "ssp.", "ssp"), "subsp.")
    InfrataxonRank<-replace(InfrataxonRank, InfrataxonRank%in%c("f", "fo", "fo."), "f.")
    InfrataxonRank<-replace(InfrataxonRank, InfrataxonRank=="var","var.")
    tab$InfrataxonRank<-as.character(InfrataxonRank)
    tab$InfrataxonName<-as.character(InfrataxonName)
    rownames(tab)[!is.na(tab$InfrataxonRank)]=paste(tab[!is.na(tab$InfrataxonRank),]$Genus, tab[!is.na(tab$InfrataxonRank),]$Species,
                                                    tab[!is.na(tab$InfrataxonRank),]$InfrataxonRank, tab[!is.na(tab$InfrataxonRank),]$InfrataxonName, sep=" ")
  }
  # Research in TreeGhatsdata the taxonomic Information
  # FoundName is the name found in the database, which can differ from the original name if there are typos
  # Research names without spelling difference
  sel <- intersect(TreeGhatsData$Name,rownames(tab));
  tab[sel,]$FoundName <- sel
  tab$Typo <- ifelse(rownames(tab)%in% sel, F, NA)
  tab$ID_TPL <-NA;tab$Status_TPL <-NA;tab$ReferenceName_TPL <-NA;tab$ReferenceAuthority_TPL <-NA;tab$Status_TBGRI=NA;tab$ReferenceName_TBGRI <-NA;tab$ReferenceAuthority_TBGRI <-NA; 
  tab$Status_proposed<-NA;tab$ReferenceName_proposed <-NA;tab$ReferenceAuthority_proposed <-NA;tab$Infrataxon<-NA;tab$Family_APGIII <-NA
  #tab$NewID_TPL<-NA;
  
  # Research names  with spelling errors maxDist
  selcor<-setdiff(rownames(tab),TreeGhatsData$Name)[!is.na(as.character(sapply(setdiff(rownames(tab),TreeGhatsData$Name),function(x) TreeGhatsData$Name[amatch(x,TreeGhatsData$Name, maxDist=max.distance)])))]
  if(length(selcor)>=1)
  {cornames<-as.character(sapply(setdiff(rownames(tab),TreeGhatsData$Name),function(x) TreeGhatsData$Name[amatch(x,TreeGhatsData$Name, maxDist=max.distance)]))
  tab[selcor,]$FoundName <- cornames[!is.na(cornames)]
  tab[selcor,]$Typo <- T
  tab$Typo[tab$FoundName=="NULL"]<-NA
  tab$FoundName[tab$FoundName=="NULL"]<-NA
  }
  sel<-!is.na(tab$Typo)&tab$FoundName!="IncompleteName"
  
  if(any(sel))
  {
    WGinfo<-NA
    WGinfo<- sapply(tab[sel,]$FoundName,function(x) TreeGhatsData[which(TreeGhatsData$Name==x),c("ID_TPL","Family_APGIII","Status_TPL","ReferenceName_TPL","ReferenceAuthority_TPL","Status_TBGRI","ReferenceName_TBGRI","ReferenceAuthority_TBGRI","Status_proposed","ReferenceName_proposed","ReferenceAuthority_proposed","Family_APGIII")]); 
    tab[sel,]$Status_TBGRI <-unlist(WGinfo["Status_TBGRI",])   
    tab[sel,]$Status_TPL <- unlist(WGinfo["Status_TPL",])
    tab[sel,]$ReferenceName_TBGRI <- unlist(WGinfo["ReferenceName_TBGRI",])
    tab[sel,]$ReferenceAuthority_TPL <- unlist(WGinfo["ReferenceAuthority_TPL",])
    tab[sel,]$ReferenceAuthority_TBGRI <- unlist(WGinfo["ReferenceAuthority_TBGRI",])
    tab[sel,]$ReferenceName_TPL <-unlist(WGinfo["ReferenceName_TPL",]) 
    tab[sel,]$ID_TPL <-unlist(WGinfo["ID_TPL",])
    tab[sel,]$Family_APGIII <-unlist(WGinfo["Family_APGIII",])
    tab[sel,]$Status_proposed <-unlist(WGinfo["Status_proposed",])
    tab[sel,]$ReferenceAuthority_proposed <-unlist(WGinfo["ReferenceAuthority_proposed",])
    tab[sel,]$ReferenceName_proposed <-unlist(WGinfo["ReferenceName_proposed",])
  } else {warning("No match in TreeGhats database")}
  #sel <- !is.na(tab$ReferenceName_proposed)
  #if (any(sel))
  #{tab[sel,]$NewID_TPL<- unlist(sapply(tolower(tab$ReferenceName_proposed[sel]),function(x) TreeGhatsData[which(TreeGhatsData$Name==x),"ID_TPL"]));} 

  ## For taxa absent from TreeGhatsData, check in PlantList ##
  taxonCheckTPL<-rownames(tab[is.na(tab$FoundName),])
  tab.plantlist <- c();
  if(length(taxonCheckTPL)>=1)
  {pb <- utils::winProgressBar(title = "progress bar", min = 0, max = length(taxonCheckTPL), width = 300)
  for(i in 1:length(taxonCheckTPL))
  {Sys.sleep(0.1);utils::setWinProgressBar(pb, i, title=paste("Check in TPL" ,round(i/length(taxonCheckTPL)*100, 0),"% done"));res=TPLck2(taxonCheckTPL[i]); tab.plantlist <- rbind(tab.plantlist,res)}
  #}
  
  rownames(tab.plantlist) <- rownames(tab[is.na(tab$FoundName),])
  tab.plantlist$NewNames<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)","",paste(tab.plantlist$New.Genus,tab.plantlist$New.Species,tab.plantlist$New.Infraspecific,sep=" "), perl=T)
  # Complete tab with infos from TPL
  sel<-is.na(tab$FoundName)
  tab[is.na(tab$FoundName),]$Typo <- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Typo"]
  tab[is.na(tab$FoundName),]$Family_APGIII<- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Family"]
  tab[is.na(tab$FoundName),]$ID_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"ID"]
  # Check several homonyms
  if(is.numeric(tab.plantlist$severalNames))
  {tab[rownames(tab.plantlist)[!is.na(tab.plantlist$severalNames)],]$Status_TPL<-"SeveralHomonyms"
  sel<-!is.na(tab$Status_TPL) & tab$Status_TPL=="SeveralHomonyms"
  tab[sel,]$FoundName<-rownames(tab[sel,])
  }  
  #Check infrataxon in tlp refnames  
  if(any(tab.plantlist$New.Infraspecific!=""&tab.plantlist$Taxonomic.status!=""))
  {
    sel<-tab.plantlist$New.Infraspecific!=""&tab.plantlist$Taxonomic.status!=""
    tab.plantlist[sel,]$NewNames=unlist(sapply(tab.plantlist[sel,]$NewNames,function(x) paste0(read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=",x,"&csv=true",sep=""), header = TRUE, sep = ",", fill = TRUE, colClasses = "character", row.names=1,as.is = TRUE)[tab.plantlist[tab.plantlist$NewNames==x,]$New.ID,c("Genus","Species","Infraspecific.rank","Infraspecific.epithet")],collapse=" ")))
    if(any(tab.plantlist[sel,]$Typo==T))
    {sel<-sel & tab.plantlist$Typo==T & tab.plantlist$Taxonomic.status!="Synonym"
    tab.plantlist[sel ,]$Infraspecific<-tab.plantlist[sel ,]$New.Infraspecific}                                                                                                                                                                
  } 
  
  # Complete tab with infos from TPL
  tab[is.na(tab$FoundName),]$Status_TPL <- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Taxonomic.status"]
  #tab[is.na(tab$FoundName),]$NewID_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"New.ID"]
  tab[is.na(tab$FoundName),]$ReferenceAuthority_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Authority"]
  tab[is.na(tab$FoundName),]$InfrataxonName<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Infraspecific"]
  tab[is.na(tab$FoundName),]$Typo <- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Typo"]
  tab[is.na(tab$FoundName),]$ReferenceName_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),]$NewNames
  sel<-tab$Status_TPL=="" &! is.na(tab$Status_TPL)
  if(any(sel)){tab[sel,]$Typo<-NA
  tab[sel,]$ReferenceName_TPL<-NA
  tab[sel,]$FoundName<-NA}
  
  ## Correct possible typos ## 
  sel <- is.na(tab$FoundName) & tab$Status_TPL %in% c("Accepted","Unresolved")
  if(any(sel)){
    tab[sel,]$FoundName <- tolower(tab.plantlist[rownames(tab)[sel],]$NewNames)}
  sel <- is.na(tab$FoundName) & tab$Status_TPL=="Synonym" & tab$Typo==F
  if(any(sel)){
    tab[sel,]$FoundName <- tolower(rownames(tab)[sel])}
  sel <-  tab$Typo==T &  tab$Status_TPL %in% c("Synonym","SeveralHomonyms") & !is.na(tab$Status_TPL)
  if(any(sel)){
    tab[sel,]$FoundName <-unlist(sapply(tab[rownames(tab)[sel],"ID_TPL"],function(x) tolower(paste0(read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=",x,"&csv=true",sep=""), header = TRUE, sep = ",", fill = TRUE, colClasses = "character", as.is = TRUE)[1,c("Genus","Species")],collapse=" "))))}
  sel<-!is.na(tab$Typo) &tab$Typo==T
  if(any(sel)){
    tab[sel,]$Species<-sapply(tab[sel,]$FoundName, function(x)  strsplit(x, " ")[[1]][2])
    tab[sel,]$Genus<-capitalize(sapply(tab[sel,]$FoundName, function(x)  strsplit(x, " ")[[1]][1]))
  }
  
  #For infrataxon names
  if(any(tab.plantlist$New.Infraspecific!="" &  tab.plantlist$Typo==T & tab.plantlist$Taxonomic.status=="Synonym"))
  {sel<-tab.plantlist$New.Infraspecific!=""& tab.plantlist$Typo==T & tab.plantlist$Taxonomic.status=="Synonym"
  IDsel<-tab.plantlist[sel,]$ID
  tab[tab$ID_TPL%in%IDsel ,]$InfrataxonName<-unlist(sapply(IDsel,function(x) 
    read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=", paste(tab[tab$ID_TPL==x & !is.na(tab$ID_TPL),]$Genus,
                                                                         tab[tab$ID_TPL==x & !is.na(tab$ID_TPL),]$Species,sep=" "), "&csv=true", sep=""), 
             header = TRUE, sep = ",", fill = TRUE, colClasses = "character", row.names=1, as.is = TRUE)[x,"Infraspecific.epithet"]))
  tab[tab$ID_TPL%in%IDsel ,]$FoundName <- paste0(tab[tab$ID_TPL%in%IDsel, c("FoundName", "InfrataxonRank", "InfrataxonName")], collapse=" ")}
  close(pb) 
  }
  
  # Check again in TreeGhatsData
  sel <- is.na(tab$Status_TBGRI) & tab$Typo==T &!is.na(tab$Typo) & !is.na(tab$Status_TPL) & tab$FoundName%in%TreeGhatsData$Name
  if (any(sel))
  {
    WGinfo<- sapply(tab$FoundName[sel],function(x) TreeGhatsData[which(TreeGhatsData$Name==x),c("Status_TBGRI","ReferenceName_TBGRI","ReferenceAuthority_TBGRI","Status_proposed")]); 
    tab[sel,]$Status_TBGRI <-unlist(WGinfo["Status_TBGRI",])   
    tab[sel,]$ReferenceName_TBGRI <- unlist(WGinfo["ReferenceName_TBGRI",])
    tab[sel,]$ReferenceAuthority_TBGRI <- unlist(WGinfo["ReferenceAuthority_TBGRI",])
  }
  
  # Infrataxon management: Count number of Infrataxa in TreeGhatsData  
  tab$Infrataxon<-NA
  sel <- !is.na(tab$Status_TBGRI) & is.na(tab$InfrataxonRank)
  if (any(sel))
  {
    InfrataxonCount<-table(paste(TreeGhatsData$Genus,TreeGhatsData$Species, sep=" "))-1
    tab$Infrataxon[sel]<-sapply(tab$FoundName[sel],function(x) InfrataxonCount[which(tolower(names(InfrataxonCount))==x)])
    tab$Infrataxon[tab$Infrataxon>1]<-"SeveralInfrataxa"
    tab$Infrataxon[tab$Infrataxon==1]<-"OneInfrataxon"
    tab$Infrataxon[tab$Infrataxon==0]<-NA
  }
  
  ## Statut proposed ##
  tab1<-tab
  tab<-tab[,-c(3:7)]
  tab$Genus<-NA;tab$Species<-NA;tab$InfrataxonRank<-NA; tab$InfrataxonName<-NA;
  sel<-!is.na(tab$FoundName) & tab$FoundName!="IncompleteName" & is.na(tab$Status_TBGRI)
  if(sum(sel)!=0)
  {
    tab[sel,]$Status_proposed<-tab[sel,]$Status_TPL
    tab[sel,]$ReferenceName_proposed<-tab[sel,]$ReferenceName_TPL
    tab[sel,]$ReferenceAuthority_proposed <-tab[sel,]$ReferenceAuthority_TPL
  }
  sel<-!is.na(tab$ReferenceName_proposed)
  if(sum(sel!=0))
  {
    tab[sel,]$Genus <- capitalize(do.call(rbind, strsplit(as.vector(tab[sel,]$ReferenceName_proposed), " "))[,1])
    tab[sel,]$Species <-  do.call(rbind, strsplit(as.vector(tab[sel,]$ReferenceName_proposed), " "))[,2]
    tab[sel,]$InfrataxonRank <-  sapply(tab[sel,]$ReferenceName_proposed, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][3], ""))
    tab[sel,]$InfrataxonName <-  sapply(tab[sel,]$ReferenceName_proposed, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][4], ""))
  }
  
  sel<-tab$FoundName=="IncompleteName" & !is.na(tab$FoundName)
  if(any(sel))
  {
    tab[sel,]$Status_proposed <- "IncompleteName"
  }
  
  # For Infrataxon the reference name proposed depend on the number of Infrataxa present in WG.
  sel <- tab$Infrataxon=="SeveralInfrataxa" & !is.na(tab$Infrataxon)
  if (any(sel))
  {
    tab$InfrataxonRank[sel]<-NA
    tab$InfrataxonName[sel]<-NA
    tab$ReferenceName_proposed[sel]<-paste(tab$Genus[sel],tab$Species[sel], sep=" ")
  }
  sel <- tab$Infrataxon=="OneInfrataxon" & !is.na(tab$Infrataxon)
  if (any(sel))
  { 
    tab$ReferenceName_proposed[sel]<-sapply(tab$ReferenceName_proposed[sel],function(x) TreeGhatsData$ReferenceName_proposed[which(paste(TreeGhatsData$Genus,TreeGhatsData$Species, sep=" ")==x & !is.na(TreeGhatsData$InfraTaxonRank))])
    tab$InfrataxonRank[sel]<-sapply(tolower(tab$ReferenceName_proposed[sel]),function(x) TreeGhatsData$InfraTaxonRank[which(TreeGhatsData$Name==x)])
    tab$InfrataxonName[sel]<-sapply(tolower(tab$ReferenceName_proposed[sel]),function(x) TreeGhatsData$InfraTaxonNames[which(TreeGhatsData$Name==x)])
   } 
  
  ## Provide some ecological information ##
  if(otherinfo & any(!is.na(tab$ReferenceName_proposed)))
  {

    tab$Origin <- NA; tab$Habit <- NA; tab$Phenology <- NA;tab$IUCN <- NA; 
    Info <-sapply(tab[!is.na(tab$ReferenceName_proposed),]$ReferenceName_proposed,function(x) TreeGhatsData[which(TreeGhatsData$Name==tolower(x)),c("Origin","Habit","Phenology","IUCN_Status")])
    Info[lengths(Info) == 0] <- NA_character_
    tab$Origin[!is.na(tab$ReferenceName_proposed)]<-unlist(Info["Origin",])
    tab$Habit[!is.na(tab$ReferenceName_proposed)]<-unlist(Info["Habit",]) 
    tab$Phenology[!is.na(tab$ReferenceName_proposed)]<-unlist(Info["Phenology",])
    tab$IUCN[!is.na(tab$ReferenceName_proposed)]<-unlist(Info["IUCN_Status",])
  }
  
  tab[tab == ""] <- NA
  if(all(is.na(tab$Infrataxon)))
  {tab<-tab[,-which(colnames(tab)=="Infrataxon")]}
  if(all(is.na(tab$InfrataxonName)))
  {tab<-tab[,-which(colnames(tab) %in% c("InfrataxonName","InfrataxonRank"))]}
  
  rownames(tab)<-orig.names
  tab <- data.frame(tab)
  
  if(!phylo)
  {
    # Return a table with original names in Rownames, and information on these taxa in other columns
    return(tab)
  } else
  {
    # Create the phylogeny corresponding to the taxa (create.phylo with default options)
    phylo <- create.phylo(tab)
    return(list(tab=tab, phylo=phylo$Scenario.3))
  }
}