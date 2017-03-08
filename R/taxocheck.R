taxocheck <- function(names, othersinfo = T, iucn=F, max.distance = 1)
{
  force(reftaxo)
  #### faut-il vraiment garder l'option spelling? on peut toujours mettre comme nom de ligne la liste initiale et les gens check par la colonne typo? ###
  # names = vector of taxa names (genus species, with space separation)
  # apg = should we provide apg3 family names
  # IMPORTANT reftaxo = base de donnees a utiliser
  # Pour l'instant la fonction ne fournit pas id tropicos
  names <- na.omit(names)
  names <- tolower(unique(str_trim(names)))
  names<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)","",names, perl=T)
  orig.names <- names;
  tab<-data.frame(FoundName=rep(NA,length(names)),Typo=rep(NA,length(names)),Genus=rep(NA,length(names)), Species=rep(NA,length(names)),InfrataxonRank=rep(NA,length(names)), InfrataxonName=rep(NA,length(names)))
  rownames(tab) <- orig.names;
  # Detect incomplete names or names with number:
  num<-c()
  for(i in 0:9){num<-rbind(num,str_detect(names, as.character(i)))}
  num<-apply(num,2,function(x) any(x))
  num<-cbind(num,sapply(rownames(tab),function(x) length(unlist(strsplit(x,split=" ")))==1),unlist(lapply(rownames(tab),function(x) strsplit(x,split=" ")[[1]][2]%in%c("sp.","sp","species"))))
  num<-apply(num,1,function(x) any(x))
  tab$FoundName<-ifelse(num==T, "incompleteName", NA)
  sel<-tab$FoundName!= "incompleteName" | is.na(tab$FoundName)
  tab[sel,]$Genus <- capitalize(do.call(rbind, strsplit(as.vector(names[sel]), " "))[,1])
  tab[sel,]$Species <-  sapply(names[sel], function(x) ifelse(length(unlist(strsplit(x, " "))) > 1, strsplit(x, " ")[[1]][2], ""))
  
  # Detect infrataxon
  vec0 <- c( "nothossp.", " nothossp ", "nothosubsp.", " nothosubsp ", "cultivar.", 
             " cultivar ",  " subfo ",  "subf."," subf ", " subproles ",  "cf.", " cf ", "aff.", " aff ",  "s.l.", "s.l ",  
             "s.str.", "s.str ", "×",  "x.", " x ", "X.", " X ",  "f.", " f ",  "fo.", " fo ", 
             " forma ", "subvar.", " subvar ",  "var.", " var ",  "subsp.", " subsp ",  
             "ssp.", " ssp ", " gama ", " grex ", "lus.", " lus ", " lusus ", " monstr ",  
             "nm.", " nm ", "prol.", " prol ", " proles ", " race ", "subvar.",  "cv.", " cv ")
  InfrataxonRank<-apply(sapply(names, function(names) sapply(vec0,function(x) ifelse(length(grep(x, names, fixed = TRUE)) > 
                                                                                       0, T, NA))),2,function(x) ifelse(all(is.na(x)),NA,names(x[!is.na(x)])))
  InfrataxonRank<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)","",InfrataxonRank, perl=T)
  
  
  if(length(unique(InfrataxonRank))>1)
  {
    for(j in 1:length(unique(InfrataxonRank[!is.na(InfrataxonRank)]))){
      names<-as.vector(sapply(names,function(x) gsub(unique(InfrataxonRank[!is.na(InfrataxonRank)])[j]," ", x, fixed = TRUE)))}
    names<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)","",names, perl=T)
    names <- ifelse(substr(names, 1, 1) == " ", substr(sp, 2, nchar(names)), names)
    InfrataxonName <- sapply(names, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][3], ""))
    InfrataxonRank<-replace(InfrataxonRank,InfrataxonRank%in%c("subsp","ssp.", "ssp"),"subsp.")
    InfrataxonRank<-replace(InfrataxonRank,InfrataxonRank%in%c("f","fo","fo."),"f.")
    InfrataxonRank<-replace(InfrataxonRank,InfrataxonRank=="var","var.")
    tab$InfrataxonRank<-as.character(InfrataxonRank)
    tab$InfrataxonName<-as.character(InfrataxonName)
    rownames(tab)[!is.na(tab$InfrataxonRank)]=paste(tab[!is.na(tab$InfrataxonRank),]$Genus,tab[!is.na(tab$InfrataxonRank),]$Species,
                                                    tab[!is.na(tab$InfrataxonRank),]$InfrataxonRank,tab[!is.na(tab$InfrataxonRank),]$InfrataxonName, sep=" ")
    
  }
  
  # FoundName is the name found in the database, which can differ from the original name if there are typos
  # Research in Western Ghats database, without spelling difference
  sel <- intersect(reftaxo$Full.name,rownames(tab));
  tab[sel,]$FoundName <- sel
  tab$Typo <- ifelse(rownames(tab)%in% sel, F, NA)
  tab$ID_TPL <-NA;tab$Status_TPL <-NA;tab$ReferenceName_TPL <-NA;tab$ReferenceAuthority_TPL <-NA;tab$Status_TBGRI=NA;tab$ReferenceName_TBGRI <-NA;tab$ReferenceAuthority_TBGRI <-NA; 
  tab$Status_proposed<-NA;tab$ReferenceName_proposed <-NA;tab$ReferenceAuthority_proposed <-NA;tab$NewID_TPL<-NA;tab$Family_APGIII <-NA
  
  # Research in Western Ghats database, reftaxo with spelling errors maxDist=2
  selcor<-setdiff(rownames(tab),reftaxo$Full.name)[!is.na(as.character(sapply(setdiff(rownames(tab),reftaxo$Full.name),function(x) reftaxo$Full.name[amatch(x,reftaxo$Full.name, maxDist=max.distance)])))]
  if(length(selcor)>=1)
  {cornames<-as.character(sapply(setdiff(rownames(tab),reftaxo$Full.name),function(x) reftaxo$Full.name[amatch(x,reftaxo$Full.name, maxDist=max.distance)]))
  tab[selcor,]$FoundName <- cornames[!is.na(cornames)]
  tab[selcor,]$Typo <- T
  sel<-c(sel,cornames[!is.na(cornames)])}
  if(length(sel)==0){warning("No match in TreeGhats database")}
  if(length(sel)>=1)
  {
    # Research in Western Ghats database, info taxonomique 
    WGinfo<-NA
    WGinfo<- sapply(tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$FoundName,function(x) reftaxo[which(reftaxo$Full.name==x),c("ID_TPL","Family_APGIII","Status_TPL","ReferenceName_TPL","ReferenceAuthority_TPL","Status_TBGRI","ReferenceName_TBGRI","ReferenceAuthority_TBGRI","Status_proposed","ReferenceName_proposed","ReferenceAuthority_proposed","Family_APGIII")]); 
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$Status_TBGRI <-unlist(WGinfo["Status_TBGRI",])   
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$Status_TPL <- unlist(WGinfo["Status_TPL",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceName_TBGRI <- unlist(WGinfo["ReferenceName_TBGRI",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceAuthority_TPL <- unlist(WGinfo["ReferenceAuthority_TPL",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceAuthority_TBGRI <- unlist(WGinfo["ReferenceAuthority_TBGRI",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceName_TPL <-unlist(WGinfo["ReferenceName_TPL",]) 
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ID_TPL <-unlist(WGinfo["ID_TPL",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$Family_APGIII <-unlist(WGinfo["Family_APGIII",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$Status_proposed <-unlist(WGinfo["Status_proposed",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceAuthority_proposed <-unlist(WGinfo["ReferenceAuthority_proposed",])
    tab[!is.na(tab$Typo)&tab$FoundName!="incompleteName",]$ReferenceName_proposed <-unlist(WGinfo["ReferenceName_proposed",])
  }
  sel <- !is.na(tab$ReferenceName_proposed)
  if (any(sel))
  {tab[sel,]$NewID_TPL<- sapply(tolower(tab$ReferenceName_proposed[sel]),function(x) reftaxo[which(reftaxo$Full.name==x),"ID_TPL"]);} 
    
  
  # For taxa absent from reftaxo, check in PlantList
  #### Ici intervient la la fonction TPLck2, modifiee pour savoir s'il y a plusieurs synomym? comment appeler une fonction si on l'inclut dans notre package? ###  
  taxonCheckTPL<-rownames(tab[is.na(tab$FoundName),])
  tab.plantlist <- c();
  if(length(taxonCheckTPL)>=1)
  {pb <- winProgressBar(title = "progress bar", min = 0,max = length(taxonCheckTPL), width = 300)
  for(i in 1:length(taxonCheckTPL))
  {Sys.sleep(0.1);setWinProgressBar(pb, i, title=paste("Check in TPL" ,round(i/length(taxonCheckTPL)*100, 0),"% done"));res=TPLck2(taxonCheckTPL[i]); tab.plantlist <- rbind(tab.plantlist,res)}
  #}
  rownames(tab.plantlist) <- rownames(tab[is.na(tab$FoundName),])
  tab.plantlist$NewNames<-gsub("(^\\s+|\\s+$|(?<=\\s)\\s)","",paste(tab.plantlist$New.Genus,tab.plantlist$New.Species,tab.plantlist$New.Infraspecific,sep=" "), perl=T)
  # Complete tab with infos from TPL
  sel<-is.na(tab$FoundName)
  tab[is.na(tab$FoundName),]$Typo <- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Typo"]
  tab[is.na(tab$FoundName),]$Family_APGIII<- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Family"]
  # Check several homonyms
  if(is.numeric(tab.plantlist$severalNames))
  {tab[rownames(tab.plantlist)[!is.na(tab.plantlist$severalNames)],]$Status_TPL<-"TPLfoundSeveralHomonyms"
  sel<-!is.na(tab$Status_TPL) & tab$Status_TPL=="TPLfoundSeveralHomonyms"
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
  tab[is.na(tab$FoundName),]$ID_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"ID"]
  tab[is.na(tab$FoundName),]$NewID_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"New.ID"]
  tab[is.na(tab$FoundName),]$ReferenceAuthority_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Authority"]
  tab[is.na(tab$FoundName),]$InfrataxonName<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Infraspecific"]
  tab[is.na(tab$FoundName),]$Typo <- tab.plantlist[rownames(tab[is.na(tab$FoundName),]),"Typo"]
  tab[is.na(tab$FoundName),]$ReferenceName_TPL<-tab.plantlist[rownames(tab[is.na(tab$FoundName),]),]$NewNames
  sel<-tab$Status_TPL=="" &! is.na(tab$Status_TPL)
  if(any(sel)){tab[sel,]$Typo<-NA
  tab[sel,]$ReferenceName_TPL<-NA
  tab[sel,]$FoundName<-NA}
  
  # Correct possible typos  
  sel <- is.na(tab$FoundName) & tab$Status_TPL %in% c("Accepted","Unresolved")
  if(any(sel)){
  tab[sel,]$FoundName <- tolower(tab.plantlist[rownames(tab)[sel],]$NewNames)}
  sel <- is.na(tab$FoundName) & tab$Status_TPL=="Synonym" & tab$Typo==F
  if(any(sel)){
  tab[sel,]$FoundName <- tolower(rownames(tab)[sel])}
  sel <-  tab$Typo==T &  tab$Status_TPL %in% c("Synonym","TPLfoundSeveralHomonyms") & !is.na(tab$Status_TPL)
  if(any(sel)){
  tab[sel,]$FoundName <-unlist(sapply(tab.plantlist[rownames(tab)[sel],"ID"],function(x) tolower(paste0(read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=",x,"&csv=true",sep=""), header = TRUE, sep = ",", fill = TRUE, colClasses = "character", as.is = TRUE)[1,c("Genus","Species")],collapse=" "))))}
  sel<-!is.na(tab$Typo) &tab$Typo==T
  if(any(sel)){
  tab[sel,]$Species<-sapply(tab[sel,]$FoundName, function(x)  strsplit(x, " ")[[1]][2])
  tab[sel,]$Genus<-capitalize(sapply(tab[sel,]$FoundName, function(x)  strsplit(x, " ")[[1]][1]))
  }
  #Correct infraTaxonNames
  if(any(tab.plantlist$New.Infraspecific!="" &  tab.plantlist$Typo==T & tab.plantlist$Taxonomic.status=="Synonym"))
  {sel<-tab.plantlist$New.Infraspecific!=""& tab.plantlist$Typo==T & tab.plantlist$Taxonomic.status=="Synonym"
  IDsel<-tab.plantlist[sel,]$ID
  tab[tab$ID_TPL%in%IDsel ,]$InfrataxonName<-unlist(sapply(IDsel,function(x) read.csv(paste("http://www.theplantlist.org/tpl1.1/search?q=",paste(tab[tab$ID_TPL==x & !is.na(tab$ID_TPL),]$Genus,tab[tab$ID_TPL==x & !is.na(tab$ID_TPL),]$Species,sep=" "),"&csv=true",sep=""), header = TRUE, sep = ",", fill = TRUE, colClasses = "character", row.names=1,as.is = TRUE)[x,"Infraspecific.epithet"]))
  tab[tab$ID_TPL%in%IDsel ,]$FoundName<-paste0(tab[tab$ID_TPL%in%IDsel,c("FoundName","InfrataxonRank","InfrataxonName")],collapse=" ")}
  close(pb) 
  }
  
  # Check again in reftaxo: is it usefull
  sel <- is.na(tab$Status_TBGRI) & tab$Typo==T &!is.na(tab$Typo) & !is.na(tab$Status_TPL) & tab$FoundName%in%reftaxo$Full.name
  if (any(sel))
  {
    WGinfo<- sapply(tab$FoundName[sel],function(x) reftaxo[which(reftaxo$Full.name==x),c("Status_TBGRI","ReferenceName_TBGRI","ReferenceAuthority_TBGRI","Status_proposed")]); 
    tab[sel,]$Status_TBGRI <-unlist(WGinfo["Status_TBGRI",])   
    tab[sel,]$ReferenceName_TBGRI <- unlist(WGinfo["ReferenceName_TBGRI",])
    tab[sel,]$ReferenceAuthority_TBGRI <- unlist(WGinfo["ReferenceAuthority_TBGRI",])
  }
  
  # Statut proposed:
  tab1<-tab
  
  tab<-tab[,-c(3:6)]
  tab$Genus<-NA;tab$Species<-NA;tab$InfrataxonRank<-NA; tab$InfrataxonName<-NA
  sel<-!is.na(tab$FoundName) & tab$FoundName!="incompleteName" & is.na(tab$Status_TBGRI)
  tab[sel,]$Status_proposed<-tab[sel,]$Status_TPL
  tab[sel,]$ReferenceName_proposed<-tab[sel,]$ReferenceName_TPL
  tab[sel,]$ReferenceAuthority_proposed <-tab[sel,]$ReferenceAuthority_TPL
  # sel<-!is.na(tab$FoundName)&tab$FoundName!="incompleteName" & !is.na(tab$Status_TBGRI) & tab$Status_TBGRI!="Absent"
  # tab[sel,]$Status_proposed<-tab[sel,]$Status_TBGRI
  # tab[sel,]$ReferenceName_proposed<-tab[sel,]$ReferenceName_TBGRI
  # tab[sel,]$ReferenceAuthority_proposed <-tab[sel,]$ReferenceAuthority_TBGRI
  # sel<-!is.na(tab$FoundName)&tab$FoundName!="incompleteName" & !is.na(tab$Status_TBGRI) & tab$Status_TBGRI=="Absent"
  # tab[sel,]$Status_proposed<-tab[sel,]$Status_TPL
  # tab[sel,]$ReferenceName_proposed<-tab[sel,]$ReferenceName_TPL
  # tab[sel,]$ReferenceAuthority_proposed <-tab[sel,]$ReferenceAuthority_TPL
  # sel<-!is.na(tab$FoundName) & tab$Status_TPL=="Absent" & tab$Status_TBGRI=="Absent" & !is.na(tab$Status_TBGRI) & !is.na(tab$Status_TPL)
  # if(any(sel))
  # {tab[sel,]$Status_proposed<-"New species"
  # tab[sel,]$ReferenceName_proposed<-capitalize(tab[sel,]$FoundName)}
  # sel<-!is.na(tab$FoundName)&tab$FoundName=="incompleteName"
  # if(any(sel))
  # {tab[sel,]$Status_proposed<-"IncompleteName"
  # tab[sel,]$FoundName<-NA}
  # sel<-tab$Status_TPL=="Absent" &!is.na(tab$Status_TPL)
  # if(any(sel))
  # {
  # tab[sel,]$Status_TPL<-NA}
  sel<-!is.na(tab$ReferenceName_proposed)
  tab[sel,]$Genus <- capitalize(do.call(rbind, strsplit(as.vector(tab[sel,]$ReferenceName_proposed), " "))[,1])
  tab[sel,]$Species <-  do.call(rbind, strsplit(as.vector(tab[sel,]$ReferenceName_proposed), " "))[,2]
  tab[sel,]$InfrataxonRank <-  sapply(tab[sel,]$ReferenceName_proposed, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][3], ""))
  tab[sel,]$InfrataxonName <-  sapply(tab[sel,]$ReferenceName_proposed, function(x) ifelse(length(unlist(strsplit(x, " "))) > 2, strsplit(x, " ")[[1]][4], ""))
  sel<-!is.na(tab$Status_TPL) & tab$Status_TPL=="TPLfoundSeveralHomonyms"
  if(any(sel))
  {
  tab[sel,]$Genus <- tab1[sel,]$Genus
  tab[sel,]$Species <-tab1[sel,]$Species  
  tab[sel,]$InfrataxonRank <- tab1[sel,]$InfrataxonRank 
  tab[sel,]$InfrataxonName<- tab1[sel,]$InfrataxonName 
  }
  
  # APG III family
  sel<-!is.na(tab$Family_APGIII) & tab$Family_APGIII==""
  if(any(sel)){tab[sel,]$Family_APGIII<-NA }
    fam.list <- c();
    fam.list=sapply(tab[!is.na(tab$Family_APGIII),]$Family_APGIII,function(x) apg_lookup(taxa =as.character(x), rank = "family")[1])
    tab[!is.na(tab$Family_APGIII),]$Family_APGIII=fam.list
  
  # OthersInfo
  if(othersinfo & any(!is.na(tab$Status_TBGRI)))
  {
    tab$Origin <- NA; tab$Habit <- NA; tab$Phenology <- NA; 
    Info <-sapply(tab[!is.na(tab$Status_TBGRI),]$FoundName,function(x) reftaxo[which(reftaxo$Full.name==x),c("Origin","Habit","Phenology")])
    tab$Origin[!is.na(tab$Status_TBGRI)]<-unlist(Info["Origin",])
    tab$Habit[!is.na(tab$Status_TBGRI)]<-unlist(Info["Habit",]) 
    tab$Phenology[!is.na(tab$Status_TBGRI)]<-unlist(Info["Phenology",]) 
  }
  
  #IUCN Status a revoir
  if(iucn)
    {
    tab$IUCN <- NA; 
    sel<-!is.na(tab$ReferenceName_proposed)
    tab$IUCN[sel]<-as.character(lets.iucn(paste(tab[sel,]$Genus,tab[sel,]$Species, sep=" ") , count = TRUE)$Status)
  }
  
  tab[tab == ""] <- NA
  # Return a table with original names in Rownames, and information on these taxa in other columns
  return(list(SpeciesTable=tab))
}

