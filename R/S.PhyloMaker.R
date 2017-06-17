S.PhyloMaker<-function (splist, tree, nodes, output.splist = T, scenarios = "S3") 
{
    splist[sapply(splist, is.factor)] <- lapply(splist[sapply(splist, 
        is.factor)], as.character)
    nodes[sapply(nodes, is.factor)] <- lapply(nodes[sapply(nodes, 
        is.factor)], as.character)
    if (any(duplicated(splist$species))) 
        warning("Duplicated species detected and removed.")
    splist <- splist[!duplicated(splist$species), ]
    Ori <- splist
    splist$species <- gsub(" ", "_", splist$species)
    splist$species <- gsub("(^[[:alpha:]])", "\\U\\1", splist$species, 
        perl = TRUE)
    splist$genus <- gsub("(^[[:alpha:]])", "\\U\\1", splist$genus, 
        perl = TRUE)
    splist$family <- gsub("(^[[:alpha:]])", "\\U\\1", splist$family, 
        perl = TRUE)
    renameNode <- data.frame(node.label = paste("N", 1:length(tree$node.label), 
        sep = ""), original.node.label = tree$node.label, stringsAsFactors = FALSE)
    tree$node.label <- paste("N", 1:length(tree$node.label), 
        sep = "")
    add.tip <- splist[which(is.na(match(splist$species, tree$tip.label))), 
        ]
    status <- rep("match(prune)", dim(splist)[1])
    status[which(is.na(match(splist$species, tree$tip.label)))] <- "match(add)"
    
    # Instead of stopping script if no new species to add, 
    # return the corresponding subtree from the megaphylogeny
    if (dim(add.tip)[1] == 0) 
        #stop("Format mismatch OR There is no species needs to be added.")
      return(drop.tip(tree,tree$tip.label[!tree$tip.label%in%splist$species]))
    
    add.tip$sort <- ""
    add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level == 
        "G", ]$genus)))] <- "G"
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == 
        "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level == 
        "F", ]$family)))] <- "F1"
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == 
        "F1", ]$genus)] <- "F2"
    a <- which(add.tip$sort == "")
    b <- as.character(add.tip$species[a])
    if (length(a) > 0) 
        print(paste("Note:", length(a), "taxa unmatched."))
    status[match(b, splist$species)] <- "unmatch"
    Ori$status <- status
    if ("S1" %in% scenarios) {
        T1 <- tree
        nodeG <- nodes[nodes$level == "G", ]
        nodeF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nodeF$family)
                num <- unique(T1$edge[, 1])[match(nodeF$node.label[n], 
                  T1$node.label)]
                len <- nodeF$node.age[n]
                T1 <- bind.tip(T1, tip.label = data$species[i], 
                  edge.length = len, where = num)
            }
        }
        data <- add.tip[add.tip$sort == "F2", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), 
                  T1$tip.label)
                if (length(n) == 1) {
                  len <- nodeF$node.age[match(data$family[i], 
                    nodeF$family)]/3
                  T1 <- bind.tip(T1, tip.label = data$species[i], 
                    where = n, edge.length = len, position = len)
                }
                if (length(n) > 1) {
                  num <- fastMRCA(T1, T1$tip.label[min(n)], 
                    T1$tip.label[max(n)])
                  len <- T1$edge.length[match(n[1], T1$edge[, 
                    2])]
                  T1 <- bind.tip(T1, tip.label = data$species[i], 
                    where = num, edge.length = len)
                }
            }
        }
        data <- add.tip[add.tip$sort == "G", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$genus[i], nodeG$genus)
                num <- unique(T1$edge[, 1])[match(nodeG$node.label[n], 
                  T1$node.label)]
                len <- nodeG$node.age[n]
                T1 <- bind.tip(T1, tip.label = data$species[i], 
                  edge.length = len, where = num)
            }
        }
        T1$edge.length <- round(T1$edge.length, 5)
        toDrop <- setdiff(1:length(T1$tip.label), which(!is.na(match(T1$tip.label, 
            splist$species))))
        T1 <- drop.tip(T1, tip = toDrop)
        Re <- which(!is.na(match(T1$node.label, renameNode$node.label)))
        noRe <- which(is.na(match(T1$node.label, renameNode$node.label)))
        T1$node.label[Re] <- renameNode$original.node.label[match(T1$node.label, 
            renameNode$node.label)[Re]]
        T1$node.label[noRe] <- ""
    }
    else {
        T1 <- NULL
    }
    if ("S2" %in% scenarios) {
        T2 <- tree
        nodeG <- nodes[nodes$level == "G", ]
        nodeF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nodeF$family)
                f <- nodeF$clade.size[n]
                if (f == 1) {
                  num <- grep(nodeF$taxa[n], T2$tip.label)
                  len <- T2$edge.length[match(num, T2$edge[, 
                    2])] * (2/3)
                  T2 <- bind.tip(T2, tip.label = data$species[i], 
                    edge.length = len, where = num, position = len)
                  nodeF$clade.size[n] <- f + 1
                  num <- grep(nodeF$taxa[n], T2$tip.label)
                  T2$node.label[match(T2$edge[match(num, 
                    T2$edge[, 2]), 1], unique(T2$edge[, 
                    1]))] <- paste("N", T2$Nnode + 1, sep = "")
                  nodeF$node.label[n] <- paste("N", T2$Nnode + 
                    1, sep = "")
                  nodeF$node.age[n] <- len
                  nodeG$node.label[match(nodeF$family[n], nodeG$family)] <- paste("N", 
                    T2$Nnode + 1, sep = "")
                  nodeG$node.age[match(nodeF$family[n], nodeG$family)] <- len
                }
                if (f > 1) {
                  nu <- which(nodeG$family == data$family[i])[sample(1:length(which(nodeG$family == 
                    data$family[i])), 1)]
                  num <- unique(T2$edge[, 1])[match(nodeG$node.label[nu], 
                    T2$node.label)]
                  len <- nodeG$node.age[nu]
                  T2 <- bind.tip(T2, tip.label = data$species[i], 
                    edge.length = len, where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort == "F2", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), 
                  T2$tip.label)
                if (length(n) == 1) {
                  len <- T2$edge.length[match(n, T2$edge[, 
                    2])]/2
                  T2 <- bind.tip(T2, tip.label = data$species[i], 
                    edge.length = len, where = n, position = len)
                }
                if (length(n) > 1) {
                  num <- fastMRCA(T2, T2$tip.label[min(n)], 
                    T2$tip.label[max(n)])
                  len <- T2$edge.length[match(n[1], T2$edge[, 
                    2])]
                  T2 <- bind.tip(T2, tip.label = data$species[i], 
                    edge.length = len, where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort == "G", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                nu <- grep(paste(data$genus[i], "_", sep = ""), 
                  T2$tip.label)
                n <- nu[sample(1:length(nu), 1)]
                num <- T2$edge[, 1][match(n, T2$edge[, 
                  2])]
                len <- T2$edge.length[match(n, T2$edge[, 
                  2])]
                T2 <- bind.tip(T2, tip.label = data$species[i], 
                  edge.length = len, where = num)
            }
        }
        T2$edge.length <- round(T2$edge.length, 5)
        toDrop <- setdiff(1:length(T2$tip.label), which(!is.na(match(T2$tip.label, 
            splist$species))))
        T2 <- drop.tip(T2, tip = toDrop)
        Re <- which(!is.na(match(T2$node.label, renameNode$node.label)))
        noRe <- which(is.na(match(T2$node.label, renameNode$node.label)))
        T2$node.label[Re] <- renameNode$original.node.label[match(T2$node.label, 
            renameNode$node.label)[Re]]
        T2$node.label[noRe] <- ""
    }
    else {
        T2 <- NULL
    }
    if ("S3" %in% scenarios) {
        T3 <- tree
        nodeG <- nodes[nodes$level == "G", ]
        nodeF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nodeF$family)
                f <- nodeF$clade.size[n]
                if (f == 1) {
                  num <- grep(nodeF$taxa[n], T3$tip.label)
                  len <- T3$edge.length[match(num, T3$edge[, 
                    2])] * (2/3)
                  T3 <- bind.tip(T3, tip.label = data$species[i], 
                    edge.length = len, where = num, position = len)
                  nodeF$clade.size[n] <- f + 1
                  num <- grep(nodeF$taxa[n], T3$tip.label)
                  T3$node.label[match(T3$edge[match(num, 
                    T3$edge[, 2]), 1], unique(T3$edge[, 
                    1]))] <- paste("N", T3$Nnode + 1, 
                    sep = "")
                  nodeF$node.label[n] <- paste("N", T3$Nnode + 
                    1, sep = "")
                  nodeF$node.age[n] <- len
                }
                if (f > 1) {
                  num <- unique(T3$edge[, 1])[match(nodeF$node.label[n], 
                    T3$node.label)]
                  len <- nodeF$node.age[n]
                  T3 <- bind.tip(T3, tip.label = data$species[i], 
                    edge.length = len, where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort != "F1", ]
        if (dim(data)[1] > 0) {
            for (i in 1:dim(data)[1]) {
                n <- grep(paste(data$genus[i], "_", sep = ""), 
                  T3$tip.label)
                if (length(n) == 1) {
                  len <- T3$edge.length[match(n, T3$edge[, 
                    2])]/2
                  T3 <- bind.tip(T3, tip.label = data$species[i], 
                    edge.length = len, where = n, position = len)
                  nodeG$clade.size[match(data$genus[i], nodeG$genus)] <- length(n) + 
                    1
                  nu <- grep(paste(data$genus[i], "_", sep = ""), 
                    T3$tip.label)
                  num <- fastMRCA(T3, T3$tip.label[nu[1]], 
                    T3$tip.label[nu[2]])
                  T3$node.label[match(num, unique(T3$edge[, 
                    1]))] <- paste("N", T3$Nnode + 1, 
                    sep = "")
                  nodeG$node.label[match(data$genus[i], nodeG$genus)] <- paste("N", 
                    T3$Nnode + 1, sep = "")
                  nodeG$node.age[match(data$genus[i], nodeG$genus)] <- len
                }
                if (length(n) > 1) {
                  num <- fastMRCA(T3, T3$tip.label[min(n)], 
                    T3$tip.label[max(n)])
                  len <- as.numeric(branching.times(T3))[match(num, 
                    unique(T3$edge[, 1]))]
                  T3 <- bind.tip(T3, tip.label = data$species[i], 
                    edge.length = len, where = num)
                }
            }
        }
        T3$edge.length <- round(T3$edge.length, 
            5)
        toDrop <- setdiff(1:length(T3$tip.label), which(!is.na(match(T3$tip.label, 
            splist$species))))
        T3 <- drop.tip(T3, tip = toDrop)
        Re <- which(!is.na(match(T3$node.label, renameNode$node.label)))
        noRe <- which(is.na(match(T3$node.label, renameNode$node.label)))
        T3$node.label[Re] <- renameNode$original.node.label[match(T3$node.label, 
            renameNode$node.label)[Re]]
        T3$node.label[noRe] <- ""
    }
    else {
        T3 <- NULL
    }
    if (output.splist == FALSE) 
        splist <- NULL
    phylo <- list(Scenario.1 = T1, Scenario.2 = T2, 
        Scenario.3 = T3, Species.list = Ori)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
}
