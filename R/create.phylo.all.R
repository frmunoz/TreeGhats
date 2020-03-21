create.phylo.all <- function(scenarios = "S3")
# Used once to generate the phylogeny of all the taxa of the database
{
  data('TreeGhatsData', package='TreeGhats', envir=environment())
  data('GBOTB.extended ', package='V.PhyloMaker', envir=environment())
  GBOTB.extended <- get("GBOTB.extended", envir=environment())
  
  TreeGhatsData <- get("TreeGhatsData", envir=environment())
  TreeGhatsData.accept <- TreeGhatsData[(!is.na(TreeGhatsData$Status_TBGRI) | !is.na(TreeGhatsData$Status_TPL)) & ((is.na(TreeGhatsData$Status_TBGRI) & TreeGhatsData$Status_TPL=="Accepted") | TreeGhatsData$Status_TBGRI=="Accepted"),]
  TreeGhatsData.phylo <- create.phylo(TreeGhatsData.accept, scenarios = scenarios)
}