create.phylo.all <- function()
# Used once to generate the phylogeny of all the taxa of the database
{
  data('TreeGhatsData', package='TreeGhats', envir=environment())
  TreeGhatsData <- get("TreeGhatsData", envir=environment())
  TreeGhatsData.accept <- TreeGhatsData[(!is.na(TreeGhatsData$Status_TBGRI) | !is.na(TreeGhatsData$Status_thePlantlist)) & ((is.na(TreeGhatsData$Status_TBGRI) & TreeGhatsData$Status_thePlantlist=="Accepted") | TreeGhatsData$Status_TBGRI=="Accepted"),]
  TreeGhatsData.phylo <- create.phylo(TreeGhatsData.accept)
}