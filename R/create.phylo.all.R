create.phylo.all <- function()
# Used once to generate the phylogeny of all the taxa of the database
{
  reftaxo.accept <- reftaxo[(!is.na(reftaxo$Status_TBGRI) | !is.na(reftaxo$Status_thePlantlist)) & ((is.na(reftaxo$Status_TBGRI) & reftaxo$Status_thePlantlist=="Accepted") | reftaxo$Status_TBGRI=="Accepted"),]
  reftaxo.phylo <- create.phylo(reftaxo.accept)
}