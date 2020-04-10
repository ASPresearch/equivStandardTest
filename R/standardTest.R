#' extractGOIDs
#'
#' This is an auxiliary function to extract GOIDs from the output of an enricgment analysis performed with clusterprofiler package
#'
#' @param enriched An object with the output of an enricgment analysis performed with clusterprofiler package
extractGOIDs <- function (enriched) {
  as.character(as.data.frame(enriched)$ID)
}

#' enrichOnto
#'
#' This function performs a standard test of equivalence between two gene lists
#'
#' @param geneList character vector containing a FIRST gene list of entrez IDs
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param orgPackage A string wih the name of the annotation package
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff A numeric value
#' @param qvalCutoff A numeric value
#' @import devtools
#' @import clusterProfiler goProfiles
#' @examples
#'  data(kidneyGeneLists)
#'  data(humanEntrezIDs)
#'  gl1 <-kidneyGeneLists[[1]]
#'  anOnto <- 'BP'
#'  enriched <- enrichOnto (geneL=gl1, geneUniverse=humanEntrezIDs, orgPackage='org.Hs.eg.db', onto=anOnto)
#'  GOIDs <- as.character(as.data.frame(enriched)$ID)
#' @export
enrichOnto <- function (geneList, 
                        geneUniverse, 
                        orgPackage='org.Hs.eg.db', 
                        onto=c("BP", "MF", "CC"),
                         pAdjustMeth="BH", 
                        pvalCutoff=0.01, 
                        qvalCutoff=0.05) {
  clusterProfiler::enrichGO(gene = geneList, universe=geneUniverse, OrgDb =orgPackage, ont=onto,
           pAdjustMethod = pAdjustMeth, pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff, readable=TRUE)
}

#' crossTabGOIDs
#'
#' This function performs a crosstabulation between two lists of enriched GOTerms
#' The lists are intended to have been obtained from enrichment analyses performed on two gene lists
#'
#' @param GO1 character vector containing a FIRST list of GO identifiers
#' @param GO2 character vector containing a SECOND gene list GO identifiers
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param listNames character vector with names of the genelists that generated the enriched GOIDs
#' @examples
#'  data(kidneyEnrichedGOIDs)
#'  GOIDs1 <-kidneyEnrichedGOIDs[[1]]
#'  GOIDs2 <-kidneyEnrichedGOIDs[[2]]
#'  names4lists<- names(kidneyEnrichedGOIDs)[1:2]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  crossTabbedGOIds <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto, 
#'  GOLevel =GOLev, listNames=names4lists)
#'@export
crossTabGOIDs <- function (GO1, GO2, onto, GOLevel, listNames=NULL)
{
  levelIDs <- goProfiles::getGOLevel('BP', level=GOLevel)
  levelInGO1 <- levelIDs %in% GO1
  levelInGO2 <- levelIDs %in% GO2
  if (is.null(listNames)) {
    dnnNames<-paste0("Enriched_in_", c("list_1", "list_2"))
  }else{
    dnnNames<-paste0("Enriched_in_", listNames)
  }
  crossedTable <- table(levelInGO1, levelInGO2, dnn=dnnNames)
}


#' StandardTest4GOIDs
#'
#' This function performs a "standard test" of equivalence between two lists of enriched GOTerms
#' The lists are intended to have been obtained from enrichment analyses performed on two gene lists
#'
#' @param GO1 character vector containing a FIRST list of GO identifiers
#' @param GO2 character vector containing a SECOND gene list GO identifiers
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param listNames character vector with names of the genelists that generated the enriched GOIDs
#' @examples
#'  data(kidneyEnrichedGOIDs)
#'  GOIDs1 <-kidneyEnrichedGOIDs[[1]]
#'  GOIDs2 <-kidneyEnrichedGOIDs[[2]]
#'  names4lists<- names(kidneyEnrichedGOIDs)[1:2]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  testedFromGOIds <- stdTest4GOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto, 
#'  GOLevel =GOLev, listNames=names4lists)
#'@export
stdTest4GOIDs <- function (GO1, GO2, onto, GOLevel, listNames=NULL)
{
  crossTabbedGOIds <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto, 
                                     GOLevel =GOLev, listNames=names4lists)
  fisher.test (crossTabbedGOIds, alt="g")
}

#' crossTab4GeneLists
#'
#' This function performs a "standard test" of equivalence between two gene lists
#'
#' @param geneList1 character vector containing a FIRST gene list of entrez IDs
#' @param geneList2 character vector containing a SECOND gene list of entrez IDs
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param orgPackage A string wih the name of the annotation package
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff A numeric value
#' @param qvalCutoff A numeric value
#' @examples
#'  data(kidneyGeneLists)
#'  data(humanEntrezIDs)
#'  gl1 <-kidneyGeneLists[[1]]
#'  gl2 <-kidneyGeneLists[[2]]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  adjMeth<- 'BH'
#'  pValCut <- 0.05
#'  qValCut <- 0.01
#'  crossTabFromGeneLists <- crossTab4GeneLists (genelist1=gl1, genelist2=gl2,
#'  geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db", onto=anOnto, GOLevel=GOLev)
#'
#'@export 
crossTabGOIDs4GeneLists <- function (genelist1, genelist2, geneUniverse, orgPackg, onto, GOLevel,
                               pAdjustMeth="BH", pvalCutoff=0.01, qvalCutoff=0.05)
{
  enriched1 <- clusterProfiler::enrichGO(gene=genelist1, universe=geneUniverse, OrgDb=orgPackg, ont=anOnto, pvalueCutoff=pvalCutoff)
  GOIDs1 <- as.character(as.data.frame(enriched1)$ID)
  enriched2 <- clusterProfiler::enrichGO(gene=genelist2, universe=geneUniverse, OrgDb=orgPackg, ont=anOnto, pvalueCutoff=pvalCutoff)
  GOIDs2 <- as.character(as.data.frame(enriched2)$ID)
  crossTabGOIDs4GeneLists <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto, GOLevel =GOLevel)
}



#' StandardTest4GeneLists
#'
#' This function performs a "standard test" of equivalence between two gene lists
#'
#' @param geneList1 character vector containing a FIRST gene list of entrez IDs
#' @param geneList2 character vector containing a SECOND gene list of entrez IDs
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param orgPackage A string wih the name of the annotation package
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff A numeric value
#' @param qvalCutoff A numeric value
#' @examples
#'  data(kidneyGeneLists)
#'  data(humanEntrezIDs)
#'  gl1 <-kidneyGeneLists[[1]]
#'  gl2 <-kidneyGeneLists[[2]]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  adjMeth<- 'BH'
#'  pValCut <- 0.05
#'  qValCut <- 0.01
#'  testedFromGeneLists <- stdTest4GeneLists (genelist1=gl1, genelist2=gl2,
#'  geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db", onto=anOnto, GOLevel=GOLev)
#'
#'@export 
stdTest4GeneLists <- function (genelist1, genelist2, geneUniverse, orgPackg, onto, GOLevel,
                               pAdjustMeth="BH", pvalCutoff=0.01, qvalCutoff=0.05)
{
  enriched1 <- enrichOnto(gene=genelist1, universe=geneUniverse, OrgDb=orgPackg, ont=anOnto, pvalueCutoff=pvalCutoff)
  GOIDs1 <- as.character(as.data.frame(enriched1)$ID)
  enriched2 <- enrichOnto(gene=genelist2, universe=geneUniverse, OrgDb=orgPackg, ont=anOnto, pvalueCutoff=pvalCutoff)
  GOIDs2 <- as.character(as.data.frame(enriched2)$ID)
  stdTest<- stdTest4GOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto, GOLevel =GOLevel)
}

