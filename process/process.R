# add marker genes from analysis package and create sample data

library(ogbox)
library(magrittr)
library(dplyr)
loadGithub('oganm/neuroExpressoAnalysis/data/mouseMarkerGenes.rda')
use_data(mouseMarkerGenes,overwrite = TRUE)


mgp_sampleProfiles = data.frame(
    Gene.Symbol = c("Gene1",'Gene2','Gene3','Gene4','Gene5','Gene6'),
    Sample01 = c(16,1,1,1,1,1),
    Sample02 = c(16,1,1,1,1,1),
    Sample03 = c(16,1,1,1,1,1),
    Sample04 = c(16,1,1,1,1,1),
    Sample05 = c(16,1,1,1,1,1),
    Sample06 = c(16,1,1,1,1,1),
    Sample07 = c(1,16,1,1,1,1),
    Sample08 = c(1,16,1,1,1,1),
    Sample09 = c(1,16,1,1,1,1),
    Sample10 = c(1,16,1,1,1,1),
    Sample11 = c(1,16,1,1,1,1),
    Sample12 = c(1,16,1,1,1,1),
    Sample13 = c(1,1,16,13,9,7),
    Sample14 = c(1,1,16,13,9,7),
    Sample15 = c(1,1,16,13,9,7),
    Sample16 = c(1,1,16,13,9,7),
    Sample17 = c(1,1,16,13,9,7),
    Sample18 = c(1,1,16,13,9,7))

use_data(mgp_sampleProfiles, overwrite = TRUE)


mgp_sampleProfilesMeta = data.frame(sampleName = paste0('Sample', formatC(1:18,width=2, flag="0")),
                                    replicate = ogbox::repIndiv(1:6,3),
                                    PMID = ogbox::repIndiv(1:6,3),
                                    CellType = ogbox::repIndiv(c('Cell A','Cell B', 'Cell C'),6),
                                    region = ogbox::repIndiv(c('Region 1', 'Region 2', 'Region 2'), 6),
                                    RegionToParent = TRUE,
                                    RegionToChildren = TRUE,
                                    stringsAsFactors = FALSE)
use_data(mgp_sampleProfilesMeta, overwrite = TRUE)


mpg_sampleRegionHiearchy = list(All = list('Region 1' ='',
                                           'Region 2' = ''))


# cropped Lesnick et al
loadGithub('oganm/neuroExpressoAnalysis/data/LesnickParkinsonsExp.rda')
loadGithub('oganm/neuroExpressoAnalysis/data/LesnickParkinsonsMeta.rda')

mgp_LesnickCroppedExpression = LesnickParkinsonsExp[LesnickParkinsonsExp$Gene.Symbol %in%
                                                        (mouseMarkerGenes$Midbrain %>%
                                                             unlist %>%
                                                             (homologene::mouse2human) %$%
                                                             humanGene),] %>%
    select(-GOTerms,-GemmaIDs)

mgp_LesnickCroppedMeta = LesnickParkinsonsMeta %>%
    mutate(disease = replaceElement(parkinson,
                                    c('TRUE' = "parkinson's", 'FALSE' = 'control'))$newVector) %>%
    select(GSM,disease)

use_data(mgp_LesnickCroppedExpression, overwrite = TRUE)
use_data(mgp_LesnickCroppedMeta, overwrite = TRUE)

library(memoise)
ogbox::sourceGithub('oganm/neuroExpressoAnalysis/R/regionize.R',lines = 5:18 )
mouseRegionHiearchy = regionHierarchy
use_data(mouseRegionHiearchy, overwrite = TRUE)
