# add marker genes from analysis package and create sample data

library(ogbox)
library(magrittr)
library(dplyr)
library(devtools)
library(viridis)
#loadGithub('oganm/neuroExpressoAnalysis/data/mouseMarkerGenes.rda')
load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenes.rda')
load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenesNCBI.rda')

load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenesPyramidalDeep.rda')
load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenesPyramidalDeepNCBI.rda')

load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenesCombined.rda')
load('../../wholeOtto/omancarci/brainGenesManuscript/data/mouseMarkerGenesCombinedNCBI.rda')

use_data(mouseMarkerGenes,overwrite = TRUE)
use_data(mouseMarkerGenesNCBI,overwrite = TRUE)

use_data(mouseMarkerGenesPyramidalDeep,overwrite = TRUE)
use_data(mouseMarkerGenesPyramidalDeepNCBI,overwrite = TRUE)

use_data(mouseMarkerGenesCombined,overwrite = TRUE)
use_data(mouseMarkerGenesCombinedNCBI,overwrite = TRUE)



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


# Lesnick et al full expression data
# loadGithub('oganm/neuroExpressoAnalysis/data/LesnickParkinsonsExp.rda')
load('../../wholeOtto/omancarci/brainGenesManuscript/lesnickPreLowExpression.rda')
loadGithub('oganm/neuroExpressoAnalysis/data/LesnickParkinsonsMeta.rda')

LesnickParkinsonsExp = annotated
mgp_LesnickParkinsonsExp =  LesnickParkinsonsExp %>% select(-GOTerms,-GemmaIDs)


mgp_LesnickParkinsonsMeta = LesnickParkinsonsMeta %>%
    mutate(disease = replaceElement(parkinson,
                                    c('TRUE' = "PD", 'FALSE' = 'Control'))$newVector) %>%
    select(GSM,disease)

use_data(mgp_LesnickParkinsonsExp, overwrite = TRUE)
use_data(mgp_LesnickParkinsonsMeta, overwrite = TRUE)


ogbox::loadGithub('oganm/neuroExpressoAnalysis/data/regionHierarchy.rda')
mouseRegionHierarchy = regionHierarchy
use_data(mouseRegionHierarchy, overwrite = TRUE)

ogbox::sourceGithub('oganm/neuroExpressoAnalysis/R/cellColors.R')
mouseCellColor = cellColors()
use_data(mouseCellColor, overwrite = TRUE)
