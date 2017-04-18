context('gene selection')

test_that('single run gene selection',{
    tempDir = tempdir()
    markerCandidates(design = mgp_sampleProfilesMeta,
                     expression = mgp_sampleProfiles,
                     outLoc = file.path(tempDir,'markers'),
                     groupNames = 'CellType',
                     regionNames = 'region',
                     PMID = 'PMID',
                     sampleName = 'sampleName',
                     replicates = 'replicate',
                     foldChangeThresh = 10,
                     minimumExpression = 8,
                     regionHierarchy = mpg_sampleRegionHiearchy,
                     geneID = 'Gene.Symbol',
                     cores = 1)
    expect_that(list.files(file.path(tempDir,'markers')),
                testthat::equals(c("All_CellType", "CellType", "Region 2_CellType")))
    markers = pickMarkersAll(file.path(tempDir,'markers'))
    expect_that(markers$All_CellType, testthat::equals(list(`Cell A` = "Gene1",
                                               `Cell B` = "Gene2",
                                               `Cell C` = c("Gene3", "Gene4", "Gene5"))))

    expect_that(markers$`Region 2_CellType`, testthat::equals(list(`Cell B` = "Gene2",
                                                                   `Cell C` = c("Gene3", "Gene4", "Gene5"))))

})
