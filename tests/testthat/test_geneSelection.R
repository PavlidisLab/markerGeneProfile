context('gene selection')

test_that('single run gene selection',{
    tempDir = tempfile()
    dir.create(tempDir)
    data(mgp_sampleProfilesMeta)
    data(mgp_sampleProfiles)
    data(mpg_sampleRegionHiearchy)
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
                     cores = 1,rotate = 0.33)
    testthat::expect_that(list.files(file.path(tempDir,'markers')),
                testthat::equals(c("All_CellType", "CellType", "Region 2_CellType")))
    markers = pickMarkersAll(file.path(tempDir,'markers'))
    testthat::expect_that(markers$All_CellType, testthat::equals(list(`Cell A` = "Gene1",
                                               `Cell B` = "Gene2",
                                               `Cell C` = c("Gene3", "Gene4", "Gene5"))))

    testthat::expect_that(markers$`Region 2_CellType`, testthat::equals(list(`Cell B` = "Gene2",
                                                                   `Cell C` = c("Gene3", "Gene4", "Gene5"))))
    unlink(tempDir)

})

test_that('gene selection without optional columns and no regions',{
    tempDir = tempfile()
    dir.create(tempDir)
    data(mgp_sampleProfilesMeta)
    data(mgp_sampleProfiles)
    data(mpg_sampleRegionHiearchy)
    mgp_sampleProfilesMeta %<>% dplyr::select(-RegionToParent,-region,-RegionToChildren)

    markerCandidates(design = mgp_sampleProfilesMeta,
                     expression = mgp_sampleProfiles,
                     outLoc = file.path(tempDir,'markers'),
                     groupNames = 'CellType',
                     regionNames = NULL,
                     PMID = 'PMID',
                     sampleName = 'sampleName',
                     replicates = 'replicate',
                     foldChangeThresh = 10,
                     minimumExpression = 8,
                     regionHierarchy = NULL,
                     geneID = 'Gene.Symbol',
                     cores = 1,rotate = 0.33)

    testthat::expect_that(list.files(file.path(tempDir,'markers')),
                          testthat::equals("CellType"))

    markers = pickMarkersAll(file.path(tempDir,'markers'))

    testthat::expect_that(markers$CellType, testthat::equals(list(`Cell A` = "Gene1",
                                                                  `Cell B` = "Gene2",
                                                                  `Cell C` = c("Gene3", "Gene4", "Gene5"))))

    unlink(tempDir)
})

test_that('file read trial',{
    markers =pickMarkers('Cortex_PyramidalDeep')
    expect_true(!"Amn" %in% markers$GabaPV)
})
