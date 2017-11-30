context('test regionize')

testthat::test_that('regionize with our region hierarchy',{
    load('n_expressoSamples.rda')
    data('mouseRegionHierarchy')
    regionized = regionize(design = n_expressoSamples,
                           regionNames = 'Region',
                           groupNames = 'PyramidalDeep',
                           regionHierarchy = mouseRegionHierarchy)

    cerebOligos = regionized$Cerebellum_PyramidalDeep %>% {.=='Oligo'} %>% which
    cortexOligos = regionized$Cortex_PyramidalDeep %>% {.=='Oligo'} %>% which
    allOligos = regionized$All_PyramidalDeep %>% {.=='Oligo'} %>% which

    testthat::expect_equal(cerebOligos,
                     c(18:23,28:30,34:36))

    testthat::expect_true(!any(cerebOligos %in% cortexOligos))
    testthat::expect_true(all(cerebOligos %in% allOligos))
    testthat::expect_true(all(cortexOligos %in% allOligos))
})
