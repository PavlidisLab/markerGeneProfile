context('estimation from sample data')

test_that("Basic estimation works",{
    temp = data.frame(Gene.Symbol = 'EpicMarker', sample1 = 4,sample2 = 8)
    genes = list(epicCellType = 'EpicMarker')
    estimates = mgpEstimate(exprData = temp,genes = genes,geneColName = 'Gene.Symbol', geneTransform = NULL)
    expect_that(estimates$estimates$epicCellType['sample2'],
                is_more_than(estimates$estimates$epicCellType['sample1']))
})

test_that("Basic estimation with tibbles",{
    temp = tibble::tibble(Gene.Symbol = 'EpicMarker', sample1 = 4,sample2 = 8)
    genes = list(epicCellType = 'EpicMarker')
    estimates = mgpEstimate(exprData = temp,genes = genes,geneColName = 'Gene.Symbol', geneTransform = NULL)
    expect_that(estimates$estimates$epicCellType['sample2'],
                is_more_than(estimates$estimates$epicCellType['sample1']))
})



test_that("Estimation in actual dataset",{
    data(mgp_LesnickParkinsonsExp)
    data(mgp_LesnickParkinsonsMeta)
    data(mouseMarkerGenes)
    medExp = mgp_LesnickParkinsonsExp %>% sepExpr() %>% {.[[2]]} %>% unlist %>% median
    mgp_LesnickParkinsonsExp = mostVariable(mgp_LesnickParkinsonsExp,threshold = medExp, threshFun= median)

    estimates = mgpEstimate(exprData = mgp_LesnickParkinsonsExp,
                genes = mouseMarkerGenes$Midbrain,
                geneColName = 'Gene.Symbol',
                groups = mgp_LesnickParkinsonsMeta$disease)
    cells = names(mouseMarkerGenes$Midbrain)

    less = sapply(1:length(estimates$estimates),function(i){
        wilcox.test(estimates$estimates[[i]][estimates$groups[[i]] %in% "PD"],
                    estimates$estimates[[i]][estimates$groups[[i]] %in% "Control"],alternative = 'less')$p.value
    })

    greater = sapply(1:length(estimates$estimates),function(i){
        wilcox.test(estimates$estimates[[i]][estimates$groups[[i]] %in% "PD"],
                    estimates$estimates[[i]][estimates$groups[[i]] %in% "Control"],alternative = 'greater')$p.value
    })

    expect_that(less[cells %in% 'Dopaminergic'], testthat::is_less_than(0.05))
    expect_that(min(less[!cells %in% 'Dopaminergic']), testthat::is_more_than(0.05))
    expect_that(min(greater), testthat::is_more_than(0.05))
})

test_that('Wrong parameters',{
    data(mgp_LesnickParkinsonsExp)
    data(mgp_LesnickParkinsonsMeta)
    data(mouseMarkerGenes)
    mockMouseGenes = mgp_LesnickParkinsonsExp$Gene.Symbol %>% homologene::human2mouse()

    mockMouseData = mgp_LesnickParkinsonsExp[match(mockMouseGenes$humanGene,mgp_LesnickParkinsonsExp$Gene.Symbol),]
    mockMouseData$Gene.Symbol = mockMouseGenes$mouseGene
    mockMouseData  = mockMouseData[!duplicated(mockMouseData$Gene.Symbol),]

    testthat::expect_warning(mgpEstimate(exprData = mockMouseData,
                            genes = mouseMarkerGenes$Midbrain,
                            geneColName = 'Gene.Symbol',
                            groups = mgp_LesnickParkinsonsMeta$disease),regexp = 'geneTransform function reduces the number of matches')

})
