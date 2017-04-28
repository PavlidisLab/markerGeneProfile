context('estimation from sample data')

test_that("Basic estimation works",{
    temp = data.frame(Gene.Symbol = 'EpicMarker', sample1 = 4,sample2 = 8)
    genes = list(epicCellType = 'EpicMarker')
    estimates = mgpEstimate(exprData = temp,genes = genes,geneColName = 'Gene.Symbol', geneTransform = NULL)
    expect_that(estimates$estimates$epicCellType['sample2'],
                is_more_than(estimates$estimates$epicCellType['sample1']))
})


test_that("Estimation in actual dataset",{
    data(mgp_LesnickCroppedExpression)
    data(mgp_LesnickCroppedMeta)
    data(mouseMarkerGenes)
    estimates = mgpEstimate(exprData = mgp_LesnickCroppedExpression,
                genes = mouseMarkerGenes$Midbrain,
                geneColName = 'Gene.Symbol',
                groups = mgp_LesnickCroppedMeta$disease)
    cells = names(mouseMarkerGenes$Midbrain)

    less = sapply(1:length(estimates$estimates),function(i){
        wilcox.test(estimates$estimates[[i]][estimates$groups[[i]] %in% "parkinson's"],
                    estimates$estimates[[i]][estimates$groups[[i]] %in% "control"],alternative = 'less')$p.value
    })

    greater = sapply(1:length(estimates$estimates),function(i){
        wilcox.test(estimates$estimates[[i]][estimates$groups[[i]] %in% "parkinson's"],
                    estimates$estimates[[i]][estimates$groups[[i]] %in% "control"],alternative = 'greater')$p.value
    })

    expect_that(less[cells %in% 'Dopaminergic'], testthat::is_less_than(0.05))
    expect_that(min(less[!cells %in% 'Dopaminergic']), testthat::is_more_than(0.05))
    expect_that(min(greater), testthat::is_more_than(0.05))
})
