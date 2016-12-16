context('estimation from sample data')

test_that("Basic estimation works",{
    temp = data.frame(Gene.Symbol = 'EpicMarker', sample1 = 4,sample2 = 8)
    genes = list(epicCellType = 'EpicMarker')
    estimates = mgpEstimate(exprData = temp,genes = genes,geneColName = 'Gene.Symbol', geneTransform = NULL)
    expect_that(estimates$estimates$epicCellType['sample2'],
                is_more_than(estimates$estimates$epicCellType['sample1']))
})


