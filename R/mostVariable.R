#' Select most variable probeset
#'
#' @description This function selects the most variable probeset each gene in an expression data. It also filters expression data of low expressed genes
#' @param x Expression dataset to filter. A data frame. Gene names must be present as a column before the expression data
#' @param genes Character. Name of the column that lists gene symbols
#' @param threshold If filtering should be done the threshold of expression
#' @param threshFun If filtering should be done, function to apply to the gene row to compare with the threshold. By default genes are removed if maximum expression is below 6
#' @export
mostVariable = function(x,genes = 'Gene.Symbol', threshold = 6, threshFun = max){
    list[,exprData]= sepExpr(x)
    rowmax = apply(exprData, 1, threshFun)
    discludeGenes = (rowmax<threshold)
    x = x[!discludeGenes,]
    exprData = exprData[!discludeGenes,]

    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    x = x[decreasingVar,]
    if (class(x)[1]=='data.table'){
        x = x[!duplicated(x[,genes, with=F]),]
    } else {
        x = x[!duplicated(x[,genes]),]

    }
    x = x[!x[,genes]=='',]
    return(x)
}
