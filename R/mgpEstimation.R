#' Full cell type profile estimation pipeline with plot outputs.
#' @description Calculates cell type profiles in whole tissue datasets using marker genes provided. Fine tunes the
#' the calculation based on experimental groups if needed. This is only useful for quickly plotting profile estimations
#' of data with discrete experimental groups
#' @param exprData data.frame. Expression data. First collumns of the expression data should include gene names in the same format
#' as the ones specified in the marker gene lists. Any other non-expression related fields must not be of type 'double'
#' @param genes a named list containing marker gene lists of each cell type
#' @param geneColName character. name of the column containing the gene names in the expression file
#' @param groups a vector stating which groups each sample belongs to
#' @param outDir output directory for plots and tables
#' @param seekConsensus logical. If TRUE any gene with negative loadings in any of the groups individually will be
#' removed. Use if there is a high likelihood of gene regulation between the groups.
#' @param groupRotations logical. should the outputs include loadings calculated for individual genes
#' @param outlierSampleRemove logical. should the outlier samples be removed from the final output
#' @param removeMinority logical. If TRUE, decides which sign is the most common and removes genes from the minority sign.
#' Note that results will always be rotated in the positive direction Setting
#' seekConsensus to TRUE makes this irrelevant.
#' @param geneTransform a function that will be applied to the gene list. the default behavior is to change mouse genes
#' to human genes. set to NULL to keep the genes as they are
#' @param comparisons a 2 x n character matrix, just 'all' or NULL. Names of groups to compare when calculating p values
#' For each column the two groups indicated in the rows will be compared by wilcox.test.
#' @param pAdjMethod character. which method to use when adjusting p values for multiple testing correction
#' @param PC integer. which principal component to use when calculating cell type profile.
#' @param estimateFile name of the file that contains the final profile estimations
#' @seealso \code{\link{mgpEstimate}}
#' @export
fullEstimate = function(exprData, # expression data
                        genes, # a list of gene lists containing marker genes for cell types
                        geneColName, # collumn name in the expression data that contains gene names
                        groups, # a vector stating which groups each sample belongs to
                        outDir, # output directory for plots and tables
                        seekConsensus=FALSE, # for trimming genes that behave differently in different groups
                        groupRotations=FALSE, # output rotations of individiual genes
                        outlierSampleRemove=FALSE, # if T outliers in each sample is removed
                        removeMinority = TRUE,
                        geneTransform = function(x){homologene::mouse2human(x)$humanGene}, # function to use when translating gene names
                        comparisons = 'all',
                        pAdjMethod = stats::p.adjust.methods, # method for multiple testing correction. defaults to holm
                        PC = 1, # which PC to use. mostly you want this to be 1
                        estimateFile = NULL
){
    estimates = mgpEstimate(exprData=exprData,
                            genes=genes,
                            geneColName=geneColName,
                            outlierSampleRemove=outlierSampleRemove,
                            groups=groups,
                            tableOut = paste0(outDir,'/',names(genes),' rotTable.tsv'),
                            indivGenePlot= paste0(outDir,'/',names(genes),' indivExp','.png'),
                            seekConsensus = seekConsensus,
                            removeMinority = removeMinority,
                            PC = PC,
                            geneTransform = geneTransform)
    estimates$estimates = ogbox::trimNAs(estimates$estimates)
    estimates$groups = ogbox::trimNAs(estimates$groups)


    if (!is.null(estimateFile) & !outlierSampleRemove){
        estimateFrame = cbind(as.data.frame(estimates$estimates), estimates$groups[[1]])
        names(estimateFrame)[ncol(estimateFrame)] = 'groups'
        utils::write.table(estimateFrame, file=estimateFile , quote=F,sep='\t')
    }

    if (groupRotations){
        groupRotations(exprData,
                       geneTransform=geneTransform,
                       genes=genes,
                       geneColName = geneColName,
                       groups=groups,
                       outDir=outDir)
    }

    plotEstimates(estimates$estimates,estimates$groups,
                  paste0(outDir,'/',
                         names(estimates$estimates),'.png'),
                  pAdjMethod = pAdjMethod,
                  comparisons = comparisons)
}






#' Plot cell type profile estimates
#' @description A convenience function to plot the estimates
#' @param estimates. a named list of cell type estimates. $estimates part of the mgpEstimate function's output.
#' @param groups A names list of groups $groups part pf the mgpEstimate function's output
#' @param plotNames a vector of strings. full address of the plots
#' @param sigTest function to be used in calculation of p values for comparisons
#' @param pAdjMethod character. which method to use when adjusting p values for multiple testing correction
#' @param comparisons a 2 by n matrix that includes group name pairs, indicating which groups should be compared when
#' calculating the p values. or "all" which will make comparison between all available groups by setting it to
#' combn(groups,2)
#' @export
plotEstimates = function(estimates,groups,plotNames= NULL, sigTest =  wilcox.test,
                         pAdjMethod = stats::p.adjust.methods,
                         comparisons = 'all' # if p value correction should happen in per plot or for the entire list of ps
){
    if(!is.null(plotNames)){
        toCreate = unique(dirname(plotNames))
    }
    sapply(toCreate,dir.create,showWarnings = F,recursive=T)


    groupNames = as.character(unique(groups[[1]]))

    if (typeof(estimates)!='list'){
        estimates = list(estimates)
    }

    estimates %<>% lapply(scale01)

    if (!is.null(comparisons)){
        if (comparisons == 'all'){
            comparisons = utils::combn(groupNames,2)
        }
    }
    # create p value lists for correction
    if (!is.null(comparisons)){
        pList = matrix(data = NA, ncol = ncol(comparisons), nrow = length(estimates))
        for (i in 1:length(estimates)) {

            for (j in 1:ncol(comparisons)){
                pList[i, j] = sigTest(estimates[[i]][groups[[i]] %in% comparisons[1,j]],
                                      estimates[[i]][groups[[i]] %in% comparisons[2,j]])$p.value
            }
        }

        # p value adjustment
        pList = matrix(stats::p.adjust(pList,pAdjMethod),nrow = length(estimates))
    }


    # plotting of things
    plots = list()
    for (i in 1:length(estimates)){
        frame = data.frame(PC1 = estimates[[i]], group = groups[[i]])
        # windowUp = max((frame$PC1)) + 1
        # windowDown = min((frame$PC1)) - 0.5

        # prepare significance text
        if (!is.null(comparisons)){
            sigText = apply(comparisons, 2, paste0, collapse = '/')

            sigText = paste0(sigText,': ',sprintf('%.5f',pList[i,]),collapse = '\n')
        }

        lePlot = ggplot2::ggplot(frame,ggplot2::aes(x=group, y = PC1)) +
            ggplot2::geom_violin( color="#C4C4C4", fill="#C4C4C4") +
            ggplot2::geom_boxplot(width=0.1,fill = 'lightblue') +
            ggplot2::geom_point(size = 3) +
            ggplot2::ggtitle(names(estimates)[i]) +
            #scale_y_continuous(limits=c(windowDown, windowUp),
            ggplot2::scale_y_continuous(limits=c(-.05, 1.05),
                                        name="Relative estimate of cell type amounts") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x  = ggplot2::element_text(size=25, angle=90),
                           axis.title.y = ggplot2::element_text(vjust=0.5, size=25),
                           axis.title.x = ggplot2::element_text(vjust=0.5, size=0) ,
                           title = ggplot2::element_text(vjust=0.5, size=25),
                           axis.text.y = ggplot2::element_text(size = 13))
        if (!is.null(comparisons)){
            lePlot = lePlot +   # annotate('text', x = 0.1 , y = windowUp ,
                ggplot2::annotate('text', x = 0.1 , y = 1.05 ,
                                  label = sigText ,
                                  hjust = 0, vjust=1, size = 4.5)
        }
        plots[i] = lePlot
        if(is.null(plotNames)){
            ggplot2::ggsave(plotNames[i],plot = lePlot, width=8,height=8)
        }
    }
    names(plots) = names(estimates)
    return(plots)
}


#' Calculate cell type profile estimations
#' @description Primary function for cell type profile estimations.
#' @param exprData data.frame. Expression data. First collumns of the expression data should include gene names in the
#' same format as the ones specified in the marker gene lists. Any other non-expression related fields must not be of
#' type 'double'
#' @param genes a named list containing marker gene lists of each cell type
#' @param geneColName character. name of the column containing the gene names in the expression file
#' @param outlierSampleRemove logical. Deprecated. If TRUE, outlier samples will be removed from the output. Outliers are calculated in the context of a group
#' @param geneTransform a function that will be applied to the gene list. the default behavior is to change mouse genes
#' to human genes. set to NULL to keep the genes as they are
#' @param groups a vector stating which groups each sample belongs to
#' @param tableOut character, directory name. If provided outputs loadings of individual genes and variance explained by
#' principal components
#' @param indivGenePlot character, directory name. If provided, plots expression of marker genes in individual groups per
#' marker gene list. Is not guaranteed to look pretty.
#' @param seekConsensus logical. If TRUE any gene with negative loadings in any of the groups individually will be
#' removed. Use if there is a high likelihood of gene regulation between the groups.
#' @param removeMinority logical. should the genes with negative loadings be removed from the estimation. Setting
#' seekConsensus to TRUE makes this irrelevant. As all negatives will be removed at that step
#' @param plotType if indivGenePlot is provided, type of plot to be saved. groupBased separates expression between groups
#' cummulative plots a single value
#' @param PC which principal component to use. For debugging purposes. Recommended value is always 1
#' @return A list.
#' \itemize{
#' \item estimates. A named vector of marker gene profiles estimated bythe function.
#' \item groups. if provided this lists the groups that the samples belong to. Unless \code{outlierSampleRemove} is set to \code{TRUE}, this will be the same as the groups for all cell types. So in most cases ignore this.
#' \item rotations. Loadings from the PCA. These values are modified to ensure the chosen PC (1st by default) has a positive direction.
#' \item trimmedPCAs. Raw PCA object used to calculate the estimates. Rotation direction is not changed
#' \item fullPCAs. Raw PCA object with all marker genes before their removal by \code{removeMinority} and \code{seekConsensus} options.
#' \item removedMarkerRatios. Ratio of markers removed by  \code{removeMinority} and \code{seekConsensus} options.
#' }
#' @export
mgpEstimate = function(exprData,
                       genes,
                       geneColName = 'Gene.Symbol',
                       outlierSampleRemove = FALSE,
                       geneTransform = function(x){homologene::mouse2human(x)$humanGene},
                       groups = NULL, # a vector designating the groups. must be defined.
                       tableOut = NULL,
                       indivGenePlot = NULL, # where should it plot individual gene expression plots.
                       seekConsensus = FALSE, # seeking concensus accross groups
                       removeMinority = TRUE,
                       plotType = c('groupBased','cummulative'), # group based plot requires groups
                       permuations = 0,
                       PC = 1){
    if(exprData[[geneColName]] %>% duplicated %>% any){
        warning('You have duplicate genes in your expression data. Function will fail if marker genes have duplicates. Please summarize your data to gene level.')
    }
    browser()

    if(!is.null(geneTransform)){
        transformedGenes = genes %>% lapply(geneTransform)
        dataGenes = exprData[[geneColName]]

        if(sum(unlist(transformedGenes) %in% dataGenes) < (sum(unlist(genes) %in% dataGenes))){
            warning('geneTransform function reduces the number of matches between marker gene list and the expression data. Default method is to transform mouse genes to human genes. If this is not what you want please set geneTransform to NULL. Please read to documentary to make sure you are doing it right.')
        }
        genes = transformedGenes

    }

    if(is.null(groups)){
        assertthat::assert_that(seekConsensus == FALSE & is.null(groups))
        list[, exp] = ogbox::sepExpr(exprData)
        groups = rep(1,ncol(exp))
    }
    if (!is.null(indivGenePlot[1])){
        if(length(indivGenePlot) == 1){
            indivGenePlot =  paste0(indivGenePlot,'/',names(genes), 'indivExp.png')
        }
        toCreate = unique(dirname(indivGenePlot))
        sapply(toCreate,dir.create,showWarnings = F,recursive=T)
    }
    if (!is.null(tableOut[1])){
        if(length(tableOut) == 1){
            tableOut =  paste0(tableOut,'/',names(genes),' rotTable.tsv')
        }
        toCreate = unique(dirname(tableOut))
        sapply(toCreate,dir.create,showWarnings = F,recursive=T)
    }
    # if a single vector is inputted, fix that
    if (typeof(genes)!='list'){
        genes = list(genes)
    }
    # if seeking consensus, check group based rotations
    if(seekConsensus){
        groupRotations = groupRotations(exprData, genes,
                                        geneColName, groups, outDir=NULL,
                                        geneTransform = NULL)
    }

    estimateOut = vector(mode = 'list', length = length(genes))
    groupsOut = vector(mode = 'list', length = length(genes))
    rotations = vector(mode = 'list', length = length(genes))

    trimmedPCAs = vector(mode = 'list', length = length(genes))
    fullPCAs = vector(mode = 'list', length = length(genes))
    removedMarkerRatios = vector(length = length(genes))

    meanUsedMarkerExpression = vector(mode = 'list', length = length(genes))
    usedMarkerExpression = vector(mode = 'list', length = length(genes))
    simpleScaledEstimation = vector(mode = 'list', length = length(genes))


    for (i in 1:length(genes)){
        # some repetition here because I want to capture the original PCA's before removing anything else/
        relevantData = exprData[exprData[, geneColName] %in% genes[[i]],]
        rownames(relevantData) = relevantData[, geneColName]
        list[,relevantExpr] = ogbox::sepExpr(relevantData)
        if (nrow(relevantData)==0){
            estimateOut[[i]]=NA
            groupsOut[[i]]=NA
            next
        }

        pca = stats::prcomp(t(relevantExpr), scale = TRUE)
        originalPCA = pca
        fullPCAs[[i]] = originalPCA
        #remove non concenting genes (based on group) if is nested because there
        #will be no group rotations to look at if seekConsensus=F some redundancy
        # exists but oh well
        if (seekConsensus){
            if (!is.na(groupRotations[i])){
                genes[[i]] = rownames(groupRotations[[i]][apply(groupRotations[[i]],1,function(x){all(x>0)}),])
            }
        }

        relevantData = exprData[exprData[, geneColName] %in% genes[[i]],]
        # what to do if no genes from the list is found
        if (nrow(relevantData)==0){
            estimateOut[[i]]=NA
            groupsOut[[i]]=NA
            next
        }

        rownames(relevantData) = relevantData[, geneColName]
        list[,relevantExpr] = ogbox::sepExpr(relevantData)


        if (!is.null(indivGenePlot[1])){
            tempExp = relevantExpr
            tempExp$gene = rownames(tempExp)
            indivGenes = reshape2::melt(tempExp)
            names(indivGenes) = c( 'gene','GSM','expression')

            indivGenes = indivGenes %>%
                dplyr::mutate(group = groups[match(GSM, colnames(tempExp))]) %>%
                # order based on the order of input
                dplyr::mutate(gene = factor(gene, levels = genes[[i]]) %>% droplevels)


            switch(plotType[1],
                   groupBased = {
                   },
                   cummulative = {
                       indivGenes$group =''
                   })

            p = ggplot2::ggplot(indivGenes,ggplot2::aes(y = expression, x = group )) +
                ggplot2::facet_wrap('gene') +
                ggplot2::geom_boxplot(fill = 'lightblue')+ ggplot2::geom_point() +
                ggplot2::theme(axis.text.x  = ggplot2::element_text( size=20,angle=90),
                               axis.title.y = ggplot2::element_text(vjust=0.5, size=20),
                               axis.title.x = ggplot2::element_text(vjust=0.5, size=0) ,
                               title = ggplot2::element_text(vjust=0.5, size=20))+
                ggplot2::scale_x_discrete(name = '')+
                ggplot2::scale_y_discrete(name = 'log2 Expression')
            (p)
            ggplot2::ggsave(filename = indivGenePlot[i], width=8,height=8)
        }

        # get rotations
        pca = stats::prcomp(t(relevantExpr), scale = TRUE)
        trimmedPCA = pca
        pca$rotation = pca$rotation * ((sum(pca$rotation[,PC])<0)*(-2)+1)

        # do not allow negative rotated genes
        if(removeMinority){
            while(any(pca$rotation[,PC]<0)){
                relevantExpr = relevantExpr[pca$rotation[,PC]>0,]
                pca = stats::prcomp(t(relevantExpr), scale = T)
                trimmedPCA = pca
                pca$rotation = pca$rotation * ((sum(pca$rotation[,PC])<0)*(-2)+1)
            }
        }

        trimmedPCAs[[i]] = trimmedPCA

        ngenes = trimmedPCA$rotation %>% nrow
        ngenesBefore = originalPCA$rotation %>% nrow

        removedMarkerRatios[[i]] = (ngenesBefore - ngenes)/ngenesBefore

        # get the final rotations
        rotations[[i]] = pca$rotation

        # output the rotation tables
        if (!is.null(tableOut[1])){

            # add variation explained as a comment on top
            file.create(tableOut[i])
            fileConn = file(tableOut[i])
            writeLines(paste0('# Variation explained: ',
                              paste(summary(pca)$importance[2,], collapse=' ')),
                       fileConn)
            close(fileConn)
            utils::write.table(pca$rotation[,PC,drop=F],
                               file = tableOut[i],
                               quote = F, row.names = T, col.names = F, sep='\t',
                               append = T)
        }


        pca$x = t(as.matrix(t(scale(t(relevantExpr))))) %*% as.matrix(pca$rotation)

        # experimental simple scaling -----
        simpleScaleWeights = ogbox::scale01(pca$rotation[,PC])/(sum(ogbox::scale01(pca$rotation[,PC])))
        simpleScale =  apply((relevantExpr * simpleScaleWeights),2,sum)
        simpleScaledEstimation[[i]] = simpleScale # this is an experimental output

        markerGeneExpression = relevantExpr


        usedMarkerExpression[[i]] = relevantExpr
        meanUsedMarkerExpression[[i]] = relevantExpr %>% apply(2,mean)

        groupsOut[[i]] = groups
        # outlier removal
        if (outlierSampleRemove){
            warning('Outlier sample removal is deprecated and will be removed. If you were using it and this caused problems for you this please mail ogan.mancarcii@gmail.com')
            groupData = sapply(unique(groups),function(x){
                pca$x[groups %in% x,PC]
            },simplify=F)
            names(groupData) = unique(groups)
            box = graphics::boxplot(groupData, plot = F)
            # because of this part, sample names are important!!! uses them
            # to match outliers
            groupsOut[[i]] = groups[!rownames(pca$x) %in% names(box$out)]
            pca$x = pca$x[!rownames(pca$x) %in% names(box$out),,drop=F]
        }
        estimateOut[[i]]=(pca$x[,PC])
    } # end of main loop over genes
    names(usedMarkerExpression) =  names(genes)
    names(meanUsedMarkerExpression) =  names(genes)

    names(estimateOut) = names(genes)
    names(groupsOut) = names(genes)
    names(rotations) = names(genes)

    names(trimmedPCAs) = names(genes)
    names(fullPCAs) = names(genes)
    names(removedMarkerRatios) = names(genes)
    names(simpleScaledEstimation) = names(genes)

    problematic = which(removedMarkerRatios > 0.4) %>% names
    if(length(problematic)>0){
        warning(
            paste0('Several cell types have a high proportion (>0.4) of their genes filtered out. Excercise caution.\nProblematic cell types are: ',
                   paste(problematic,collapse = ', ')))
    }

    output = list(estimates=estimateOut,
                  groups=groupsOut,
                  rotations  = rotations,
                  trimmedPCAs = trimmedPCAs ,
                  fullPCAs = fullPCAs,
                  removedMarkerRatios = removedMarkerRatios,
                  usedMarkerExpression = usedMarkerExpression,
                  meanUsedMarkerExpression = meanUsedMarkerExpression,
                  simpleScaledEstimation = simpleScaledEstimation)

    return(output)
}


matchingGenes = function(genes,exprData, geneColName = 'Gene.Symbol',fun = median, tolerance=0.05){
    medExpression = exprData %>% ogbox::sepExpr() %>% {.[[2]]} %>%  apply(1,fun)
    names(medExpression) = exprData[[geneColName]]

    ogbox::pickRandom(labels = genes, allValues = medExpression, tolerance = tolerance)

}





#' Calculates rotations based on each group
#' @param exprData data.frame. Expression data. First collumns of the expression data should include gene names in the
#' same format as the ones specified in the marker gene lists. Any other non-expression related fields must not be of
#' type 'double'
#' @param genes a named list containing marker gene lists of each cell type
#' @param geneColName character. name of the column containing the gene names in the expression file
#' @param groups a vector stating which groups each sample belongs to
#' @param outDir if provided group rotations will be saved there
#' @param geneTransform a function that will be applied to the gene list. the default behavior is to change mouse genes
#' to human genes. set to NULL to keep the genes as they are
#' @export
groupRotations = function(exprData, genes,geneColName, groups, outDir,
                          geneTransform = function(x){homologene::mouse2human(x)$humanGene})
{
    if (typeof(genes)!='list'){
        genes = list(genes)
    }

    allRotations = vector(mode = 'list', length=length(unique(genes)))
    for (i in 1:length(genes)){

        rotations = vector(mode = 'list', length=length(unique(groups)))

        if(!is.null(geneTransform)){
            genes[[i]] = geneTransform(genes[[i]])
        }

        relevantData = exprData[exprData[, geneColName] %in% genes[[i]],]
        if (nrow(relevantData)==0){
            allRotations[[i]]=NA
            next
        }

        rownames(relevantData) = relevantData[, geneColName]
        list[,relevantExpr] = ogbox::sepExpr(relevantData)

        for (j in 1:length(unique(groups))){
            pca = stats::prcomp(t(relevantExpr[groups %in% unique(groups)[j]]), scale = T)
            pca$rotation = pca$rotation * ((sum(pca$rotation[,1])<0)*(-2)+1)
            rotations[[j]] = pca$rotation[,1]
        }

        rotations = as.data.frame(rotations)
        names(rotations) = unique(groups)
        allRotations[[i]] = rotations
        if (!is.null(outDir)){
            utils::write.table(rotations[order(apply(rotations,1,sum),decreasing=T),],
                               file = paste0(outDir,'/',names(genes)[i], ' groupRots'), quote=F,sep = '\t')
        }
    }
    names(allRotations) = names(genes)
    invisible(allRotations)
}
