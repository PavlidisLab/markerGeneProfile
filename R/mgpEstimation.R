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
#' @param removeNegatives logical. should the genes with negative loadings be removed from the estimation. Setting
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
                        removeNegatives = TRUE,
                        geneTransform = function(x){homologene::mouse2human(x)$humanGene}, # function to use when translating gene names
                        comparisons = 'all',
                        pAdjMethod = p.adjust.methods, # method for multiple testing correction. defaults to holm
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
                                 removeNegatives = removeNegatives,
                                 PC = PC,
                                 geneTransform = geneTransform)
    estimates$estimates = ogbox::trimNAs(estimates$estimates)
    estimates$groups = ogbox::trimNAs(estimates$groups)


    if (!is.null(estimateFile) & !outlierSampleRemove){
        estimateFrame = cbind(as.data.frame(estimates$estimates), estimates$groups[[1]])
        names(estimateFrame)[ncol(estimateFrame)] = 'groups'
        write.table(estimateFrame, file=estimateFile , quote=F,sep='\t')
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
plotEstimates = function(estimates,groups,plotNames, sigTest =  wilcox.test,
                         pAdjMethod = p.adjust.methods,
                         comparisons = 'all' # if p value correction should happen in per plot or for the entire list of ps
){
    toCreate = unique(dirname(plotNames))
    sapply(toCreate,dir.create,showWarnings = F,recursive=T)


    groupNames = as.character(unique(groups[[1]]))

    if (typeof(estimates)!='list'){
        estimates = list(estimates)
    }

    estimates %<>% lapply(scale01)

    if (!is.null(comparisons)){
        if (comparisons == 'all'){
            comparisons = combn(groupNames,2)
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
        pList = matrix(p.adjust(pList,pAdjMethod),nrow = length(estimates))
    }


    # plotting of things
    for (i in 1:length(estimates)){
        frame = data.frame(PC1 = estimates[[i]], group = groups[[i]])
        # windowUp = max((frame$PC1)) + 1
        # windowDown = min((frame$PC1)) - 0.5

        # prepare significance text
        if (!is.null(comparisons)){
            sigText = apply(comparisons, 2, paste0, collapse = '/')

            sigText = paste0(sigText,': ',sprintf('%.5f',pList[i,]),collapse = '\n')
        }

        lePlot = ggplot2::ggplot(frame,aes(x=group, y = PC1)) +
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
                           title = element_text(vjust=0.5, size=25),
                           axis.text.y = element_text(size = 13))
        if (!is.null(comparisons)){
            lePlot = lePlot +   # annotate('text', x = 0.1 , y = windowUp ,
                ggplot2::annotate('text', x = 0.1 , y = 1.05 ,
                                  label = sigText ,
                                  hjust = 0, vjust=1, size = 4.5)
        }
        (lePlot)
        ggplot2::ggsave(plotNames[i],width=8,height=8)

    }

}


#' Calculate cell type profile estimations
#' @description Primary function for cell type profile estimations.
#' @param exprData data.frame. Expression data. First collumns of the expression data should include gene names in the
#' same format as the ones specified in the marker gene lists. Any other non-expression related fields must not be of
#' type 'double'
#' @param genes a named list containing marker gene lists of each cell type
#' @param geneColName character. name of the column containing the gene names in the expression file
#' @param outlierSampleRemove logical. If TRUE, outlier samples will be removed from the output. Outliers are calculated in the context of a group
#' @param synonymTaxID Taxonomy identifier of the source of cell type markers. If provided, synonyms of the genes will
#' be added as markers, not recommended since unrelated genes can share names
#' @param geneTransform a function that will be applied to the gene list. the default behavior is to change mouse genes
#' to human genes. set to NULL to keep the genes as they are
#' @param groups a vector stating which groups each sample belongs to
#' @param tableOut character, filename. If provided outputs loadings of individual genes and variance explained by
#' principal components
#' @param indivGenePlot a character vector. If provided, plots expression of marker genes in individual groups per
#' marker gene list. Is not guaranteed to look pretty.
#' @param seekConsensus logical. If TRUE any gene with negative loadings in any of the groups individually will be
#' removed. Use if there is a high likelihood of gene regulation between the groups.
#' @param removeNegatives logical. should the genes with negative loadings be removed from the estimation. Setting
#' seekConsensus to TRUE makes this irrelevant. As all negatives will be removed at that step
#' @param plotType if indivGenePlot is provided, type of plot to be saved. groupBased separates expression between groups
#' cummulative plots a single value
#' @param PC which principal component to use. For debugging purposes. Recommended value is always 1
#' @export
mgpEstimate = function(exprData,
                            genes,
                            geneColName = 'Gene.Symbol',
                            outlierSampleRemove = FALSE,
                            synonymTaxID = NULL, # do you want to add synonyms? no you don't. don't touch this
                            geneTransform = function(x){homologene::mouse2human(x)$humanGene},
                            groups = NULL, # a vector designating the groups. must be defined.
                            tableOut = NULL,
                            indivGenePlot = NULL, # where should it plot individual gene expression plots.
                            seekConsensus = F, # seeking concensus accross groups
                            removeNegatives = TRUE,
                            plotType = c('groupBased','cummulative'), # group based plot requires groups
                            PC = 1){
    if(is.null(groups)){
        assertthat::assert_that(seekConsensus==F)
        list[, exp] = sepExpr(exprData)
        groups = rep(1,ncol(exp))
    }

    if (!is.null(indivGenePlot[1])){
        toCreate = unique(dirname(indivGenePlot))
        sapply(toCreate,dir.create,showWarnings = F,recursive=T)
    }
    if (!is.null(tableOut[1])){
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
                                        geneTransform = geneTransform,
                                        synonymTaxID = synonymTaxID)
    }

    estimateOut = vector(mode = 'list', length = length(genes))
    groupsOut = vector(mode = 'list', length = length(genes))
    rotations = vector(mode = 'list', length = length(genes))
    for (i in 1:length(genes)){
        if(!is.null(geneTransform)){
            genes[[i]] = geneTransform(genes[[i]])
        }
        # add gene synonyms to the list as well. you don't want to do this with gemma annotations
        if (!is.null(synonymTaxID)){
            genes[[i]] == unlist(geneSynonym::geneSynonym(genes=genes[[i]],tax=synonymTaxID))
        }


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

            p = ggplot2::ggplot(indivGenes,aes(y = expression, x = group )) +
                ggplot2::facet_wrap('gene') +
                ggplot2::geom_boxplot(fill = 'lightblue')+ geom_point() +
                ggplot2::theme(axis.text.x  = element_text( size=20,angle=90),
                               axis.title.y = element_text(vjust=0.5, size=20),
                               axis.title.x = element_text(vjust=0.5, size=0) ,
                               title = element_text(vjust=0.5, size=20))+
                ggplot2::scale_x_discrete(name = '')+
                ggplot2::scale_y_discrete(name = 'log2 Expression')
            (p)
            ggplot2::ggsave(filename = indivGenePlot[i], width=8,height=8)
        }

        # get rotations
        pca = prcomp(t(relevantExpr), scale = T)
        pca$rotation = pca$rotation * ((sum(pca$rotation[,PC])<0)*(-2)+1)

        # do not allow negative rotated genes
        if(removeNegatives){
            while(any(pca$rotation[,PC]<0)){
                relevantExpr = relevantExpr[pca$rotation[,PC]>0,]
                pca = prcomp(t(relevantExpr), scale = T)
                pca$rotation = pca$rotation * ((sum(pca$rotation[,PC])<0)*(-2)+1)
            }
        }


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
            write.table(pca$rotation[,PC,drop=F],
                        file = tableOut[i],
                        quote = F, row.names = T, col.names = F, sep='\t',
                        append = T)
        }


        pca$x = t(as.matrix(t(scale(t(relevantExpr))))) %*% as.matrix(pca$rotation)
        groupsOut[[i]] = groups
        # outlier removal
        if (outlierSampleRemove){
            groupData = sapply(unique(groups),function(x){
                pca$x[groups %in% x,PC]
            },simplify=F)
            names(groupData) = unique(groups)
            box = boxplot(groupData, plot = F)
            # because of this part, sample names are important!!! uses them
            # to match outliers
            groupsOut[[i]] = groups[!rownames(pca$x) %in% names(box$out)]
            pca$x = pca$x[!rownames(pca$x) %in% names(box$out),,drop=F]
        }
        estimateOut[[i]]=(pca$x[,PC])
    } # end of main loop over genes
    names(estimateOut) = names(genes)
    names(groupsOut) = names(genes)
    names(rotations) = names(genes)
    output = list(estimates=estimateOut,groups=groupsOut, rotations  = rotations)

    return(output)
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
#' @param synonymTaxID Taxonomy identifier of the source of cell type markers. If provided, synonyms of the genes will
#' be added as markers, not recommended since unrelated genes can share names
#' @export
groupRotations = function(exprData, genes,geneColName, groups, outDir,
                          geneTransform = function(x){mouse2human(x)$humanGene},
                          synonymTaxID = NULL)
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
        if (!is.null(synonymTaxID)){
            genes[[i]] == unlist(geneSynonym(genes=genes[[i]],tax=synonymTaxID))
        }
        relevantData = exprData[exprData[, geneColName] %in% genes[[i]],]
        if (nrow(relevantData)==0){
            allRotations[[i]]=NA
            next
        }

        rownames(relevantData) = relevantData[, geneColName]
        list[,relevantExpr] = ogbox::sepExpr(relevantData)

        for (j in 1:length(unique(groups))){
            pca = prcomp(t(relevantExpr[groups %in% unique(groups)[j]]), scale = T)
            pca$rotation = pca$rotation * ((sum(pca$rotation[,1])<0)*(-2)+1)
            rotations[[j]] = pca$rotation[,1]
        }

        rotations = as.data.frame(rotations)
        names(rotations) = unique(groups)
        allRotations[[i]] = rotations
        if (!is.null(outDir)){
            write.table(rotations[order(apply(rotations,1,sum),decreasing=T),],
                        file = paste0(outDir,'/',names(genes)[i], ' groupRots'), quote=F,sep = '\t')
        }
    }
    names(allRotations) = names(genes)
    invisible(allRotations)
}
