#' Write candidate marker genes to a file
#'
#' Writes candidate marker genes to a file along with other required information (fold change and silhouette coefficient)
#' Selects candidates based on fold change to the median expression of other samples and a minimum expression in the
#' cell type
#'
#' @param design data.frame. Metadata for the samples.
#' @param expression data.frame. Expression data for the samples. Gene names should be included as a column. Any other non expression
#' data not be of type \code{double}
#' @param outLoc Directory to save the candidate genes
#' @param groupNames column of the \code{design} with cell type names. If multiple columns are provided, selection will
#' be performed for each independently
#' @param regionNames column of the \code{design} with region names. Multiple regions can be listed separated by commas
#' @param PMID column of the \code{design} with pubmed identifiers. This is required to identify the samples coming
#' from the same study
#' @param rotate double. percentage of samples to be removed. 0.33 is used for the study
#' @param cores Number of cores to use for paralellization
#' @param debug Nothing to see here
#' @param sampleName column of the \code{design} with sample names matching the column names of the \code{expression} data
#' @param replicates column of the \code{design} that groups replicates
#' @param foldChangeThresh minimum fold change required for selection
#' @param minimumExpression minimum level of expression for a marker gene in its cell type
#' @param background level of expression that should be considered as background level
#' @param regionHiearchy hiearchy of regions.
#' @param seed seed for random generation. if NULL will be set to random
#' @export
markerCandidates = function(design,
                            expression,
                            outLoc,
                            groupNames,
                            regionNames=NULL,
                            PMID = 'PMID',
                            rotate = NULL,
                            cores = 1,
                            debug=NULL,
                            sampleName = 'sampleName',
                            replicates = 'originalIndex',
                            foldChangeThresh = 10,
                            minimumExpression = 8,
                            background = 6,
                            regionHierarchy = NULL,
                            geneID = 'Gene.Symbol',
                            seed = NULL){
    # source('R/regionHierarchy.R')
    # so that I wont fry my laptop
    if (parallel::detectCores()<cores){
        cores = parallel::detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }

    cl<-parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    #gene selector, outputs selected genes and their fold changes
    foldChange = function (group1, group2, f = 10){


        groupAverage1 = group1



        groupAverage2 = tryCatch({apply(group2, 2, median)},
                                 error = function(cond){
                                     print('fuu')
                                     return(group2)
                                 })

        g19 = groupAverage1 < (log(f,base=2) + background) & groupAverage1 > minimumExpression
        g16 = groupAverage1  < background
        g29 = groupAverage2 < (log(f,base=2) +background) & groupAverage2 > minimumExpression
        g26 = groupAverage2 < background
        # this is a late addition preventing anything that is below 8 from being
        # selected. ends up removing the the differentially underexpressed stuff as well
        gMinTresh = groupAverage1 > minimumExpression


        tempGroupAv2 = vector(length = length(groupAverage2))

        tempGroupAv2[g26 & g19] =apply(group2[, g26 & g19,drop=F], 2, max)
        # legacy
        tempGroupAv2[g16 & g29] =apply(group2[, g16 & g29,drop=F], 2, min)



        add1 = g19 & g26 & groupAverage1>tempGroupAv2
        # add2 = g29 & g16 & tempGroupAv2>groupAverage1


        fold = groupAverage1 - groupAverage2
        # take everything below 6 as the same when selecting
        # fold =  sapply(groupAverage1,max,6) - sapply(groupAverage2,max,6)
        chosen =  which(({(fold >= (log(f)/log(2))) & !(g19 & g26) } | add1 )&gMinTresh)
        return(
            data.frame(index = chosen, foldChange = fold[chosen])
        )
    }

    giveSilhouette = function(daGeneIndex,groupInfo1,groupInfo2){
        compareGroups  = c(groupInfo2, list(all = unlist(groupInfo2)))

        silos = sapply(compareGroups, function(groupInfo2){
            clustering = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo1)+1)
            clustering = clustering[1:nrow(design) %in% c(groupInfo1, groupInfo2)]
            data = (exprData[ (1:nrow(design) %in% c(groupInfo1, groupInfo2)),  daGeneIndex])
            cluster = list(clustering = clustering, data = data)
            silo = cluster::silhouette(cluster,stats::dist(data))
            return(mean(silo[,3]))
        })

        mainSil = silos[length(silos)]

        minSil = min(silos[-length(silos)])

        return(c(mainSil,minSil))
    }
    # data prep. you transpose exprData -----
    #design = read.design(designLoc)

    #expression = read.csv(exprLoc, header = T)
    list[geneData, exprData] = ogbox::sepExpr(expression)

    if (!all(colnames(exprData) %in% design[[sampleName]])){
        if(is.null(rotate)){
            print('Unless you are rotating samples, something has gone terribly wrong!')
        }
        exprData = exprData[,colnames(exprData) %in% design[[sampleName]]]
    }

    design = design[match(colnames(exprData),design[[sampleName]]),]

    exprData = t(exprData)
    noReg = F
    if (is.null(regionNames)){
        regionNames = 'dummy'
        design[,regionNames] = 'dummy'
        noReg = T
    }


    regionGroups = memoReg(design,regionNames,groupNames,regionHierarchy)
    # concatanate new region based groups to design and to groupNames so they'll be processed normally
    if (!noReg){
        design = cbind(design,regionGroups)
        groupNamesEn = c(groupNames, names(regionGroups))
    } else {
        groupNamesEn = groupNames
    }

    # generate nameGroups to loop around -----
    nameGroups = vector(mode = 'list', length = length(groupNamesEn))


    names(nameGroups) = c(groupNamesEn)

    for (i in 1:length(groupNamesEn)){
        nameGroups[[i]] = design[,groupNamesEn[i]]
    }
    nameGroups = nameGroups[unlist(lapply(lapply(lapply(nameGroups,unique),ogbox::trimNAs),length)) > 1]
    #debug exclude
    if (!is.null(debug)){
        nameGroups = nameGroups[names(nameGroups) %in% debug]
        groupNamesEn = groupNamesEn[groupNamesEn %in% debug]
    }
    groupNamesEn = names(nameGroups)

    # the main loop around groups ------
    if (!is.null(seed)){
        doRNG::registerDoRNG(seed)
    } else {
        doRNG::registerDoRNG()
    }
    #foreach::foreach (i = 1:length(nameGroups)) %dorng% {
        for (i in 1:length(nameGroups)){
        #debub point for groups
        typeNames = ogbox::trimNAs(unique(nameGroups[[i]]))
        realGroups = vector(mode = 'list', length = length(typeNames))
        names(realGroups) = typeNames
        for (j in 1:length(typeNames)){
            realGroups[[j]] = which(nameGroups[[i]] == typeNames[j])
        }


        if (!is.null(rotate)){
            # this part equalizes representation from individual studies when rotating.
            print('yayay')
            realGroups2 = lapply(realGroups, function(x){
                articles = design[x,PMID]
                minRepresentation = articles %>%
                    table(useNA = 'ifany') %>%
                    min
                lapply (1:length(unique(articles)),function(j){
                    # this turned into a list because if it is not a list, single length vectors behave differently
                    # in sample.
                    if (length( x[articles %in% unique(articles)[j]]) ==1){
                        return(x[articles %in% unique(articles)[j]])
                    }
                    x[articles %in% unique(articles)[j]] %>%
                        sample(size=minRepresentation,replace=FALSE) %>% unlist #%>% #if you decide to remove samples per study comment this part in, delete the part below
                    #sample(.,size = length(.)-round(length(.)*rotate), replace= FALSE)
                }) %>% unlist
            })
            removed = unlist(realGroups)[!unlist(realGroups) %in% unlist(realGroups2)]
            realGroups = realGroups2

            # if rotation is checked, get a subset of the samples. result is rounded. so too low numbers can make it irrelevant
            realGroups2 = lapply(realGroups,function(x){
                if(length(x)==1){
                    warning('Samples with single replicates. Bad brenna! bad!')
                    return(x)
                }
                sort(sample(x,length(x)-round(length(x)*rotate)) %>% unlist)
            })
            removed = c(removed, unlist(realGroups)[!unlist(realGroups) %in% unlist(realGroups2)])
            realGroups = realGroups2
        }
        tempExpr = exprData[unlist(realGroups),]
        tempDesign = design[unlist(realGroups),]




        # replicateMeans ------
        # inefficient if not rotating but if you are not rotating you are only doing it once anyway

        indexes = unique(tempDesign[[replicates]])



        repMeanExpr = sapply(1:length(indexes), function(j){
            tryCatch({
                apply(tempExpr[tempDesign[[replicates]] == indexes[j],], 2,mean)},
                error= function(e){
                    if (is.null(rotate)){
                        print('unless you are rotating its not nice that you have single replicate groups')
                        print('you must be ashamed!')
                        print(j)
                    }
                    tempExpr[tempDesign[[replicates]] == indexes[j],]
                })
        })
        repMeanExpr = t(repMeanExpr)
        repMeanDesign = tempDesign[match(indexes,tempDesign[[replicates]]),]

        # since realGroups is storing the original locations required for
        # silhouette store the new locations to be used with repMeanExpr here
        # use the old typeNames since that cannot change
        realGroupsRepMean =  vector(mode = 'list', length = length(typeNames))
        print(names(nameGroups)[i])
        for (j in 1:length(typeNames)){
            realGroupsRepMean[[j]] = which(repMeanDesign[,groupNamesEn[i]] == typeNames[j])
        }
        names(realGroupsRepMean) = typeNames

        # groupMeans ----
        #take average of every group

        groupAverages = sapply(realGroupsRepMean, function(j){
            groupAverage = apply(repMeanExpr[j,,drop=F], 2, mean)

        })
        groupAverages = t(groupAverages)

        # creation of output directories ----
        # dir.create(paste0(outLoc ,'/Marker/' , names(nameGroups)[i] , '/'), showWarnings = F,recursive = T)
        # dir.create(paste0(outLoc , '/Relax/' , names(nameGroups)[i] , '/'), showWarnings = F, recursive =T)
        dir.create(paste0(outLoc ,'/' , names(nameGroups)[i] , '/'), showWarnings = F, recursive =T)
        if (!is.null(rotate)){
            utils::write.table(removed,
                               file = paste0(outLoc,'/',names(nameGroups)[i] , '/removed'),
                               col.names=F)
        }

        # for loop around groupAverages
        for (j in 1:nrow(groupAverages)){
            # cell type specific debug point
            #if (names(realGroups)[j]=='GabaOxtr'){
            #  print('loyloy')
            #}
            fileName = paste0(outLoc ,'/', names(nameGroups)[i], '/',  names(realGroups)[j])
            # fileName2 = paste0(outLoc , '/Marker/' , names(nameGroups)[i] , '/' , names(realGroups)[j])

            # find markers. larger than 10 fold change to every other group
            #             isMarker = apply(groupAverages,2,function(x){
            #                 all(x[-j] + log(10, base=2) < x[j])
            #             })
            #
            #             fMarker = data.frame(geneData$Gene.Symbol[isMarker], groupAverages[j,isMarker], apply(groupAverages[-j,isMarker,drop=F],2,max), apply(groupAverages[-j,isMarker,drop=F],2,min))
            fChange = foldChange(groupAverages[j, ], groupAverages[-j,,drop=F] ,foldChangeThresh)
            fChangePrint = data.frame(geneNames = geneData[[geneID]][fChange$index], geneFoldChange= fChange$foldChange )
            fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]

            #silhouette. selects group members based on the original data matrix
            # puts them into two clusters to calculate silhouette coefficient
            groupInfo1 = realGroups[[j]]
            groupInfo2 = realGroups[-j]

            # silo = vector(length = nrow(fChangePrint))
            if (!nrow(fChangePrint) == 0){
                silo = 1:nrow(fChangePrint) %>% sapply(function(t){
                    giveSilhouette(which(geneData[[geneID]] == fChangePrint$geneNames[t]),
                                   groupInfo1,
                                   groupInfo2)
                }) %>% t

                fChangePrint = cbind(fChangePrint, silo)
            } else {
                fChangePrint = data.frame(fChangePrint, silo=numeric(0))
            }


            print(fileName)
            # print(nameGroups[[i]])
            utils::write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
            # write.table(fMarker, quote = F, row.names = F, col.names = F, fileName2)

        }# end of for around groupAverages

    } # end of foreach loop around groups
    parallel::stopCluster(cl)
} # end of function



#' @export
pickMarkersAll = function(genesLoc,lilah=F,regex='*',...){
    allGenLocs = list.dirs(genesLoc)
    allGenLocs = allGenLocs[-1]
    allGenLocs = grep(regex,allGenLocs,value=T)
    geneLists = lapply(allGenLocs, pickMarkers,...)
    names(geneLists) = basename(allGenLocs)
    return(geneLists)
}

#' Pick marker genes out of candidate lists
#'
#' Picks marker genes out of candidate lists. Behaves differently based on the number of columns the files have.
#' @export
pickMarkers = function(geneLoc, rotationThresh = 0.95,silhouette = 0.5,foldChange = 10,lilah = F){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    filenames = filenames[!filenames %in% 'removed']
    fileContents = lapply(paste0(geneLoc,'/', filenames), function(x){
        tryCatch(
            utils::read.table(x,stringsAsFactors=FALSE),
            error = function(e){
                NULL
            })
    })
    names(fileContents) = filenames

    # empty first file protection
    lengths = sapply(fileContents,length)
    destinedToBeFirst = which.max(lengths>0)

    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst

    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = names(fileContents)

    if (ncol(fileContents[[1]])==3 & lilah == F){
        # this if for a combination of fold change and silhouette coefficient
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in%
                                                unlist(sapply((1:length(fileContents))[-i], function(x){
                                                    fileContents[[x]][,1]
                                                })),]

            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))>silhouette
                                                        & as.numeric(as.character(tempContent$V2))>log(foldChange,base=2)])
        }
    }else if (ncol(fileContents[[1]])==3 & lilah == T){
        # this if for lilah's selection method
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in%
                                                unlist(sapply((1:length(fileContents))[-i], function(x){
                                                    fileContents[[x]][,1]
                                                })),]
            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))*
                                                            as.numeric(as.character(tempContent$V2))>2])
        }
    } else if (ncol(fileContents[[1]])==1){
        # this is for a mere gene list
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else if(ncol(fileContents[[1]])==2){
        # this is for selection of percentages from confidence output
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>rotationThresh])
        }
    } else if(ncol(fileContents[[1]])==4){
        # this is for 10 fold changed markers. none of the other collumns matter. just get the genes dammit
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else {
        stop('What kind of gibberish is this')
    }

    puristList = vector(mode = 'list', length = length(geneList))

    # if the file only has a single column, this means it is the final list. consider putting this to a new
    # place later since the point of this function was supposed to be filtering the genes as well
    if (ncol(fileContents[[1]])!=1){
        for (i in 1:length(geneList)){
            puristList[[i]] = ogbox::trimElement(geneList[[i]], unlist(geneList[-i]))
        }
    } else {
        puristList = geneList
    }
    names(puristList) = names(geneList)
    puristList = lapply(puristList, as.character)

    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst


    return(puristList)
}

# lapply(c('Eno2','Mog'), findGene, list)
#' Find if a marg
#' @export
findGene = function(gene,list){
    out = lapply(list, function(x){
        ogbox::findInList(gene,x)
    })
    matches = out[lapply(out,length)>0]
    if (length(matches)<1){
        return(NULL)
    }
    matches = sapply(1:length(matches), function(i){
        paste0(names(matches[i]),'_', names(list[[names(matches[i])]][matches[[i]]]))
    })
    return(matches)
}


#' @export
rotateSelect = function(rotationOut,rotSelOut,cores=4, lilah=F, ...){

    # so that I wont fry my laptop
    if(!is.na(parallel::detectCores())){
        if (parallel::detectCores()<cores){
            cores = parallel::detectCores()
            print('max cores exceeded')
            print(paste('set core no to',cores))
        }
    }
    doMC::registerDoMC(cores)

    dirFols = list.dirs(rotationOut, recursive = F)

    loopAround = list.dirs(dirFols[1],full.names = F)

    # in server version full.names input of list.dirs do not work. This fixes it. Might add this to ogbox as an overwrite.
    loopAround = gsub(paste0(dirFols[1],'/'),'',loopAround)
    loopAround = loopAround[!grepl(dirFols[1],loopAround)]

    loopAround = loopAround[!loopAround %in% '']
    #loopAround = loopAround [-which(loopAround %in% c('Relax','Marker',''))]
    #loopAroundRel = loopAround[grepl('Relax',loopAround)]
    #loopAroundMar = loopAround[grepl('Marker',loopAround)]
    dir.create(paste0(rotSelOut), showWarnings = F)

    # for relaxed selection. forces unique selection with matching criteria in puristOut
    # for (i in loopAroundRel){
    foreach::foreach (i  = loopAround) %dopar% {
        print(i)
        dir.create(paste0(rotSelOut,'/',i),recursive = T, showWarnings = F)
        files = list.files(paste0(dirFols[1],'/',i))
        # remove the list of removed samples from the mix
        files = files[!files %in% 'removed']

        pureConfidence = vector(mode = 'list', length =length(files))
        for (j in dirFols){
            #print(paste0(j,'/',i))
            pureSample = pickMarkers(paste0(j,'/',i), lilah,
                                     ...
            )
            pureConfidence = mapply(c,pureSample,pureConfidence,SIMPLIFY=FALSE)
            #print(names(pureConfidence))
            if(length(pureConfidence)>length(files)){
                stop('dafaq man')
            }
        }
        confidence = lapply(pureConfidence,function(x){table(x)/length(dirFols)})

        #         if (any(grepl('removed',names(confidence)))){
        #             confidence = confidence[-which(grepl('removed',names(confidence)))]
        #         }

        for (j in 1:length(confidence)){
            # genes = names(confidence[[j]])[confidence[[j]]>0.95]
            write.table(confidence[[j]],row.names=F,quote=F,col.names=F,
                        file = paste0(rotSelOut,'/',i,'/',
                                      names(confidence)[j]))
        }
        #return(invisible())
    }
}

