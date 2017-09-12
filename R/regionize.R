#' Create cell type annotations per region.
#' @description For each brain region, based on a defined hierarchy, it choses which samples to take in
#' @param design A data.frame containing the design information
#' @param regionNames a string naming the column of the design table that indicates the region names
#' @param groupNames a string vector naming  the columns of the design table that indicates cell type naming schemes.
#' @param regionHierarchy a list of lists that represents the region hierarchy.
#' @return A list containing character vectors that annotates samples. If the sample will not be considered for the region
#' it will be represented as an NA, otherwise it'll have the cell type name indicated in the region
#' @export
regionize = function(design,regionNames,groupNames, regionHierarchy = NULL){
    # if a hierarchy is not provided, create one with a single layer
    if (is.null(regionHierarchy)){
        regionHierarchy = vector(mode='list',
                                 length = design[,regionNames] %>% unique %>% length)
        regionHierarchy = lapply(regionHierarchy,function(x){''})
        names(regionHierarchy)  = design[,regionNames] %>% unique
    }

    # if region to parent and to children is not specificed, take the default as true for all
    if(is.null(design$RegionToParent)[1]){
        design$RegionToParent = TRUE
    }
    if(is.null(design$RegionToChildren)[1]){
        design$RegionToChildren = TRUE
    }

    regionsTree = regionHierarchy %>% unlist %>% names %>% strsplit(split='\\.') %>% unlist %>% unique

    # some non used samples have undefined regions, trim those since you don't use them anyway
    regionsData = design[design[,groupNames, drop=F]%>% apply(1,function(x){all(!is.na(x))}),
                         regionNames] %>% as.character %>% strsplit(',') %>% unlist %>% (ogbox::trimNAs) %>% unique

    # ensure that every region in the dataset is somewhere in the tree
    assertthat::assert_that(length(regionsData[!(regionsData %in%  regionsTree)]) == 0)


    # get the region names with their parents so children can inherit their parents sins
    childrenRegions =  regionHierarchy %>%
        unlist %>%
        names

    regions = childrenRegions %>% strsplit(split='\\.') %>% lapply(function(x){
        sapply(1:length(x), function(i){
            paste(x[1:i],collapse='.')
        })

    }) %>% unlist %>% unique


    regionBased = expand.grid(groupNames, regions)
    names(regionBased) = c('groupName','region')
    regionBased %<>% dplyr::mutate(groupName = as.character(groupName)) %>% dplyr::mutate(region = as.character(region))

    regionGroups = vector(mode = 'list', length = nrow(regionBased))
    names(regionGroups) = paste0(regionBased$region,'_',regionBased$groupName)

    # have translation to short names
    regions = data.frame(region = regions)  %>%
        dplyr::mutate(shortNames = region %>% as.character %>% strsplit(split='\\.') %>% lapply(function(x){x[length(x)]}) %>% unlist)

    # get a full name for all listed regions per sample (samples can have multiple regions if they are specifically
    # extracted from a lower level region but we want to use them for others)
    regionList = design[,regionNames] %>% as.character %>% strsplit(',') %>% lapply(function(x){
        as.character(regions$region[match(x,regions$shortNames)])
    })

    # i is the region groups to assign
    # j rotates over samples and decides whether or not to assign them to the region
    for (i in 1:nrow(regionBased)){
        regionGroups[[i]] = design[,as.character(regionBased$groupName[i])]

        # decide which samples should belong to the region
        inRegion = sapply(1:length(regionList), function(j){
            # get the ones that are directly listed in the region
            inReg = regionBased$region[i] %in%  regionList[[j]]

            # get ones that are in children nodes this is done because a larger brain region will include cell types
            # taken from more specific brain regions
            if (design[j,'RegionToParent']){
                inReg = inReg | (
                    regionList[[j]] %>%
                        grepl(paste0('^',regionBased$region[i]),  .) %>% any)
            }

            # get ones that are in parent nodes. this is done because a smaller brain region will include cell types that
            # are isolated from the larger region. do not add it if we have the specific cell type from the region
            if(design[j,'RegionToChildren']){
                # take a quick peek to make sure the cell type doesn't exist in a region closer to the core region
                otherSamples = regionList[design[,regionBased$groupName[i]] == design[j,regionBased$groupName[i]]] %>% .[sapply(.,function(x){!is.null(x)})]
                otherSamples %<>% unlist
                override  = regionList[[j]] %>% sapply(function(x){
                    gsub(paste0('\\Q',x,'\\E'),'', otherSamples[grepl(paste0('\\Q',x,'\\E'), otherSamples)]) %>%
                        nchar %>% {.>0} %>% any
                }) %>% all
                if (!override){
                    toChild = regionBased$region[i] %>%
                        grepl(paste0('^',ogbox::regexMerge( regionList[[j]], exact=TRUE)), .)
                    inReg = inReg | toChild
                }



            }
            return(inReg)

        })

        regionGroups[[i]][!inRegion] =  NA
    }
    names(regionGroups) = gsub(pattern='^.*\\.','',names(regionGroups))
    return(regionGroups)
}

#' Memoised version of regionize function
#'
#' This function is a \code{\link[memoise]{memoise}}d version of the \code{\link{regionize}} function
#'@export
memoReg = memoise::memoise(regionize)
