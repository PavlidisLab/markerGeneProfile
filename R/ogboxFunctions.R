# allows directly assigning list outputs to variables
# http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
#' @name list
#'
#' @title Get individual components of a list output
#'
#' @description Separates a list into components. Usually used to split outputs
#' from functions that returns list. Taken from stackoverflow question http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
#'
#' @examples
#' list[a,b] = list(a=c(1,2,3), b= c(2,3,4))
#' list[c,] = list(c=c(1,2,3), d= c(2,3,4))
#' @export
list <- structure(NA,class="result")


#' @title Get individual components of a list output
#'
#' @description See \code{link{list}}
#' @export
"[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
        a <- args[[i]]
        if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
}



#' Seprate expression data from gene annotations
#'
#' Separetes numberic component of gene expression data from non-numeric gene
#' annotations. Detects the first non-double column in the data frame and
#' separates the data into two parts
#'
#' @param data Expression data with gene annotations. First columns must inlcude
#' the gene annotations and they must not have data of type "double".
#'
#' @return  A list of length 2. 1st element includes the gene annoation data, the
#' second element includes the gene expression data.
#'
#' @export
sepExpr = function(data){
    if (class(data)[1] %in% c('data.frame','tbl_df')){
        # unlist is there to enable tbl_df to work
        for (i in 1:ncol(data)){
            if ('double'==typeof(data[,i] %>% unlist)){
                expBound = i
                break
            }
        }
        geneData = data[,1:(expBound-1),drop=F]
        exprData = data[,expBound:ncol(data),drop=F]
        return(list(geneData,exprData))
    } else if (class(data)[1] =='data.table'){
        for (i in 1:ncol(data)){
            if ('double'==typeof(data[[i]])){
                expBound = i
                break
            }
        }
        geneData = data[,1:(expBound-1),with=F]
        exprData = data[,expBound:ncol(data), with = F]
        return(list(geneData,exprData))
    }
}


trimNAs = function(aVector) {
    return(aVector[!is.na(aVector)])
}


# removes given elements in a vector by shortening it
trimElement = function (aVector,e){
    return(aVector[!(aVector %in% e)])
}

# seeks for a given object in a single layered list
findInList = function(object, aList){
    indices = vector()
    for (i in 1:length(aList)){
        if (object %in% aList[[i]]){
            indices = c(indices, i)
        }
    }
    return(indices)
}


# merges regex's with an or clause. search for multiple regexes
regexMerge = function(regexList, exact = FALSE){
    assertthat::assert_that(is.logical(exact))
    exact = as.character(exact)
    out = switch (exact,
                  'FALSE' = {paste0('(',paste0(regexList,collapse=')|('),')')},
                  'TRUE' =  {paste0('(\\Q',paste0(regexList,collapse='\\E)|(\\Q'),'\\E)')}
    )
    out = paste0('(',out,')')
}



#' Scales vector into a 0-1 interval
#' @export
scale01 = function(x){
    scaleToInt(x,1,0)
}

#' Scales vector to arbitrary interval
#' @export
scaleToInt = function(x, max,min){
    scaleFun = scaleIntervals(max(x,na.rm = TRUE),min(x, na.rm=TRUE),max,min)
    scaleFun(x)
}

scaleIntervals = function(max,min,maxOut,minOut){
    a = (max - min)/(maxOut - minOut)
    b = max - maxOut*a
    if(a != 0){
        return(teval(paste0("function(x){(x - ",b,")/",a,'}')))
    }else{
        mean = (maxOut - minOut)/2
        return(teval(paste0("function(x){out = rep(",mean,",length(x));names(out)=names(x);return(out)}")))
    }
}



#direct text eval
teval = function(daString,...){
    eval(parse(text=daString),...)
}

# http://stackoverflow.com/questions/18122548/display-names-of-column-of-recursive-list-as-tree
# displays a list as a tree by their names
#' Displays a list as a tree
#' @export
nametree <- function(X, prefix1 = "", prefix2 = "", prefix3 = "", prefix4 = "")
    if( is.list(X) )
        for( i in seq_along(X) ) {
            cat( if(i<length(X)) prefix1 else prefix3, names(X)[i], "\n", sep="" )
            prefix <- if( i<length(X) ) prefix2 else prefix4
            nametree(
                X[[i]],
                paste0(prefix, "├──"),
                paste0(prefix, "│  "),
                paste0(prefix, "└──"),
                paste0(prefix, "   ")
            )
        }


#' A violin plot overlayed by a boxplot
#'
#' This function returns a list of ggplot elements that overlays a boxplot over
#' a violin plot.
#'
#' @param data Default dataset to use for plot. If not already a data.frame,
#' will be converted to one by \link[ggplot2]{fortify}. If not specified, must
#' be suppled in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot. If not
#' specified, must be suppled in each layer added to the plot.
#'
#'
#' @export
geom_boxvio = function(data=NULL, mapping = NULL){
    list(ggplot2::geom_violin(color="#C4C4C4", fill="#C4C4C4",data=data, mapping = mapping),
         ggplot2::geom_boxplot(width=0.1,fill = 'lightblue',data=data, mapping = mapping),
         cowplot::theme_cowplot(),
         ggplot2::theme(axis.text.x  = ggplot2::element_text(size=25),
                        axis.title.y = ggplot2::element_text(vjust=0.5, size=25),
                        axis.title.x = ggplot2::element_text(vjust=0.5, size=0) ,
                        title = ggplot2::element_text(vjust=0.5, size=25),
                        axis.text.y = ggplot2::element_text(size = 13)))
}


#' Turn every member of Vector to a color from the palette
#' @param vector a vector
#' @param dictionary a named vector
#' @param NAcolor that to put instead of NAs
#' @export
toColor = function(vector, palette = NULL,NAcolor = 'white'){
    if(class(vector) == "factor"){
        vector = as.character(vector)
    }
    if(!is.null(palette) & !is.null(names(palette))){
        assertthat::assert_that(all(vector %in%names(palette)))
    }
    if(is.null(names(palette)) & !is.null(palette)){
        names(palette) = unique(vector)
    }
    out = replaceElement(vector,dictionary = palette,NAreplace = NAcolor)
    names(out) = c('cols','palette')
    return(out)
}

#' Replaces elements of a vector based on the dictonary provided
#' @param vector a vector
#' @param dictionary a named vector or just a vector if labels are provided
#' @param labels a character vector
#' @param NAreplace that to put instead of NAs
#' @return A list with $newVector and $dictionary.
replaceElement = function(vector, dictionary = NULL,labels = NULL, NAreplace = NA){
    if(class(vector) == "factor"){
        vector = as.character(vector)
    }
    #vector = as.factor(vector)
    uniq = unique(vector) %>% trimNAs
    if (is.null(dictionary[1])){
        dictionary = rainbow(length(uniq))
    }

    if(is.null(labels) & !is.null(names(dictionary))){
        labels = names(dictionary)
    } else if(is.null(labels) & is.null(names(dictionary))){
        labels = uniq
    }

    cols = vector(length = length(vector))
    #to match palette names to uniq names so that custom naming is possible
    # dictionary = dictionary[trimNAs(match(uniq,labels))]
    #names(dictionary) = uniq


    for (i in 1:length(dictionary)){
        vector[vector == labels[i]]= dictionary[i]
    }

    vector[is.na(vector)] = NAreplace

    if(length(dictionary)>0){
        names(dictionary) = labels
    }
    out = list(newVector = vector , dictionary = dictionary)
    return(out)
}
