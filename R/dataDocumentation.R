#' Brain cell type specific genes
#'
#' A list of brain cell type specific genes selected from mouse pooled cell type
#' microarray data. Published in Mancarci et al. 2016.
'mouseMarkerGenes'

#' Mock dataset for cell type profiles
#'
#' This is a simplified mock dataset that represents the minimal information
#' required for marker gene selection. Gene symbols or other form of gene
#' identifier to be outputted should be present in a column before the
#' expression data. If there are additional columns that are not the expression
#' data from samples they must not be of type double. In the Mancarci et al.
#' study, expression data from samples is microarray data from purified brain cell types.
#' Marker gene selection pipeline is designed to accomodate such microarray data that is on
#' log2 scale. To use any other form of data, it must be converted to log scale and parameters
#' for marker gene selection needs to be adjusted accordingly.
#' @seealso \code{\link{mpg_sampleProfilesMeta}}, \code{\link{mpg_sampleRegionHiearchy}}, input of \code{\link{mpg_markerCandidates}}
'mgp_sampleProfiles'


#' Mock metadata for the mock dataset for cell type profiles
#'
#' This is a simplified mock dataset that represents the minimal information
#' required for marker gene selection. The metadata columns must include:
#' \itemize{
#' \item sample names that match the column names of profile
#' expression matrix ('sampleName' in this example)
#'  \item a column that indicates which samples are replicates of each other
#'  ('replicate' in this example)
#' \item a column that indicates which samples are from the same study
#' ('PMID' in this example)
#'  \item a column that indicates which cell type the
#' sample represents. } Additionally a column can be added to represent the
#' brain region that the sample is extracted from ('region' in this example). If
#' multiple regions apply they can be added to the same column and be separated
#' with commas without spaces. Gene selection functions will consider regions
#' seperately and output results for each of them. If desired one can create a
#' region hierarchy (see \code{\link{mpg_sampleRegionHiearchy}}) in order to
#' indicate regions that are subsets of each other.
#' @seealso \code{\link{mpg_sampleProfiles}},
#'   \code{\link{mpg_sampleRegionHiearchy}}, input of
#'   \code{\link{mpg_markerCandidates}}
'mgp_sampleProfilesMeta'


#' Mock hierarchy of regions
#'
#' This is a simplified mock hiearchy of regions. All represents whole tissue
#' while Region 1 and Region 2 are two subregions. If provided to gene selection
#' functions genes will be selected for All, Region 1 and Region 2 regions
#' seperately
'mpg_sampleRegionHiearchy'


#' Crop of Lesnick et al. dataset of Parkinson's disease patients and controls
#'
#' This dataset is provided to be used as an example for marker gene estimation.
#' @seealso \code{\link{mgp_LesnickCroppedMeta}}
'mgp_LesnickParkinsonsExp'

#' Crop of Lesnick et al. dataset of Parkinson's disease patients and controls
#'
#' This dataset is provided to be used as an example for marker gene estimation. Most of the metadata
#' columns are left behind. Full dataset is available in the analyiss pakcage from Mancarci et al.
#' @seealso \code{\link{mgp_LesnickCroppedExpression}}
'mgp_LesnickParkinsonsMeta'
