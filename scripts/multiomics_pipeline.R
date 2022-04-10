#!/usr/bin/Rscript
library(argparser, quietly=TRUE)
library(igraph)
library(mixOmics)

#' Load and parse data
#'
#' Load in omics data into a diablo-compatible format: infile_path -> dataframe
#' @param infile_path path to tab separated input file.
#' @param offset add offset to data. Defaults to 0 (no offset).
#' @param missing_as set missing values to this. Defaults to NA.
#' @param rmna remove NA values if a feature has only NA. Defaults to TRUE.
#' @export
# @examples
# parse_data("infile_path")
parse_data <- function(infile_path, offset=0, missing_as=NA, rmna=TRUE) {
  # load in omics data into a diablo-compatible format: infile_path -> dataframe
  print("Parsing file:")
  print(infile_path)
  data <- read.table(infile_path, sep="\t", header=TRUE, row.names=1) + offset
  if (is.na(missing_as)) {
    data[which(data == 0, arr.ind=TRUE)] <- NA
  } else {
    data[which(data == 0, arr.ind=TRUE)] = missing_as
  }
  if (rmna == TRUE) {
    print("Dropping features where all values are NA")
    data <- data[, unlist(lapply(data, function(x) !all(is.na(x))))]
  }
  return(data)
}

#' Load and parse mappings
#'
#' Load in id:name data into a mapping table, format: infile_path -> dataframe
#' Feature names can only reach a certain length and must be unique.
#' So you can use this if you want to remap identifiers out to something else.
#' @param infile_path path to tab separated input file.
#' @seealso [multiomics::remap_data()]
#' @export
# @examples
# parse_mappings("infile_path")
parse_mappings <- function(infile_path) {
  # load in id:name data into a mapping table, format: infile_path -> dataframe
  print("Parsing mappings...")
  print(infile_path)
  data <- read.table(infile_path, sep="\t", header=TRUE, row.names=2)
  data["X"] <- NULL
  return(data)
}

#' Map feature names to feature id
#'
#' Remap feature names to feature id, format: dataframe, dataframe -> dataframe
#' Feature names can only reach a certain length and must be unique.
#' So you can use this if you want to remap identifiers out to something else.
#' @param data path to dataframe with data
#' @param mapping path to dataframe with mapping
#' @seealso [multiomics::parse_mappings()]
#' @export
# @examples
# remap_data(dataframe_with_data, dataframe_with_mapping)
remap_data <- function(data, mapping) {
  # map feature names to feature id, format: dataframe, dataframe -> dataframe
  # TODO: doesnt handle duplicates in row ids (some feature id map to same name)
  print("Remapping data with map files (duplicate row names will be renamed!)")
  data <- merge(t(data), mapping, by=0, all.x=TRUE)
  row.names(data) <- make.names(data$val, unique=TRUE)
  data["Row.names"] <- NULL
  data["val"] <- NULL
  return(t(data))
}

#' Show proportion of NA values in data
#'
#' Show proportion of NA values in data: dataframe -> plot
#' If your NA values are represented as 0 or something else, convert them first.
#' @param data dataframe with na values.
#' @param name name the plot. Defaults to "".
#' @export
# @examples
# show_na_prop(dataframe_with_na, "some data")
show_na_prop <- function(data, name="") {
  # first, sum the number of missing values per variable
  sum_na_per_var <- apply(data, 2, function(x) {sum(is.na(x))})
  prop_na_per_var <- sum_na_per_var/nrow(data)
  na_prop_thresholds <- seq(0,1,0.05)
  sum_na_all_vars <- sapply(na_prop_thresholds, function(x) sum(prop_na_per_var >= x))
  prop_na_all_vars <- sum_na_all_vars/ncol(data)

  # show proportion of features with more than x proportion of missing values
  # for different values of x
  plot(x = na_prop_thresholds , y = prop_na_all_vars, xlab = "Proportion of missing values",
       ylab = "Proportion of features", main= paste(name, 'missing value proportions for unfiltered data'))
}

#' Remove NA values in data by proportion
#'
#' Remove NA values from data and count: dataframe with na -> dataframe clean
#' If your NA values are represented as 0 or something else, convert them first.
#' @param data list of dataframes with na values.
#' @param na_prop drop column if proportion of NA >= na_prop. Defaults to 0.3.
#' @param zero_as_na treat zero values as NA. Defaults to TRUE.
#' @seealso [multiomics::remove_na_class()], [multiomics::remove_na_prop()], [multiomics::remove_novar()]
#' @export
# @examples
# remove_na_prop(dataframe_with_data_na, na_prop=0.3, zero_as_na=TRUE)
remove_na_prop <- function(data, na_prop=0.3, zero_as_na=TRUE) {
  # mapply(function(x) remove_na_prop_(x, na_prop), data)
  data_new <- list()
  for(i in 1:length(data)){
    data_tmp <- remove_na_prop_(data[[i]], na_prop)
    data_new <- append(data_new, list(data_tmp))
  }
  names(data_new) <- names(data)
  return(data_new)
}

remove_na_prop_ <- function(data, na_prop=0.3, zero_as_na=TRUE) {
  # we want to compare dimensions later to see how many features were dropped
  data_na <- data
  if (zero_as_na == TRUE) {data_na[data_na == 0] <- NA}
  # first, sum the number of missing values per variable
  sum_na_per_var <- apply(data_na, 2, function(x){sum(is.na(x))})
  # these variables could be removed
  remove_var <- which(sum_na_per_var/ncol(data_na) >= na_prop)
  if (length(remove_var) > 0) {
    print("Removing this many columns from data:")
    print(length(remove_var))
    data_na = data_na[, -c(remove_var)]
  } else {
    print("No columns removed from data")
  }
  print("Data dimensions (original):")
  print(dim(data))
  print("Data dimensions (filtered):")
  print(dim(data_na))
  print("Proportion of missing values (original):")
  print(sum(is.na(data)) / (nrow(data) * ncol(data)))
  print("Proportion of missing values (filtered):")
  print(sum(is.na(data_na)) / (nrow(data_na) * ncol(data_na)))
  return(data_na)
}

#' Remove NA values in data by class
#'
#' Remove NA values from data and count: dataframe with na -> dataframe clean
#' If your NA values are represented as 0 or something else, convert them first.
#' @param data dataframe with na values.
#' @param classes 1-column dataframe with labelled samples matching order.
#' @param missing_as this value to represent missing values. Defaults to NA.
#' @seealso [multiomics::remove_na_class()], [multiomics::remove_na_prop()], [multiomics::remove_novar()]
#' @export
# @examples
# remove_na_class(dataframe_with_data_na, classes, missing_as=NA)
remove_na_class <- function(data, classes, missing_as=NA) {
  # where there is NA for a whole class of features, remove feature column
  if (is.na(missing_as)) {
    data[is.na(data)] <- 0
  } else {
    data[which(data == missing_as, arr.ind=TRUE)] <- 0
  }
  print(dim(data))
  print("Dropping features where at least one class is NA")
  uniq <- unique(classes)
  subsets <- list()
  for (i in uniq) {subsets[[i]] <- data[grep(i, rownames(data)), ]}
  subsets <- lapply(lapply(lapply(subsets, colSums), data.frame), t)
  subsets <- do.call("rbind", subsets)
  rownames(subsets) = uniq

  if (is.na(missing_as)) {
    subsets[which(subsets == 0, arr.ind=TRUE)] <- NA
  } else {
    subsets[which(subsets == 0, arr.ind=TRUE)] <- missing_as
  }

  subsets = t(na.omit(t(subsets)))
  data = data[, c(colnames(subsets))]

  if (is.na(missing_as)) {
    data[which(data == 0, arr.ind=TRUE)] <- NA
  } else {
    data[which(data == 0, arr.ind=TRUE)] <- missing_as
  }
  return(data)
}

#' Remove 0 variance columns in data
#'
#' Remove 0 variance columns in data: dataframe with na -> dataframe clean
#' @param data dataframe with na values.
#' @seealso [multiomics::remove_na_class()], [multiomics::remove_na_prop()], [multiomics::remove_novar()]
#' @export
# @examples
# remove_novar(dataframe_with_data_na)
remove_novar <- function(data) {
  # samples with zero variance are meaningless for PCA: dataframe -> dataframe
  # print("Dimensions before removing invariant columns:")
  # print(dim(data))
  data <- data[, which(apply(data, 2, var, na.rm=TRUE) != 0)]
  # print("Dimensions after removing invariant columns:")
  # print(dim(data))
  return(data)
}

#' Show distribution of feature variance in data
#'
#' Show distribution of feature variance in data: dataframe -> plot of variances
#' @param data dataframe with na values.
#' @param do_log2 perform log2 transform for visualisation. Defaults to TRUE.
#' @seealso [multiomics::remove_na_class()], [multiomics::remove_na_prop()], [multiomics::remove_novar()]
#' @export
# @examples
# show_variance(dataframe_with_data)
show_variance <- function(data, do_log2=TRUE) {
  variances <- lapply(
    lapply(data, data.frame),
    function(y) sapply(y, function(x) var(x, na.rm=TRUE))
  )
  variances <- lapply(variances, unname)
  if (do_log2) {variances <- lapply(variances, log2)}
  names(variances) <- names(data)
  warning(
    "If variance is too low, downstream steps may not work! You can inspect \
    the distribution of variance per feature with the box plots, and filter \
    or drop variables as needed."
  )
  boxplot(
    variances, xlab="Omics", ylab="Feature variance (log2)",
    main="Feature variances (log2) per omics"
  )
}

#' Load in class data
#'
#' Load in class data for data: infile_path -> vector (of strings)
#' @param infile_path path to tab separated input classes file with header
#' @export
# @examples
# parse_classes(infile_path)
parse_classes <- function(infile_path) {
  # load in class data for diablo: infile_path -> vector (of strings)
  data <- read.table(infile_path, sep="\t", header=TRUE, row.names=1)
  return(unlist(as.vector(t(data))))
}

#' Create the design matrix
#'
#' Create the design matrix: list of named dataframes -> matrix
#' If your NA values are represented as 0 or something else, convert them first.
#' @param data list of named dataframes
#' @param link linkage, (discrimination) 0 <-> 1 (correlation). Defaults to 0.1.
#' @export
# @examples
# create_design(infile_path, link=0.1)
create_design <- function(data, link=0.1) {
  # create design matrix from data: dataframe, link -> matrix (design)
  design <- matrix(link, ncol = length(data), nrow = length(data),
                   dimnames = list(names(data), names(data)))
  diag(design) <- 0
  return(design)
}

#' Turn zero into NA
#'
#' Turn zero into NA: dataframe with 0 -> dataframe with NA
#' @param data dataframe
#' @export
# @examples
# zero_to_na(data)
zero_to_na <- function(data) {
  data[which(data == 0, arr.ind=TRUE)] <- NA
  return(data)
}

#' Count missing values
#'
#' Count missing values: list of named dataframes -> list of dataframes
#' @param data dataframe
#' @export
# @examples
# count_missing(data)
count_missing <- function(data) {
  # count missing values in data: dataframe -> list(dataframe, double)
  # return 0 values to NA
  ids_na <- is.na(data)
  pct_na <- sum(is.na(data)) / (nrow(data) * ncol(data))
  print("Proportion of missing values in data:")
  print(pct_na)
  return(list(ids_na=ids_na, pct_na=pct_na))
}

impute_missing_ <- function(data, ncomp=10, block_name="", outdir="./") {
  # impute missing values with nipals: dataframe (with NA) -> dataframe
  print("Number of components for imputation:")
  print(ncomp)
  # data <- data[!is.na(apply(data, 1, function(x) var(x, na.rm = TRUE))),
  #              !is.na(apply(data, 2, function(x) var(x, na.rm = TRUE)))]
  nipals_tune <- nipals(data, ncomp=ncomp)
  barplot(nipals_tune$eig, main=paste(block_name, "Screeplot (nipals imputed)"),
    xlab="Number of components", ylab="Explained variance"
  )
  nipals_impute <- impute.nipals(data, ncomp = ncomp)
  outfile_path <- paste(outdir, "/", "data_", block_name, "_imputed.tsv", sep="")
  write.table(
    as.data.frame(nipals_impute), file=outfile_path, quote=FALSE, sep="\t"
  )
  return(nipals_impute)
}

#' Impute missing values
#'
#' Impute missing values: list of named dataframes with NA -> list of dataframes
#' @param data dataframe
#' @param ncomps number of components
#' @param outdir write imputed data to outfile path
#' @seealso [multiomics::impute_missing()], [multiomics::replace_missing()], [multiomics::cor_imputed_unimputed()]
#' @export
# @examples
# impute_missing(data, 10, "./")
impute_missing <- function(data, ncomps, outdir) {
  # mapply(function(x, y, z)
  #   impute_missing_(x, y, z, outdir), data, ncomps, names(data), SIMPLIFY=FALSE
  # )
  data_new <- list()
  for(i in 1:length(data)){
    data_tmp <- impute_missing_(data[[i]], ncomps[[i]], names(data)[[i]])
    data_new <- append(data_new, list(data_tmp))
  }
  names(data_new) <- names(data)
  return(data_new)
}

replace_missing_ <- function(data, imputed) {
  mask <- is.na(data)
  imputed[!mask] <- data[!mask]
  return(imputed)
}

#' Replace missing values
#'
#' Replace missing values: list of named dataframes with NA -> list of dataframes
#' This replaces missing values in the original dataframe with the corresponding imputed value
#' @param data dataframe with original data
#' @param imputed dataframe with imputed data
#' @seealso [multiomics::impute_missing()], [multiomics::replace_missing()], [multiomics::cor_imputed_unimputed()]
#' @export
# @examples
# replace_missing(data, imputed)
replace_missing <- function(data, imputed) {
  # mapply(function(x, y) replace_missing_(x, y), data, imputed)
  data_new <- list()
  for(i in 1:length(data)){
    data_tmp <- replace_missing_(data[[i]], imputed[[i]])
    data_new <- append(data_new, list(data_tmp))
  }
  names(data_new) <- names(data)
  return(data_new)
}

#' Plots PCA without performing multilevel decomposition
#'
#' Plots PCA for non-repeated measurements
#' @param data dataframe with original data
#' @param classes 1-column dataframe with classes
#' @param pch 1-column dataframe with longitudinal measurement info. Defaults to NA.
#' @param title plot title. Defaults to ""
#' @param ncomp number of components to plot, if 0 estimate automatically. Defaults to 0.
#' @param show show the plot in the R graphics output, does not save. Defaults to FALSE.
#' @seealso [multiomics::plot_pca_single()], [multiomics::plot_pca_multilevel()]
#' @export
# @examples
# plot_pca_single(data, classes, pch=NA, title="", ncomp=0, show=FALSE)
plot_pca_single <- function(data, classes, pch=NA, title="", ncomp=0, show=FALSE) {
  # do pca on individual classes: dataframe, vector, vector -> outfile_path.pdf
  names <- names(data)
  print("Removing 0 variance columns from data...")
  data <- lapply(data, remove_novar)
  if (ncomp == 0) {ncomp <- dim(classes)[1]}

  print("Performing PCA on data...")
  data_pca <- lapply(data, pca, ncomp=ncomp, center=TRUE, scale=TRUE)
  if (show == TRUE) {
    for(i in 1:length(data_pca)) {
      print(sprintf("Summary of PCA on %s...", names(data_pca)[i]))
      print(data_pca[[i]])
    }
  }

  print("Plotting explained variance of PCA components...")
  for(i in 1:length(data_pca)) {
    plot(data_pca[[i]], main=paste(names[[i]], "screeplot"))
  }
  if (length(pch) > 1) {
    print("Plotting PCA by groups...")
    # mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=FALSE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 1/2"), pch=pch), data_pca, names)
    # mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=FALSE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 1/3"), pch=pch), data_pca, names)
    # mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=FALSE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 2/3"), pch=pch), data_pca, names)
    for(i in 1:length(data_pca)){
      plotIndiv(
        data_pca[[i]], comp=c(1,2), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (1,2)"), pch=pch, legend.title = "Class", legend.title.pch = "Secondary Class",
      )
      plotIndiv(
        data_pca[[i]], comp=c(2,3), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (2,3)"), pch=pch, legend.title = "Class", legend.title.pch = "Secondary Class",
      )
      plotIndiv(
        data_pca[[i]], comp=c(1,3), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (1,3)"), pch=pch, legend.title = "Class", legend.title.pch = "Secondary Class",
      )
    }

  } else {
    print("Plotting PCA by groups...")
    # mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 1/2")), data_pca, names)
    # mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 1/3"), pch=pch), data_pca, names)
    # mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
    #   group=classes, legend=TRUE, ncomp=ncomp,
    #   title=paste(title, y, "PCA 2/3"), pch=pch), data_pca, names)
    for(i in 1:length(data_pca)) {
      plotIndiv(
        data_pca[[i]], comp=c(1,2), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (1,2)"), legend.title = "Class",
      )
      plotIndiv(
        data_pca[[i]], comp=c(2,3), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (2,3)"), legend.title = "Class",
      )
      plotIndiv(
        data_pca[[i]], comp=c(1,3), ind.names=FALSE, group=classes, legend=TRUE,
        ncomp=ncomp, title=paste(title, names[[i]], "PCs (1,3)"), legend.title = "Class",
      )
    }
  }
  return(data_pca)
}

#' Plot correlation circle plots and biplots
#'
#' Plot correlation circle plots and biplots:
#'   list of dataframes, list of pc, classes dataframe, names vector -> plots
#' @param data list of dataframes with original data
#' @param data_pca list of dataframes with principal component data
#' @param classes 1-column dataframe with classes
#' @param names vector of names for each dataframe in list
#' @export
# @examples
# plot_additional(data, data_pca, classes, names)
plot_additional <- function(data, data_pca, classes, names) {
  # correlation circle and biplots: dataframe, list, vector -> outfile_path.pdf
  print("Plotting correlation circle plots...")
  # mapply(function(x, y) plotVar(x, comp=c(1, 2),
  #   title=paste(y, "PCA 1/2"), var.names=FALSE), data_pca, names)
  for(i in 1:length(data_pca)){
    plotVar(
      data_pca[[i]], comp=c(1, 2), title=paste(names[[i]], "PCA 1/2"),
      var.names=FALSE
    )
  }

  print("Plotting biplots...")
  # mapply(function(x, y, z) biplot(y, cex=0.7, xlabs=paste(classes, 1:nrow(x)),
  #   main=paste(z, "Biplot")), data, data_pca, names)
  for(i in 1:length(data)){
    biplot(
      data_pca[[i]], cex=0.7, xlabs=paste(classes, 1:nrow(data[[i]])),
      main=paste(names[[i]], "PCA 1/2")
    )
  }
}

#' Plots PCA with multilevel decomposition
#'
#' Plots PCA with multilevel decomposition to account for repeated measurements
#' @param data list of dataframes with original data
#' @param classes 1-column dataframe with classes
#' @param pch 1-column dataframe with longitudinal measurement info. Defaults to NA.
#' @param title plot title. Defaults to ""
#' @param ncomp number of components to plot, if 0 estimate automatically. Defaults to 0.
#' @param show show the plot in the R graphics output, does not save. Defaults to FALSE.
#' @seealso [multiomics::plot_pca_single()], [multiomics::plot_pca_multilevel()]
#' @export
# @examples
# plot_pca_multilevel(data, classes, pch=NA, title="", ncomp=0, show=FALSE)
plot_pca_multilevel <- function(data, classes, pch, title="", ncomp=0, show=FALSE) {
  names <- names(data)

  print("Removing 0 variance columns from data...")
  data <- lapply(data, remove_novar)

  if (ncomp == 0) {ncomp = dim(classes)[1]}

  print("Performing PCA on data...")
  data_pca <- lapply(data, pca, ncomp=ncomp, center=TRUE, scale=TRUE,
    multilevel=pch)
  if (show == TRUE) {
    print("Showing PCA multilevel component contribution...")
    print(data_pca)
  }

  print("Plotting PCA multilevel component contribution...")
  # mapply(function(x, y) plot(x, main=paste(y, "Screeplot multilevel")),
  #   data_pca, names)
  for(i in 1:length(data_pca)){
    plot(data_pca[[i]], main=paste(names[[i]], "Screeplot multilevel"))
  }

  print("Plotting PCA multilevel...")
  # mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=FALSE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PCA M 1/2"), pch=pch), data_pca, names)
  # mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=FALSE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PCA M 1/3"), pch=pch), data_pca, names)
  # mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=FALSE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PCA M 2/3"), pch=pch), data_pca, names)
  for(i in 1:length(data_pca)){
    plotIndiv(
      data_pca[[i]], comp=c(1,2), ind.names=FALSE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PCA M 1/2"), pch=pch
    )
    plotIndiv(
      data_pca[[i]], comp=c(1,3), ind.names=FALSE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PCA M 1/3"), pch=pch
    )
    plotIndiv(
      data_pca[[i]], comp=c(2,3), ind.names=FALSE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PCA M 2/3"), pch=pch
    )
  }

  return(data_pca)
}

#' Perform PLSDA
#'
#' Perform PLSDA (can account for repeated measurments if specified):
#' @param data list of dataframes with original data
#' @param classes 1-column dataframe with classes
#' @param pch 1-column dataframe with repeated measurments. Defaults to NA.
#' @param title string to name plots with. Defaults to "".
#' @param ncomp integer assigning number of components for PLSDA. Defaults to 0.
#' @param contrib string to determine loading weights "max", "min". Defaults to "max".
#' @param outdir string to outfile directory. Defaults to "./"
#' @param mappings dataframe of mappings for features. Defaults to NULL
#' @param dist string describing distance metric "centroids.dist", "max.dist", "mahalanobis.dist". Defaults to "centroids.dist"
#' @param bg boolean describing if background should be plotted. Defaults to TRUE
#' @param validation character specifying "loo" or "M-fold" cross-validation. Defaults to "loo".
#' @param folds if M-fold validation, number of folds. Defaults to 10. No effect for "loo"
#' @param nrepeat if M-fold validation, number of repeats. Defaults to 10. No effect for "loo"
#' @param near_zero_var if variance approaches zero, enable this option
#' @seealso [multiomics::classify_plsda()], [multiomics::classify_splsda()]
#' @export
# @examples
# classify_plsda(data, classes, pch=NA, title="", ncomp=0, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE)
classify_plsda <- function(data, classes, pch=NA, title="", ncomp=0,
  contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE,
  validation="loo", nrepeat=10, folds=10, near_zero_var=FALSE) {
  # mapply(function(x, y) classify_plsda_(
  #   x, classes, pch, y, ncomp, contrib, outdir, mappings, dist, bg), data, title,
  #   SIMPLIFY=FALSE)
  data_new <- list()
  for(i in 1:length(data)){
    data_tmp <- classify_plsda_(
      data[[i]], classes, pch, title[[i]], ncomp, contrib, outdir, mappings,
      dist, bg, validation=validation, nrepeat=nrepeat, folds=folds,
      near_zero_var
    )
    data_new <- append(data_new, list(data_tmp))
  }
  names(data_new) <- names(data)
  return(data_new)
}

classify_plsda_ <- function(data, classes, pch=NA, title="", ncomp=0,
  contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE,
  validation="loo", nrepeat=10, folds=10, near_zero_var=FALSE) {
  # discriminate samples: list, vector, bool, integer -> list
  # single or multilevel PLS-DA
  if (length(pch) > 1) {
    print("Plotting multi level partial least squares discriminant analysis")
    pch <- c(as.factor(pch))
    title_plt <- paste(title, "PLSDA multi")
    data_plsda <- plsda(
      data, Y=classes, multilevel=c(as.factor(pch)), ncomp=ncomp,
      near.zero.var=near_zero_var
    )

    if (ncomp > 1) {
      if (!is.na(bg)) {
        bg = background.predict(data_plsda,comp.predicted=2,dist=dist)
        plotIndiv(data_plsda, ind.names=FALSE, group=classes,
          legend=TRUE, pch=pch, title=paste(title_plt, "1/2"), comp=c(1,2),
          ellipse=TRUE, background=bg
        )
      }
    }

    if (ncomp > 1) {
      plotIndiv(data_plsda, ind.names=FALSE, group=classes,
        legend=TRUE, pch=pch, title=paste(title_plt, "1/2"), comp=c(1,2),
        ellipse=TRUE,
      )
    }

    if (ncomp > 2) {
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        pch=pch, title=paste(title_plt, "1/3"), comp=c(1,3), ellipse=TRUE
      )
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        pch=pch, title=paste(title_plt, "2/3"), comp=c(2,3), ellipse=TRUE
      )
    }

  } else {
    print("Plotting multi level partial least squares discriminant analysis")
    title_plt <- paste(title, "PLSDA single")
    data_plsda <- plsda(data, Y=classes, ncomp=ncomp)

    if (ncomp > 1) {
      if (!is.na(bg)) {
        bg <- background.predict(data_plsda,comp.predicted=2,dist=dist)
        plotIndiv(data_plsda, ind.names=FALSE, group=classes,
          legend=TRUE, title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
          background=bg
        )
      }
    }

    if (ncomp > 1) {
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
      )
    }
    if (ncomp > 2) {
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        title=paste(title_plt, "1/3"), comp=c(1,3), ellipse=TRUE
      )
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        title=paste(title_plt, "2/3"), comp=c(2,3), ellipse=TRUE
      )
    }
  }

  print("Getting performance metrics")
  print("Plotting error rates...")
  metrics <- mixOmics::perf(
    data_plsda, progressBar=TRUE, auc=TRUE,
    validation=validation, nrepeat=nrepeat, folds=folds
  )
  print(metrics$error.rate)
  plot(metrics, main="Error rate PLSDA", col=color.mixo(5:7), sd=TRUE)

  # supporess roc printing automatically to stdout
  sink("/dev/null")
  # roc <- mapply(function(x) auroc(data_plsda, roc.comp=x), seq(ncomp))
  roc <- list()
  for(i in 1:ncomp){
    data_tmp <- auroc(data_plsda, roc.comp=i)
    roc <- append(roc, list(data_tmp))
  }
  sink()

  # consider using plotvar for some datasets
  if (ncomp > 1) {
    print("Plotting arrow plot...")
    plotArrow(data_plsda, ind.names=FALSE, legend=TRUE, title="PLSDA")
  }

  print("Getting loadings and plotting clustered image maps")

  # setup colour map for clustered image plots
  colours_class <- color.mixo(
    1:length(unique(classes))
  )[as.numeric(as.factor(classes))]

  if (length(pch) > 1) {
    colours_pch <- color.mixo(
      1:length(unique(pch))
    )[as.numeric(as.factor(pch))]
    colours_cim <- cbind(colours_class, colours_pch)
  } else {colours_cim <- data.frame(colours_class)}

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data_plsda$names$colnames$X, nchar))) > 8) {
    show_cols = FALSE
  } else {
    show_cols = TRUE
  }

  cim(data_plsda, title="PLSDA", row.sideColors=colours_cim,
    legend=list(title="Status"), col.names=show_cols
  )
  for (comp in seq(ncomp)) {
    cim(data_plsda, comp=comp, title=paste("PLSDA Component", comp),
      row.sideColors=colours_cim, legend=list(title="Status"),
      col.names=show_cols
    )
    plotLoadings(data_plsda, contrib="max", comp=comp, max.name.length=16,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "PLSDA max loadings"))
    plotLoadings(data_plsda, contrib="min", comp=comp, max.name.length=16,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "PLSDA min loadings"))
    loading_max <- plotLoadings(data_plsda, contrib="max", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    loading_min <- plotLoadings(data_plsda, contrib="min", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    title <- gsub(" ", "_", title)
    path_max <- paste(outdir, "/", title, "_", comp, "_PLSDA_max.txt", sep="")
    path_min <- paste(outdir, "/", title, "_", comp, "_PLSDA_min.txt", sep="")
    print("Writing PLSDA loadings to:")
    print(path_max)
    print(path_min)
    write.table(as.data.frame(loading_max), file=path_max, quote=FALSE, sep="\t")
    write.table(as.data.frame(loading_min), file=path_min, quote=FALSE, sep="\t")
  }
  return(list(data_plsda=data_plsda, perf_plsda=metrics))
}

#' Tune sPLSDA
#'
#' Tune sPLSDA optimal components
#' @param data list of dataframes with original data.
#' @param classes 1-column dataframe with classes.
#' @param names vector of names for dataframes.
#' @param multilevel 1-column dataframe with repeated measurments. Defaults to NULL.
#' @param ncomp integer assigning number of components for PLSDA. Defaults to 0.
#' @param nrepeat integer of repetitions for cross-validation. Defaults to 10.
#' @param logratio string to outfile directory. Defaults to "./"
#' @param test_keepX vector of integers for feature selection. Defaults to c(5, 50, 100).
#' @param validation string describing validation "loo", "M-fold". Defaults to "loo"
#' @param folds integer describing folds for M-fold validation only. Defaults to 10.
#' @param dist string describing distance metric "centroids.dist", "max.dist", "mahalanobis.dist". Defaults to "centroids.dist"
#' @param cpus integer number of cpus for parallel processing. Defaults to 2.
#' @param progressBar boolean showing progress bar. Defaults to TRUE.
#' @param near_zero_var if variance approaches zero, enable this option
#' @seealso [multiomics::classify_splsda()]
#' @export
# @examples
# tune_splsda(data, classes, names, multilevel=NULL, ncomp=3, nrepeat=10, logratio="none", test_keepX=c(5, 50, 100), validation="loo", folds=10, dist="centroids.dist", cpus=2, progressBar=TRUE)
tune_splsda <- function(data, classes, names, multilevel=NULL, ncomp=3, nrepeat=10,
  logratio="none", test_keepX=c(5, 50, 100), validation="loo", folds=10,
  dist="centroids.dist", cpus=2, progressBar=TRUE, near_zero_var=FALSE) {
    # mapply(function(x, y) tune_splsda_(x, classes, names, multilevel, ncomp,
    #   nrepeat, logratio, test_keepX, validation, folds, dist, cpus, progressBar),
    #   data, names, SIMPLIFY=FALSE)
    data_new <- list()
    for(i in 1:length(names)){
      data_tmp <- tune_splsda_(
        data[[i]], classes, names[[i]], multilevel, ncomp, nrepeat, logratio,
        test_keepX, validation, folds, dist, cpus, progressBar, near_zero_var
      )
      data_new <- append(data_new, list(data_tmp))
    }
    names(data_new) <- names(data)
    return(data_new)
  }

tune_splsda_ <- function(data, classes, names, multilevel=NULL, ncomp=0, nrepeat=10,
  logratio="none", test_keepX=c(5, 50, 100), validation="loo", folds=10,
  dist="centroids.dist", cpus=2, progressBar=TRUE, near_zero_var=FALSE) {
  if (ncomp == 0) {ncomp <- (length(test_keepX))}
  # tune splsda components
  tuned <- tune.splsda(data, Y=classes, multilevel=multilevel,
    ncomp=ncomp, nrepeat=nrepeat, logratio=logratio, test.keepX=test_keepX,
    validation=validation, folds=folds, dist=dist, cpus=cpus,
    progressBar=progressBar, near.zero.var=near_zero_var, measure="BER"
  )
  if (nrepeat == 1 | validation == "loo") { sd <- FALSE } else { sd <- TRUE }
  print(mixOmics:::plot.tune.splsda(tuned, sd=sd, main=names))
  return(tuned)
}

#' Perform sPLSDA
#'
#' Perform sPLSDA (can account for repeated measurments if specified):
#' @param data list of dataframes with original data
#' @param classes 1-column dataframe with classes
#' @param pch 1-column dataframe with repeated measurments. Defaults to NA.
#' @param title string to name plots with. Defaults to "".
#' @param ncomp integer assigning number of components for PLSDA. Defaults to 0.
#' @param keepX vector of integers for feature selection. Defaults to NULL.
#' @param contrib string to determine loading weights "max", "min". Defaults to "max".
#' @param outdir string to outfile directory. Defaults to "./"
#' @param mappings dataframe of mappings for features. Defaults to NULL
#' @param dist string describing distance metric "centroids.dist", "max.dist", "mahalanobis.dist". Defaults to "centroids.dist"
#' @param bg boolean describing if background should be plotted. Defaults to TRUE
#' @param validation character specifying "loo" or "M-fold" cross-validation. Defaults to "loo".
#' @param folds if M-fold validation, number of folds. Defaults to 10. No effect for "loo"
#' @param nrepeat if M-fold validation, number of repeats. Defaults to 10. No effect for "loo"
#' @param near_zero_var if variance approaches zero, enable this option
#' @seealso [multiomics::classify_plsda()], [multiomics::classify_splsda()], [multiomics::tune_splsda()]
#' @export
# @examples
# classify_splsda(data, classes, pch=NA, title="", ncomp=0, keepX=NULL, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE)
classify_splsda <- function(data, classes, pch=NA, title="", ncomp=NULL,
  keepX=NULL, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist",
  bg=TRUE, validation="loo", nrepeat=10, folds=10, near_zero_var=FALSE) {
  # mapply(function(x, y, c, k) classify_splsda_(
  #   x, classes, pch, y, c, k, contrib, outdir
  # ), data, title, ncomp, keepX, SIMPLIFY=FALSE)
  print(keepX)
  print(ncomp)
  data_new <- list()
  for(i in 1:length(data)){
    if (length(ncomp) > 1) {ncomp_tmp <- ncomp[[i]]} else {ncomp_tmp <- ncomp}
    data_tmp <- classify_splsda_(
      data[[i]], classes, pch, title[[i]], ncomp_tmp, keepX[[i]],
      validation=validation, nrepeat=nrepeat, folds=folds, near_zero_var
    )
    data_new <- append(data_new, list(data_tmp))
  }
  names(data_new) <- names(data)
  return(data_new)
}

classify_splsda_ <- function(data, classes, pch=NA, title="", ncomp=NULL,
  keepX=NULL, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist",
  bg=TRUE, validation="loo", nrepeat=10, folds=10, near_zero_var=FALSE) {
  # discriminate samples: list, vector, bool, integer, vector -> list
  # single or multilevel sPLS-DA
  if (is.null(keepX)) {
    print("The number of variables selected on each component is not selected!")
    q(status=1)
  }
  if (is.null(ncomp)) {
    print("Invalid number of components, inferring from keepX")
    ncomp = length(keepX)
    print(ncomp)
  }
  print("splsda components:")
  print(ncomp)
  print("number of variables on each component:")
  print(keepX)
  classes <- as.factor(classes)

  if (length(pch) > 1) {
    data_splsda <- splsda(
      data,Y=classes, multilevel=pch, ncomp=ncomp, keepX=keepX,
      near.zero.var=near_zero_var
    )
  } else {
    data_splsda <- splsda(
      data, Y=classes, ncomp=ncomp, keepX=keepX, near.zero.var=near_zero_var
    )
  }

  if (ncomp > 1) {
    if (!is.na(bg)) {
      bg <- background.predict(data_splsda, comp.predicted=2, dist=dist)
      plotIndiv(data_splsda, ind.names=FALSE, group=classes,
        legend=TRUE, pch=pch, title=paste(title, "sPLSDA multi 1/2"),
        comp=c(1,2), ellipse=TRUE, background=bg
      )
    } else {
      plotIndiv(data_splsda, ind.names=FALSE, group=classes,
        legend=TRUE, pch=pch, title=paste(title, "sPLSDA multi 1/2"),
        comp=c(1,2), ellipse=TRUE, background=NULL
      )
    }
  }

  if (ncomp > 1) {
    plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "sPLSDA multi 1/2"), comp=c(1,2), ellipse=TRUE
    )
  }

  if (ncomp > 2) {
    plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "sPLSDA multi 1/3"), comp=c(1,3), ellipse=TRUE
    )
    plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "sPLSDA multi 2/3"), comp=c(2,3), ellipse=TRUE
    )
  }

  print("Getting performance metrics")
  print("Plotting error rates...")
  metrics <- mixOmics::perf(
    data_splsda, validation=validation, folds=folds, nrepeat=nrepeat,
    progressBar=TRUE, auc=TRUE, near_zero_var=low_var
  )
  print(metrics$error.rate)
  plot(metrics, main="Error rate sPLSDA", col=color.mixo(5:7), sd=TRUE)
  print("Plotting stability of sPLSDA...")
  plot(metrics$features$stable[[1]], type="h", main="Comp 1", las=2,
    ylab="Stability", xlab="Features", xaxt='n'
  )
  if (ncomp > 1) {
    plot(metrics$features$stable[[2]], type="h", main="Comp 2", las=2,
      ylab="Stability", xlab="Features", xaxt='n'
    )
  }
  if (ncomp > 2) {
    plot(metrics$features$stable[[3]], type="h", main="Comp 3", las=2,
      ylab="Stability", xlab="Features", xaxt='n'
    )
  }
  sink("/dev/null")
  # roc <- mapply(function(x) auroc(data_splsda, roc.comp=x), seq(ncomp))
  roc <- list()
  for(i in 1:ncomp){
    data_tmp <- auroc(data_splsda, roc.comp=i)
    roc <- append(roc, list(data_tmp))
  }
  sink()

  if (ncomp > 1) {
    print("Plotting arrow plot...")
    plotArrow(data_splsda, ind.names=FALSE, legend=TRUE, title="sPLSDA")
  }
  print("Getting loadings and plotting clustered image maps")

  # setup colour map for clustered image plots
  colours_class <- color.mixo(
    1:length(unique(classes))
  )[as.numeric(as.factor(classes))]

  if (length(pch) > 1) {
    colours_pch <- color.mixo(
      1:length(unique(pch))
    )[as.numeric(as.factor(pch))]
    colours_cim <- cbind(colours_class, colours_pch)
  } else {colours_cim <- data.frame(colours_class)}

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data_splsda$names$colnames$X, nchar))) > 8) {
    show_cols <- FALSE
  } else {
    show_cols <- TRUE
  }
  cim(data_splsda, title="sPLSDA", row.sideColors=colours_cim,
    legend=list(title="Status"), col.names=show_cols
  )

  short <- make.names(sapply(colnames(data), strtrim, 6, USE.NAMES=FALSE), unique=TRUE)
  for (comp in seq(ncomp)) {
    cim(data_splsda, comp=comp, title=paste("sPLSDA Component", comp),
      row.sideColors=colours_cim, legend=list(title="Status"),
      col.names=show_cols
    )
    plotLoadings(data_splsda, contrib="max", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=short, size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "sPLSDA max loadings"))
    plotLoadings(data_splsda, contrib="min", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=short, size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "sPLSDA min loadings"))
    loading_max <- plotLoadings(data_splsda, contrib="max", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    loading_min <- plotLoadings(data_splsda, contrib="min", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    title <- gsub(" ", "_", title)
    path_max <- paste(outdir, "/", title, "_", comp, "_sPLSDA_max.txt", sep="")
    path_min <- paste(outdir, "/", title, "_", comp, "_sPLSDA_min.txt", sep="")
    print("Writing sPLSDA loadings to:")
    print(path_max)
    print(path_min)
    write.table(as.data.frame(loading_max), file=path_max, quote=FALSE, sep="\t")
    write.table(as.data.frame(loading_min), file=path_min, quote=FALSE, sep="\t")
  }
  return(list(data_splsda=data_splsda, perf_splsda=metrics))
}

#' Plot PLSDA output
#'
#' Plot PLSDA output (PCA-like plots)
#' @param data list of dataframes with original data
#' @param classes 1-column dataframe with classes
#' @param pch 1-column dataframe with repeated measurments. Defaults to NA.
#' @param title string to name plots with. Defaults to "".
#' @param ncomp integer assigning number of components for PLSDA. Defaults to 0.
#' @export
# @examples
# plot_plsda(data, classes, pch, title="", ncomp=0)
plot_plsda <- function(data, classes, pch, title="", ncomp=0) {
  names <- names(data)
  print("Plotting PLSDA component contribution...")
  # mapply(function(x, y) plot(x, main=paste(y, "Screeplot multilevel")),
  #   data, names)
  for(i in 1:length(names)){
    plot(data[[i]], main=paste(names[[i]], "Screeplot multilevel"))
  }

  print("Plotting plsda...")
  for(i in 1:length(names)){
    plotIndiv(
      data_pca[[i]], comp=c(1,2), ind.names=TRUE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PLSDA 1/2"), pch=pch
    )
    plotIndiv(
      data_pca[[i]], comp=c(1,3), ind.names=TRUE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PLSDA 1/3"), pch=pch
    )
    plotIndiv(
      data_pca[[i]], comp=c(2,3), ind.names=TRUE, group=classes, legend=TRUE,
      ncomp=ncomp, title=paste(title, names[[i]], "PLSDA 2/3"), pch=pch
    )
  }
  # mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PLSDA 1/2"), pch=pch), data_pca, names)
  # mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PLSDA 1/3"), pch=pch), data_pca, names)
  # mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
  #   group=classes, legend=TRUE, ncomp=ncomp,
  #   title=paste(title, y, "PLSDA 2/3"), pch=pch), data_pca, names)
  return(data_pca)
}

#' Plot PLSDA loadings
#'
#' Plot PLSDA loadings (weight of variable contribution to category)
#' @param data list of dataframes with original data
#' @param contrib string to determine loading weights "max", "min". Defaults to "max".
#' @param ncomp integer assigning number of components for plot. Defaults to 2.
#' @param method string describing method of computation. Defaults to "median"
#' @param ndisplay integer number of features to display. Defaults to 50.
#' @param title string to name plots with. Defaults to "Loadings".
#' @export
# @examples
# plot_loadings(data, contrib="max", ncomp=2, method="median", ndisplay=50, title="Loadings")
plot_loadings <- function(data, contrib="max", ncomp=2, method="median",
  ndisplay=50, title="Loadings") {
  # mapply(function(x, y) plot_loadings_(x, contrib, ncomp, method, ndisplay,
  #   y), data, title)
  for (i in length(data)) {
    plot_loadings_(data[[i]], contrib, ncomp, method, ndisplay, title[[i]])
  }
}

plot_loadings_ <- function(data, contrib="max", ncomp=2, method="median",
  ndisplay=50, title="Loadings") {
  name_var <- names(data)
  loadings <- plotLoadings(
    data, contrib, ncomp, method, ndisplay, name_var, title
  )
  print(loadings)
  path <- paste(title, "txt", sep=".")
  print(path)
  write.table(as.data.frame(loadings), file=path, quote=FALSE, sep="\t")
  return(loadings)
}

#' Plot correlation between unimputed and imputed data
#'
#' Quality control plot to check if data has mutated significantly
#' @param pca_withna list of dataframes with original data
#' @param pca_impute list of dataframes with imputed data
#' @param names vector of names for each data category
#' @seealso [multiomics::impute_missing()], [multiomics::replace_missing()], [multiomics::cor_imputed_unimputed()]
#' @export
# @examples
# cor_imputed_unimputed(pca_withna, pca_impute, names)
cor_imputed_unimputed <- function(pca_withna, pca_impute, names) {
  # plots a heatmap of correlations: -> list of df, list of df, vector of names
  print("Plotting correlation between unimputed and imputed components")
  # mapply(function(x, y, z) print(
  #   ggplot2::ggplot(reshape2::melt(cor(x$variates$X, y$variates$X)),
  #   ggplot2::aes(Var1, Var2, fill=value)) +
  #   ggplot2::ggtitle(paste(z, "Correlation between imputed and unimputed data")) +
  #   ggplot2::geom_tile() +
  #   ggplot2::scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0,
  #     limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
  #   ggplot2::theme_minimal()),
  # pca_withna, pca_impute, names)
  for (i in 1:length(names)) {
    print(
      ggplot2::ggplot(reshape2::melt(
        cor(pca_withna[[i]]$variates$X, pca_impute[[i]]$variates$X)
      ), ggplot2::aes(Var1, Var2, fill=value)) +
      ggplot2::ggtitle(
        paste(names[[i]], "Correlation between imputed and unimputed data")
      ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low="blue", high="red", mid="white", midpoint=0, limit=c(-1,1),
        space="Lab", name="Pearson\nCorrelation"
      ) +
      ggplot2::theme_minimal()
    )
  }
}

#' Tune sPLSDA multi-block (DIABLO) number of components
#'
#' Tune sPLSDA multi-block optimal components
#' @param data list of dataframes with original data.
#' @param classes 1-column dataframe with classes.
#' @param design matrix with linkage
#' @param ncomp integer assigning number of components for multi-block sPLSDA. Defaults to 0.
#' @param cpus integer number of cpus for parallel processing. Defaults to 2.
#' @param near_zero_var if variance approaches zero, enable this option
#' @seealso [multiomics::create_design()], [multiomics::tune_diablo_ncomp()], [multiomics::tune_diablo_keepx()], [multiomics::run_diablo()]
#' @export
# @examples
# tune_diablo_ncomp(data, classes, ncomp=0, design, cpus=2)
tune_diablo_ncomp <- function(data, classes, design, ncomp=0, cpus=1,
  near_zero_var=FALSE, nrepeat=1, validation="loo", folds=1) {
  # First, we fit a DIABLO model without variable selection to assess the global
  # performance and choose the number of components for the final DIABLO model.
  # The function perf is run with 10-fold cross validation repeated 10 times.
  print("Finding optimal number of components for DIABLO...")
  if (ncomp == 0) {
    ncomp <- length(unique(classes)) - 1
    ncomp <- min(2, ncomp)
  }
  sgccda_res <- block.splsda(
    X=data, Y=classes, ncomp=ncomp, design=design, near.zero.var=near_zero_var
  )

  # this code takes a couple of min to run
  perf_diablo <- perf(
    sgccda_res, validation=validation, folds=folds, nrepeat=nrepeat, cpus=cpus
  )
  print(perf_diablo$error.rate)
  plot(perf_diablo, main="DIABLO optimal components")
  print(perf_diablo$choice.ncomp)

  print("Plotting stability of DIABLO components...")
  for (i in perf_diablo$features$stable$nrep1) {
    for (j in names(perf_diablo$features$stable$nrep1)) {
      if (any(is.infinite(i$comp1)) == FALSE) {
        plot(i$comp1, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 1"), xaxt='n'
        )
        print(paste(j, "Comp 1"))
      }
      if (any(is.infinite(i$comp2)) == FALSE) {
        plot(i$comp2, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 2"), xaxt='n'
        )
        print(paste(j, "Comp 2"))
      }
      if (any(is.infinite(i$comp3)) == FALSE) {
        plot(i$comp1, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 3"), xaxt='n'
        )
        print(paste(j, "Comp 3"))
      }
    }
  }
  return(perf_diablo)
}

#' Tune sPLSDA multi-block (DIABLO) features to keep
#'
#' Tune sPLSDA multi-block features to keep
#' @param data list of dataframes with original data.
#' @param classes 1-column dataframe with classes.
#' @param ncomp integer assigning number of components for multi-block sPLSDA. Defaults to 0.
#' @param design matrix with linkage
#' @param test_keepX vector of integers for feature selection. Defaults to c(5, 50, 100).
#' @param dist string describing distance metric "centroids.dist", "max.dist", "mahalanobis.dist". Defaults to "centroids.dist"
#' @param cpus integer number of cpus for parallel processing. Defaults to 2.
#' @param progressBar boolean showing progress bar. Defaults to TRUE.
#' @param validation character specifying "loo" or "M-fold" cross-validation. Defaults to "loo".
#' @param folds if M-fold validation, number of folds. Defaults to 10. No effect for "loo"
#' @param nrepeat if M-fold validation, number of repeats. Defaults to 10. No effect for "loo"
#' @param near_zero_var if variance approaches zero, enable this option
#' @seealso [multiomics::create_design()], [multiomics::tune_diablo_ncomp()], [multiomics::tune_diablo_keepx()], [multiomics::run_diablo()]
#' @export
# @examples
# tune_diablo_keepx(data, classes, ncomp, design, test_keepX=c(5,50,100), cpus=2, dist="centroids.dist", progressBar=TRUE)
tune_diablo_keepx <- function(data, classes, ncomp, design,
  test_keepX=c(5,50,100), cpus=2, dist="centroids.dist", progressBar=TRUE,
  validation="loo", folds=10, nrepeat=10, near_zero_var=FALSE) {
  # This tuning function should be used to tune the keepX parameters in the
  #   block.splsda function.
  # We choose the optimal number of variables to select in each data set using
  # the tune function, for a grid of keepX values. Note that the function has
  # been set to favor the small-ish signature while allowing to obtain a
  # sufficient number of variables for downstream validation / interpretation.
  # See ?tune.block.splsda.
  print("Tuning keepX parameter...")
  print(test_keepX)
  test_keepX <- mapply(function(name, dims) list(name=dims), names(data),
    rep(list(test_keepX))
  )
  print(test_keepX)
  # test_keepX <- list()
  # for(i in 1:length(rep(list(test_keepX)))){
  #   data_tmp <- list(names(data)[[i]]=rep(list(test_keepX))[[i]])
  #   test_keepX <- append(test_keepX, list(data_tmp))
  # }

  cpus <- NULL
  tune_data <- tune.block.splsda(
    X=data, Y=classes, ncomp=ncomp, test.keepX=test_keepX, design=design,
    validation=validation, folds=folds, nrepeat=nrepeat, dist=dist,
    progressBar=progressBar, near.zero.var=near_zero_var)
  list_keepX <- tune_data$choice.keepX
  return(list_keepX)
}

#' Force unique blocks
#'
#' Force each feature in each block of omics data to have unique names
#' @param data dataframe with original data.
#' @export
# @examples
# force_unique_blocks(data)
force_unique_blocks <- function(data) {
  mapply(data, names(data), FUN = function(blockmat, blockname){
    colnames(blockmat) <- paste0(blockname, '_', colnames(blockmat) )
    blockmat
  }, SIMPLIFY = FALSE)
}

#' Run the multi-omic pipeline
#'
#' Run the multi-block splsda on the multi omics data
#' @param data list of dataframes with quantitative multi-omics data.
#' @param classes 1-column dataframe with classes.
#' @param ncomp integer assigning number of components for multi-block sPLSDA.
#' @param design matrix with linkage.
#' @param keepx vector of integers for feature selection. Defaults to NULL.
#' @param near_zero_var if variance approaches zero, enable this option
#' @export
# @examples
# run_diablo(data, classes, ncomp, design, keepx=NULL)
run_diablo <- function(data, classes, ncomp, design, keepx=NULL,
  near_zero_var=FALSE) {
  # this is the actual part where diablo is run
  print("Running DIABLO...")
  for (i in names(keepx)) {
    if (ncomp < length(keepx[[i]])) {
      keepx[[i]] <- head(keepx[[i]], n=ncomp)
    }
  }

  data <- block.splsda(
    X=data, Y=classes, ncomp=ncomp, keepX=keepx, design=design,
    near.zero.var=near_zero_var
  )
  return(data)
}

#' Make plots for diablo
#'
#' Run the multi-block splsda on the multi omics data
#' @param data list of dataframes with quantitative multi-omics data.
#' @param ncomp integer assigning number of components for multi-block sPLSDA. Defaults to 0.
#' @param outdir string to outfile directory. Defaults to "./"
#' @param data_names names of individual omics data blocks. Defaults to NA.
#' @param keepvar name of kept variable. Defaults to "".
#' @param cutoff display correlations that exceed this threshold
#' @export
# @examples
# plot_diablo(data, ncomp=0, outdir="./", data_names=NA, keepvar="")
plot_diablo <- function(data, ncomp=0, outdir="./", data_names=NA, keepvar="",
  cutoff=0.95) {
  # plot the diablo data with a series of diagnostic plots

  # need to make a function to squeeze sample names automatically and remap
  trim_names_ <- function(data, trim=6) {
    all_names <- data
    long_names <- which(sapply(data, nchar, USE.NAMES=FALSE) > trim)
    if (length(long_names) == 0) {return()}
    original_names <- data[long_names]
    for (i in long_names) { data[i] <- as.character(i) }
    # later map these back
    maptable <- data.frame(from=long_names, to=original_names)
    return(list(data=data, all_names=all_names, maptable=maptable))
  }
  trimmed_names <- lapply(data$names$colnames, trim_names_)
  block_to_trim <- names(trimmed_names[lapply(trimmed_names, length) > 0])

  # replace names in all associated columns for visualisation only
  replace_names_ <- function(data, trim=6) {
    all_names <- head(data$names$colnames, n=-1)
    split_ <- function(block, all_names) {
      all_names <- gsub("__FEATUREID", "", all_names[[block]])
      all_names <- gsub(paste("_", block, sep=""), "", all_names)
      return(all_names)
    }
    splitted <- mapply(
      function(x) split_(x, all_names), names(all_names), SIMPLIFY=FALSE
    )
    # splitted <- list()
    # for(i in 1:length(names(all_names))) {
    #   data_tmp <- split_(names(all_names)[[i]], all_names)
    #   splitted <- append(splitted, data_tmp)
    # }
    truncate_ <- function(names, trim) {
      ifelse(nchar(names) > trim, paste0(strtrim(names, trim), ''), names)
    }
    truncated <- lapply(splitted, truncate_, trim)
    # make a copy, dont want to overwrite
    data_vis <- data
    # map truncated values to all locations
    for (i in names(all_names)) {
      row.names(data_vis$loadings[[i]]) <- make.unique(
        sapply(truncated[[i]], toString), sep="__"
      )
      data_vis$names$colnames[[i]] <- make.unique(
        sapply(truncated[[i]], toString), sep="__"
      )
      colnames(data_vis$X[[i]]) <- make.unique(
        sapply(truncated[[i]], toString), sep="__"
      )
    }
    return(list(data_vis=data_vis, truncated=truncated))
  }
  data_vis_names <- replace_names_(data, trim=6)
  data_vis <- data_vis_names$data_vis
  truncated <- data_vis_names$truncated
  print("Plotting correlation between components...")
  sink("/dev/null")
  # roc <- mapply(function(x) auroc(data, roc.comp=x), seq(ncomp))
  roc <- list()
  for(i in 1:ncomp){
    data_tmp <- auroc(data, roc.comp=i)
    roc <- append(roc, list(data_tmp))
  }
  sink()
  # mapply(function(x) plotDiablo(data, ncomp=x), seq(ncomp))
  for(i in 1:ncomp) {plotDiablo(data, ncomp=i)}
  if (ncomp > 1) {
    print("Plotting individual samples into space spanned by block components...")
    plotIndiv(data_vis, ind.names=FALSE, legend=TRUE, title='DIABLO',
      ellipse=TRUE
    )
    print("Plotting arrow plot...")
    plotArrow(data_vis, ind.names=FALSE, legend=TRUE, title='DIABLO')
    print("Plotting correlation circle plot...")
    plotVar(data_vis, style='graphics', legend=TRUE, comp=c(1,2),
      title="DIABLO 1/2", var.names=FALSE
    )
  }
  if (ncomp > 2) {
    plotVar(data_vis, style='graphics', legend=TRUE, comp=c(1,3),
      title="DIABLO 1/3", var.names=FALSE
    )
    plotVar(data_vis, style='graphics', legend=TRUE, comp=c(2,3),
      title="DIABLO 2/3", var.names=FALSE
    )
  }
  print("Plotting circos from similarity matrix...")
  # cant remove feature labels, need to make label size 0.001 or lower
  corr_diablo <- circosPlot(
    data, cutoff=cutoff, line=TRUE, size.legend=0.5, size.variables=0.001,
    var.names=truncated
  )
  corr_out <- paste(outdir,"/DIABLO_var_",keepvar,"_correlations.txt",sep="")
  write.table(corr_diablo, file=corr_out, sep="\t", quote=FALSE)
  print("Plotting relevance network from similarity matrix...")
  cyto <- network(
    data, blocks=c(1,2), color.node=c('darkorchid','lightgreen'), cutoff=cutoff
  )
  cyto_out <- paste(outdir, "/DIABLO_var_", keepvar, "_network.graphml", sep="")
  igraph::write.graph(cyto$gR, cyto_out, format="graphml")
  print("Plotting overall heatmap...")

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data$names$colnames, nchar))) > 6) {
    show_cols = FALSE
  } else {
    show_cols = TRUE
  }
  cimDiablo(data, size.legend=0.5, col.names=show_cols)
  block_to_trim <- names(trimmed_names[lapply(trimmed_names, length) > 0])
  print("Plotting loading weight of selected variables on each component...")
  for (comp in seq(ncomp)) {
    for (i in block_to_trim) {
      data$names$colnames[[i]] = trimmed_names[[i]][["data"]]
    }
    cimDiablo(data, comp=comp, size.legend=0.5, col.names=show_cols)
    plotLoadings(data, contrib="max", comp=comp, max.name.length=6,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(comp, "DIABLO max loadings"))
    plotLoadings(data, contrib="min", comp=comp, max.name.length=6,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(comp, "DIABLO min loadings"))
    for (i in block_to_trim) {
      data$names$colnames[[i]] = trimmed_names[[i]][["all_names"]]
    }

    for (j in data_names) {
      for (i in block_to_trim) {
        data$names$colnames[[i]] <- trimmed_names[[i]][["data"]]
      }
      plotLoadings(data, contrib="max", comp=comp, block=j,
        max.name.length=6, method='median', ndisplay=20,
        name.var=colnames(data), plot=TRUE,
        title=paste(comp, j, "DIABLO max loadings"), size.name=0.6
      )
      plotLoadings(data, contrib="min", comp=comp, block=j,
        max.name.length=6, method='median', ndisplay=20,
        name.var=colnames(data), plot=TRUE,
        title=paste(comp, j, "DIABLO min loadings"), size.name=0.6
      )
      for (i in block_to_trim) {
        data$names$colnames[[i]] <- trimmed_names[[i]][["all_names"]]
      }

      loading_max <- plotLoadings(
        data, contrib="max", comp=comp, block=j, method='median', ndisplay=NULL,
        name.var=colnames(data), plot=FALSE
      )
      loading_min <- plotLoadings(
        data, contrib="min", comp=comp, block=j, method='median', ndisplay=NULL,
        name.var=colnames(data), plot=FALSE
      )

      path_max <- paste(
        outdir, "/", j, "_", comp, "_DIABLO_var_", keepvar, "_max.txt", sep=""
      )
      path_min <- paste(
        outdir, "/", j, "_", comp, "_DIABLO_var_", keepvar, "_min.txt", sep=""
      )
      print("Writing DIABLO loadings to:")
      print(path_max)
      print(path_min)
      write.table(as.data.frame(loading_max),file=path_max,quote=FALSE,sep="\t")
      write.table(as.data.frame(loading_min),file=path_min,quote=FALSE,sep="\t")
    }
  }
}

#' Assess multi-block sPLSDA (DIABLO) performance
#'
#' Assess multi-block sPLSDA (DIABLO) performance
#' @param data list of dataframes with quantitative multi-omics data.
#' @param dist string describing distance metric "centroids.dist", "max.dist", "mahalanobis.dist". Defaults to "centroids.dist"
#' @param ncomp integer assigning number of components for multi-block sPLSDA. Defaults to 0.
#' @export
#' @keywords Internal
assess_performance <- function(data, dist, ncomp) {
  # review performance of diablo
  # remember to use the same distance metric which had the max value!
  # print("Assessing performance...")
  # perf_diablo = perf(data, validation='loo', M=10, nrepeat=10, dist=dist)
  # perf.diablo  # lists the different outputs

  # Performance with Majority vote
  # print(perf_diablo$MajorityVote.error.rate)

  # ROC and AUC criteria are not particularly insightful in relation to the
  # performance evaluation of our methods, but can complement the analysis.
  print("Plotting ROC...")
  # mapply(function(x) auroc(data, x), seq(ncomp))
  for(i in 1:ncomp) {auroc(data, roc.comp=i)}
  return(perf_diablo)
}

#' Predict with weights
#'
#' Predict diablo with weights
#' @param data dataframe with quantitative multi-omics data.
#' @param test dataframe with test data
#' @param classes 1-column dataframe of classes
#' @export
#' @keywords Internal
predict_diablo <- function(data, test, classes) {
  # prepare test set data: here one block (proteins) is missing
  print("Predicting data on an external test set...")
  predict.diablo <- predict(data, newdata = test)
  # the warning message will inform us that one block is missing
  #predict.diablo # list the different outputs
  print("Getting confusion matrix...")
  confusion_matrix <- get.confusion_matrix(
    truth=classes, predicted=predict.diablo$WeightedVote$max.dist[,2])
  print(confusion_matrix)
  print(get.BER(confusion_matrix))
}
