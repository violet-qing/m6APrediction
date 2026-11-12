#' @title Encode DNA Sequences
#' @description Convert DNA sequences into nucleotide-position factors.
#' @param dna_strings A character vector of equal-length DNA sequences
#'   containing uppercase A, T, C and G.
#' @return A data frame with one factor column for each nucleotide position.
#' @examples
#' m6APrediction:::dna_encoding(c("ATCGA", "GCTAT"))
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])

  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )

  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)

  seq_df[] <- lapply(
    seq_df, factor,
    levels = c("A", "T", "C", "G")
  )

  return(seq_df)
}
#' @title Predict m6A Sites for Multiple Samples
#' @description Predict m6A probabilities and classes for multiple samples
#'   using a trained random forest model.
#' @param ml_fit A trained random forest classification model.
#' @param feature_df A data frame containing gc_content, RNA_type,
#'   RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer.
#' @param positive_threshold A probability cutoff. Default is 0.5.
#'   Probabilities above this cutoff are classified as Positive.
#' @return The input data frame with predicted_m6A_prob and
#'   predicted_m6A_status added.
#' @examples
#' model <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' example_data <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' ))
#' prediction_multiple(model, example_data)
#' @export
#' @import randomForest
prediction_multiple <- function(
    ml_fit, feature_df, positive_threshold = 0.5
) {
  stopifnot(all(c(
    "gc_content", "RNA_type", "RNA_region",
    "exon_length", "distance_to_junction",
    "evolutionary_conservation", "DNA_5mer"
  ) %in% colnames(feature_df)))

  prediction_df <- feature_df

  prediction_df$RNA_type <- factor(
    prediction_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )

  prediction_df$RNA_region <- factor(
    prediction_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  prediction_df <- cbind(
    prediction_df,
    dna_encoding(as.character(prediction_df$DNA_5mer))
  )

  probabilities <- stats::predict(
    ml_fit,
    newdata = prediction_df,
    type = "prob"
  )

  feature_df$predicted_m6A_prob <- probabilities[, "Positive"]

  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df)
}
#' @title Predict an m6A Site for One Sample
#' @description Predict the m6A probability and class for a single sample.
#' @param ml_fit A trained random forest classification model.
#' @param gc_content GC content of the sample.
#' @param RNA_type RNA type: mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region RNA region: CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Exon length.
#' @param distance_to_junction Distance to the splice junction.
#' @param evolutionary_conservation Evolutionary conservation score.
#' @param DNA_5mer A five-nucleotide DNA sequence in uppercase.
#' @param positive_threshold A probability cutoff. Default is 0.5.
#'   Probabilities above this cutoff are classified as Positive.
#' @return A character vector containing predicted_m6A_prob and
#'   predicted_m6A_status.
#' @examples
#' model <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' example_data <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' ))
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = example_data$gc_content[1],
#'   RNA_type = example_data$RNA_type[1],
#'   RNA_region = example_data$RNA_region[1],
#'   exon_length = example_data$exon_length[1],
#'   distance_to_junction = example_data$distance_to_junction[1],
#'   evolutionary_conservation =
#'     example_data$evolutionary_conservation[1],
#'   DNA_5mer = example_data$DNA_5mer[1]
#' )
#' @export
prediction_single <- function(
    ml_fit, gc_content, RNA_type, RNA_region,
    exon_length, distance_to_junction,
    evolutionary_conservation, DNA_5mer,
    positive_threshold = 0.5
) {
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  result <- prediction_multiple(
    ml_fit, feature_df, positive_threshold
  )

  return(c(
    predicted_m6A_prob = result$predicted_m6A_prob[1],
    predicted_m6A_status = result$predicted_m6A_status[1]
  ))
}
