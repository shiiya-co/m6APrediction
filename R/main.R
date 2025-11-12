#' DNA Sequence Encoding (Internal Function)
#'
#' Converts DNA/RNA sequence strings into a data frame with factor levels.
#' It automatically converts "U" to "T".
#' This is an internal helper function for the prediction functions.
#'
#' @param dna_strings A character vector of DNA/RNA sequences.
#' @return A data frame where each column represents a nucleotide position with factor values.
#' @keywords internal
dna_encoding <- function(dna_strings){
  dna_strings_T <- toupper(gsub("U", "T", dna_strings))

  nn <- nchar(dna_strings_T[1])
  seq_m <- matrix(unlist(strsplit(dna_strings_T, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Entries Based on Multiple Features
#'
#' Uses a trained random forest model to predict m6A status from a data frame
#' containing multiple samples and features.
#'
#' @import randomForest
#' @importFrom stats predict
#' @param feature_df A data frame that must contain the following columns: "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer".
#' @param model A trained random forest model object.
#' @param positive_threshold A numeric value between 0 and 1, used as the cutoff to classify predictions as Positive or Negative. Default is 0.5.
#' @return Returns the input data frame with two additional columns: "predicted_m6A_prob" (the probability of being positive) and "predicted_m6A_status" (the predicted status).
#' @export
#' @examples
#' \dontrun{
#'   rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'   example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'   predictions <- prediction_multiple(feature_df = example_df, model = rf_model)
#'   print(head(predictions))
#' }
prediction_multiple <- function(feature_df, model, positive_threshold = 0.5){
  required_cols <- c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  if (!all(required_cols %in% colnames(feature_df))) {
    stop("Input data frame is missing required columns.")
  }

  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  encoded_dna <- dna_encoding(feature_df$DNA_5mer)

  prediction_input_df <- cbind(feature_df, encoded_dna)

  predicted_probs <- predict(model, newdata = prediction_input_df, type = "prob")

  feature_df$predicted_m6A_prob <- predicted_probs[, 2]
  feature_df$predicted_m6A_status <- ifelse(feature_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")

  return(feature_df)
}

#' Predict a Single m6A Site Based on Multiple Features
#'
#' Uses a trained random forest model to predict the status of a single m6A site
#' based on a set of input features.
#'
#' @import randomForest
#' @param model A trained random forest model object.
#' @param gc_content GC content (numeric).
#' @param RNA_type RNA type (character, e.g., "mRNA").
#' @param RNA_region RNA region (character, e.g., "CDS").
#' @param exon_length Exon length (numeric).
#' @param distance_to_junction Distance to junction (numeric).
#' @param evolutionary_conservation Evolutionary conservation (numeric).
#' @param DNA_5mer A 5-character DNA sequence (e.g., "GGACA").
#' @param positive_threshold The cutoff for classifying as positive. Default is 0.5.
#' @return A named vector containing "predicted_m6A_prob" and "predicted_m6A_status".
#' @export
#' @examples
#' \dontrun{
#'   rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'   single_pred <- prediction_single(
#'     model = rf_model,
#'     gc_content = 0.6,
#'     RNA_type = "mRNA",
#'     RNA_region = "CDS",
#'     exon_length = 12.5,
#'     distance_to_junction = 5,
#'     evolutionary_conservation = 0.8,
#'     DNA_5mer = "GGACU"
#'   )
#'   print(single_pred)
#' }
prediction_single <- function(model, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  input_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = as.character(RNA_type),
    RNA_region = as.character(RNA_region),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = as.character(DNA_5mer),
    stringsAsFactors = FALSE
  )

  result_df <- prediction_multiple(feature_df = input_df, model = model, positive_threshold = positive_threshold)

  prob_value <- result_df$predicted_m6A_prob[1]
  status_value <- result_df$predicted_m6A_status[1]

  returned_vector <- c(prob_value, status_value)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")

  return(returned_vector)
}
