import::here(R6, R6Class)
import::here(magrittr, "%>%")
import::here(purrr, map, reduce, imap)
import::here(dplyr, .all = TRUE)
import::here(readr, .all = TRUE)
import::here(jsonlite, .all = TRUE)
import::here(stringr, .all = TRUE)

#' @export
computed_matrix <- R6Class("computed_matrix", list(
    metadata = NULL,
    sample_matrix_list = list(),
    sample_mean_coverage = list(),
    region_total_coverage = list(),
    initialize = function(metadata, sample_matrix) {
        self$metadata <- metadata
        self$sample_matrix_list <- sample_matrix
    },
    calculate_mean_coverage_each_sample = function() {
        dflist <- self$sample_matrix_list
        sample_mean_coverage <- dflist %>%
            imap(function(df, name) {
                splited_by_group_label <- df %>%
                    group_by(group_label) %>%
                    group_split() %>%
                    setNames(df$group_label %>% unique())

                mean_coverage_per_group <- splited_by_group_label %>%
                    imap(function(df, group_label) {
                        df <- df %>% select(!group_label)
                        mtx <- df %>% as.matrix()
                        colMeans <- mtx %>% colMeans()
                        return(colMeans)
                    })
                return(mean_coverage_per_group)
            })
        self$sample_mean_coverage <- sample_mean_coverage
    },
    calculate_total_coverage_each_region = function() {
        total_coverage_each_region <- self$sample_matrix_list %>%
            imap(function(df, name) {
                df %>%
                    group_by(group_label) %>%
                    group_split() %>%
                    setNames(df$group_label %>% unique()) %>%
                    imap(function(df, group_label) {
                        df %>%
                            select(!group_label) %>%
                            as.matrix() %>%
                            rowSums()
                    })
            })
        self$region_total_coverage <- total_coverage_each_region
    }
))

#' @export
read_computed_matrix <- function(matrix_path, offset = 6, df_obj = NULL) {
    if (is.null(df_obj)) {
        df <- read_tsv(
            skip = 1,
            file = matrix_path,
            col_names = F,
            col_types = cols(
                X1 = "c",
                X2 = "i",
                X3 = "i",
                X4 = "c",
                X5 = "d",
                X6 = "c",
                .default = "d"
            )
        ) %>%
            mutate_if(
                is.numeric,
                ~ replace(., is.na(.), 0)
            )
    } else {
        df <- df_obj
    }
    metadata <- readLines(gzfile(matrix_path), 1) %>%
        str_remove("^@") %>%
        fromJSON()
    df_splited <- list()
    df_splited <- map(
        seq_len(length(metadata$sample_labels)),
        function(i) {
            sample_left_bond <- offset + 1 + metadata$sample_boundaries[i]
            sample_right_bond <- offset + metadata$sample_boundaries[i + 1]
            df_sample <- df[, sample_left_bond:sample_right_bond]
            if ("group_labels" %in% names(metadata)) {
                group_bonds <- metadata$group_boundaries
                group_labels <- metadata$group_labels
                group_label <- map(seq_len(length(metadata$group_labels)), function(j) {
                    group_ <- rep_len(
                        group_labels[j],
                        group_bonds[j + 1] - group_bonds[j]
                    )
                    return(group_)
                }) %>%
                    reduce(~ c(.x, .y))
                df_sample$group_label <- group_label
            }
            return(df_sample)
        }
    ) %>%
        setNames(metadata$sample_labels)
    res <- computed_matrix$new(
        metadata = metadata,
        sample_matrix = df_splited
    )
    return(res)
}