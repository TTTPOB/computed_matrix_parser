suppressPackageStartupMessages({
    library(tidyverse)
    library(jsonlite)
})

#' @export
read_computed_matrix <- function(matrix_path, offset = 6, df_obj = NULL) {
    if (is.null(df_obj)) {
        df <- read_tsv(
            skip = 1,
            file = df_path,
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
    metadata <- readLines(gzfile(df_path), 1) %>%
        str_remove("^@") %>%
        fromJSON(json_str = .)
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
    )
    res <- c(df_splited, list(metadata))
    names(res) <- c(metadata$sample_labels, "metadata")
    return(res)
}