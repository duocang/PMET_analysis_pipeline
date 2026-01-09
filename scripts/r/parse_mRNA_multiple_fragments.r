suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})


# 从命令行获取参数
args <- commandArgs(trailingOnly = TRUE)

# 确保提供了足够的参数
if (length(args) < 3) {
  cat("Usage: script_name.R binomial_thresholds_file fimohits_dir output_dir\n")
  quit(status = 1)
}

# 获取参数值
binomial_thresholds_file <- args[1]
fimohits_dir             <- args[2]
output_dir               <- args[3]

dir.create(output_dir, showWarnings = FALSE)

binomial_thresholds <- fread(binomial_thresholds_file) %>% as.data.frame()

for (i in 1:nrow(binomial_thresholds)) {
  motif <- binomial_thresholds[i, "V1"]
  binomial_threshold <- binomial_thresholds[i, "V2"]
  a <- file.path(fimohits_dir, paste0(motif, ".txt")) %>%
    fread() %>% as.data.frame()

  a$V2 <- str_remove(a$V2, "__\\d+") %>% str_remove("__")
  a <- a %>%
    arrange(V2, V7) %>%
    group_by(V2) %>%
    arrange(V7) %>%
    slice_head(n = 5) %>%
    filter(V7 < binomial_threshold)

  write.table(a, file.path(output_dir, paste0(motif, ".txt")),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t")
}