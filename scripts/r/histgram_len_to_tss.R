#!/usr/bin/env Rscript

# 加载ggplot2库
library(ggplot2)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  # 第一个参数是文件路径
output_file <- ifelse(length(args) > 1, args[2], paste(dirname(input_file), "/histogram_distance_tss.png", sep=""))

# 读取数据
data <- read.table(input_file, header = FALSE, col.names = c("Gene", "Distance"))
# 过滤掉大于10000的距离值
data <- data[data$Distance <= 10000, ]
data <- data[data$Distance > 0, ]
# 绘制直方图
p <- ggplot(data, aes(x = Distance)) +
  geom_histogram(binwidth = 100, fill = "#1ba784", color = "#1ba784") +
  labs(title = "Distance to TSS Histogram", x = "Distance", y = "Count") +
  theme_minimal()

# 保存图像
ggsave(output_file, plot = p, width = 8, height = 5)

# 打印保存图像的路径
cat("          Histogram saved to: histogram_distance_tss.png")
