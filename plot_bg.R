#!/bin/env Rscript

# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Bna background evaluate plot")

# Add command line arguments
p <- add_argument(p, "--chrLen", help = "a file contain chromosome file", type = "character")
p <- add_argument(p, "--marker", help = "a file contain marker position, created by bg.R", type = "character")
p <- add_argument(p, "--region", help = "a file contain introgression fragment", type = "character")
p <- add_argument(p, "--output", help = "output file prefix", type = "character")

# Parse the command line arguments
argv <- parse_args(p)

library(tidyverse)
library(plot3D)

## 读取数据
# 染色体长度
chr <- read_tsv(argv$chrLen)

# 标记位置
marker <- read_tsv(argv$marker) %>%
  select(CHROM = chr, POS = mid, VALUE = sum)

# 导入片段位置
fragment_file <- argv$region
if (exists("fragment_file")) {
  fragment <- read_tsv(fragment_file)
}

test <- FALSE
if (test) {
  chr <- read_tsv("./chromosome_length.txt")
  marker <- read_tsv("./WS100.bg_estimate.txt") %>%
    select(CHROM = chr, POS = mid, VALUE = sum)
  fragment <- read_tsv("./WS100.region.txt")
  argv$output <- "WS100.region.txt"
}

# 定义标记颜色梯度
color <- colorRampPalette(c("blue", "red"))(max(marker$VALUE)-min(marker$VALUE)+1)

## 打开画板
png(filename = paste(argv$output, "png", sep = "."), width = 8, height = 7, units = "in", res = 500)
# 产生面板
par(mar = c(5, 5, 1, 6))
plot(0, 0, xlim = c(0, max(chr$LEN)), ylim = c (-nrow(chr)*3+1, 0), type = "n", 
     axes = F, xlab = "Position (Mb)", ylab = "Chromosome", mar = c(0, 0, 0, 1))
  # 坐标轴
axis(side = 1, at = seq(0, ceiling(max(chr$LEN)/10000000)*10000000, by = 10000000), 
     labels = paste(seq(0, ceiling(max(chr$LEN)/1000000), by = 10), "M", sep = ""))
axis(side = 2, at = -seq(1, 3*nrow(chr), by = 3), labels = chr$CHROM, las="2")

for(i in 1:nrow(chr)){
  # 画染色体
  rect(xleft = 0, ybottom = -3*i+1, xright = chr$LEN[i], ytop = -3*i+3 , col = "grey", lwd = 2)

  if (exists("fragment")) {
    # 画标记
    marker_sub <- marker[marker$CHROM==chr$CHROM[i], ]
    for (j in 1:nrow(marker_sub)) {
      lines(x = c(marker_sub$POS[j], c(marker_sub$POS[j])), y = c(-3*i+1, -3*i+2), col = color[marker_sub$VALUE[j]+1], lwd = 0.2)
    }
    # 画导入片段
    fragment_sub <- fragment[fragment$CHROM==chr$CHROM[i], ]
    if (nrow(fragment_sub)) {
      for (j in 1:nrow(fragment_sub)) {
        rect(xleft = fragment_sub$START[j], ybottom = -i*3+2, xright = fragment_sub$END[j], ytop = -i*3+3, col = fragment_sub$COLOR[j], border = NA)
      }
    }
  }else{
    # 画标记
    marker_sub <- marker[marker$CHROM==chr$CHROM[i], ]
    for (j in 1:nrow(marker_sub)) {
      lines(x = c(marker_sub$POS[j], c(marker_sub$POS[j])), y = c(-3*i+1, -3*i+3), col = color[marker_sub$VALUE[j]+1], lwd = 0.2)
    }
  }
}
# 绘制图例
colkey(col=color, clim=range(marker$VALUE), clab = "", add=TRUE, length = 0.75, side = 4, tick = 0.01, width = 0.75)
dev.off()

