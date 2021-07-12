#!/bin/env Rscript
#cut -f1,2,4,5,6,42,13,73 Bnabg.filter.SNPs.txt | grep -vP '^scaffold\d' | grep -v '\./\.' | grep -v '0/0.*0/0.*0/0' | grep -v '1/1.*1/1.*1/1' | grep -v '0/1.*0/1.*0/1' > QS_Q10.txt &
#Rscript bg.R -v QS_Q01.txt -r s4270 -d QS_T01 -s QS_Q01 &
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Bna background evaluate")

# Add command line arguments
p <- add_argument(p, "--variation", help = "variation info file", type = "character")
p <- add_argument(p, "--recurrent", help = "recurrent parent name", type = "character")
p <- add_argument(p, "--donor", help = "donor parent name", type = "character")
p <- add_argument(p, "--sample", help = "sample name", type = "character")
p <- add_argument(p, "--minQ", help = "the min mapQ", type = "numeric", default = 300)

#p <- add_argument(p, "--bam", help="input: bam file", type="character")
#p <- add_argument(p, "--gtf", help="input: gtf file", type="character")
#p <- add_argument(p, "--output", help="output prefix", type="character")
#p <- add_argument(p, "--nthread", help="thread number", type="integer", default="1")
#p <- add_argument(p, "--attrType", help="calculate expression for gene_id or transcript_id, default: gene_id", type="character", default="gene_id")
#p <- add_argument(p, "--strandSpecific", help="if strand-specific read counting should be performed, it should be one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded), default: 0", type="integer", default="0")
#p <- add_argument(p, "--featureType", help="the feature type used to select rows in the GTF annotation which will be used for read summarization", type="character", default="exon")

# Parse the command line arguments
argv <- parse_args(p)


library(tidyverse)
library(ggsci)
library(cowplot)

filename <- argv$variation
recurrentP <- argv$recurrent
donorP <- argv$donor
Sample <- argv$sample
minQ <- argv$minQ
#filename <- "./FCC_F11.txt"
#recurrentP <- "s8804"
#donorP <- "s6672"
#Sample <- "FCC_F11"
#minQ <- 300

df <- read_tsv(file = filename) %>%
  select(CHROM, POS, REF, ALT, QUAL, "recurrentP" = all_of(recurrentP), "donorP" = all_of(donorP), "Sample" = all_of(Sample)) %>%
  mutate(CHROM = str_remove(CHROM, "scaffold"))


dd<-df %>% filter(QUAL > minQ)

l1 <- nchar(dd$REF)	#求REF字符串长度
l2 <- nchar(dd$ALT)	#求ALT字符串长度

#过滤indel行，保留SNP信息
lse <- l1 == 1 & l2 == 1
dd <- dd[lse, ]
colnames(dd)
dd1 <- dd %>%
  separate(recurrentP, c("recurrentP.geno", "recurrentP", "recurrentP.depth"), sep = ":", convert = T) %>%
  separate(donorP, c("donorP.geno", "donorP", "donorP.depth"), sep = ":", convert = T) %>%
  separate(Sample, c("Sample.geno", "Sample", "Sample.depth"), sep = ":", convert = T)
pdf(file = paste(Sample, "depth.pdf", sep = "_"))
par(mfrow = c(3, 1))
hist(dd1 %>% pull(recurrentP.depth), breaks = 1000, xlim = c(0, 50), xlab = NULL, 
     main = paste(recurrentP, "depth frequency", sep = " "))
hist(dd1 %>% pull(donorP.depth), breaks = 1000, xlim = c(0, 50), xlab = NULL,
     main = paste(donorP, "depth frequency", sep = " "))
hist(dd1 %>% pull(Sample.depth), breaks = 1000, xlim = c(0, 50), xlab = NULL,
     main = paste(Sample, "depth frequency", sep = " "))
plot(density(dd1 %>% pull(recurrentP.depth), width = 1.5), 
     xlim = c(0, 50), main = paste(recurrentP, "depth density", sep = " "))
plot(density(dd1 %>% pull(donorP.depth), width = 1.5), 
     xlim = c(0, 50), main = paste(donorP, "depth density", sep = " "))
plot(density(dd1 %>% pull(Sample.depth), width = 1.5), 
     xlim = c(0, 50), main = paste(Sample, "depth density", sep = " "))
maxRD <- 100
minRD <- 5
maxDD <- 100
minDD <- 4
maxSD <- 50
minSD <- 4
dev.off()

dd2 <- dd1 %>%
  filter(recurrentP.depth > minRD, recurrentP.depth < maxRD,
         donorP.depth > minDD, donorP.depth < maxDD,
         Sample.depth > minSD, Sample.depth < maxSD) %>%
  filter(recurrentP.geno == "0/0" | recurrentP.geno == "1/1") %>%
  filter(donorP.geno == "0/0" | donorP.geno == "1/1") %>%
  filter(Sample.geno == "0/0" | Sample.geno == "1/1") %>%
  filter(recurrentP.geno != donorP.geno)

dd2 <- dd2 %>% mutate(genotype = if_else(Sample.geno == "0/1", 1,
                                  if_else(Sample.geno == recurrentP.geno, 0, 2)))

chrs <- unique(dd2$CHROM)
bg <- tibble(chr = NULL, start = NULL, end = NULL, sum = NULL)
for (chr in chrs) {
  df <- filter(dd2, CHROM == chr)
  g1 <- df %>% pull(genotype)
  no.na.id <- which(!is.na(g1))
  wind.sum <- c(); wind.start <- c(); wind.end <- c(); wind.chr <- c()
  for (i in 1:(length(no.na.id) - 14)) {
    wind.sum[i] <- sum(g1[no.na.id[i]:no.na.id[i+14]], na.rm = T)
    wind.start[i] <- df$POS[no.na.id[i]]
    wind.end[i] <- df$POS[no.na.id[i+14]]
    wind.chr[i] <- chr
  }
  bg_tmp <- tibble(chr = wind.chr, start = wind.start, end = wind.end, sum = wind.sum)
  bg <- rbind(bg, bg_tmp)
}

bg <- bg %>% mutate(mid = (start + end) / 2)

write_tsv(x = bg, path = paste(Sample, "bg_estimate.txt", sep = "."))

p <- ggplot(bg, aes(x = mid, y = sum)) +
  geom_point(aes(color = sum)) +
  scale_x_continuous(breaks = c(0, 10000000, 20000000, 30000000, 40000000, 
                                50000000, 60000000, 70000000, 80000000),
                     labels = c("0M", "10M", "20M", "30M", "40M", 
                                "50M", "60M", "70M", "80M")) +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = NULL, y = NULL, title = paste("Sample: ", Sample, ";  Recurrent Parent: ", recurrentP, ";  Donor Parent: ", donorP, sep = "")) +
  theme_half_open() +
  theme(legend.position = "none") +
  facet_grid(chr ~ .)
ggsave(p, filename = paste(Sample, "png", sep = "."), height = 15, width = 10, dpi = 500)

#list.files()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#dd1 <- dd %>% 
#  separate(recurrentP, c(paste(recurrentP, "geno", sep = "."), recurrentP, paste(recurrentP, "depth", sep = ".")), sep = ":", convert = T) %>%
#  separate(donorP, c(paste(donorP, "geno", sep = "."), donorP, paste(donorP, "depth", sep = ".")), sep = ":", convert = T) %>%
#  separate(Sample, c(paste(Sample, "geno", sep = "."), Sample, paste(Sample, "depth", sep = ".")), sep = ":", convert = T)
#par(mfrow = c(3, 1))
#hist(dd1[, paste(recurrentP, "depth", sep = ".")][[1]], breaks = 1000, xlim = c(0, 50),
#     xlab = NULL, main = paste(recurrentP, "depth frequency", sep = " "))
#hist(dd1[, paste(donorP, "depth", sep = ".")][[1]], breaks = 1000, xlim = c(0, 50),
#     xlab = NULL, main = paste(donorP, "depth frequency", sep = " "))
#hist(dd1[, paste(Sample, "depth", sep = ".")][[1]], breaks = 1000, xlim = c(0, 50),
#     xlab = NULL, main = paste(Sample, "depth frequency", sep = " "))
#
#
#plot(density(dd1[, paste(recurrentP, "depth", sep = ".")][[1]], width = 1.5), 
#     xlim = c(0, 50), main = paste(recurrentP, "depth density", sep = " "))
#plot(density(dd1[, paste(donorP, "depth", sep = ".")][[1]], width = 1.5), 
#     xlim = c(0, 50), main = paste(donorP, "depth density", sep = " "))
#plot(density(dd1[, paste(Sample, "depth", sep = ".")][[1]], width = 1.5), 
#     xlim = c(0, 50), main = paste(Sample, "depth density", sep = " "))
#dev.off()
#
#
#dd2 <- dd1 %>% 
#  filter(dd1[, paste(recurrentP, "depth", sep = ".")] > minRD &
#           dd1[, paste(recurrentP, "depth", sep = ".")] < maxRD &
#           dd1[, paste(donorP, "depth", sep = ".")] > minDD &
#           dd1[, paste(donorP, "depth", sep = ".")] < maxDD &
#           dd1[, paste(Sample, "depth", sep = ".")] > minSD &
#           dd1[, paste(Sample, "depth", sep = ".")] < maxSD)
#dd2 <- dd2 %>% filter(dd2[, paste(recurrentP, "geno", sep = ".")] == "0/0" | dd2[, paste(recurrentP, "geno", sep = ".")] == "1/1")
#dd2 <- dd2 %>% filter(dd2[, paste(donorP, "geno", sep = ".")] == "0/0" | dd2[, paste(donorP, "geno", sep = ".")] == "1/1")
#dd2 <- dd2 %>% filter(dd2[, paste(recurrentP, "geno", sep = ".")] != dd2[, paste(donorP, "geno", sep = ".")])
#
#
#dd2 %>% mutate(genotype = if_else(dd2[]))
#
#
#
#
#hh <- hh %>% mutate(genotype = if_else(SYX_A9.geno == "0/1", 1,
#                                 if_else(SYX_A9.geno == Y633.geno, 0, 2)))
#
#
#sample <- "genotype"
#chrs <- unique(hh$CHROM)
#bg <- data.frame(chr = NULL, start = NULL, end = NULL, sum = NULL)
#for (chr in chrs) {
#  df <- filter(hh, CHROM == chr)
#  g1 <- df[, sample][[1]]
#  no.na.id <- which(!is.na(g1))
#  wind.sum <- c(); wind.start <- c(); wind.end <- c(); wind.chr <- c()
#  for (i in 1:(length(no.na.id) - 14)) {
#    wind.sum[i] <- sum(g1[no.na.id[i]:no.na.id[i+14]], na.rm = T)
#    wind.start[i] <- df$POS[no.na.id[i]]
#    wind.end[i] <- df$POS[no.na.id[i+14]]
#    wind.chr[i] <- chr
#  }
#  bg_tmp <- data.frame(chr = wind.chr, start = wind.start, end = wind.end, sum = wind.sum)
#  bg <- rbind(bg, bg_tmp)
#}
#
#
#bg <- bg %>% mutate(mid = start + end)
#
#
#
#
#
#
#p <- ggplot(bg, aes(x = mid, y = sum)) +
#  geom_point(aes(color = sum)) +
#  scale_color_gradient(low = "blue", high = "red") +
#  theme_half_open() +
#  facet_grid(chr ~ .)
#ggsave(p, filename = "bg.png", height = 15, width = 10, dpi = 500)
