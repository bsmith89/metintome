# Interaction analysis of ling microbiome from non-smokers.
# Values in data set are relative abundance of the given OTU from
# a sample size of 1000.
source("metintome.r")

data = read.csv("bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000

top20 = sort(rel_abund(data), decreasing = T)[1:20]
data_top20only = data[, names(top20)]