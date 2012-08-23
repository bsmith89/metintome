# Interaction analysis of ling microbiome from non-smokers.
# Values in data set are relative abundance of the given OTU from
# a sample size of 1000.
source("metintome.r")

data = read.csv("bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000

top20 = sort(rel_abunds(data), decreasing = T)[1:20]
sum(top20)

data_top20only = data[, names(top20)]

rel_abunds(data_top20only)
sum(rel_abunds(data_top20only))

totals = apply(data_top20only, c(1), sum)

results = analyze(data_top20only, n = 1000, assign_intra = F)
percentile_heatmap(results$percentile$spearman, cutoff = 0.001)

percentile_scores$sim_data$cov[1,1,]
