# Interaction analysis of ling microbiome from non-smokers.
# Values in data set are relative abundance of the given OTU from
# a sample size of 1000.
source("metintome.r")

data = read.csv("bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000

data_20 = data[,1:20]

totals = apply(data, 2, sum)
data_all = data[, which(totals != 0)]

# If I cut out all but the top 20 OTUs
obs_rxs20 = data_20
method20 = "spearman"
reps20 = dim(obs_rxs)[1]
sample_sizes20 = apply(obs_rxs20, 1, sum)
obs_rel_abunds20 = rel_abund(obs_rxs20)
neut_rxs_sample20 = sim_neut(1000, reps20, sample_sizes20, obs_rel_abunds20)
obs_inter_mat20 = inter_mat_1trial(obs_rxs20, method = method20)
neut_inter_mat_sample20 = inter_mat(neut_rxs_sample20, method = method20)
comparison20 = percentile(obs_inter_mat20, neut_inter_mat_sample20)
percentile_heatmap(comparison20[1:20, 1:20],
                   cutoff1 = 0.005, cutoff2 = 0.001)

obs_rxs = data_all
method = "spearman"
reps = dim(obs_rxs)[1]
sample_sizes = apply(obs_rxs, 1, sum)
obs_rel_abunds = rel_abund(obs_rxs)
neut_rxs_sample = sim_neut(1000, reps, sample_sizes, obs_rel_abunds)
obs_inter_mat = inter_mat_1trial(obs_rxs, method = method)
neut_inter_mat_sample = inter_mat(neut_rxs_sample, method = method)
comparison = percentile(obs_inter_mat, neut_inter_mat_sample)
percentile_heatmap(comparison[1:20, 1:20],
                   cutoff1 = 0.005, cutoff2 = 0.001)
