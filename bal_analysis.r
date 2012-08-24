# Interaction analysis of ling microbiome from non-smokers.
# Values in data set are relative abundance of the given OTU from
# a sample size of 1000.
source("metintome.r")

data = read.csv("bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000

totals = apply(data, 2, sum)
data_all = data[, which(totals != 0)]

# 
# inter_analysis = function(obs_rxs, cutoff = 0.01){
#   reps = dim(obs_rxs)[1]
#   sample_sizes = apply(obs_rxs, 1, sum)
#   
#   obs_rel_abunds = rel_abund(obs_rxs)
#   neutral_sim_results = sim_neut(100, reps = reps,
#                                  sample_size = sample_sizes,
#                                  rel_abunds = obs_rel_abunds)
#   for (method in c("spearman", "kendall", "covar")){
#     obs_inter = inter_mat_1trial(obs, method = method)
#     sim_inters = inter_mat(neutral_sim_results, method = method)
#     perc_mat = percentile(obs_inter, sim_inters)
#     percentile_heatmap(perc_mat, cutoff = cutoff)
#     perc_mat = percentile(sim_inter, sim_inters)
#     percentile_heatmap(perc_mat, cutoff = cutoff)
#   }
# }

# obs_rxs = data_top20
obs_rxs = data_all
method = "spearman"
reps = dim(obs_rxs)[1]
sample_sizes = apply(obs_rxs, 1, sum)
obs_rel_abunds = rel_abund(obs_rxs)
neut_rxs_sample = sim_neut(1000, reps, sample_sizes, obs_rel_abunds)
obs_inter_mat = inter_mat_1trial(obs_rxs, method = method)
neut_inter_mat_sample = inter_mat(neut_rxs_sample, method = method)
# why are all of the entires positive?
comparison = percentile(obs_inter_mat, neut_inter_mat_sample)
percentile_heatmap(comparison[1:50, 1:50], cutoff = 0.001)