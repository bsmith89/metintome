source("metintome.r")

num_top_species = 50
method = "spearman"

data_lung = read.csv("data/lung.csv", header = T, sep = ",")
data_lung = data_lung[, -c(1:2)]
data_lung = data_lung * 1000
species_totals_lung = apply(data_lung, 2, sum)
common_species_lung = names(sort(species_totals_lung,
                                 decreasing = T))[1:num_top_species]
data_oral = read.csv("data/oral.csv", header = T, sep = ",")
data_oral = data_oral[, -c(1:2)]
data_oral = data_oral * 1000
species_totals_oral = apply(data_oral, 2, sum)
common_species_oral = names(sort(species_totals_oral,
                                 decreasing = T))[1:num_top_species]

shared_common_species = intersect(common_species_lung,
                                  common_species_oral)
both_common_species = union(common_species_lung, common_species_oral)

# Generate Lung Interaction Scores
obs_lung = data_lung[, which(species_totals_lung != 0)]
obs_rabund_lung = rabund(obs_lung)
obs_inter_lung = inter_mat_1trial(obs_lung, method = method)
reps_lung = dim(obs_lung)[1]
sample_sizes_lung = apply(obs_lung, 1, sum)
sims_lung = sim_neut(1000, reps_lung,
                     sample_sizes_lung, obs_rabund_lung)
sim_inter_lung = inter_mat(sims_lung, method = method)
comp_lung_to_sim = percentile(obs_inter_lung, sim_inter_lung)

# Generate Oral Interaction Scores
obs_oral = data_oral[, which(species_totals_oral != 0)]
obs_rabund_oral = rabund(obs_oral)
obs_inter_oral = inter_mat_1trial(obs_oral, method = method)
reps_oral = dim(obs_oral)[1]
sample_sizes_oral = apply(obs_oral, 1, sum)
sims_oral = sim_neut(1000, reps_oral,
                     sample_sizes_oral, obs_rabund_oral)
sim_inter_oral = inter_mat(sims_oral, method = method)
comp_oral_to_sim = percentile(obs_inter_oral, sim_inter_oral)

# Correlation percentile heatmaps for common species in Lung and Oral
percentile_heatmap(comp_lung_to_sim[common_species_lung,
                                    common_species_lung],
                   cutoff1 = 0.005, cutoff2 = 0.001,
                   main = "Lung, lung common")

percentile_heatmap(comp_oral_to_sim[common_species_oral,
                                    common_species_oral],
                   cutoff1 = 0.005, cutoff2 = 0.001,
                   main = "Oral, oral common")

# Correlation percentile heatmaps for species commonly found in either
# environment

percentile_heatmap(comp_lung_to_sim[shared_common_species,
                                    shared_common_species],
                   cutoff1 = 0.005, cutoff2 = 0.001, cluster = F, 
                   main = "Lung, shared common")

percentile_heatmap(comp_oral_to_sim[shared_common_species,
                                    shared_common_species],
                   cutoff1 = 0.005, cutoff2 = 0.001, cluster = F, 
                   main = "Oral, shared common")

# What fraction of the interactions which are significant are shared?
# Change sign?

cutoff = 0.005
comp_lung_shared = lower_left(comp_lung_to_sim[shared_common_species, shared_common_species])
comp_oral_shared = lower_left(comp_oral_to_sim[shared_common_species, shared_common_species])
pos_signif_lung = comp_lung_shared > 1 - cutoff
pos_signif_oral = comp_oral_shared > 1 - cutoff
neg_signif_lung = comp_lung_shared < cutoff
neg_signif_oral = comp_oral_shared < cutoff
signif_lung = pos_signif_lung | neg_signif_lung
signif_oral = pos_signif_oral | neg_signif_oral

(((length(shared_common_species) ** 2) - 
  length(shared_common_species))) / 2
# Number of possible interactions: 741 for 39 OTUs

length(which(pos_signif_lung))
length(which(pos_signif_oral))
length(which(neg_signif_lung))
length(which(neg_signif_oral))
length(which(signif_lung & signif_oral))
length(which(signif_lung | signif_oral))
length(which(signif_lung & ! signif_oral))
length(which(signif_oral & ! signif_lung))
length(which(neg_signif_lung & pos_signif_oral))
length(which(pos_signif_oral & neg_signif_lung))
length(which(neg_signif_lung & ! signif_oral))
length(which(pos_signif_lung & ! signif_oral))
length(which(neg_signif_oral & ! signif_lung))
length(which(pos_signif_oral & ! signif_lung))
# We don't find any swapping between postive and negative correlations
# in either direction (Good, 'cause how would we explain it?). We do 
# find a variety of switches between significant and insignificant in 
# both directions for both environments.  There is clearly, however, a
# bias towards significance and more negative significance in the
# lung.


# What would a correlation heatmap look like if we are actually
# getting some replicates from "different" communities?

pareto_location = 1
pareto_shape = 1
num_species = 500
abund1 = rpareto(num_species, pareto_location, pareto_shape)
abund2 = rpareto(num_species, pareto_location, pareto_shape)
abund3 = rpareto(num_species, pareto_location, pareto_shape)
abund4 = rpareto(num_species, pareto_location, pareto_shape)
abund1 = namespecies(abund1)
abund2 = namespecies(abund2)
abund3 = namespecies(abund3)
abund4 = namespecies(abund4)
rabund1 = abund1/sum(abund1)
rabund2 = abund2/sum(abund2)
rabund3 = abund3/sum(abund3)
rabund4 = abund4/sum(abund4)
obs_reps1 = sim_neut_1trial(10, 1000, rabund1)
obs_reps2 = sim_neut_1trial(10, 1000, rabund2)
obs_reps3 = sim_neut_1trial(10, 1000, rabund3)
obs_reps4 = sim_neut_1trial(10, 1000, rabund4)
data_multi = rbind(obs_reps1, obs_reps2, obs_reps3, obs_reps4)

species_totals_multi = apply(data_multi, 2, sum)
common_species_multi = names(sort(species_totals_multi,
                                  decreasing = T))[1:num_top_species]

# Generate "Multi Environment" Interaction Scores
obs_multi = data_multi[, which(species_totals_multi != 0)]
obs_rabund_multi = rabund(obs_multi)
obs_inter_multi = inter_mat_1trial(obs_multi, method = method)
reps_multi = dim(obs_multi)[1]
sample_sizes_multi = apply(obs_multi, 1, sum)
sims_multi = sim_neut(1000, reps_multi,
                     sample_sizes_multi, obs_rabund_multi)
sim_inter_multi = inter_mat(sims_multi, method = method)
comp_multi_to_sim = percentile(obs_inter_multi, sim_inter_multi)

# Correlation percentile heatmaps for common species in "Multi
# Environment"
percentile_heatmap(comp_multi_to_sim[common_species_multi,
                                     common_species_multi],
                   cutoff1 = 0.005, cutoff2 = 0.001,
                   main = "Multi, multi common")
