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
obs_lung = data_lung[which(species_totals_lung != 0)]
obs_rabund_lung = rabund(obs_lung)
obs_inter_lung = inter_mat_1trial(obs_lung, method = method)
reps_lung = dim(obs_lung)[1]
sample_sizes_lung = apply(obs_lung, 1, sum)
sims_lung = sim_neut(1000, reps_lung,
                     sample_sizes_lung, obs_rabund_lung)
sim_inter_lung = inter_mat(sims_lung, method = method)
comp_lung_to_sim = percentile(obs_inter_lung, sim_inter_lung)

# Generate Oral Interaction Scores
obs_oral = data_oral[which(species_totals_oral != 0)]
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

length(which(signif_lung & signif_oral))
length(which(signif_lung & ! signif_oral))
length(which(! signif_lung & signif_oral))
length(which(neg_signif_lung & pos_signif_oral))
length(which(pos_signif_oral & neg_signif_lung))
length(which(neg_signif_lung & ! signif_oral))
length(which(pos_signif_lung & ! signif_oral))
length(which(neg_signif_oral & ! signif_lung))
length(which(pos_signif_oral & ! signif_lung))

# What would a correlation heatmap look like if we are actually
# getting some replicates from "different" communities?



