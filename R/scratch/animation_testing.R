a = readRDS("/Users/matthewcheng/Desktop/PostDoc/Spatial Assessments and Sablefish/sablefish_cie_sims_2026/outputs/mse_results/spatial_noblock_scenarios/three_region_test.RDS")
library(ggplot2)
library(gganimate)
anim_df <- do.call(rbind, lapply(65:95, function(yr) {
  mod_yr <- yr - 65 + 1

  # SSB
  ssb_em <- colSums(a$models[[mod_yr]][[1]]$rep$SSB)
  ssb_om <- colSums(a$SSB[, 1:yr, 1])
  ssb_re <- (ssb_em[1:yr] - ssb_om) / ssb_om

  # Rec
  rec_em <- colSums(a$models[[mod_yr]][[1]]$rep$Rec)
  rec_om <- colSums(a$Rec[, 1:yr, 1])
  rec_re <- (rec_em[1:yr] - rec_om) / rec_om

  # totoal boimass
  totbiom_em <- colSums(a$models[[mod_yr]][[1]]$rep$Total_Biom)
  totbiom_om <- colSums(a$Total_Biom[, 1:yr, 1])
  totbiom_re <- (totbiom_em[1:yr] - totbiom_om) / totbiom_om

  # q
  q_em <- as.vector(a$models[[mod_yr]][[1]]$rep$srv_q[,,1])[1:yr]
  q_om <- as.vector(a$srv_q[, 1:yr, 1, 1])
  q_re <- (q_em - q_om) / q_om

  rbind(
    data.frame(year = 1:length(ssb_em), value = ssb_em, type = 'EM', panel = 'SSB',           frame_yr = yr),
    data.frame(year = 1:yr,             value = ssb_om,  type = 'OM', panel = 'SSB',           frame_yr = yr),
    data.frame(year = 1:yr,             value = ssb_re,  type = 'RE', panel = 'SSB Rel. Error', frame_yr = yr),
    data.frame(year = 1:length(totbiom_em), value = totbiom_em, type = 'EM', panel = 'Tot Biom',           frame_yr = yr),
    data.frame(year = 1:yr,             value = totbiom_om,  type = 'OM', panel = 'Tot Biom',           frame_yr = yr),
    data.frame(year = 1:yr,             value = totbiom_re,  type = 'RE', panel = 'Tot Biom Rel. Error', frame_yr = yr),
    data.frame(year = 1:length(rec_em), value = rec_em, type = 'EM', panel = 'Rec',           frame_yr = yr),
    data.frame(year = 1:yr,             value = rec_om,  type = 'OM', panel = 'Rec',           frame_yr = yr),
    data.frame(year = 1:yr,             value = rec_re,  type = 'RE', panel = 'Rec Rel. Error', frame_yr = yr),
    data.frame(year = 1:yr,             value = q_em,    type = 'EM', panel = 'q',             frame_yr = yr),
    data.frame(year = 1:yr,             value = q_om,    type = 'OM', panel = 'q',             frame_yr = yr),
    data.frame(year = 1:yr,             value = q_re,    type = 'RE', panel = 'q Rel. Error',  frame_yr = yr)
  )
}))

p <- ggplot(anim_df, aes(x = year, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_hline(data = data.frame(panel = c('SSB Rel. Error', 'q Rel. Error'), yint = 0),
             aes(yintercept = yint), linetype = 'dashed', color = 'black') +
  scale_color_manual(values = c('EM' = 'red', 'OM' = 'black', 'RE' = 'blue')) +
  facet_wrap(~panel, ncol = 1, scales = 'free_y') +
  labs(x = 'Year', y = NULL, title = 'Terminal year: {closest_state}') +
  theme_bw() +
  transition_states(frame_yr, transition_length = 1, state_length = 2) +
  ease_aes('linear')

animate(p, nframes = 100, fps = 10, width = 800, height = 900)
anim_save('re_animation.gif')
