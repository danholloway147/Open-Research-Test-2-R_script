##############################################################################
# Halstock Wood SSSI Habitat Survey Report - 3D NMDS, envfit, PERMANOVA, betadisper, etc...
# Purpose: To Explore Understorey Plant Community Composition Along a Riparian Distance Gradient in Halstock Wood SSSI
# Author: Dan A. Holloway
# Contact: Dan2.holloway@live.uwe.ac.uk
# Date: 10/08/2026
# License: MIT License
##############################################################################

install.packages(c("vegan", "ggplot2", "ggrepel", "vegan3d", "patchwork"))
library(vegan)
library(ggplot2)
library(ggrepel)
library(vegan3d)
library(patchwork)

##############################################################################
# 1. Load and prepare data
##############################################################################

dframe1 <- read.csv("Species_Data.csv")
species_data <- dframe1[, !(names(dframe1) %in% c("plot_id", "transect", "quadrat"))]

dframe2 <- read.csv("Env_Factors.csv", stringsAsFactors = TRUE)

dframe2$distance_group <- factor(dframe2$distance_from_river_m,
                                 levels = sort(unique(dframe2$distance_from_river_m)))

##############################################################################
# 2. Run 3D NMDS
##############################################################################

set.seed(1) 
model1 <- metaMDS(species_data, k = 3, trymax = 100, autotransform = TRUE)
model1
stress_val <- round(model1$stress, 3)
stressplot(model1)

##############################################################################
# 2b. Justify dimensionality - stress vs k scree check
##############################################################################

stress_by_k <- sapply(2:5, function(k) {
  m <- metaMDS(species_data, k = k, trymax = 100, autotransform = TRUE, trace = FALSE)
  m$stress
})
scree_df <- data.frame(k = 2:5, stress = stress_by_k)
print(scree_df)

fig_scree <- ggplot(scree_df, aes(x = k, y = stress)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "#08519c") +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "red") +
  annotate("text", x = 2.2, y = 0.21, label = "stress = 0.2 threshold", colour = "red", size = 3.2, hjust = 0) +
  scale_x_continuous(breaks = 2:5) +
  labs(title = "NMDS stress by number of dimensions",
       x = "Number of dimensions (k)", y = "Stress") +
  theme_bw(base_size = 13)

print(fig_scree)
ggsave("Figure0_stress_scree.png", fig_scree, width = 6, height = 5, dpi = 400)

# Stress diagnostic (Shepard plot) 

# 1. View on screen first 
stressplot(model1)
title(main = paste0("Shepard plot (stress = ", stress_val, ")"))

# 2. Save to file 
png("stressplot.png", width = 1600, height = 1600, res = 300)
result <- try({
  stressplot(model1)
  title(main = paste0("Shepard plot (stress = ", stress_val, ")"))
})
dev.off()
if (inherits(result, "try-error")) {
  warning("stressplot() failed - check model1 was created successfully before this step")
}

##############################################################################
# 3. Pull NMDS scores into data frames for ggplot
##############################################################################

site_scores <- as.data.frame(scores(model1, display = "sites"))
site_scores$plot_id <- dframe1$plot_id
site_scores <- merge(site_scores, dframe2, by = "plot_id")

species_scores <- as.data.frame(scores(model1, display = "species"))
species_scores$species <- rownames(species_scores)

##############################################################################
# 4. Fit environmental variables (envfit)
##############################################################################

habmodel <- envfit(model1 ~ distance_from_river_m + soil_moisture_mg_g +
                     soil_ph + canopy_cover_pct,
                   data = dframe2, perm = 999, choices = 1:3)
habmodel

# Extract vectors for ggplot
env_vectors <- as.data.frame(scores(habmodel, display = "vectors"))
env_vectors$var <- rownames(env_vectors)
env_vectors$r2 <- habmodel$vectors$r
env_vectors$pval <- habmodel$vectors$pvals

# Nice looking labels :)
env_vectors$label <- c(
  distance_from_river_m = "Distance from river",
  soil_moisture_mg_g     = "Soil moisture",
  soil_ph                = "Soil pH",
  canopy_cover_pct        = "Canopy cover"
)[env_vectors$var]

arrow_mult <- 2 # scales arrows  - adjust if its  needed

##############################################################################
# 5. ggplot theme
##############################################################################

theme_nmds <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey30", size = 11)
  )

distance_pal <- c("#08519c", "#3182bd", "#6baed6", "#fd8d3c", "#e6550d")

##############################################################################
# 6. Figure 1 - NMDS Axes 1 & 2
##############################################################################

hulls_12 <- do.call(rbind, lapply(split(site_scores, site_scores$distance_group), function(d) {
  d[chull(d$NMDS1, d$NMDS2), ]
}))

fig1 <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2, colour = distance_group, fill = distance_group)) +
  geom_polygon(data = hulls_12, alpha = 0.15, colour = NA) +
  geom_point(size = 3, shape = 21, colour = "black", stroke = 0.3) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1 * arrow_mult, yend = NMDS2 * arrow_mult),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey20", linewidth = 0.6) +
  geom_text_repel(data = env_vectors,
                  aes(x = NMDS1 * arrow_mult, y = NMDS2 * arrow_mult, label = label),
                  inherit.aes = FALSE, colour = "grey10", fontface = "bold", size = 3.8) +
  scale_colour_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  scale_fill_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  labs(title = "NMDS ordination of understorey plant community composition",
       subtitle = paste0("Axes 1-2 | 3D solution, stress = ", stress_val),
       x = "NMDS Axis 1", y = "NMDS Axis 2") +
  coord_equal() +
  theme_nmds

print(fig1)
ggsave("Figure1_NMDS_axis1_2.png", fig1, width = 8, height = 6.5, dpi = 400)

##############################################################################
# 7. Figure 2 - NMDS Axes 1 & 3
##############################################################################

hulls_13 <- do.call(rbind, lapply(split(site_scores, site_scores$distance_group), function(d) {
  d[chull(d$NMDS1, d$NMDS3), ]
}))

fig2 <- ggplot(site_scores, aes(x = NMDS1, y = NMDS3, colour = distance_group, fill = distance_group)) +
  geom_polygon(data = hulls_13, alpha = 0.15, colour = NA) +
  geom_point(size = 3, shape = 21, colour = "black", stroke = 0.3) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1 * arrow_mult, yend = NMDS3 * arrow_mult),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey20", linewidth = 0.6) +
  geom_text_repel(data = env_vectors,
                  aes(x = NMDS1 * arrow_mult, y = NMDS3 * arrow_mult, label = label),
                  inherit.aes = FALSE, colour = "grey10", fontface = "bold", size = 3.8) +
  scale_colour_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  scale_fill_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  labs(title = "NMDS ordination of understorey plant community composition",
       subtitle = paste0("Axes 1-3 | 3D solution, stress = ", stress_val),
       x = "NMDS Axis 1", y = "NMDS Axis 3") +
  coord_equal() +
  theme_nmds

print(fig2)
ggsave("Figure2_NMDS_axis1_3.png", fig2, width = 8, height = 6.5, dpi = 400)

##############################################################################
# 8. Figure 3 - NMDS Axes 2 & 3
##############################################################################

hulls_23 <- do.call(rbind, lapply(split(site_scores, site_scores$distance_group), function(d) {
  d[chull(d$NMDS2, d$NMDS3), ]
}))

fig3 <- ggplot(site_scores, aes(x = NMDS2, y = NMDS3, colour = distance_group, fill = distance_group)) +
  geom_polygon(data = hulls_23, alpha = 0.15, colour = NA) +
  geom_point(size = 3, shape = 21, colour = "black", stroke = 0.3) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS2 * arrow_mult, yend = NMDS3 * arrow_mult),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey20", linewidth = 0.6) +
  geom_text_repel(data = env_vectors,
                  aes(x = NMDS2 * arrow_mult, y = NMDS3 * arrow_mult, label = label),
                  inherit.aes = FALSE, colour = "grey10", fontface = "bold", size = 3.8) +
  scale_colour_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  scale_fill_manual(values = distance_pal, name = "Distance from\nriver (m)") +
  labs(title = "NMDS ordination of understorey plant community composition",
       subtitle = paste0("Axes 2-3 | 3D solution, stress = ", stress_val),
       x = "NMDS Axis 2", y = "NMDS Axis 3") +
  coord_equal() +
  theme_nmds

print(fig3)
ggsave("Figure3_NMDS_axis2_3.png", fig3, width = 8, height = 6.5, dpi = 400)

##############################################################################
# 9. Combined 3-D figure - looks cool, can't use it tho cos word docs don't allow it
##############################################################################

fig_combined <- (fig1 + theme(legend.position = "none")) +
  (fig2 + theme(legend.position = "none")) +
  fig3 +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(title = "3D NMDS ordination - all axis combinations") &
  theme(legend.position = "right")

ggsave("Figure_combined_all_axes.png", fig_combined, width = 16, height = 6, dpi = 400)

##############################################################################
# 10. Interactive/rotatable 3D plot - looks cool, but not useale again...
##############################################################################

ord3d <- ordiplot3d(model1, display = "sites", pch = 16,
                    col = distance_pal[as.numeric(site_scores$distance_group)],
                    angle = 40, main = paste0("3D NMDS (stress = ", stress_val, ")"))
text(ord3d$xyz.convert(scores(model1, display = "sites")),
     labels = site_scores$plot_id, cex = 0.6, col = "grey20", pos = 3)
legend("topleft", legend = levels(site_scores$distance_group),
       col = distance_pal, pch = 16, title = "Distance (m)", bty = "n")

install.packages("plotly")
library(plotly)

fig_interactive <- plot_ly(
  data = site_scores,
  x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
  color = ~distance_group, colors = distance_pal,
  text = ~plot_id,
  hoverinfo = "text",
  type = "scatter3d", mode = "markers",
  marker = list(size = 5)
) %>%
  layout(title = paste0("3D NMDS (stress = ", stress_val, ") - hover to see plot ID"),
         scene = list(xaxis = list(title = "NMDS1"),
                      yaxis = list(title = "NMDS2"),
                      zaxis = list(title = "NMDS3")))

fig_interactive 

# Save as HTML file 
install.packages("htmlwidgets")
htmlwidgets::saveWidget(fig_interactive, "NMDS_3D_interactive.html")


##############################################################################
# 11. PERMANOVA 
##############################################################################

permanova1 <- adonis2(species_data ~ distance_group, data = dframe2,
                      permutations = 999, method = "bray")
permanova1

##############################################################################
# 12. Betadisper 
##############################################################################

dist_matrix <- vegdist(species_data, method = "bray")
disp1 <- betadisper(dist_matrix, dframe2$distance_group)
disp1
anova(disp1)
permutest(disp1, permutations = 999)

#  boxplot
disp_df <- data.frame(distance_group = dframe2$distance_group,
                      distance_to_centroid = disp1$distances)

fig_disp <- ggplot(disp_df, aes(x = distance_group, y = distance_to_centroid, fill = distance_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.6) +
  scale_fill_manual(values = distance_pal, guide = "none") +
  labs(title = "Multivariate dispersion by distance group",
       subtitle = "Distance to group centroid (Bray-Curtis dissimilarity)",
       x = "Distance from river (m)", y = "Distance to centroid") +
  theme_nmds

print(fig_disp)
ggsave("Figure4_betadisper_boxplot.png", fig_disp, width = 7, height = 5.5, dpi = 400)

##############################################################################
# 13. Ordisurf 
##############################################################################

png("ordisurf_distance.png", width = 2000, height = 1600, res = 300)
plot(model1, choices = c(1,2), dis = "site", type = "n",
     main = "Ordisurf: distance from river fitted to NMDS axes 1-2")
ordisurf(model1, dframe2$distance_from_river_m, choices = c(1,2), add = TRUE, col = "forestgreen")
points(model1, display = "sites", choices = c(1,2), pch = 21,
       bg = distance_pal[as.numeric(site_scores$distance_group)], cex = 1.2)
legend("topright", legend = levels(site_scores$distance_group),
       pt.bg = distance_pal, pch = 21, title = "Distance (m)", bty = "n")
dev.off()

##############################################################################
# 14. Summary table 
##############################################################################

envfit_summary <- data.frame(
  Variable = env_vectors$label,
  r2 = round(env_vectors$r2, 3),
  p_value = env_vectors$pval
)
print(envfit_summary)

cat("\n--- Key results ---\n")
cat("NMDS stress (k=3):", stress_val, "\n")
cat("PERMANOVA (distance_group): R2 =", round(permanova1$R2[1], 3),
    ", F =", round(permanova1$F[1], 3),
    ", p =", permanova1$`Pr(>F)`[1], "\n")
cat("Betadisper ANOVA p-value:", round(anova(disp1)$`Pr(>F)`[1], 3), "\n")
anova(disp1)$`Pr(>F)`[1]
permutest(disp1, permutations = 999)

TukeyHSD(disp1)
