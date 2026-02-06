####################################################################################
####### The biodiversity we ignore: a global meta-analysis on taxonomic bias #######
####################################################################################

####################################################################################
# R script to generate the analysis
# Analysis performed with R (v. 4.4.1) and Rstudio (v. 2025.09.0+387)
# Created by Stefano Mammola and Martino Adamo
# Initiated 16/12/2025
#####################################################################################

####################################################################################
# Preparation ----------------------------------------------------------------------
####################################################################################

# Loading R packages ------------------------------------------------------

if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("dplyr", "gghalves", "ggplot2","ggh4x", "ggimage", "gridExtra", "metafor", "rnaturalearth", "rnaturalearthdata", "patchwork","stringr","readxl","tidyr")

# Function ----------------------------------------------------------------

#Function to model with moderators and extract results
extract_moderator_results <- function(data, moderator, i,
                                      min_n = 5,
                                      min_ID = 3) {
  
  # Guard: insufficient data 
  if (nrow(data) < min_n || length(unique(data$ID)) < min_ID) {
    return(
      data.frame(
        label = i,
        Moderator_group = moderator,
        Moderator = NA, 
        N_study = NA,
        N_est   = NA,
        b     = NA,
        ci.lb = NA,
        ci.ub = NA,
        ES    = NA,
        L     = NA,
        U     = NA,
        status = "skipped_small_sample"
      )
    )
  }
  
  mods_formula <- as.formula(paste0("~ 0 + ", moderator))
  
  # Try fitting the model
  model <- tryCatch(
    metafor::rma.mv(
      yi, vi,
      mods   = mods_formula,
      random = list(~1 | ID),
      data   = data,
      control = list(optimizer = "optim", optmethod = "BFGS")
    ),
    error = function(e) NULL
  )
  
  # Guard: Non-convergence / failure
  if (is.null(model)) {
    return(
      data.frame(
        label = i,
        Moderator_group = moderator,
        Moderator = NA, 
        N_study = NA,
        N_est   = NA,
        b     = NA,
        ci.lb = NA,
        ci.ub = NA,
        ES    = NA,
        L     = NA,
        U     = NA,
        status = "model_failed"
      )
    )
  }
  
  # Extract moderator levels
  levs <- sub(paste0("^", moderator), "", rownames(model$beta))
  
  # Compute number of studies and estimates per level
  N_study <- sapply(levs, function(l) length(unique(data$ID[data[[moderator]] == l])))
  N_est   <- sapply(levs, function(l) sum(data[[moderator]] == l, na.rm = TRUE))
  
  # Build output
  out <- data.frame(
    label = i,
    Moderator_group = moderator,
    Moderator = levs,
    N_study = N_study,
    N_est   = N_est,
    b     = model$b,
    ci.lb = model$ci.lb,
    ci.ub = model$ci.ub,
    status = "ok"
  ) |> dplyr::mutate_if(is.character, as.factor)
  
  out$ES <- (exp(out$b)     - 1) / (exp(out$b)     + 1)
  out$L  <- (exp(out$ci.lb) - 1) / (exp(out$ci.lb) + 1)
  out$U  <- (exp(out$ci.ub) - 1) / (exp(out$ci.ub) + 1)
  
  out
}

# Settings ----------------------------------------------------------------

# Plotting parameters
theme_set(theme_bw())
theme_update(
  legend.background = element_blank(), #No background (legend)
  plot.background = element_blank(), #No background
  panel.grid = element_blank(), #No gridlines
  axis.text  = element_text(size = 10, colour = "grey10"), #Size and color of text
  axis.title = element_text(size = 12, colour = "grey10") #Size and color of text
)

#Colors
my_colors <- c("#66a61e","#d95f02", "#7570b3", "#e7298a", "#1b9e77")

# Load data ---------------------------------------------------------------
db <- readxl::read_excel("Data/Database_meta_analysis_bias.xlsx", na = "NA") |>
  dplyr::mutate_if(is.character, as.factor) |> 
  dplyr::arrange(Pearson_r) 
str(db)

#Distinct db
db2 <- db |> 
  dplyr::distinct(ID, .keep_all = TRUE) |> 
  droplevels()

####################################################################################
# Summary stats --------------------------------------------------------------------
####################################################################################

#Number of papers
nrow(db2) #78

#Number of estimates
nrow(db) #905

#Average number of publication per study
mean(table(db$ID))
sd(table(db$ID))

#Number of journals
nlevels(db$Journal) #54

#Range of years 
range(db$Year)

table(db$Taxon_macro)/sum(table(db$Taxon_macro))

#Summary of estimates in each combination of response and predictor
(table_tot <- table(db$Dependent_zoomed,db$Independent_zoomed))

####################################################################################
# Summary Figure  1 ----------------------------------------------------------------
####################################################################################

# studies by year --------------------------------------------------------------

A <- db |> 
      dplyr::select(ID, Year) |> 
      dplyr::group_by(Year) |> 
      dplyr::distinct() |> 
      dplyr::mutate(freq = n_distinct(ID)) |>
  ggplot(aes(x = Year)) +
  geom_bar(stat = "count", fill = "#Aec6cf", color = "grey70") +
  labs(title = "A", x = NULL, y = "Number of studies") +
  scale_x_continuous(breaks = seq(min(db$Year), max(db$Year), by = 5)) +
  theme_minimal(base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_text(margin = margin(t=10), color = "grey40"),
        axis.title.y = element_text(color = "grey40"),
        axis.text.x =  element_text(color = "grey40", angle = 90, hjust = -5),
        plot.title = element_text(size = 16, color = "#1a1a1a"))

# studies by country -----------------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

geo <- read.delim("geo.tsv") # manually edited version of the file created at line 227 (check of country names match)

world_map <- world |> 
  dplyr::left_join(geo, by =  "sovereignt") |> 
  dplyr::filter(sovereignt != "Antarctica")

B <-
  ggplot(world_map) +
  geom_sf(aes(fill = count-29), color = "grey40", linewidth = 0.1) +
  scale_fill_gradient(low = "#7570b3", high = "grey10", name = "") + #"#d95f02"
  ggtitle("B") +
  theme_void() + theme(plot.title = element_text(size = 16, color = "#1a1a1a", hjust = 0.05))

# studies by taxon -------------------------------------------------------------

tax <- db |> 
  dplyr::select(ID,Taxon_macro,Taxon,Taxon_verbatim) |>
  dplyr::distinct() |>
  dplyr::mutate(Taxon = str_replace_all(Taxon, "Invertebrates", "Arthropods")) |>
  dplyr::mutate(Taxon = str_replace_all(Taxon, "Vertebrates", "Multiple")) |>
  dplyr::mutate(group = case_when(Taxon == "Multiple" ~ "Multiple", TRUE ~ "Other taxa"))

img_df <- tax |> 
  dplyr::count(Taxon) |>
  dplyr::mutate(
    image = case_when(
      Taxon == "Plants"        ~ "img/plants.png",
      Taxon == "Birds"         ~ "img/birds.png",
      Taxon == "Mammals"       ~ "img/mammal.png",
      Taxon == "Fish"          ~ "img/fish.png",
      Taxon == "Arthropods"    ~ "img/spid.png",
      Taxon == "Reptiles"      ~ "img/reptile.png",
      Taxon == "Amphibians"    ~ "img/frogs.png",
      TRUE                     ~ NA_character_)) |>
  dplyr::mutate(group = case_when(Taxon == "Multiple" ~ "Multiple", TRUE ~ "Other taxa"))

C <- 
  ggplot2::ggplot(img_df, aes(x = reorder(Taxon,-n), y = n)) +
  geom_col(aes(fill = Taxon), fill = "#Aec6cf", color = "grey70") +
  geom_image(
    data = img_df,
    aes(x = Taxon, y = n + 0.5, image = image),
    size = 0.2,
    inherit.aes = T,
    na.rm = TRUE) +
  ylim(0,35) +
  labs(title = "C", x = NULL, y = "Number of studies") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x", drop = TRUE) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10), color = "grey40"),
    axis.title.y = element_text(color = "grey40"),
    axis.text.x  = element_text(color = "grey40", angle = 90, hjust = 1),
    plot.title   = element_text(size = 16, color = "#1a1a1a"),
    legend.position = "none",
    strip.text = element_blank())

# variables --------------------------------------------------------------------

df_counts <- 
  db |> 
  dplyr::select(ID, Independent_macro) |> 
  dplyr::distinct() |>
  dplyr::count(Independent_macro, name = "count")

units <- df_counts |> 
  tidyr::uncount(count) |> 
  dplyr::mutate(id = row_number())

ncol <- 10  

units <- units  |> dplyr::mutate(x = (id - 1) %% ncol + 1, y = (id - 1) %/% ncol + 1) 

#Sort
units$Independent_macro <- factor(units$Independent_macro,
                               levels = c(levels(units$Independent_macro)[1:4],
                                          levels(units$Independent_macro)[6],
                                          levels(units$Independent_macro)[5])) #put level "Other" last


D <- 
  ggplot(units, aes(x, y, fill = Independent_macro)) +
  geom_tile(color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = c(my_colors,"grey70"), drop = FALSE) +
  coord_equal() + labs(title = "D") +
  scale_y_reverse() +
  labs(fill = "Explanatory variable") +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, color = "#1a1a1a"))

# compose the image panel ------------------------------------------------------

lay <- rbind(c(1,2,2),c(3,3,4))

panel <- grid.arrange(A, B, C, D, layout_matrix = lay)

#ggsave(plot = panel, filename = "figure1.tiff", device = "tiff", dpi = 320) #example

####################################################################################
# Data analysis -----------------------------------------------------------------
####################################################################################

# Derive Fischer's Z and its variance -------------------------------------

db$Pearson_r <- pmin(pmax(db$Pearson_r, -0.999), 0.999)
db <- metafor::escalc(measure = "ZCOR", ri = Pearson_r, ni = N, data = db)

#dropping general estimates in the category Other (too few observation)
db_subset <- db
db_subset <- db_subset[db_subset$Independent_macro != "Other",] |> droplevels()

########################@###################################################
# Publication bias --------------------------------------------------------
########################@###################################################

#Funnel plot (aggregating to one estimate per study)
dat_study <- aggregate(
  cbind(yi, vi) ~ ID,
  data = db,
  FUN = mean
)

pdf(file = "Figures/Figure_S2.pdf", width = 5, height = 4)
metafor::funnel(metafor::rma(yi, vi, data = dat_study))
dev.off()

# Egger regression
db$SE <- sqrt(db$vi)

model_egger <- rma.mv(
  yi,
  vi,
  mods = ~ SE,
  random = list(~1 | ID),
  data = db
)
summary(model_egger) #A significant SE slope --> small-study effects

# Time-lag bias
model_time <- rma.mv(
  yi,
  vi,
  mods = ~ scale(Year),
  random = list(~1 | ID),
  data = db
)
summary(model_time)

########################@###################################################
# Generating all models in a loop -----------------------------------------
########################@###################################################

# Moderators
moderators <- c(
  "Study_type",
  "Taxon",
  "Geography_macro",
  "Dependent_zoomed",
  "Independent_zoomed"
)

#Setting baselines
db_subset$Taxon_macro     <- relevel(db_subset$Taxon_macro, "Multiple") 
db_subset$Geography_macro <- relevel(db_subset$Geography_macro, "Global") 

result_for_plot <- list() 
k <- 1

for (i in levels(db_subset$Independent_macro)) {
  
  # subset
  data_i <- db_subset[db_subset$Independent_macro == i, ] |> droplevels()
  
  # Base model
  model_base <- metafor::rma.mv(
    yi, vi,
    random = list(~1 | ID),
    data = data_i,
    control = list(optimizer = "optim", optmethod = "BFGS")
  )
  
  result_base <- data.frame(
    label = i,
    Moderator_group = "Base model",
    Moderator = "Overall",
    N_study = length(unique(data_i$ID)),
    N_est   = nrow(data_i),
    b     = model_base$b,
    ci.lb = model_base$ci.lb,
    ci.ub = model_base$ci.ub,
    status = "ok"
  ) |> dplyr::mutate_if(is.character, as.factor)
  
  result_base$ES <- (exp(result_base$b)     - 1) / (exp(result_base$b)     + 1)
  result_base$L  <- (exp(result_base$ci.lb) - 1) / (exp(result_base$ci.lb) + 1)
  result_base$U  <- (exp(result_base$ci.ub) - 1) / (exp(result_base$ci.ub) + 1)
  
  # Moderator models (see custom funtion at the beginning of the script)
  results_moderators <- lapply(
    moderators,
    extract_moderator_results,
    data = data_i,
    i    = i
  )
  
  results_moderators <- do.call(rbind, results_moderators)
  rownames(results_moderators) <- NULL
  
  # Store results
  result_for_plot[[k]] <- rbind(result_base, results_moderators)
  k <- k + 1
}

# Final database
result_for_plot <- do.call(rbind, result_for_plot)
rownames(result_for_plot) <- NULL

# plot base model ---------------------------------------------------------

#subsetting the database for base plot estimates
result_tot <- result_for_plot[result_for_plot$Moderator_group == "Base model",] |> droplevels()
result_tot$label_plot <- result_tot$label

#rename label in the main db
db_subset$label_plot <- db_subset$Independent_macro
levels(db_subset$label_plot) <- levels(result_tot$label_plot)

#sort
result_tot$label_plot <- factor(result_tot$label_plot, levels = rev(result_tot$label_plot))
db_subset$label_plot  <- factor(db_subset$label_plot, levels = rev(result_tot$label_plot))

#plot
(plot_base_model <- ggplot() +
  
  gghalves::geom_half_violin(
    data = db_subset,
    aes(x = label_plot, y = Pearson_r, fill = label_plot),
    side = "r",
    #fill = "grey70",
    color = NA,
    width = 0.6,
    trim = FALSE,
    position = position_nudge(x = 0.1)
  ) +
  
  geom_jitter(
    data = db_subset,
    aes(x = label_plot, y = Pearson_r, col = label_plot),
    size = 0.3,
    #color = "grey70",
    width = 0.05
  ) +
  
  geom_pointrange(
    data = result_tot,
    aes(x = label_plot, y = ES, ymin = L, ymax = U),
    size = 1
  ) +
    
  geom_text(data = result_tot, aes( y = -Inf, x = label_plot, label = paste0(N_study, " | ", N_est)), hjust = -0.05, size = 3,color = "grey30") +
    
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  coord_flip() +
  scale_fill_manual(values = rev(my_colors))+
  scale_color_manual(values = rev(my_colors))+
  ylim(-1, 1) +
  xlab("") +
  ylab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval")))+
  theme(legend.position = "none")
)

#subset for type
result_ind <- result_for_plot[result_for_plot$Moderator_group == "Independent_zoomed",] |> droplevels()

result_ind$Moderator2 <- paste(result_ind$label, result_ind$Moderator ,sep = " - ")

result_ind$Moderator2 <- factor(result_ind$Moderator2,
                               levels = rev(result_ind$Moderator2))

(plot_ind <- ggplot(
  result_ind,
  aes(y = Moderator2, x = ES, xmin = L, xmax = U, col = label)
) +
    geom_pointrange(size = 0.9) +
    geom_text(aes( x = -Inf, label = paste0(N_study, " | ", N_est)), hjust = -0.2, size = 3,color = "grey30") +
    geom_vline(xintercept = 0, lty = 2, color = "grey50") +
    coord_cartesian(xlim = c(-1, 1)) +
    xlab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval"))) +
    ylab("")+
    scale_fill_manual(values = my_colors)+
    scale_color_manual(values = my_colors)+
    theme(legend.position = "none"))

pdf(file = "Figures/Figure_2.pdf", width = 13, height =7)

(plot_base_model + plot_ind) +
  plot_annotation(tag_levels = list(c("A","B")))

dev.off()

# Now plot of each moderator ---------------------------------------------------

x_limits <- c(-0.8, 0.8)
strip <- strip_themed(background_x = elem_list_rect(fill = my_colors))

#subset for type
result_type <- result_for_plot[result_for_plot$Moderator_group == "Study_type",]

(plot_type <- ggplot(
  result_type,
  aes(y = Moderator, x = ES, xmin = L, xmax = U)
) +
  geom_pointrange(size = 0.9) +
  geom_text(aes( x = -Inf, label = paste0(N_study, " | ", N_est)), hjust = -0.1, size = 2,color = "grey30") +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  ggh4x::facet_wrap2(~ label, nrow = 1, strip = strip) +
  coord_cartesian(xlim = x_limits) +
  xlab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval"))) +
  ylab("Methodology"))

#subset for dependent
result_dep <- result_for_plot[result_for_plot$Moderator_group == "Dependent_zoomed",]

result_dep$Moderator <- factor(result_dep$Moderator,
                             levels = levels(result_dep$Moderator))

(plot_dep <- ggplot(
  result_dep,
  aes(y = Moderator, x = ES, xmin = L, xmax = U)
) +
    geom_pointrange(size = 0.9) +
    geom_text(aes( x = -Inf, label = paste0(N_study, " | ", N_est)), hjust = -0.1, size = 2,color = "grey30") +
    geom_vline(xintercept = 0, lty = 2, color = "grey50") +
    ggh4x::facet_wrap2(~ label, nrow = 1, strip = strip) +
    coord_cartesian(xlim = x_limits) +
    xlab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval"))) +
    ylab("Dependent variable"))

#subset for geography
result_geo <- result_for_plot[result_for_plot$Moderator_group == "Geography_macro",]

#Sort
result_geo$Moderator <- factor(result_geo$Moderator,
                           levels = rev(c("Global",
                                          "Africa",
                                          "Asia",
                                          "Europe",
                                          "N-America",
                                          "CS-America",
                                          "Oceania")))

(plot_geo <- ggplot(
  result_geo,
  aes(y = Moderator, x = ES, xmin = L, xmax = U)
) +
    geom_text(aes( x = -Inf, label = paste0(N_study, " | ", N_est)), hjust = -0.1, size = 2,color = "grey30") +
    geom_pointrange(size = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey50") +
    ggh4x::facet_wrap2(~ label, nrow = 1, strip = strip) +
    coord_cartesian(xlim = x_limits) +
    xlab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval"))) +
    ylab("Geography"))

#subset for taxonomy
result_taxa <- result_for_plot[result_for_plot$Moderator_group == "Taxon",]

result_taxa$Moderator <- factor(result_taxa$Moderator,
                            levels = rev(c("Multiple",
                                           "Plants",
                                           "Vertebrates",
                                           "Invertebrates",
                                           "Amphibians",
                                           "Birds",
                                           "Fish",
                                           "Mammals",
                                           "Reptiles")))

(plot_taxa <- ggplot(
  result_taxa,
  aes(y = Moderator, x = ES, xmin = L, xmax = U)
) +
    geom_text(aes( x = -Inf, label = paste0(N_study, " | ", N_est)), hjust = -0.1, size = 2,color = "grey30") +
    geom_pointrange(size = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey50") +
    ggh4x::facet_wrap2(~ label, nrow = 1, strip = strip) +
    coord_cartesian(xlim = x_limits) +
    xlab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval"))) +
    ylab("Taxon (versus control)"))

#Final plot

#remove strip & xaxis as needed from plots
plot_dep  <- plot_dep  + theme(strip.text = element_blank())
plot_geo  <- plot_geo  + theme(strip.text = element_blank())
plot_taxa <- plot_taxa + theme(strip.text = element_blank())

plot_type <- plot_type + theme(
  axis.text.x  = element_blank(),
  axis.title.x = element_blank()
)

plot_dep <- plot_dep + theme(
  axis.text.x  = element_blank(),
  axis.title.x = element_blank()
)

plot_geo <- plot_geo + theme(
  axis.text.x  = element_blank(),
  axis.title.x = element_blank()
)

#Assemble
pdf(file = "Figures/Figure_3.pdf", width = 10, height = 10)

(final_plot <- 
  plot_type +
  plot_dep +
  plot_geo +
  plot_taxa +
  plot_layout(
    ncol = 1,
    heights = c(.3, .5, .8, 1)  # tweak depending on number of rows
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major.y = element_blank()
  )+
    plot_annotation(tag_levels = list(c("A","B","C","D"))))

dev.off()

sessionInfo()
#end
