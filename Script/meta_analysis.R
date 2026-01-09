####################################################################################
# TITLE #
####################################################################################

##################
# R script to generate the analysis
# Analysis performed with R (v. XXX) and Rstudio (v. XXX)
# Created by Stefano Mammola
# Initiated 16/12/2025
##################

####################################################################################
# Preparation ----------------------------------------------------------------------
####################################################################################

# Loading R packages ------------------------------------------------------

stringr, tidyr

if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("dplyr", "gghalves", "ggplot2","metafor", "patchwork", "readxl")

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
# Data exploration -----------------------------------------------------------------
####################################################################################

#Number of papers
nrow(db2) #78

#Number of estimates
nrow(db) #905

#Number of journals
nlevels(db$Journal) #54

#Publication by year
db2 |>
  ggplot2::ggplot(aes(Year))+
  geom_bar(fill = "#2c7fb8", color = "black")+
  labs(
    x = "",
    y = "Number of publications"
  ) +
  scale_x_continuous(breaks = seq(min(db2$Year), max(db2$Year), by = 2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Studies by country
table(db2$Geography_macro) |>
  data.frame() |>
  ggplot2::ggplot(aes(x = Freq, y = Var1))+
  geom_col(fill = "#2c7fb8")+
  labs(
    y = "",
    x = "Number of publications"
  )

#Studies by country
table(db2$Study_type) |>
  data.frame() |>
  ggplot2::ggplot(aes(x = Freq, y = Var1))+
  geom_col(fill = "#2c7fb8")+
  labs(
    y = "",
    yx = "Number of publications"
  )

#Summary of estimates in each combination of response and predictor
(table_tot <- table(db$Dependent_zoomed,db$Independent_zoomed))

####################################################################################
# Data analysis -----------------------------------------------------------------
####################################################################################

# Derive Fischer's Z and its variance -------------------------------------

db$Pearson_r <- pmin(pmax(db$Pearson_r, -0.999), 0.999)
db <- metafor::escalc(measure = "ZCOR", ri = Pearson_r, ni = N, data = db)

# Base model - All estimates ----------------------------------------------

model_tot <- metafor::rma.mv(yi, vi, random = list(~1 | ID), data = db)

result_tot <- data.frame(label = paste0("Intercept (", nrow(db2),", ", nrow(db),")"),
                                b     = model_tot$b,
                                ci.lb = model_tot$ci.lb,
                                ci.ub = model_tot$ci.ub,
                                ES    = ((exp(model_tot$b)-1))/((exp(model_tot$b)+1)),
                                L     = ((exp(model_tot$ci.lb)-1)/(exp(model_tot$ci.lb)+1)),
                                U     = ((exp(model_tot$ci.ub)-1)/(exp(model_tot$ci.ub)+1)))

(plot_base_model <- ggplot(data = result_tot) +
    xlab("")+
    ylab("")+
    gghalves::geom_half_violin(
    data = db,
    aes(x = result_tot$label, y = Pearson_r),
    side = "r",     
    fill = "grey70",
    color = NA,
    width = 0.6,
    trim = FALSE,
    position = position_nudge(x = 0.1))+
  
    geom_jitter(data = db, aes(x = result_tot$label, y = Pearson_r), 
                size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Methodological moderators ------------------------------------------------
model_type <- metafor::rma.mv(yi, vi, 
                              mods = ~ 0 + Study_type, 
                              random = list(~1 | ID), 
                              data = db)

#Extract estimates
levs <- levels(db$Study_type)

result_type <- data.frame(
  Label = levs,
  b     = model_type$b,
  ci.lb = model_type$ci.lb,
  ci.ub = model_type$ci.ub,
  ES    = ((exp(model_type$b)-1))/((exp(model_type$b)+1)),
  L     = ((exp(model_type$ci.lb)-1)/(exp(model_type$ci.lb)+1)),
  U     = ((exp(model_type$ci.ub)-1)/(exp(model_type$ci.ub)+1))
)

#rename the label
result_type$Label <- paste0(
  levs, " (",
  sapply(levs, function(l) length(unique(db[db$Study_type == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Study_type == l,])),
  ")"
)


label_map <- setNames(result_type$Label, result_type$Label)

db$Label <- factor(
  result_type$Label[match(db$Study_type, levs)],
  levels = result_type$Label
)

(plot_type <- ggplot(data = result_type) +
    xlab("")+
    ylab("")+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r),
                 size = 0.3, color = "grey70", width = 0.05)+
    
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Geographic moderators ------------------------------------------------
db$Geography_macro <- relevel(db$Geography_macro, "Global") #setting baseline

model_geo <- metafor::rma.mv(yi, vi, mods = ~ 0 + Geography_macro, random = list(~1 | ID), data = db)

#Extract estimates
levs <- levels(db$Geography_macro)

result_geo <- data.frame(
  Label = levs,
  b     = model_geo$b,
  ci.lb = model_geo$ci.lb,
  ci.ub = model_geo$ci.ub,
  ES    = ((exp(model_geo$b)-1))/((exp(model_geo$b)+1)),
  L     = ((exp(model_geo$ci.lb)-1)/(exp(model_geo$ci.lb)+1)),
  U     = ((exp(model_geo$ci.ub)-1)/(exp(model_geo$ci.ub)+1))
)

#rename the label
result_geo$Label <- paste0(
  levs, " (",
  sapply(levs, function(l) length(unique(db[db$Geography_macro == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Geography_macro == l,])),
  ")"
)

label_map <- setNames(result_geo$Label, result_geo$Label)

db$Label <- factor(
  result_geo$Label[match(db$Geography_macro, levs)],
  levels = result_geo$Label
)

(plot_geo <- ggplot(data = result_geo) +
    xlab("")+
    ylab("")+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r),
                size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Taxon moderators ------------------------------------------------
db$Taxon_macro <- relevel(db$Taxon_macro, "Multiple") #setting baseline

model_taxa <- metafor::rma.mv(yi, vi, mods = ~ 0 + Taxon_macro, random = list(~1 | ID), data = db)

#Extract estimates
levs <- levels(db$Taxon_macro)

result_taxa <- data.frame(
  Label = levs,
  b     = model_taxa$b,
  ci.lb = model_taxa$ci.lb,
  ci.ub = model_taxa$ci.ub,
  ES    = ((exp(model_taxa$b)-1))/((exp(model_taxa$b)+1)),
  L     = ((exp(model_taxa$ci.lb)-1)/(exp(model_taxa$ci.lb)+1)),
  U     = ((exp(model_taxa$ci.ub)-1)/(exp(model_taxa$ci.ub)+1))
)

#rename the label
result_taxa$Label <- paste0(
  levs, " (",
  sapply(levs, function(l) length(unique(db[db$Taxon_macro == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Taxon_macro == l,])),
  ")"
)

label_map <- setNames(result_taxa$Label, result_taxa$Label)

db$Label <- factor(
  result_taxa$Label[match(db$Taxon_macro, levs)],
  levels = result_taxa$Label
)

(plot_taxa <- ggplot(data = result_taxa) +
    xlab("")+
    ylab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval")))+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r),
                size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Trait moderators ------------------------------------------------
db$Independent <- paste(db$Independent_macro, db$Independent_zoomed, sep = " - ") |> as.factor()

model_trait <- metafor::rma.mv(yi, vi, mods = ~ 0 + Independent, random = list(~1 | ID), data = db)

#Extract estimates
levs <- levels(db$Independent)

result_trait <- data.frame(
  Label = levs,
  b     = model_trait$b,
  ci.lb = model_trait$ci.lb,
  ci.ub = model_trait$ci.ub,
  ES    = ((exp(model_trait$b)-1))/((exp(model_trait$b)+1)),
  L     = ((exp(model_trait$ci.lb)-1)/(exp(model_trait$ci.lb)+1)),
  U     = ((exp(model_trait$ci.ub)-1)/(exp(model_trait$ci.ub)+1))
)

#rename the label
result_trait$Label <- paste0(
  levs, " (",
  sapply(levs, function(l) length(unique(db[db$Independent == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Independent == l,])),
  ")"
)

result_trait$Category <- trimws(sub(" -.*$", "", result_trait$Label))

label_map <- setNames(result_trait$Label, result_trait$Label)

db$Label <- factor(
  result_trait$Label[match(db$Independent, levs)],
  levels = result_trait$Label
)

(plot_trait <- ggplot(data = result_trait) +
    xlab("")+
    ylab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval")))+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r, fill = Independent_macro, color = Independent_macro),
                size = 0.3, width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U, fill = Category, color = Category), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1)+theme(legend.position = "none"))

# Dependent moderators ---------------------------------------------
model_dep <- metafor::rma.mv(yi, vi, mods = ~ 0 + Dependent_zoomed, random = list(~1 | ID), data = db)

#Extract estimates
levs <- levels(db$Dependent_zoomed)

result_dep <- data.frame(
  Label = levs,
  b     = model_dep$b,
  ci.lb = model_dep$ci.lb,
  ci.ub = model_dep$ci.ub,
  ES    = ((exp(model_dep$b)-1))/((exp(model_dep$b)+1)),
  L     = ((exp(model_dep$ci.lb)-1)/(exp(model_dep$ci.lb)+1)),
  U     = ((exp(model_dep$ci.ub)-1)/(exp(model_dep$ci.ub)+1))
)

#rename the label
result_dep$Label <- paste0(
  levs, " (",
  sapply(levs, function(l) length(unique(db[db$Dependent_zoomed == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Dependent_zoomed == l,])),
  ")"
)

label_map <- setNames(result_dep$Label, result_dep$Label)

db$Label <- factor(
  result_dep$Label[match(db$Dependent_zoomed, levs)],
  levels = result_dep$Label
)

(plot_dep <- ggplot(data = result_dep) +
    xlab("")+
    ylab("")+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r),
                size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# General figure ----------------------------------------------------------

left  <- plot_base_model

middle <- wrap_plots(
  plot_dep,
  plot_trait,
  nrow = 2
)

right <- wrap_plots(
  plot_type,
  plot_geo,
  plot_taxa,
  nrow = 3
)

design <- "
  12
  32
  32
  32
  32
"
(left + middle + right) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A","C","","B")))

pdf(file = "Figures/Figure_1.pdf", width = 13, height = 9)

(left + middle + right) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A","C","","B")))

dev.off()

####################################
###### *** TESTING STUFF **** ######
####################################

## Final interaction????

# All possible combinations
all_combinations <- expand.grid(
  Ind = unique(db$Independent_macro),
  Dep = unique(db$Dependent_zoomed),
  stringsAsFactors = FALSE
)

# Count how many studies / estimates exist for each combination
combination_counts <- all_combinations %>%
  rowwise() %>%
  mutate(
    n_studies = length(unique(db$ID[db$Independent_macro == Ind & db$Dependent_zoomed == Dep])),
    n_estimates = sum(db$Independent_macro == Ind & db$Dependent_zoomed == Dep),
    present = n_studies > 0  # TRUE if this combination exists in your data
  )


db$Macro <- db$Independent_macro
db$Zoomed <- db$Independent_zoomed
db$Dep <- db$Dependent_zoomed

db$Ind <- interaction(db$Macro, db$Zoomed, sep = " - ", drop = TRUE)

model_interaction <- metafor::rma.mv(yi, vi, mods = ~ 0 + Ind * Dep, random = list(~1 | ID), data = db)

result_interaction <- unique(data.frame(
  Label = rownames(model_interaction$b),
  b        = model_interaction$b,
  ci.lb    = model_interaction$ci.lb,
  ci.ub    = model_interaction$ci.ub
))

result_interaction$ES <- (exp(result_interaction$b)    - 1) / (exp(result_interaction$b)    + 1)
result_interaction$L  <- (exp(result_interaction$ci.lb) - 1) / (exp(result_interaction$ci.lb) + 1)
result_interaction$U  <- (exp(result_interaction$ci.ub) - 1) / (exp(result_interaction$ci.ub) + 1)

result_interaction <- result_interaction |>
   dplyr::filter(str_detect(Label, ":")) |>
  tidyr::separate(Label, into = c("Ind", "Dep"), sep = ":") |>
  dplyr::mutate(
    Ind = str_remove(Ind, "^(Ind|Dep)"),
    Dep = str_remove(Dep, "^(Ind|Dep)")
  )

# Final labels
result_interaction$n_studies <- mapply(
  function(m, s)
    length(unique(db$ID[db$Ind == m & db$Dep == s])),
  result_interaction$Ind, result_interaction$Dep
)

result_interaction$n_estimates <- mapply(
  function(m, s)
    sum(db$Ind == m & db$Dep == s),
  result_interaction$Ind, result_interaction$Dep
)

result_interaction$Label <- paste0(result_interaction$Ind, " (", result_interaction$n_studies, ", ", result_interaction$n_estimates, ")")
result_interaction$Label <- factor(result_interaction$Label, levels = sort(as.character(result_interaction$Label)))

ggplot(result_interaction) +
  facet_grid(. ~ Dep, scales = "fixed", space = "fixed") +
  geom_pointrange(
    aes(x = Ind, y = ES, ymin = L, ymax = U),
    size = 1
  ) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  ylim(-.6, .5) +
  xlab("") +
  ylab(expression(
    paste("Effect size [", italic("r"), "] ± 95% CI")
  )) +
  theme(legend.position = "none")


library(metafor)
library(dplyr)
library(ggplot2)

## ---------------------------
## 1. Prepare factors
## ---------------------------

db <- db %>%
  mutate(
    Ind = interaction(
      Independent_macro,
      Independent_zoomed,
      sep = " - ",
      drop = TRUE
    ),
    Dep = factor(Dependent_zoomed)
  )

ind_levels <- levels(db$Ind)
dep_levels <- levels(db$Dep)

## ---------------------------
## 2. Fit model
## ---------------------------

m <- rma.mv(
  yi,
  vi,
  mods   = ~ Ind * Dep,
  random = ~ 1 | ID,
  data   = db
)

## ---------------------------
## 3. Build full grid
## ---------------------------

newdat <- expand.grid(
  Ind = ind_levels,
  Dep = dep_levels,
  stringsAsFactors = FALSE
)

newdat$Ind <- factor(newdat$Ind, levels = ind_levels)
newdat$Dep <- factor(newdat$Dep, levels = dep_levels)

## ---------------------------
## 4. Build model matrix and KEEP estimable rows
## ---------------------------

X <- model.matrix(~ Ind * Dep, data = newdat)

## Columns actually estimated in the model
coef_names <- colnames(m$X)

## Keep only columns common to both
X_common <- X[, colnames(X) %in% coef_names, drop = FALSE]

## Identify rows that are fully estimable
estimable <- apply(X_common, 1, function(r) all(is.finite(r)))

newdat <- newdat[estimable, , drop = FALSE]
X_common <- X_common[estimable, , drop = FALSE]

## ---------------------------
## 5. Predict
## ---------------------------

pred <- predict(m, newmods = X_common)

newdat <- newdat %>%
  mutate(
    b     = pred$pred,
    ci.lb = pred$ci.lb,
    ci.ub = pred$ci.ub
  )

## Fisher-z to r
newdat <- newdat %>%
  mutate(
    ES = (exp(b)     - 1) / (exp(b)     + 1),
    L  = (exp(ci.lb) - 1) / (exp(ci.lb) + 1),
    U  = (exp(ci.ub) - 1) / (exp(ci.ub) + 1)
  )

## ---------------------------
## 6. Add sample sizes
## ---------------------------

newdat <- newdat %>%
  rowwise() %>%
  mutate(
    n_studies   = n_distinct(db$ID[db$Ind == Ind & db$Dep == Dep]),
    n_estimates = sum(db$Ind == Ind & db$Dep == Dep)
  ) %>%
  ungroup()

## ---------------------------
## 7. Plot
## ---------------------------

ggplot(newdat) +
  facet_grid(. ~ Dep) +
  geom_pointrange(
    aes(x = Ind, y = ES, ymin = L, ymax = U),
    size = 1
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  ylim(-0.6, 0.5) +
  xlab("") +
  ylab(expression(
    paste("Effect size [", italic("r"), "] ± 95% CI")
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )



library(metafor)
library(dplyr)
library(ggplot2)

## ---------------------------
## 1. Prepare factors
## ---------------------------


db$Ind <- factor(db$Independent_macro)
db$DepM <- factor(db$Dependent_zoomed)

ind_levels  <- levels(db$Ind)
depm_levels <- levels(db$DepM)

## ---------------------------
## 2. Fit interaction model
## ---------------------------

m <- rma.mv(
  yi,
  vi,
  mods   = ~ Ind * DepM,
  random = ~ 1 | ID,
  data   = db
)

## ---------------------------
## 3. Build full Ind × DepM grid
## ---------------------------

newdat <- expand.grid(
  Ind  = ind_levels,
  DepM = depm_levels,
  stringsAsFactors = FALSE
)

newdat$Ind  <- factor(newdat$Ind,  levels = ind_levels)
newdat$DepM <- factor(newdat$DepM, levels = depm_levels)

## ---------------------------
## 4. Model matrix & estimable cells
## ---------------------------

X <- model.matrix(~ Ind * DepM, data = newdat)

## keep only columns that exist in fitted model
X_common <- X[, colnames(X) %in% colnames(m$X), drop = FALSE]

## identify estimable rows
estimable <- apply(X_common, 1, function(r) all(is.finite(r)))

newdat   <- newdat[estimable, , drop = FALSE]
X_common <- X_common[estimable, , drop = FALSE]

## ---------------------------
## 5. Predict marginal cell effects
## ---------------------------

pred <- predict(m, newmods = X_common)

newdat <- newdat %>%
  mutate(
    b     = pred$pred,
    ci.lb = pred$ci.lb,
    ci.ub = pred$ci.ub
  )

## Fisher-z to r
newdat <- newdat %>%
  mutate(
    ES = (exp(b)     - 1) / (exp(b)     + 1),
    L  = (exp(ci.lb) - 1) / (exp(ci.lb) + 1),
    U  = (exp(ci.ub) - 1) / (exp(ci.ub) + 1)
  )

## ---------------------------
## 6. Sample size information
## ---------------------------

newdat <- newdat %>%
  rowwise() %>%
  mutate(
    n_studies   = n_distinct(db$ID[db$Ind == Ind & db$DepM == DepM]),
    n_estimates = sum(db$Ind == Ind & db$DepM == DepM)
  ) %>%
  ungroup()

## optional: drop empty cells
newdat <- newdat %>% filter(n_estimates > 0)

## ---------------------------
## 7. Plot
## ---------------------------

ggplot(newdat) +
  facet_grid(. ~ DepM) +
  geom_pointrange(
    aes(x = Ind, y = ES, ymin = L, ymax = U),
    size = 1
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  ylim(-0.6, 0.5) +
  xlab("") +
  ylab(expression(
    paste("Effect size [", italic("r"), "] ± 95% CI")
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )


## --------------------------------------------------
##  Hierarchical moderator
## --------------------------------------------------

db$Macro    <- factor(db$Independent_macro)
db$Specific <- factor(db$Independent_zoomed)

db$Ind_full <- interaction(db$Macro, db$Specific, sep = " - ", drop = TRUE)

model_trait_hier <- rma.mv(yi, vi, mods   = ~ 0 + Macro / Specific, random = list(~1 | ID), data = db)

X <- model.matrix(model_trait_hier)
pred <- predict(model_trait_hier, newmods = X)

result_trait <- unique(data.frame(
  Macro    = db$Macro,
  Specific = db$Specific,
  Ind_full = db$Ind_full,
  b        = pred$pred,
  ci.lb    = pred$ci.lb,
  ci.ub    = pred$ci.ub
))

result_trait$ES <- (exp(result_trait$b)    - 1) / (exp(result_trait$b)    + 1)
result_trait$L  <- (exp(result_trait$ci.lb) - 1) / (exp(result_trait$ci.lb) + 1)
result_trait$U  <- (exp(result_trait$ci.ub) - 1) / (exp(result_trait$ci.ub) + 1)

# Final labels
result_trait$n_studies <- mapply(
  function(m, s)
    length(unique(db$ID[db$Macro == m & db$Specific == s])),
  result_trait$Macro, result_trait$Specific
)

result_trait$n_estimates <- mapply(
  function(m, s)
    sum(db$Macro == m & db$Specific == s),
  result_trait$Macro, result_trait$Specific
)

result_trait$Label <- paste0(result_trait$Macro, " - ", result_trait$Specific," (", result_trait$n_studies, ", ", result_trait$n_estimates, ")")
result_trait$Label <- factor(result_trait$Label, levels = sort(as.character(result_trait$Label)))

db$Label <- factor(
  result_trait$Label[match(db$Ind_full, result_trait$Ind_full)],
  levels = result_trait$Label
)

ggplot(result_trait) +
  xlab("") +
  ylab(expression(paste(
    "Effect size [", italic("r"), "]" %+-% "95% Confidence interval"
  ))) +
  
  geom_pointrange(
    aes(x = Label, y = ES, ymin = L, ymax = U,
        color = Macro, fill = Macro),
    size = 1
  ) +
  
  geom_jitter(
    data = db,
    aes(x = Label, y = Pearson_r,
        color = Macro, fill = Macro),
    size = 0.3, width = 0.05
  ) +
  
  geom_pointrange(
    aes(x = Label, y = ES, ymin = L, ymax = U,
        color = Macro, fill = Macro),
    size = 1
  ) +
  
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  coord_flip() +
  ylim(-1, 1) +
  theme(legend.position = "none")