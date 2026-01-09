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
db <- read_excel("Data/Database_meta_analysis_bias.xlsx", na = "NA") |>
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

#### plotting

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


