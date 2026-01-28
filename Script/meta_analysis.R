####################################################################################
# TITLE #
####################################################################################

##################
# R script to generate the analysis
# Analysis performed with R (v. 4.4.1) and Rstudio (v. XXX)
# Created by Stefano Mammola
# Initiated 16/12/2025
##################

####################################################################################
# Preparation ----------------------------------------------------------------------
####################################################################################

# Loading R packages ------------------------------------------------------

if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("dplyr", "gghalves", "ggplot2","metafor", "patchwork","stringr","readxl","tidyr")

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

# db <- db[db$Independent_macro != "Taxonomic bias",]
# db <- droplevels(db)

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

#Average number of publication per study
mean(table(db$ID))
sd(table(db$ID))

#Number of journals
nlevels(db$Journal) #54

#Range of years 
range(db$Year)

table(db$Taxon_macro)/sum(table(db$Taxon_macro))

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

# Publication bias --------------------------------------------------------

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
    # geom_jitter(data = db, aes(x = Label, y = Pearson_r),
    #              size = 0.3, color = "grey70", width = 0.05)+
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

#Sort
result_geo$Label <- factor(result_geo$Label,
                   levels = rev(c("Global (32, 377)",
                              "Africa (3, 53)",
                              "Asia (3, 56)",
                              "Europe (30, 189)",
                              "N-America (9, 98)",
                              "CS-America (5, 71)",
                              "Oceania (4, 61)")))

#Plot
(plot_geo <- ggplot(data = result_geo) +
    xlab("")+
    ylab("")+
    # geom_jitter(data = db, aes(x = Label, y = Pearson_r),
    #             size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Taxon moderators ------------------------------------------------
db$Taxon <- relevel(db$Taxon, "Multiple") #setting baseline

model_taxa <- metafor::rma.mv(yi, vi, mods = ~ 0 + Taxon, random = list(~1 | ID), data = db)

#Extract estimates
levs <- levels(db$Taxon)

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
  sapply(levs, function(l) length(unique(db[db$Taxon == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db[db$Taxon == l,])),
  ")"
)

label_map <- setNames(result_taxa$Label, result_taxa$Label)

db$Label <- factor(
  result_taxa$Label[match(db$Taxon, levs)],
  levels = result_taxa$Label
)

#Sort
result_taxa$Label <- factor(result_taxa$Label,
                           levels = rev(c("Multiple (20, 104)",
                                          "Plants (21, 141)",
                                          "Vertebrates (5, 29)",
                                          "Invertebrates (6, 40)",
                                          "Amphibians (4, 87)",
                                          "Birds (18, 216)",
                                          "Fish (4, 29)",
                                          "Mammals (14, 168)",
                                          "Reptiles (4, 91)")))

#plot
(plot_taxa <- ggplot(data = result_taxa) +
    xlab("")+
    ylab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval")))+
    # geom_jitter(data = db, aes(x = Label, y = Pearson_r),
    #             size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-1,1))

# Trait moderators ------------------------------------------------

#dropping general estimates on taxonomic bias without specific traits
db_subset <- db
db_subset <- db_subset[db_subset$Independent_macro != "Taxonomic bias",]
db_subset <- db_subset[db_subset$Independent_macro != "Other",] |> droplevels()

db_subset$Independent <- paste(db_subset$Independent_macro, db_subset$Independent_zoomed, sep = " - ") |> as.factor()

model_trait <- metafor::rma.mv(yi, vi, mods = ~ 0 + Independent, random = list(~ 1 | ID), data = db_subset)

#Extract estimates
levs <- levels(db_subset$Independent)

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
  sapply(levs, function(l) length(unique(db_subset[db_subset$Independent == l,]$ID))),
  ", ",
  sapply(levs, function(l) nrow(db_subset[db_subset$Independent == l,])),
  ")"
)

result_trait$Category <- trimws(sub(" -.*$", "", result_trait$Label))

label_map <- setNames(result_trait$Label, result_trait$Label)

db_subset$Label <- factor(
  result_trait$Label[match(db_subset$Independent, levs)],
  levels = result_trait$Label
)

#Sort
result_trait$Label <- factor(result_trait$Label,
                            levels = rev(result_trait$Label))

#plot
(plot_trait <- ggplot(data = result_trait) +
    xlab("")+
    ylab(expression(paste("Effect size [",italic("r"),"]" %+-% "95% Confidence interval")))+
    # geom_jitter(data = db_subset, aes(x = Label, y = Pearson_r, fill = Independent_macro, color = Independent_macro),
    #             size = 0.3, width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U, fill = Category, color = Category), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-.3,.5)+theme(legend.position = "none"))

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

#Sort
result_trait$Label <- factor(result_trait$Label,
                             levels = rev(result_trait$Label))

#plot
(plot_dep <- ggplot(data = result_dep) +
    xlab("")+
    ylab("")+
    geom_jitter(data = db, aes(x = Label, y = Pearson_r),
                size = 0.3, color = "grey70", width = 0.05)+
    geom_pointrange(aes(x = Label, y = ES, ymin = L, ymax = U), size = 1) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +  
    coord_flip()+
    ylim(-.3,1))

# General figure ----------------------------------------------------------

left  <- plot_base_model

right <- wrap_plots(
  plot_dep,
  plot_trait,
  nrow = 2
)

left_up <- plot_type

left_low <- wrap_plots(
  plot_geo,
  plot_taxa,
  nrow = 2
)

design <- "
  12
  32
  42
  42
  42
  42
"
(left + right + left_up + left_low) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A","D","","B","C")))

pdf(file = "Figures/Figure_1.pdf", width = 13, height = 9)

(left + right + left_up + left_low) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A","D","","B","C")))

dev.off()

# test --------------------------------------------------------------------

#dropping general estimates on taxonomic bias without specific traits
db_subset <- db
db_subset <- db_subset[db_subset$Independent_macro != "Taxonomic bias",]
db_subset <- db_subset[db_subset$Independent_macro != "Other",] |> droplevels()

db_subset$Independent <- paste(db_subset$Independent_macro, db_subset$Independent_zoomed, sep = " - ") |> as.factor()

SUBSET    <- list()
MODEL     <- list()

result_for_plot <- data.frame(label = NULL,
                              b     = NULL,
                              ci.lb = NULL,
                              ci.ub = NULL,
                              ES    = NULL,
                              L     = NULL,
                              U     = NULL)

predictors_to_analyse <- levels(db_subset$Independent)

for (i in 1:length(predictors_to_analyse)){  
  
  #subset the predictor
  data_i  <- db_subset[db_subset$Independent == predictors_to_analyse[i], ]
  
  #fitting the model
  model_i <- rma.mv(yi, vi, random =  ~ 1 | ID, data = data_i) 
  
  #extracting coefficients
  result_for_plot_i <- data.frame(label = paste(predictors_to_analyse[i]," (" ,
                                                length(unique(data_i$ID)),", ",
                                                nrow(data_i),")",sep=''),
                                  b     = model_i$b,
                                  ci.lb = model_i$ci.lb,
                                  ci.ub = model_i$ci.ub,
                                  ES    = ((exp(model_i$b)-1))/((exp(model_i$b)+1)),
                                  L     = ((exp(model_i$ci.lb)-1)/(exp(model_i$ci.lb)+1)),
                                  U     = ((exp(model_i$ci.ub)-1)/(exp(model_i$ci.ub)+1)))
  
  #store the data 
  SUBSET[[i]]     <- data_i
  MODEL[[i]]      <- model_i
  result_for_plot <- rbind(result_for_plot,result_for_plot_i)
  
}


(plot1 <- ggplot(data= result_for_plot, aes(x=label, y=ES, ymin=L, ymax=U)) +
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    #geom_hline(yintercept=0.373, lty=3, col="grey70") +  
    geom_pointrange(size= 1) + ylim(-.3,0.5) + coord_flip())
