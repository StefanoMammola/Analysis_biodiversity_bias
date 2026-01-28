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


# Function ----------------------------------------------------------------

#Fucntion to model with moderators and extract results
extract_moderator_results <- function(data, moderator, i,
                                      min_n = 5,
                                      min_ID = 3) {
  
  # Guard: insufficient data 
  if (nrow(data) < min_n || length(unique(data$ID)) < min_ID) {
    return(
      data.frame(
        label = paste0(
          i, " (",
          length(unique(data$ID)), ", ",
          nrow(data), ")"
        ),
        Moderator_group = moderator,
        Moderator = NA,
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
        label = paste0(
          i, " (",
          length(unique(data$ID)), ", ",
          nrow(data), ")"
        ),
        Moderator_group = moderator,
        Moderator = NA,
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
  
  out <- data.frame(
    label = paste0(
      i, " (",
      length(unique(data$ID)), ", ",
      nrow(data), ")"
    ),
    Moderator_group = moderator,
    Moderator = sub(
      paste0("^", moderator),
      "",
      rownames(model$beta)
    ),
    b     = model$b,
    ci.lb = model$ci.lb,
    ci.ub = model$ci.ub,
    status = "ok"
  )
  
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

#dropping general estimates in the category Other (too few observation)
db_subset <- db
db_subset <- db_subset[db_subset$Independent_macro != "Other",] |> droplevels()

########################@###################################################
# Generating all models in a loop -----------------------------------------
########################@###################################################

# Moderators
moderators <- c(
  "Study_type",
  "Geography_macro",
  "Taxon_macro",
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
    label = paste0(
      i, " (",
      length(unique(data_i$ID)), ", ",
      nrow(data_i), ")"
    ),
    Moderator_group = "Base model",
    Moderator = "Overall",
    b     = model_base$b,
    ci.lb = model_base$ci.lb,
    ci.ub = model_base$ci.ub,
    status = "ok"
  )
  
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

result_tot <- result_for_plot[result_for_plot$Moderator_group == "Base model",]

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

(plot1 <- ggplot(data= result_for_plot, aes(x=label, y=ES, ymin=L, ymax=U, col= Discipline)) +
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    #geom_hline(yintercept=0.373, lty=3, col="grey70") +  
    geom_pointrange(size= 1) + ylim(-.3,0.5) + coord_flip())

