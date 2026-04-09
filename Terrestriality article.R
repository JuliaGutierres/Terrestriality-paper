# Loading packages
library(janitor)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(dunn.test)
library(purrr)
library(tidyverse)
library(ggpubr)
library(DHARMa)
library(rstatix)
library(multcompView)

getwd()
setwd("C:/Users/jgutierres/Downloads")

# Load and clean data
data_tb <- read_delim(
  "sem_scans_tool.csv",
  delim = ";",
  col_types = cols(
    `display sexual` = col_double(),
    inter_vert = col_double(),
    vocalizacao = col_double(),
    copula_sexual = col_double()
  )
) %>%
  clean_names()

glimpse(data_tb)
n_distinct(data_tb$scan)


# Long format and expand counts
data_tb_long <- data_tb %>%
  pivot_longer(
    cols = desloc:copula_sexual,
    names_to = "activity",
    values_to = "n"
  ) %>%
  filter(n > 0)

data_tb_long <- data_tb_long %>%
  uncount(weights = n) %>%
  mutate(terrestrial = ifelse(stratum == "ground", 1, 0))

# Calculating proportion of individuals per estrato per scan
data_tb_prop_scan <- data_tb_long %>%
  group_by(scan, stratum, group, data) %>%
  summarise(n_estrato = n(), .groups = "drop") %>%
  group_by(scan) %>%
  mutate(prop_estrato = n_estrato / sum(n_estrato)) %>%
  ungroup()


# Ensuring all scan x statrum combinations exist
scans_group <- data_tb_long %>%
  distinct(scan, group, data)

all_scan_strata <- expand.grid(scan = unique(data_tb_long$scan),
                               stratum = unique(data_tb_long$stratum)) %>%
  left_join(scans_group, by = "scan")


prop_scan_stratum <- all_scan_strata %>%
  left_join(data_tb_prop_scan, by = c("scan", "stratum", "group", "data")) %>%
  mutate(prop_stratum = ifelse(is.na(prop_estrato), 0, prop_estrato))

# Mean proportion per day (method of proportions) by group
meangroup_day <- prop_scan_stratum %>%
  group_by(data, group, stratum) %>%
  summarise(
    prop_day = mean(prop_stratum, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

# Mean proportion per day (general, all bandos)
mean_day <- prop_scan_stratum %>%
  group_by(data, stratum) %>%
  summarise(
    prop_day = mean(prop_stratum, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  ) %>%
  mutate(group = "GERAL")

# Overall mean by stratum and group (across days)
mean_group_stratum <- meangroup_day %>%
  group_by(stratum, group) %>%
  summarise(
    proportion = mean(prop_day, na.rm = TRUE),
    sd = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = sd / sqrt(n_days),
    .groups = "drop"
  )

# Overall mean by stratum (general)
mean_stratum <- mean_day %>%
  group_by(stratum) %>%
  summarise(
    proportion = mean(prop_day, na.rm = TRUE),
    sd = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = sd / sqrt(n_days),
    .groups = "drop"
  ) %>%
  mutate(group = "General")

# Combining results
bind_rows(mean_group_stratum, mean_stratum)


# Wilcoxon test: comparing groups within solo stratum
wilcox.test(
  prop_day ~ group,
  data = meangroup_day %>% filter(stratum == "ground"),
  exact = FALSE
)

##########################################################################
################################# SEASON #################################
##########################################################################

# Proportion by scan 
prop_season_scan <- data_tb_long %>%
  group_by(data, period, scan, group) %>%
  summarise(
    total = n(),
    ground = sum(terrestrial == 1, na.rm = TRUE),
    prop_ground = ground / total,
    .groups = "drop"
  )

# Mean by day (method of proportions)
mean_season <- prop_season_scan %>%
  group_by(data, period) %>%
  summarise(
    prop_day = mean(prop_ground, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

# General (method of proportions)
mean_season %>%
  group_by(period) %>%
  summarise(
    mean_prop = mean(prop_day, na.rm = TRUE),
    se = sd(prop_day, na.rm = TRUE),
    n_dias = n(),
    sd = se / sqrt(n_dias),
    .groups = "drop"
  )


# By group
season_group_day <- prop_season_scan %>%
  group_by(data, period, group) %>%
  summarise(
    prop_day = mean(prop_ground, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

mean_season_group <- season_group_day %>%
  group_by(period, group) %>%
  summarise(
    mean_prop = mean(prop_day, na.rm = TRUE),
    dp = sd(prop_day, na.rm = TRUE),
    n_dias = n(),
    se = dp / sqrt(n_dias),
    .groups = "drop"
  )

print(mean_season_group)

mean_season_group %>%
  group_by(period, group) %>%
  summarise(n_dias = n())

# Wilcoxon test: groups in different periods
season_group_day %>%
  group_by(period) %>%
  group_map(~wilcox.test(prop_day ~ group, data = .x), .keep = TRUE)


# Wilcoxon test: season in different groups
season_group_day %>%
  group_by(group) %>%
  group_map(~wilcox.test(prop_day ~ period, data = .x), .keep = TRUE)

# Graphic
plot1 <- mean_season_group %>%
  mutate(
    period = recode(period,
                    "rainy" = "Rainy",
                    "dry" = "Dry"),
    group = recode(group,
                   "two" = "Group 2",
                   "one" = "Group 1"),
    period = factor(period, levels = c("Rainy", "Dry")),
    group = factor(group, levels = c("Group 1", "Group 2"))
  )

ggplot(plot1, aes(x = group, y = mean_prop * 100, fill = group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = (mean_prop - se) * 100,
                    ymax = (mean_prop + se) * 100),
                width = 0.2) +
  facet_wrap(~ period) +
  scale_fill_manual(
    name = " ",
    values = c("Group 1" = "#00BFC4",
               "Group 2" = "#F8766D")
  ) +
  scale_x_discrete(labels = NULL) +
  labs(
    x = NULL,
    y = "Ground use (%)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )


###########################################################################
############################## GENERAL AREAS ##############################
###########################################################################

# Proportion by scan 
prop_area_scan <- data_tb_long %>%
  group_by(data, general_local, scan, group) %>%
  summarise(
    total = n(),
    ground = sum(terrestrial == 1, na.rm = TRUE),
    prop_ground = ground / total,
    .groups = "drop"
  )

# Mean by day (method of proportions)
prop_area_day <- prop_area_scan %>%
  group_by(data, general_local) %>%
  summarise(
    prop_day = mean(prop_ground, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

# General (method of proportions)
general_area <- prop_area_day %>%
  group_by(general_local) %>%
  summarise(
    media_proporcao = mean(prop_day, na.rm = TRUE),
    se = sd(prop_day, na.rm = TRUE),
    n_dias = n(),
    sd = se / sqrt(n_dias),
    .groups = "drop"
  )

print(general_area)


# By group
mean_day_group <- prop_area_scan %>%
  group_by(data, general_local, group) %>%
  summarise(
    prop_day = mean(prop_ground, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

# General (method of proportions)
general_group <- mean_day_group %>%
  group_by(general_local, group) %>%
  summarise(
    mean_prop = mean(prop_day, na.rm = TRUE),
    dp = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = dp / sqrt(n_days),
    .groups = "drop"
  )

print(general_group)

# Wilcoxon test: by groups in different areas 
mean_day_group %>%
  group_by(general_local) %>%
  group_map(~wilcox.test(prop_day ~ group, data = .x), .keep = TRUE)


# Wilcoxon test: by areas in different groups
mean_day_group %>%
  group_by(group) %>%
  group_map(~wilcox.test(prop_day ~ general_local, data = .x), .keep = TRUE)


# Recode and reorder group and general_local
general_group <- general_group %>%
  mutate(
    # Sobrescrevemos a coluna group com os novos nomes
    group = recode(group, "two" = "Group 2", "one" = "Group 1"),
    # Transformamos em fator para garantir a ordem Group 1 -> Group 2
    group = factor(group, levels = c("Group 1", "Group 2")),
    # Ajustamos a coluna de local
    general_local = recode(general_local, "anthropic" = "Anthropic", "natural" = "Natural")
  )

ggplot(general_group, aes(x = group, y = mean_prop * 100, fill = group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = (mean_prop - se) * 100,
                    ymax = (mean_prop + se) * 100),
                width = 0.2) +
  facet_wrap(~ general_local) +
  scale_fill_manual(
    values = c("Group 1" = "#00BFC4", "Group 2" = "#F8766D"),
    name = NULL
  ) +
  scale_x_discrete(labels = NULL) +
  labs(
    x = NULL,
    y = "Ground use (%)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

###########################################################################
############################## SPECIFIC AREAS #############################
###########################################################################

# Proportion by scan
prop_sparea_scan <- data_tb_long %>%
  group_by(data, specific_local, scan, group) %>%
  summarise(
    total = n(),
    ground = sum(terrestrial == 1, na.rm = TRUE),
    prop_ground = ground / total,
    .groups = "drop"
  )

# Mean by day (method of proportions)
mean_sparea_day <- prop_sparea_scan %>%
  group_by(data, specific_local) %>%
  summarise(
    prop_day = mean(prop_ground, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  )

# General (method of proportions)
general_sparea <- mean_sparea_day %>%
  group_by(specific_local) %>%
  summarise(
    media_proporcao = mean(prop_day, na.rm = TRUE),
    dp = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = dp / sqrt(n_days),
    .groups = "drop"
  )

print(general_sparea)


# Kruskal-Wallis test 
kruskal.test(prop_day~ specific_local, data = mean_sparea_day)


##################################################################################
#################################### GLM  ########################################
##################################################################################

# Including provision date 
data_model <- data_tb_long %>%
  mutate(
    data_format = parse_date_time(data, orders = c("dmy", "ymd")), 
    supplementation = if_else(data_format <= as.Date("2024-07-01"), "yes", "no"),
    supplementation = as.factor(supplementation)
  )


# Filtering out with supplemented food
ali_supplemented <- c(
  "alim_folha", "alim_fruto", "alim_artro", "alim_flor", "alim_caule",
  "ali_semente", "ali_antro", "alim_lixo", "alim_raiz", "alim_leite"
)

# Items
supplemented <- c(
  "melzinho", "coca cola", "mexerica", "pera", "mamão", "mandioca",
  "banana", "bolo", "pão", "pao", "cana", "tomate supl", "inhame",
  "mandioca/jaca", "cana/inhame", "SUPLEMENT/pera", "ovo",
  "garrafa de agua", "garrafa", "bolacha"
)


data_filter <- data_model %>%
  filter(!(activity %in% ali_supplemented & outros %in% supplemented))


data_filter$sex_age <- relevel(factor(data_filter$sex_age), ref = "AM")

data_filter$period <- relevel(factor(data_filter$period), ref = "rainy")

data_filter$supplementation <- relevel(factor(data_filter$supplementation),
                                       ref = "no")
data_filter$general_local <- relevel(factor(data_filter$general_local), 
                                     ref = "natural")


model<- glm(terrestrial ~ sex_age + period + supplementation + general_local, 
                family = binomial, data = data_filter)

summary(model)

# GLM Diagnosis 
simulationOutput <- simulateResiduals(fittedModel = mod2_ref, n=1000, plot = TRUE)
plot(simulationOutput)

# Plot model
library(sjPlot)

plot_model(model, 
           show.values = TRUE, 
           value.offset = .3,
           axis.labels = c(
             "general_localanthropic" = "Anthropic Area",
             "supplementationyes"     = "Food provisioning",
             "perioddry"           = "Dry Season ",
             "sex_ageJU"             = "Juveniles",
             "sex_ageAF"             = "Adult Females"
           ),
           title = " ",
           vline.color = "grey20") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"))

ggsave("model_coefficients.png", width = 8, height = 5, dpi = 300)


library(ggeffects)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)


color_main <- "#3B6FB6"

theme_pub <- function() {
  theme_classic(base_size = 14) +
    theme(
      text = element_text(family = "serif"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, color = "black"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.tag = element_text(size = 16, face = "bold"),
      plot.margin = margin(15, 15, 15, 15)
    )
}


plot_marginal <- function(model, term, xlab) {
  
  pred <- ggpredict(model, terms = term) %>%
    as.data.frame()
  
  ggplot(pred, aes(x = x, y = predicted)) +
    
    geom_point(size = 4, color = color_main) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.15,
      linewidth = 0.7,
      color = color_main
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    coord_cartesian(ylim = c(0, 0.55)) +
    scale_x_discrete(expand = expansion(add = 0.3)) +
    labs(
      x = xlab,
      y = "Predicted probability"
    ) +
    
    theme_pub()
}

p1 <- plot_marginal(model, "sex_age", "Sex-age class")
p2 <- plot_marginal(model, "period", "Season")
p3 <- plot_marginal(model, "general_local", "Area")
p4 <- plot_marginal(model, "supplementation", "Food supplementation")

final_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(15, 15, 15, 15)
  )

print(final_plot)

ggsave(
  "Predicted.tiff",
  final_plot,
  width = 8,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


###########################################################################
######### SEX-AGE #########################################################
###########################################################################

# Proportion of terrestrial use per scan within each sex-age class
prop_sex_scan <- data_tb_long %>%
  group_by(sex_age, data, scan) %>% 
  summarise(
    total = n(),
    ground = sum(terrestrial == 1, na.rm = TRUE),
    prop_scan = ifelse(total > 0, ground / total, NA_real_)
  ) %>%
  ungroup()

# Proportion per day 
prop_sex_day <- prop_sex_scan %>%
  group_by(sex_age, data) %>%
  summarise(prop_day = mean(prop_scan, na.rm = TRUE)) %>%
  ungroup()

# General statistics by sex-age class
general_sex <- prop_sex_day %>%
  group_by(sex_age) %>%
  summarise(
    mean_prop = mean(prop_day, na.rm = TRUE),
    dp = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = dp / sqrt(n_days)
  ) %>%
  ungroup()

print(general_sex)

# Kruskal Wallis Test
kruskal.test(prop_day ~ sex_age, data = prop_sex_day)


# Sex-age within groups

# Proportion per scan by group and sex-age class
prop_sex_group <- data_tb_long %>%
  group_by(group, sex_age, data, scan) %>%
  summarise(
    total = n(),
    ground = sum(terrestrial == 1, na.rm = TRUE),
    prop_scan = ifelse(total > 0, ground / total, NA_real_)
  ) %>%
  ungroup()

# Proportion per day
prop_sex_group_day <- prop_sex_group %>%
  group_by(group, sex_age, data) %>%
  summarise(prop_day = mean(prop_scan, na.rm = TRUE)) %>%
  ungroup()

# General statistics
general_group_sex <- prop_sex_group_day %>%
  group_by(group, sex_age) %>%
  summarise(
    mean_prop = mean(prop_day, na.rm = TRUE),
    dp = sd(prop_day, na.rm = TRUE),
    n_dias = n(),
    se = dp / sqrt(n_dias)
  ) %>%
  ungroup()

print(general_group_sex)


# Kruskal-Wallis 
prop_sex_group_day %>%
  group_by(group) %>%
  group_map(~kruskal.test(prop_day ~ sex_age, data = .x), .keep = TRUE)

##################################################################################
####################### Proportion of tool use and all activities ################
##################################################################################

# Load and clean data
tool <- read_delim("com_tool.csv", 
                        delim = ";",
                        col_types = cols(
                          # Forçamos as colunas de comportamento a serem números
                          urina = col_double(),
                          vocalizacao = col_double(),
                          copula_sexual = col_double(),
                          # Se houver outras lgl que deveriam ser dbl, adicione aqui
                          .default = col_guess() 
                        )) %>% 
  clean_names() %>% 
  # Remove as colunas que começam com 'x' ou 'dot' (geralmente as vazias do Excel)
  select(-starts_with("x"), -starts_with("dot"), -starts_with("..."))

# Transforming in long data
data_long_tool <- tool %>%
  pivot_longer(cols = tool_use:copula_sexual, 
               names_to = "activity", 
               values_to = "n") %>% 
  filter(n > 0)

data_expanded_tool <- data_long_tool %>%
  uncount(weights = n) %>% 
  mutate(terrestrial = ifelse(stratum == "ground", 1, 0))

# Classifing behaviors 
tool_allbehav <- data_expanded_tool %>%
  filter(activity != "urina") %>%
  mutate(activity = case_when(
    str_starts(activity, "ali_agua") ~ "Drinking",
    str_starts(activity, "ali") ~ "Feeding",
    str_starts(activity, "for") ~ "Foraging",
    str_starts(activity, "desloc") ~ "Travelling",
    str_starts(activity, "desc") ~ "Resting",
    str_starts(activity, "tool") ~ "Tool use",
    str_starts(activity, "manip")~ "Manipulation",
    str_starts(activity, "inter")~ "Interaction",
    activity %in% c("obs", "auto_atividade", "vigilancia") ~ "Self-activitie",
    activity %in% c("copula_sexual", "copula não", "soli_copula", 
                     "display sexual", "catacao", "lipsmacking", 
                     "brincadeira", "vocalizacao","amamentando", 
                     "agonismo") ~ "Social",
    TRUE ~ activity
  ))

unique(tool_allbehav$activity)


tool_all_combos <- expand.grid(
  scan = unique(tool_allbehav$scan),
  activity = unique(tool_allbehav$activity)
)

# Proportion by scan
tool_prop_scan <- tool_allbehav %>%
  group_by(scan) %>%
  mutate(total_scan = n()) %>%
  group_by(scan, activity) %>%
  summarise(
    prop_activity = n() / unique(total_scan),
    .groups = "drop"
  )

tool_prop_scan_full <- tool_all_combos %>%
  left_join(tool_prop_scan, by = c("scan", "activity")) %>%
  mutate(prop_activity = ifelse(is.na(prop_activity), 0, prop_activity))

# Mean proportion
general_activity <- tool_prop_scan_full %>%
  group_by(activity) %>%
  summarise(
    mean_prop = mean(prop_activity),
    sd = sd(prop_activity),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

sum(general_activity$mean_prop)


#############################################################################
# Ground use in each activity ###############################################
#############################################################################


# Proportion of terrestrial use by activity per scan 
prop_activity_scan <- tool_allbehav %>%
  group_by(data, scan, activity) %>%
  summarise(
    total = n(),
    terrestrial = sum(terrestrial == 1, na.rm = TRUE),
    prop_terrestrial = terrestrial / total,
    .groups = "drop"
  )

# Mean proportion of terrestrial use per day 
prop_activity_day <- prop_activity_scan %>%
  group_by(data, activity) %>%
  summarise(
    prop_day = mean(prop_terrestrial, na.rm = TRUE),
    .groups = "drop"
  )

# Mean proportion across all days 
general_ground_activity <- prop_activity_day %>%
  group_by(activity) %>%
  summarise(
    mean_prop_terrestrial = mean(prop_day, na.rm = TRUE),
    sd = sd(prop_day, na.rm = TRUE),
    n_days = n(),
    se = sd / sqrt(n_days),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prop_terrestrial))

# Kruskal-Wallis Test 
kruskal.test(prop_day ~ activity, data = prop_activity_day)

# Dunn test
tool_allbehav$activity <- gsub("-", ".", tool_allbehav$activity)
prop_activity_day$activity <- gsub("-", ".", prop_activity_day$activity)

dunn_res <- dunn.test(tool_allbehav$terrestrial, tool_allbehav$activity, method = "bonferroni", list = TRUE)
pvals <- dunn_res$P.adjusted
names(pvals) <- gsub(" - ", "-", dunn_res$comparisons)
letters <- multcompLetters(pvals)$Letters
letters_df <- data.frame(activity = names(letters), letters = letters)

# Plot
prop_activity_day <- prop_activity_day %>%
  left_join(letters_df, by = "activity")

prop_activity_day$legend <- gsub("_", " ", prop_activity_day$activity)


plot_data <- general_ground_activity %>%
  mutate(activity = gsub("-", ".", activity)) %>%
  left_join(letters_df, by = "activity") %>%
  mutate(legend = gsub("\\.", "-", activity),
         legend = gsub("_", " ", legend))

ggplot(plot_data, aes(x = reorder(legend, -mean_prop_terrestrial), 
                      y = mean_prop_terrestrial * 100)) +
  geom_col(fill = "gray40", width = 0.7) +
  geom_errorbar(aes(ymin = (mean_prop_terrestrial - se) * 100,
                    ymax = (mean_prop_terrestrial + se) * 100),
                width = 0.2, color = "gray30") +
   geom_text(aes(label = paste0(round(mean_prop_terrestrial * 100, 1), "%")),
            vjust = -1.2, size = 5, fontface = "bold") +
  geom_text(aes(y = (mean_prop_terrestrial + se) * 100 + 4, label = letters),
            size = 5, fontface = "bold") +
  labs(
    x = "Activities",
    y = "Ground use (%)"
  ) +
  theme_minimal(base_size = 19) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

###################### TOOL USE #####################

# Loading data 
seasons <- read_excel("C:/Users/jgutierres/Downloads/tool_use.xlsx", sheet = "c3")


# Comparing season
tabela <- table(seasons$epoca, seasons$manip)
chisq.test(tabela)        

# Season balanced   
rainy_events <- 37
dry_events <- 83
rainy_days <- 331
dry_days <- 271


success <- c(rainy_events, dry_events)      # c(37, 83)
trials  <- c(rainy_days,   dry_days)        # c(331, 271)

# Proportion test 
prop.test(success, trials)



### Food provision ####

data_local <- seasons %>%
  mutate(provisao = if_else(as.Date(data) <= as.Date("2024-07-31"), "yes", "no"))

table(data_local$provisao)

# Chi-squared test for tool use and provisioning
teste_qui <- chisq.test(table(data_local$provisao))
print(teste_qui)


# Food provision balanced 

resumo_provisao <- seasons %>%
  mutate(
    data     = as_date(data),
    provisao = if_else(data <= as_date("2024-07-31"), "yes", "no")
  ) %>%
  group_by(provisao) %>%
  summarise(
    eventos = n(),
    dias = as.integer(max(data) - min(data) + 1),
    taxa_dia = eventos / dias,
    taxa_100dias = taxa_dia * 100,
    .groups = "drop"
  )

print(resumo_provisao)

prop.test(
  x = resumo_provisao$eventos,
  n = resumo_provisao$dias
)


# Preparing data
rate_season <- data.frame(
  Condition = c("Rainy", "Dry"),
  Rate = c(37/331, 83/271) * 100, 
  Category = "Season"
)


rate_supp <- resumo_provisao %>%
  mutate(
    Condition = ifelse(provisao == "yes", "Yes", "No"),
    Rate = taxa_100dias, 
    Category = "Supplementation"
  ) %>%
  select(Condition, Rate, Category)

# Plot
plot_data <- bind_rows(rate_season, rate_supp) %>%
  mutate(
    Category = factor(Category, levels = c("Season", "Supplementation")),
    Condition = factor(Condition, levels = c("Dry", "Rainy", "No", "Yes"))
  )

### Plot  ###
ggplot(plot_data, aes(x = Condition, y = Rate)) +
  geom_bar(fill = "gray40",
           color = "black",
           stat = "identity",
           width = 0.4) +
  facet_wrap(~ Category, 
             scales = "free_x", 
             nrow = 1, 
             strip.position = "bottom") +
  labs(
    x = "",
    y = "Tool-use rate"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y =  element_blank(),
    axis.text.x = element_text(face = "bold")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


ggsave("tool_use_rate.png", width = 10, height = 4, dpi = 300, bg = "white")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  


