library(tidyverse)
library(taxize)

SpeciesPerPlotID <- readRDS("SpeciesPerPlotID.rds")

SpeciesPerPlotID$LatArt <- stringr::str_conv(SpeciesPerPlotID$LatArt, "UTF-8")

UniqueLatart <- data.frame(LatArt = unique(SpeciesPerPlotID$LatArt))

UniqueLatartTax <- list()

for(i in 1:nrow(UniqueLatart)){
  try({UniqueLatartTax[[i]] <-gnr_resolve(UniqueLatart$LatArt[i], data_source_ids = "11", canonical = T, highestscore = T, best_match_only = T)})
  if((i %% 100) == 0){
    message(paste(i, "of", nrow(UniqueLatart), Sys.time()))
  }
}

UniqueLatartTax <- UniqueLatartTax %>%
  purrr::reduce(bind_rows) %>%
  dplyr::rename(LatArt = user_supplied_name) %>%
  dplyr::select(LatArt, score, matched_name2)


SpeciesPerPlotID <- SpeciesPerPlotID %>%
  left_join(UniqueLatartTax) %>%
  dplyr::filter(!is.na(matched_name2))

SpeciesPerPlotIDCorrected <-  SpeciesPerPlotID %>%
  dplyr::select(-LatArt, -score) %>%
  rename(LatArt = matched_name2)

saveRDS(SpeciesPerPlotIDCorrected, "SpeciesPerPlotIDCorrected.rds")


Test <- readRDS("CRS_Final.rds")

Test$Species <- stringr::str_conv(Test$Species, "UTF-8")

UniqueCrsSpecies <- data.frame(Species = unique(Test$Species))

UniqueCrsSpeciesTax <- list()

for(i in 1:nrow(UniqueCrsSpecies)){
  try({UniqueCrsSpeciesTax[[i]] <-gnr_resolve(UniqueCrsSpecies$Species[i], data_source_ids = "11", canonical = T, highestscore = T, best_match_only = T)})
  if((i %% 100) == 0){
    message(paste(i, "of", nrow(UniqueLatart), Sys.time()))
  }
}

UniqueCrsSpeciesTax <- UniqueCrsSpeciesTax %>%
  purrr::reduce(bind_rows) %>%
  dplyr::rename(Species = user_supplied_name) %>%
  dplyr::select(Species, score, matched_name2)


Test <- Test %>%
  left_join(UniqueCrsSpeciesTax) %>%
  dplyr::filter(!is.na(matched_name2))

CorrectedCRS <-  Test %>%
  dplyr::select(-Species) %>%
  rename(Species = matched_name2)

Corrected <- list()

Taxa <- unique(CorrectedCRS$Species)

for(i in 1:length(Taxa)){
  Corrected[[i]] <- CorrectedCRS %>%
    dplyr::filter(Species == Taxa[i]) %>%
    dplyr::filter(score == max(score)) %>%
    sample_n(size = 1)
}


Corrected <- Corrected %>%
  purrr::reduce(bind_rows)

PlotCRS <- Corrected %>%
  rename(LatArt = Species) %>%
  right_join(SpeciesPerPlotIDCorrected) %>%
  distinct() %>%
  dplyr::filter(!is.na(score))


#### Calculate CRS per site

CRS <- function(leaf_area, leaf_dry_mass, leaf_fresh_mass, Species){
  Strategy <- tibble::tribble(
    ~Strategy,  ~C,  ~S,  ~R,
    "C", 90L,  5L,  5L,
    "C/CR", 73L,  5L, 23L,
    "C/CS", 73L, 23L,  5L,
    "CR", 48L,  5L, 48L,
    "C/CSR", 54L, 23L, 23L,
    "CS", 48L, 48L,  5L,
    "CR/CSR", 42L, 17L, 42L,
    "CS/CSR", 42L, 42L, 17L,
    "R/CR", 23L,  5L, 73L,
    "CSR", 33L, 33L, 33L,
    "S/CS", 23L, 73L,  5L,
    "R/CSR", 23L, 23L, 54L,
    "S/CSR", 23L, 54L, 23L,
    "R",  5L,  5L, 90L,
    "SR/CSR", 17L, 42L, 42L,
    "S",  5L, 90L,  5L,
    "R/SR",  5L, 23L, 73L,
    "S/SR",  5L, 73L, 23L,
    "SR",  5L, 48L, 48L
  )
  ## Calculate leaf suculent index
  LSI <- (leaf_fresh_mass - leaf_dry_mass)/(leaf_area/10)
  ## Calculate leaf water content
  LWC <- ((leaf_fresh_mass-leaf_dry_mass)/leaf_fresh_mass)*100
  ## leaf mass per area
  LMA <- leaf_dry_mass/(leaf_area/1000)
  ## leaf_dry_matter_content
  LDMC <- (leaf_dry_mass/leaf_fresh_mass)*100
  ## Specific leaf Area
  SLA <- leaf_area/leaf_dry_mass
  ## Transformations
  #SQRT_max_LA
  SQRT_max_LA <- sqrt(leaf_area/894205)*100
  LDW2 <- leaf_area/SLA
  LFW2 <- (100/LDMC)*LDW2
  LSI2 <- (LFW2-LDW2)/(leaf_area/10)
  LDMC2 <- ifelse(LSI2 > 5, ((100-(LDW2*100)/LFW2)), ((LDW2*100)/LFW2))
  Logit_LDM <- log((LDMC2/100)/(1-(LDMC2/100)))
  Log_SLA <- log(SLA)
  PCA2_C <- -0.8678 + 1.6464*SQRT_max_LA
  PCA1_S <- 1.3369+0.000010019*(1-exp(-2.2303E-12*Logit_LDM))+4.5835*(1-exp(-0.2328*Logit_LDM))
  PCA1_R <- -57.5924 + 62.6802*exp(-0.0288*Log_SLA)
  Min_C <- 0
  Min_S <- -0.756451214853076
  Min_R <- -11.3467682227961
  Negative_Outlier_c_C <- ifelse(PCA2_C < Min_C, Min_C, PCA2_C)
  Negative_Outlier_c_S <- ifelse(PCA1_S < Min_S, Min_S, PCA1_S)
  Negative_Outlier_c_R <- ifelse(PCA1_R < Min_R, Min_R, PCA1_R)

  Max_C <- 57.3756711966087
  Max_S <- 5.79158377609218
  Max_R <- 1.10795515716546
  Positive_Outlier_c_C <- ifelse(PCA2_C > Max_C, Max_C, PCA2_C)
  Positive_Outlier_c_S <- ifelse(PCA1_S > Max_S, Max_S, PCA1_S)
  Positive_Outlier_c_R <- ifelse(PCA1_R > Max_R, Max_R, PCA1_R)

  Positive_Translation_coef_C <- abs(Min_C)
  Positive_Translation_coef_S <- abs(Min_S)
  Positive_Translation_coef_R <- abs(Min_R)

  C_positive_translation <- Positive_Outlier_c_C + Positive_Translation_coef_C

  S_positive_translation <- Positive_Outlier_c_S + Positive_Translation_coef_S

  R_positive_translation <- Positive_Outlier_c_R + Positive_Translation_coef_R

  Range_PCA2_C_positive_translation <- Max_C + (abs(Min_C))
  Range_PCA2_S_positive_translation <- Max_S + (abs(Min_S))
  Range_PCA2_R_positive_translation <- Max_R + (abs(Min_R))

  Proportion_of_total_variability_C <- (C_positive_translation/Range_PCA2_C_positive_translation)*100
  Proportion_of_total_variability_S <- (S_positive_translation/Range_PCA2_S_positive_translation)*100
  Proportion_of_total_variability_R <- 100-((R_positive_translation/Range_PCA2_R_positive_translation)*100)

  Percentage_Conversion_Coefficient <- 100/(Proportion_of_total_variability_C + Proportion_of_total_variability_S + Proportion_of_total_variability_R)

  Mat <- data.frame(C =c(90, 73, 73, 48, 54, 48, 42, 42, 23, 33, 23, 23, 23, 5,17, 5, 5, 5, 5),
                    S = c(5,	5,	23,	5,	23,	48,	17,	42,	5,	33,	73,	23,	54,	5, 42,	90,	23,	73,	48),
                    R =c(5,	23,	5,	48,	23,	5,	42,	17,	73,	33,	5,	54,	23,	90,	42,	5,	73,	23,	48))

  #Variances <- ($M5-BN5)^2 + ($N5-CG5)^2 + ($O5-CZ5)^2

  C <- Proportion_of_total_variability_C*Percentage_Conversion_Coefficient
  S <- Proportion_of_total_variability_S*Percentage_Conversion_Coefficient
  R <- Proportion_of_total_variability_R*Percentage_Conversion_Coefficient

  Red <- round((255/100*C),0)
  Green <- round((255/100*S),0)
  Blue <- round((255/100*R),0)

  Results <- data.frame(Species = Species,
                        leaf_area = leaf_area,
                        leaf_dry_mass = leaf_dry_mass,
                        leaf_fresh_mass = leaf_fresh_mass,
                        leaf_succulence_index = LSI,
                        leaf_water_content = LWC,
                        leaf_mass_per_area = LMA,
                        C = C,
                        R = R,
                        S = S,
                        Red = Red,
                        Green = Green,
                        Blue = Blue)

  Results <- Results[complete.cases(Results),]
  Results <- Results %>% dplyr::filter(C >= 0 & C <= 100,
                                       S >= 0 & S <= 100,
                                       R >= 0 & R <= 100)

  Results$RGB <- rgb(red = Results$Red, green = Results$Green, blue = Results$Blue, maxColorValue = 255)



  Results$Strategy <- NA

  for(i in 1:nrow(Results)){
    Sol <- dist(Strategy[,2:4] %>% bind_rows(Results[i,8:10])) %>% as.matrix()
    Results$Strategy[i] <- Strategy$Strategy[c(1:19)[min(Sol[-20,20]) == Sol[-20,20]]]
  }



  return(Results)
}


ForMedianSite <- PlotCRS %>%
  dplyr::select("plot", "leaf_area", "leaf_dry_mass", "leaf_fresh_mass") %>%
  dplyr::group_by(plot) %>%
  summarise(leaf_area = median(leaf_area), leaf_dry_mass = median(leaf_dry_mass), leaf_fresh_mass = median(leaf_fresh_mass), n = n()) %>%
  dplyr::filter(n > 3)

CRS_Per_Site <- CRS(leaf_area = ForMedianSite$leaf_area,
            leaf_dry_mass = ForMedianSite$leaf_dry_mass,
            leaf_fresh_mass = ForMedianSite$leaf_fresh_mass,
            Species = ForMedianSite$plot)


library(plotly)

CRS_Per_Site <- CRS_Per_Site %>% mutate(Label = paste(paste("plot:",CRS_Per_Site$Species), paste("Strategy:", CRS_Per_Site$Strategy), sep = "\n "))

# axis layout

axis <- function(title) {

  list(

    title = title,

    titlefont = list(

      size = 20

    ),

    tickfont = list(

      size = 15

    ),

    tickcolor = 'rgba(0,0,0,0)',

    ticklen = 5

  )

}


fig <- CRS_Per_Site %>% plot_ly()

fig <- fig %>% add_trace(

  type = 'scatterternary',

  mode = 'markers',

  a = ~C,

  b = ~R,

  c = ~S,

  text = ~Label,

  marker = list(

    symbol = 100,

    color = CRS_Per_Site$RGB,

    size = 14,

    line = list('width' = 2)

  )

)

fig <- fig %>% layout(

  title = "Simple Ternary Plot with Markers",

  ternary = list(

    sum = 100,

    aaxis = axis('C'),

    baxis = axis('R'),

    caxis = axis('S')

  )

)


fig


saveRDS(CRS_Per_Site, "CRS_Per_Site.rds")


PlotCRSBySpecies <- PlotCRS %>% dplyr::filter(plot %in% CRS_Per_Site$Species)

saveRDS(PlotCRSBySpecies, "PlotCRSBySpecies.rds")
