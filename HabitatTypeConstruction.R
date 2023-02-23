library(Artscore)
library(tidyverse)

Reference <- Artscore::Habitat_List |>
  dplyr::rename(habtype = Code)

Spatial <- readRDS("~/Documents/Danish-Flora-and-Vegetation.github.io/MapAndnmds/Spatial.rds") |>  as.data.frame() |>
  dplyr::select("habtype", "MajorHab") |>
  dplyr::mutate(habtype = as.character(habtype)) |>
  dplyr::distinct() |>
  left_join(Reference) |>
  dplyr::mutate(MajorHabName = NA) |>
  dplyr::arrange(habtype)


write_csv(Spatial, "Habitat_Types.csv")


DF <- structure(list(habtype = c("7220", "7230", "2130", "7110", "7140",
                           "6210", "6200", "6410", "1330", "2190", "1300", "1340", "6400",
                           "7200", NA, "2140", "2170", "2250", "4010", "4030", "9100", "4000",
                           "7100", "6230", "7150", "3100", "9999", "2120", "2110", "2180",
                           "6120", "1220", "7120", "2320", "9998", "7210", "2160", "6430",
                           "1100", "9120", "1310", "1210", "5130", "9150", "9110", "2100",
                           "2310", "9130", "9190", "3130", "7000", "1200", "3260", "2330",
                           "2300", "1320", "1230", "8220"), MajorHab = c("72", "72", "21",
                                                                         "71", "71", "62", "62", "64", "13", "21", "13", "13", "64", "72",
                                                                         NA, "21", "21", "22", "40", "40", "91", "40", "71", "62", "71",
                                                                         "31", "99", "21", "21", "21", "61", "12", "71", "23", "99", "72",
                                                                         "21", "64", "11", "91", "13", "12", "51", "91", "91", "21", "23",
                                                                         "91", "91", "31", "70", "12", "32", "23", "23", "13", "12", "82"
                           ), habitat_name = c("Kildevæld", "Rigkær", "Grå/grøn klit",
                                               "Højmose", "Hængesæk", "Kalkoverdrev", NA, "Tidvis våd eng",
                                               "Strandeng", "Klitlavning", NA, "Indlandssalteng", NA, NA, NA,
                                               "Klithede", NA, "Enebærklit", "Våd hede", "Tør hede", NA,
                                               NA, NA, "Surt overdrev", "Tørvelavning", NA, NA, NA, NA, "Skovklit",
                                               "Tørt kalkoverdrev", NA, "Nedbrudt højmose", "Revling-indlandsklit",
                                               NA, "Avneknippemose", NA, NA, NA, "Bøg på mor med kristtorn",
                                               NA, NA, "Enekrat", "Bøg på kalk", "Bøg på mor", NA, "Visse-indlandsklit",
                                               "Bøg på muld", "Stilkege-krat", NA, NA, NA, NA, "Græs-indlandsklit",
                                               NA, NA, NA, NA), MajorHabName = c(NA, NA, NA, NA, NA, NA, NA,
                                                                                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                                                 NA, NA, NA)), row.names = c(NA, -58L), class = "data.frame")


MajorHabitat <- Spatial |>
  dplyr::select(MajorHab, MajorHabName) |>
  dplyr::distinct() |>
  arrange(MajorHab)

MajorHabitat$MajorHabName <- c(NA, NA, "Strandeng", "Klit", "Klit", "Klit", )

  structure(list(MajorHab = c("11", "12", "13", "21", "22", "23",
                              "31", "32", "40", "51", "61", "62", "64", "70", "71", "72", "82",
                              "91", "99", NA), MajorHabName = c(NA, NA, NA, NA, NA, NA, NA,
                                                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)), row.names = c(NA,
                                                                                                                                    -20L), class = "data.frame")
