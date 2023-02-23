
PlotCommunity <- readRDS("~/Documents/Danish-Flora-and-Vegetation.github.io/MapAndnmds/PlotCommunity.rds")

Spatial <- readRDS("~/Documents/Danish-Flora-and-Vegetation.github.io/MapAndnmds/Spatial.rds")

ForModel <- Spatial |>
  left_join(PlotCommunity) |>
  nrow()
