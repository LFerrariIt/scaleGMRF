# code to obtain `nwEngland_sp` and nwEngland_adj_mat`

require(R2BayesX)
require(spBayesSurv)
require(sp)
require(spdep)

nw_england <- R2BayesX::read.bnd(
  system.file(
    "otherdata/nwengland.bnd",
    package = "spBayesSurv"
  )
)
nwEngland_sp <- R2BayesX::bnd2sp(nw_england)
nwEngland_adj_mat <- unname(spdep::nb2mat(spdep::poly2nb(nwEngland_sp), style = "B"))

usethis::use_data(nwEngland_sp, overwrite = TRUE)
usethis::use_data(nwEngland_adj_mat, overwrite = TRUE)
