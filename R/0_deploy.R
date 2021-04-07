# library(devtools)
# library(here)

### Run once:
#create_package(here())
#use_mit_license()




#Set up all functions as package functions
# source(here("R/0_gam_functions.R"))
# ls()

#"fit_RE_gam"       "gam_simul_CI"     "Mode"             "plot_gam_diff"    "predict_gam_diff" "predict_gam_int"
# use_r("fit_RE_gam")
# use_test()
#
# use_r("gam_simul_CI")
# use_test()
#
# use_r("Mode")
# use_test()
#
# use_r("plot_gam_diff")
# use_test()
#
# use_r("predict_gam_diff")
# use_test()
#
# use_r("predict_gam_int")
# use_test()
# use_r("predict_gam_emm")
# use_test()

# use_r("fit_coxph")
# use_test()

# use_r("predict_gam_HR")
# use_test()

# use_r("plot_sig_heatmap")
# use_test()



#Load functions (use in console)
#devtools::load_all()

#After running Code/insert_roxygen for all functions, run this:
#document()


#Declare package dependencies
# use_package("washb")
# use_package("dplyr")
# use_package("mgcv")
# use_package("MASS")
# use_package("scam")





#To install the package locally:
#install()


#To install off github:
#devtools::install_github("washb-eed-substudies/washbgam"")




