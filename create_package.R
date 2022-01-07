library(tidyverse)
setwd("~/Dropbox (ASU)/tmp_code")

devtools::create("dotPlot")
devtools::document(pkg="./dotPlot")
<<<<<<< HEAD
devtools::load_all(path="./dotPlot",export_all=F)
=======
devtools::load_all(path="./dotPlot", export_all = FALSE)
>>>>>>> 74404e0ba7aa714287d2d02a4c0ec656c89aa184

#make bibary data available
x = 1:10
usethis::use_data(x)


#make raw data available
system.file('extdata','sample.json',package='dotPlot')
system.file('extdata','samples',package='dotPlot')

#make vignettes
#usethis::use_vignette("introduction")

#testing funcs
input1 = "dotPlot/inst/extdata/sample.json"
input2 = "dotPlot/inst/extdata/samples"

dotPlot::plot_dot_matrix(input1)
dotPlot::plot_dot_matrix(input2,TRUE)

dotPlot::generate_dot_matrix(input1)



