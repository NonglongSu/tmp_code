devtools::create("dotPlot")
devtools::document(pkg="./dotPlot")
devtools::load_all(path="./dotPlot")

#make bibary data available
x = 1:10
usethis::use_data(x)


#make raw data available
system.file('extdata','sample.json',package='dotPlot')
system.file('extdata','samples',package='dotPlot')

#make vignettes
#usethis::use_vignette("introduction")


