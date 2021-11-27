jfile       = system.file('extdata','sample.json',package='dotPlot')
sample_json = jsonlite::fromJSON(jfile)
usethis::use_data(sample_json)