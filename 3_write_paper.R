##################################################################################
source("F:/_SM/SMAO/Felles/autostat/Libraries/r_tools/r_tools.R")
req(c("sfsmisc","deSolve","data.table","boot","rmarkdown","knitr"))

setwd("F:/_SM/_SM-Felles/Utbrudd/Ebola 2014/Ebola_models_updated")

rmarkdown::render("test.Rmd","word_document")