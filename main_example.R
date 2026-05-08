library(plm)
library(lmtest)
library(harmonicmeanp)

source("granger_utils_PEGranger.R")

#load example data
data <- read.csv("example_data.csv")
data <- pdata.frame(example_data, index = c("id", "time"))

#run the power-enhanced panel Granger causality test
result <- PE_PGCT(
  data = data,
  order = 1,
  index = c("id", "time")
)

print(result)
