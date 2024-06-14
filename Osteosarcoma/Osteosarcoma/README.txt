Supplementary file for the paper “Ultra-high dimensional confounder selection algorithm comparison, with application
  to radiomics data” by I. Baldé, and D. Ghosh, for publication in Statistics in Medicine.

There are Rmarkdown documents, three R examples of code (SIS + GOAL, SIS + OAL and CBS), and osteosarcoma data, to help the reader reproduce similar results to the ones obtained in the paper with osteosarcoma study. 

Note that R packages "devtools", "lqa", "Ball", "glmnet", "grpreg", "MASS" and "cytominer" have to be installed in order to use this code. 



```{r}
# install.packages("devtools") # install "devtool" package first
library(devtools)
# install_github("cran/lqa")  # install "lqa" package from GitHub
library(lqa)
# install.packages("Ball")
library(Ball)
# install.packages("glmnet")
library(glmnet)
# install.packages("grpreg")
library(grpreg)
# install.packages("MASS")
library(MASS) # version 3.3.1
# install.packages("cytominer")
library(cytominer)
```
