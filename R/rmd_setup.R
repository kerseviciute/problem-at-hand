library(reticulate)
library(knitr)
reticulate::use_condaenv('problem-at-hand-mne')
knitr::knit_engines$set(python = reticulate::eng_python)
setOption <- options
setOption(knitr.graphics.rel_path = FALSE)
