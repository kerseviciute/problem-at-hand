library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(limma)

snakemake <- readRDS('.features.Rmd.RDS')

features <- foreach(file = 'output/problem-at-hand/S2/run1/features.csv', .combine = rbind) %do% {
  sid <- gsub('.*/.*/.*/(.+)/features.csv', '\\1', file)
  fread(file) %>%
    .[ , Run := sid ]
} %>%
  .[ , FeatureID := paste(Channel, Frequency, sep = '_') ] %>%
  .[ , EventID := paste(EventID, Run, sep = '_') ] %>%
  .[ , Transformed := scale(log(Feature)) ] %>%
  .[ , Event := paste0('X', Event) ] %>%
  .[ Event == 'X0', Event := 'Control' ] %>%
  .[ Event == 'X2', Event := 'Thumb' ] %>%
  .[ Event == 'X3', Event := 'Index' ] %>%
  .[ Event == 'X4', Event := 'Middle' ] %>%
  .[ Event == 'X5', Event := 'Ring' ] %>%
  .[ Event == 'X6', Event := 'Pinky' ] %>%
  .[ , Event := factor(Event, levels = c('Control', 'Thumb', 'Index', 'Middle', 'Ring', 'Pinky')) ]

# Create feature matrix
x <- features %>% reshape2::dcast(FeatureID ~ EventID, value.var = 'Transformed')
featureIds <- x[ , 1 ]
x <- x[ , -1 ]
x <- apply(x, 2, as.numeric)
rownames(x) <- featureIds

x[ 1:5, 1:5 ]

# Model: Channel + Frequency (feature) ~ Event
uniqueEvents <- features[ , list(EventID, Event, Run) ] %>%
  unique()
modelData <- as.data.frame(uniqueEvents)
rownames(modelData) <- uniqueEvents[ , EventID ]

formula <- as.formula('~ 0 + Event')
design <- model.matrix(formula, data = droplevels(modelData))
head(design)

contrasts <- makeContrasts(contrasts = list(
  Finger1 = 'EventThumb - (EventControl + EventIndex + EventMiddle + EventRing + EventPinky) / 5'
  # Finger1 = 'EventIndex - (EventControl + EventThumb + EventMiddle + EventRing + EventPinky) / 5'
), levels = design)

colnames(contrasts) <- 'EventThumb'

fit <- lmFit(x, design, method = 'robust', maxit = 100)
fitc <- contrasts.fit(fit, contrasts)
fitc <- eBayes(fitc)

res <- foreach(contrast = colnames(contrasts), .combine = rbind) %dopar% {
  dt <- topTable(fitc, coef = contrast, number = Inf, sort.by = 'none')
  return(data.table(ID = rownames(x), as.data.frame(dt), Contrast = contrast))
}

id <- res[ order(P.Value) ][1][ , ID ]

res[ order(P.Value) ]
res[ adj.P.Val < 0.05 ]

features[ FeatureID == id ] %>%
  ggplot(aes(x = Event, y = Transformed)) +
  geom_boxplot() +
  ggtitle(id)
