
# Interaction BRT --------------------------------------------------------------
set.seed(4) # For reproducible results
int.mod <- gbm(
  formula = outcome ~ .,
  distribution = "bernoulli",
  data = int.set[,c("outcome", "animal.group", trait.cols, nd.cols)],
  shrinkage = 0.01,
  n.trees = 12000,
  interaction.depth = 35,
  n.minobsinnode = 1,
  cv.folds = 10,
  n.cores = NULL, # default to all  cores
  verbose = T
)



# Gut Passage Time BRT ---------------------------------------------------------

set.seed(4)
gpt.mod <- gbm(
  formula = log.gpt ~ .,
  distribution = "gaussian",
  data = gpt.set[,c("log.gpt", "animal.group", trait.cols)],
  shrinkage = 0.01,
  n.minobsinnode = 2,
  n.trees = 1300,
  interaction.depth = 20,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)



# Movement BRT -----------------------------------------------------------------

set.seed(4)
move.mod <- gbm(
  formula = log.speed ~ .,
  distribution = "gaussian",
  data = move.set[,c("log.speed", "animal.group", move.trait.cols)],
  shrinkage = 0.005,
  n.trees = 1000,
  n.minobsinnode = 4,
  interaction.depth = 20,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)



# Correction BRT ---------------------------------------------------------------

# NOTE that the displacement correction values are contingent on the
# gpt.mod and move.mod analyses. Correction values are calculated using
# the following code (assuming gpt.mod and move.mod hadn't been run before).
# For reproducibility, I've already done these calculations for the version 
# of the displacement.set that is loaded in datasets.for.gbm.RData, so we don't
# need to run the following code.

# displacement.set$gpt.pred <- predict(object = gpt.mod,
#                                      newdata = displacement.set,
#                                      n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
#                                      type = "response") %>% exp()
# 
# displacement.set$move.pred <- predict(object = move.mod,
#                                       newdata = displacement.set,
#                                       n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
#                                       type = "response") %>% exp() * 1000 / 60
# 
# displacement.set$corr <- displacement.set$displacement / (displacement.set$gpt.pred * displacement.set$move.pred)
# 
# displacement.set$log.corr <- log(displacement.set$displacement) / log(displacement.set$gpt.pred * displacement.set$move.pred)
# displacement.set$log.corr <- log(displacement.set$displacement / (displacement.set$gpt.pred * displacement.set$move.pred))


set.seed(4)
corr.mod <- gbm(
  formula = log.corr ~ .,
  distribution = "gaussian",
  data = displacement.set[,c("log.corr", "animal.group", trait.cols)],
  shrinkage = 0.001,
  n.trees = 6000,
  interaction.depth = 15,
  n.minobsinnode = 1,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)



# Germination BRT --------------------------------------------------------------

set.seed(4)
germ.mod <- gbm(
  formula = yi ~ .,
  distribution = "gaussian",
  data = germ.set[,c("yi", "animal.group", trait.cols, germ.method.cols)],
  shrinkage = 0.02,
  n.trees = 1500,
  interaction.depth = 25,
  n.minobsinnode = 5,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)




# Introduced species predictions -----------------------------------------------

set.seed(4)
introd.mod <- gbm(
  formula = outcome ~ .,
  distribution = "bernoulli",
  data = introd.set[,c("outcome", "animal.group", trait.cols, nd.cols)], 
  shrinkage = 0.01, 
  n.minobsinnode = 10,
  n.trees = 12000,
  interaction.depth = 40, 
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)




# Participation BRT ------------------------------------------------------------

set.seed(4) # Plotting consistent results
particip.mod <- gbm(
  formula = outcome ~ .,
  distribution = "bernoulli",
  data = animal.particip.set[,c("outcome", "animal.group", animal.trait.cols)],
  shrinkage = 0.1,
  n.trees = 1200,
  interaction.depth = 40, # 30
  n.minobsinnode = 1,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)



# Normalized degree BRT --------------------------------------------------------

set.seed(4)
animal.nd.mod <- gbm(
  formula = animal.nd ~ .,
  distribution = "gaussian",
  data = animal.nd.set.no0[,c("animal.nd", "animal.group", animal.trait.cols)],
  shrinkage = 0.001,
  n.trees = 9000,
  interaction.depth = 25,
  n.minobsinnode = 1,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)



# Quantity z score BRT ---------------------------------------------------------

set.seed(4)
quant.mod <- gbm(
  formula = z.quant ~ .,
  distribution = "gaussian",
  data = quant.set[,c("z.quant", "animal.group", animal.trait.cols)],
  shrinkage = 0.005,
  n.trees = 2000,
  interaction.depth = 40,
  n.minobsinnode = 1,
  cv.folds = 10,
  n.cores = NULL,
  verbose = T
)


# Save model output
save(int.mod, gpt.mod, move.mod, corr.mod, germ.mod, 
     introd.mod,
     file = "data-and-model-outputs/gbm.mods1.RData")

save(particip.mod, animal.nd.mod, quant.mod,
     file = "data-and-model-outputs/gbm.mods2.RData")
