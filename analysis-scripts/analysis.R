# Functions, load packages, and pull in data ----------------------------

# Source from functions.R
source("./analysis-scripts/functions.R")

packages <- c("tidyverse", "gbm", "caret", "BIOMOD", "plotrix",
              "base2grob", "maps", "grid", "RColorBrewer",
              "rasterVis")

ipak(packages)


# Source data (and model outputs for long-running analyses)
# This takes a little while

lapply(list.files("./data-and-model-outputs", full.names = T),
       load,
       .GlobalEnv)

int.set$outcome <- as.numeric(as.character(int.set$outcome))
animal.particip.set$outcome <- as.numeric(as.character(animal.particip.set$outcome))


# For analyses focused on interactions involving only native species 
# or including a non-native species, will get subsets of only these
introd.test.set <- int.set %>% 
  filter(plant.native.status == "introduced" | animal.native.status == "introduced")
introd.set <- int.set %>% 
  filter(plant.native.status == "native" & animal.native.status == "native")




# Table S1 dispersal process model results -------------------------------------

# Run boosted regression tree models. Because these may take days to 
# run, the model outputs are loaded in the code above from 
# the file "./data-and-model-outputs/gbm.mods.RData". Uncomment to run these.

# source("./analysis-scripts/run_gbm_models.R")


# First get a sample size matrix for each of these models
table.s1.n <- matrix(NA, ncol = 3, nrow = 6) %>% as.data.frame()
rownames(table.s1.n) <- c("Interaction", "Novel interaction", "Gut passage time",
                          "Movement speed", "Displacement correction",
                          "Germination effect")
colnames(table.s1.n) <- c("N", "N animal spp.", "N plant spp.")

n.fun <- function(set){
  c(set %>% dim() %>% .[1],
    set$animal.phylo.id %>% l.unique(),
    set$plant.phylo.id %>% l.unique())
}

table.s1.n[1:6, 1:3] <- sapply(list(int.set, introd.set, gpt.set, move.set, 
                                    displacement.set, germ.set), n.fun) %>% t()


### Interaction

# Int mod
int.preds <- predict(object = int.mod,
                     newdata = int.set,
                     n.trees = gbm.perf(int.mod, method = "cv", plot.it = F),
                     type = "response")

int.conf.mat <- confusionMatrix(factor(ifelse(int.preds>0.5,1,0)), 
                                factor(int.set[,int.mod$response.name]))

int.cv.conf.mat <- confusionMatrix(factor(ifelse(plogis(int.mod$cv.fitted)>0.5,1,0)), 
                                   factor(int.set[,int.mod$response.name]))

# Int based on mass
int.mass.mod <- glm(int.set[,int.mod$response.name] ~ int.set$body_mass_median * int.set$seed.dry.mass,
                    family = "binomial")

int.mass.preds <- predict(object = int.mass.mod,
                          newdata = int.set,
                          type = "response")
int.mass.conf.mat <- confusionMatrix(factor(ifelse(int.mass.preds>0.5,1,0)), 
                                     factor(int.set[,int.mod$response.name]))


# Interaction results for table
table.s1.int <- matrix(NA, ncol = 3, nrow = 4) %>% as.data.frame()
rownames(table.s1.int) <- c("AUC", "Accuracy", "Kappa", "TSS")
colnames(table.s1.int) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.int[,1] <- c(gbm.roc.area(int.set[,int.mod$response.name], int.preds),
                      int.conf.mat$overall[c("Accuracy", "Kappa")],
                      BIOMOD::TSS.Stat(int.conf.mat$table))

table.s1.int[,2] <- c(gbm.roc.area(int.set[,int.mod$response.name], plogis(int.mod$cv.fitted)),
                      int.cv.conf.mat$overall[c("Accuracy", "Kappa")],
                      BIOMOD::TSS.Stat(int.cv.conf.mat$table))

table.s1.int[,3] <- c(gbm.roc.area(int.set[,int.mod$response.name], int.mass.preds),
                      int.mass.conf.mat$overall[c("Accuracy", "Kappa")],
                      BIOMOD::TSS.Stat(int.mass.conf.mat$table))
table.s1.int[1:4, 1:3] <- sapply(table.s1.int, function(x) round(x, 3))





### Novel interactions

# Introd mod
introd.preds <- predict(object = introd.mod,
                        newdata = introd.set,
                        n.trees = gbm.perf(introd.mod, method = "cv", plot.it = F),
                        type = "response")

introd.test.preds <- predict(object = introd.mod,
                             newdata = introd.test.set,
                             n.trees = gbm.perf(introd.mod, method = "cv", plot.it = F),
                             type = "response")

introd.conf.mat <- confusionMatrix(factor(ifelse(introd.preds>0.5,1,0)), 
                                   factor(introd.set[,introd.mod$response.name]))

introd.test.conf.mat <- confusionMatrix(factor(ifelse(introd.test.preds>0.5,1,0)), 
                                        factor(introd.test.set[,introd.mod$response.name]))

# Introd based on mass
introd.mass.mod <- glm(introd.set[,introd.mod$response.name] ~ introd.set$body_mass_median * introd.set$seed.dry.mass,
                       family = "binomial")

introd.mass.preds <- predict(object = introd.mass.mod,
                             newdata = introd.set,
                             type = "response")
introd.mass.conf.mat <- confusionMatrix(factor(ifelse(introd.mass.preds>0.5,1,0)), 
                                        factor(introd.set[,introd.mod$response.name]))


# Interaction results for table
table.s1.introd <- matrix(NA, ncol = 3, nrow = 4) %>% as.data.frame()
rownames(table.s1.introd) <- c("AUC", "Accuracy", "Kappa", "TSS")
colnames(table.s1.introd) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.introd[,1] <- c(gbm.roc.area(introd.set[,introd.mod$response.name], introd.preds),
                         introd.conf.mat$overall[c("Accuracy", "Kappa")],
                         BIOMOD::TSS.Stat(introd.conf.mat$table))

table.s1.introd[,2] <- c(gbm.roc.area(introd.set[,introd.mod$response.name], plogis(introd.mod$cv.fitted)),
                         introd.test.conf.mat$overall[c("Accuracy", "Kappa")],
                         BIOMOD::TSS.Stat(introd.test.conf.mat$table))

table.s1.introd[,3] <- c(gbm.roc.area(introd.set[,introd.mod$response.name], introd.mass.preds),
                         introd.mass.conf.mat$overall[c("Accuracy", "Kappa")],
                         BIOMOD::TSS.Stat(introd.mass.conf.mat$table))
table.s1.introd[1:4, 1:3] <- sapply(table.s1.introd, function(x) round(x, 3))





### Gut passage time

# GPT mod
gpt.preds <- predict(object = gpt.mod,
                     newdata = gpt.set,
                     n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                     type = "response")

# GPT based on mass
gpt.mass.mod <- lm(gpt.set[,gpt.mod$response.name] ~ gpt.set$body_mass_median * gpt.set$seed.dry.mass)

gpt.mass.preds <- predict(object = gpt.mass.mod,
                          newdata = gpt.set,
                          type = "response")


# GPT results for table
table.s1.gpt <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s1.gpt) <- c("R2", "NRMSE")
colnames(table.s1.gpt) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.gpt[,1] <- c(cor(gpt.set[,gpt.mod$response.name], gpt.preds)^2,
                      hydroGOF::nrmse(gpt.preds, gpt.set[,gpt.mod$response.name]))

table.s1.gpt[,2] <- c(cor(gpt.set[,gpt.mod$response.name], gpt.mod$cv.fitted)^2,
                      hydroGOF::nrmse(gpt.mod$cv.fitted, gpt.set[,gpt.mod$response.name]))

table.s1.gpt[,3] <- c(cor(gpt.set[,gpt.mod$response.name], gpt.mass.preds)^2,
                      hydroGOF::nrmse(gpt.mass.preds, gpt.set[,gpt.mod$response.name]))
table.s1.gpt[1, 1:3] <- round(table.s1.gpt[1,], 3)





### Movement speed

# move mod
move.preds <- predict(object = move.mod,
                      newdata = move.set,
                      n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                      type = "response")

# move based on mass
move.mass.mod <- lm(move.set[,move.mod$response.name] ~ move.set$body_mass_median)

move.mass.preds <- predict(object = move.mass.mod,
                           newdata = move.set,
                           type = "response")


# move results for table
table.s1.move <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s1.move) <- c("R2", "NRMSE")
colnames(table.s1.move) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.move[,1] <- c(cor(move.set[,move.mod$response.name], move.preds)^2,
                       hydroGOF::nrmse(move.preds, move.set[,move.mod$response.name]))

table.s1.move[,2] <- c(cor(move.set[,move.mod$response.name], move.mod$cv.fitted)^2,
                       hydroGOF::nrmse(move.mod$cv.fitted, move.set[,move.mod$response.name]))

table.s1.move[,3] <- c(cor(move.set[,move.mod$response.name], move.mass.preds)^2,
                       hydroGOF::nrmse(move.mass.preds, move.set[,move.mod$response.name]))
table.s1.move[1, 1:3] <- round(table.s1.move[1,], 3)





### Displacement correction

# corr mod
corr.preds <- predict(object = corr.mod,
                      newdata = displacement.set,
                      n.trees = gbm.perf(corr.mod, method = "cv", plot.it = F),
                      type = "response")

# corr based on mass - first need to do some calculations if the
# run_gbm_models.R script isn't executed above
displacement.set$gpt.pred <- predict(object = gpt.mod,
                                     newdata = displacement.set,
                                     n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                                     type = "response") %>% exp()

displacement.set$move.pred <- predict(object = move.mod,
                                      newdata = displacement.set,
                                      n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                                      type = "response") %>% exp() * 1000 / 60

displacement.set$log.corr <- log(displacement.set$displacement / (displacement.set$gpt.pred * displacement.set$move.pred))
corr.mass.mod <- lm(displacement.set[,corr.mod$response.name] ~ displacement.set$body_mass_median * displacement.set$seed.dry.mass)

corr.mass.preds <- predict(object = corr.mass.mod,
                           newdata = displacement.set,
                           type = "response")


# corr results for table
table.s1.corr <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s1.corr) <- c("R2", "NRMSE")
colnames(table.s1.corr) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.corr[,1] <- c(cor(displacement.set[,corr.mod$response.name], corr.preds)^2,
                       hydroGOF::nrmse(corr.preds, displacement.set[,corr.mod$response.name]))

table.s1.corr[,2] <- c(cor(displacement.set[,corr.mod$response.name], corr.mod$cv.fitted)^2,
                       hydroGOF::nrmse(corr.mod$cv.fitted, displacement.set[,corr.mod$response.name]))

table.s1.corr[,3] <- c(cor(displacement.set[,corr.mod$response.name], corr.mass.preds)^2,
                       hydroGOF::nrmse(corr.mass.preds, displacement.set[,corr.mod$response.name]))
table.s1.corr[1, 1:3] <- round(table.s1.corr[1,], 3)





### Germination effect

# germ mod
germ.preds <- predict(object = germ.mod,
                      newdata = germ.set,
                      n.trees = gbm.perf(germ.mod, method = "cv", plot.it = F),
                      type = "response")

# germ based on mass
germ.mass.mod <- lm(germ.set[,germ.mod$response.name] ~ germ.set$body_mass_median * germ.set$seed.dry.mass)

germ.mass.preds <- predict(object = germ.mass.mod,
                           newdata = germ.set,
                           type = "response")


# germ results for table
table.s1.germ <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s1.germ) <- c("R2", "NRMSE")
colnames(table.s1.germ) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s1.germ[,1] <- c(cor(germ.set[,germ.mod$response.name], germ.preds)^2,
                       hydroGOF::nrmse(germ.preds, germ.set[,germ.mod$response.name]))

table.s1.germ[,2] <- c(cor(germ.set[,germ.mod$response.name], germ.mod$cv.fitted)^2,
                       hydroGOF::nrmse(germ.mod$cv.fitted, germ.set[,germ.mod$response.name]))

table.s1.germ[,3] <- c(cor(germ.set[,germ.mod$response.name], germ.mass.preds)^2,
                       hydroGOF::nrmse(germ.mass.preds, germ.set[,germ.mod$response.name]))
table.s1.germ[1, 1:3] <- round(table.s1.germ[1,], 3)



### Now put together this whole table

table.s1 <- matrix(NA, nrow =(4*2 + 2*4), ncol = 7) %>% as.data.frame()
colnames(table.s1) <- c(colnames(table.s1.n), "Metric", colnames(table.s1.int))

table.s1[c(1,5,9,11,13,15), 1:3] <- table.s1.n

table.s1[,5:7] <- rbind(table.s1.int, table.s1.introd,
                        table.s1.gpt, table.s1.move,
                        table.s1.corr, table.s1.germ)

table.s1[,4] <- list(table.s1.int, table.s1.introd,
                     table.s1.gpt, table.s1.move,
                     table.s1.corr, table.s1.germ) %>% lapply(rownames) %>% unlist()
table.s1[is.na(table.s1)] <- ""

table.s1 %>% write.csv("./figures/table.s1.csv")





# Figure S2 prediction summary -------------------------------------------------

# To visualize model results, will use two approaches. First, just plotting
# observed values versus predictions. Second, an example of these from an 
# example network (M_SD_005) in the form of a network diagram to give a sense
# for how these individual processes can be predicted for pairwise interactions

ex.set <- subset(int.set, net.id == "M_SD_005")
germ.added.cols <- c("whole fruit", "greenhouse/nursery soil", "yes") %>% t() %>% as.data.frame()
colnames(germ.added.cols) <- germ.method.cols

ex.set <- cbind(ex.set, germ.added.cols)

ex.set$compare.against <- factor(ex.set$compare.against, levels = levels(germ.set$compare.against)) # "mechanically cleaned" "whole fruit" 
ex.set$planting.location <- factor(ex.set$planting.location, levels = levels(germ.set$planting.location)) # "field soil" "greenhouse/nursery soil" "other" "petri dish (or similar)"
ex.set$feedtrial <- factor(ex.set$feedtrial, levels = levels(germ.set$feedtrial)) # "both" "no" "yes"


int.ex <- gpt.ex <- move.ex <- corr.ex <- germ.ex <- ex.set

int.ex$value <- predict(object = int.mod,
                        newdata = ex.set,
                        n.trees = gbm.perf(int.mod, method = "cv", plot.it = F),
                        type = "response")

gpt.ex$value <- predict(object = gpt.mod,
                        newdata = ex.set,
                        n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                        type = "response")

move.ex$value <- predict(object = move.mod,
                         newdata = ex.set,
                         n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                         type = "response")

corr.ex$value <- predict(object = corr.mod,
                         newdata = ex.set,
                         n.trees = gbm.perf(corr.mod, method = "cv", plot.it = F),
                         type = "response")

germ.ex$value <- predict(object = germ.mod,
                         newdata = ex.set,
                         n.trees = gbm.perf(germ.mod, method = "cv", plot.it = F),
                         type = "response")

order.net <- function(net.to.order, ref.net){
  
  rs <- rowSums(ref.net > 0)
  cs <- colSums(ref.net > 0)
  
  row.ordered <- rownames(ref.net)[order(rs, decreasing = T)]
  col.ordered <- colnames(ref.net)[order(cs, decreasing = T)]
  
  return(net.to.order[row.ordered, col.ordered])
  
}

par(mfrow=c(1,4))

ex.net <- net.spread(split.by = "net.id", 
                     split.vals = "M_SD_005",
                     tax.type = "phylo.id", 
                     data.type = "quant",
                     long.df = ex.set)[[1]]
ex.net <- order.net(ex.net, ex.net)

int.net <- net.spread(split.by = "net.id", 
                      split.vals = "M_SD_005",
                      tax.type = "phylo.id", 
                      data.type = "quant",
                      long.df = int.ex)[[1]]
int.net <- order.net(int.net, ex.net)

gpt.net <- net.spread(split.by = "net.id", 
                      split.vals = "M_SD_005",
                      tax.type = "phylo.id", 
                      data.type = "quant",
                      long.df = gpt.ex)[[1]] 
gpt.net <- order.net(gpt.net, ex.net)

move.net <- net.spread(split.by = "net.id", 
                       split.vals = "M_SD_005",
                       tax.type = "phylo.id", 
                       data.type = "quant",
                       long.df = move.ex)[[1]] 
move.net <- order.net(move.net, ex.net)

corr.net <- net.spread(split.by = "net.id", 
                       split.vals = "M_SD_005",
                       tax.type = "phylo.id", 
                       data.type = "quant",
                       long.df = corr.ex)[[1]] 
corr.net <- order.net(corr.net, ex.net)

germ.net <- net.spread(split.by = "net.id", 
                       split.vals = "M_SD_005",
                       tax.type = "phylo.id", 
                       data.type = "quant",
                       long.df = germ.ex)[[1]] 
germ.net <- order.net(germ.net, ex.net)



dev.off()
pdf("./figures/prediction.summary.pdf", width = 7.25, height = 3.5)
opmar <- par()$mar
par(mfrow = c(2,5))

mar.factor <- 1.3

obs.pred.cex <- 0.6

xlab.line.factor <- 1.25

r2.cex <- 0.5
r2.adj <- 0.1

xlab.line <- 2.6
ylab.line <- 3.25

par(mar = c(6.1,4.1, 4.1, 2.1)/mar.factor, cex = 0.5)


plant.rgb <- rgb(90, 180, 172, maxColorValue = 255)
animal.rgb <- rgb(216, 179, 101, maxColorValue = 255)

germ.rgb <- function(alpha) rgb(27, 158, 119, alpha*255, maxColorValue = 255)
corr.rgb <- function(alpha) rgb(217, 95, 2, alpha*255, maxColorValue = 255)
gpt.rgb <- function(alpha) rgb(117, 112, 179, alpha*255, maxColorValue = 255)
move.rgb <- function(alpha) rgb(231, 41, 138, alpha*255, maxColorValue = 255)

# Panel a
int.mat <- confusionMatrix(factor(ifelse(int.preds>0.5,"yes","no")),
                           factor(ifelse(int.set$outcome > 0, "yes", "no")))
my.fourfoldplot(int.mat$table, margin.text = F, 
                col = c(rgb(0.9,0.9,0.9),
                        rgb(0.5,0.5,0.5)))
ff.cex <- 0.9
text(0.45,0.45, "False\nnegative", cex = ff.cex)
text(-0.45,0.45, "True\nnegative", cex = ff.cex, col = "white")
text(0.45,-0.45, "True\npositive", cex = ff.cex, col = "white")
text(-0.45,-0.45, "False\npositive", cex = ff.cex)
put.fig.letter("a", cex = 1.5, font = 2, pos = 2)
#mtext("Predicted interactions", adj = 0, side = 3, line = lab.line)
put.fig.letter("      Interaction\n", pos = 4, cex = 1.25)

# Panel b
par(mar=c(5.1,6.1,4.1,1.1)/mar.factor)
plot(gpt.preds ~ gpt.set$log.gpt, las = 1,
     frame = F,
     pch = 16,
     #col = rgb(0.5,0.5,0.5),
     col = gpt.rgb(1),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
mtext("Observed gut passage\ntime (min.)",
      side = 1,
      line = xlab.line * xlab.line.factor,
      cex = obs.pred.cex)
mtext("Predicted gut passage\ntime (min.)",
      side = 2,
      line = ylab.line * 0.7,
      cex = obs.pred.cex)

axis(1, las = 1, at = log(c(10,100,1000,10000)), labels = c(10,parse(text='10^2'),parse(text='10^3'),parse(text='10^4')))
axis(2, las = 1, at = log(c(10,100,1000,10000)), labels = c(10,parse(text='10^2'),parse(text='10^3'),parse(text='10^4')))
gpt.lm <- lm(gpt.preds ~ gpt.set$log.gpt)

mtext(side = 3,
      line = -2,
      cex = r2.cex,
      adj = r2.adj,
      bquote(R^2 == .(round(summary(gpt.lm)$r.squared, 3))))
mtext(side = 3,
      line = -3,
      cex = r2.cex,
      adj = r2.adj,
      paste("n =", nrow(gpt.set)))

put.fig.letter("b", cex = 1.5, font = 2, pos = 2)
put.fig.letter("      Gut passage\n      time", pos = 4, cex = 1.25)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)


# Panel c
ax.vals <- c(6,36,216)
plot(move.preds ~ move.set$log.speed, las = 1,
     frame = F,
     xlim = log(c(4,216)),
     ylim = log(c(4,216)),
     pch = 16,
     #col = rgb(0.5,0.5,0.5),
     col = move.rgb(1),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
mtext("Observed speed (km/h)",
      side = 1,
      line = xlab.line,
      cex = obs.pred.cex)
mtext("Predicted speed (km/h)",
      side = 2,
      line = ylab.line,
      cex = obs.pred.cex)
move.lm <- lm(move.preds ~ move.set$log.speed)

mtext(side = 3,
      line = -2,
      cex = r2.cex,
      adj = r2.adj,
      bquote(R^2 == .(round(summary(move.lm)$r.squared, 3))))
mtext(side = 3,
      line = -3,
      cex = r2.cex,
      adj = r2.adj,
      paste("n =", nrow(move.set)))

axis(1, las = 1, at = log(ax.vals), labels = ax.vals)
axis(2, las = 1, at = log(ax.vals), labels = ax.vals)
put.fig.letter("c", cex = 1.5, font = 2, pos = 2)
put.fig.letter("      Movement\n", pos = 4, cex = 1.25)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)


# Panel d
ax.vals <- c(0.00005, 0.001, 0.02)
plot(corr.preds ~ displacement.set$log.corr, las = 1,
     frame = F,
     xlim = log(range(ax.vals)),
     ylim = log(c(0.00007, 0.01)),
     pch = 16,
     #col = rgb(0.5,0.5,0.5),
     col = corr.rgb(1),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
mtext("Observed value",
      side = 1,
      line = xlab.line,
      cex = obs.pred.cex)
mtext("Predicted value",
      side = 2,
      line = ylab.line,
      cex = obs.pred.cex)
axis(1, las = 1, at = log(ax.vals), labels = c("0.00005", "0.001", "0.02"))
axis(2, las = 1, at = log(c(0.0001, 0.001, 0.01)), labels = c("0.0001", "0.001", "0.01"))

corr.lm <- lm(corr.preds ~ displacement.set$log.corr)

mtext(side = 3,
      line = -2,
      cex = r2.cex,
      adj = r2.adj,
      bquote(R^2 == .(round(summary(corr.lm)$r.squared, 3))))
mtext(side = 3,
      line = -3,
      cex = r2.cex,
      adj = r2.adj,
      paste("n =", nrow(displacement.set)))

axis(1, las = 1, at = log(c(6,36,216)), labels = c(6,36,216))
axis(2, las = 1, at = log(c(6,36,216)), labels = c(6,36,216))
put.fig.letter("d", cex = 1.5, font = 2, pos = 2)
put.fig.letter("      Displacement\n      correction", pos = 4, cex = 1.25)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)


# Panel e
plot(germ.preds ~ germ.set$yi, las = 1,
     ylim = c(-10,10),
     frame = F,
     pch = 16,
     #col = rgb(0.5,0.5,0.5),
     col = germ.rgb(1),
     xlab = "",
     ylab = "",)
mtext("Observed germination\neffect (log-odds)",
      side = 1,
      line = xlab.line * xlab.line.factor,
      cex = obs.pred.cex)
mtext("Predicted germination\neffect (log-odds)",
      side = 2,
      line = ylab.line * 0.7,
      cex = obs.pred.cex)
germ.lm <- lm(germ.preds ~ germ.set$yi)

mtext(side = 3,
      line = -2,
      cex = r2.cex,
      adj = r2.adj,
      bquote(R^2 == .(round(summary(germ.lm)$r.squared, 3))))
mtext(side = 3,
      line = -3,
      cex = r2.cex,
      adj = r2.adj,
      paste("n =", nrow(germ.set)))

put.fig.letter("e", cex = 1.5, font = 2, pos = 2)
put.fig.letter("      Germination\n      effect", pos = 4, cex = 1.25)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)



# Panel f
par(mar = c(0.1,6.1,4.1,1.5)/mar.factor)

net.plot(int.net, grey.breaks = 6, show.net = int.net, row.text = T, 
         obs.rect.net = ex.net, leg.title = "Probability") #"Interaction\nprobability"

points(10.6,-11.5,
       pch = 15,
       cex = 1,
       col = "lightblue")
text(10.8,-11.8,
     cex = 0.8,
     "Intx.\nobs.",
     pos = 4)

lab.line <- 1.5

mtext("Example from", 
      side = 3, line = 3.25,
      adj = -0.8, cex = 0.55, font = 3, xpd = T)
mtext("Cordillera, Puerto Rico:", 
      side = 3, line = 2,
      adj = 1.5, cex = 0.55, font = 3, xpd = T)
put.fig.letter("f", cex = 1.5, font = 2, pos = 2)

# Panel g
net.plot(gpt.net, grey.breaks = 6, show.net = int.net, row.text = T, leg.title = "Minutes", exp.leg = T) #"Gut\npassage\ntime"
put.fig.letter("g", cex = 1.5, font = 2, pos = 2)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)

# Panel h
net.plot(move.net, grey.breaks = 6, show.net = int.net, row.text = T, leg.title = "km/h", exp.leg = T) #"Gut\npassage\ntime"
put.fig.letter("h", cex = 1.5, font = 2, pos = 2)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)

# Panel i
net.plot(exp(corr.net), grey.breaks = 6, show.net = int.net, row.text = T, leg.title = "Value", exp.leg = F, sig.digits = 1) #"Gut\npassage\ntime"
put.fig.letter("i", cex = 1.5, font = 2, pos = 2)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)

# Panel j
net.plot(germ.net, grey.breaks = 6, show.net = int.net, row.text = T, leg.title = "Effect size") #"Germination\neffect"
put.fig.letter("j", cex = 1.5, font = 2, pos = 2)
abline(v = grconvertX(0, from = "nfc", to = "user"), xpd = T)

par(mar = opmar)

dev.off()





# Figure 1 functional outcomes examples ----------------------------------------

# First, make a subset of int.set showing all observed pairwise species 
# interactions. Then use the gbm models to predict dispersal distance
# and the impact of gut passge for germination for each of these cases.

fig1.set <- int.set[which(int.set$value > 0),]
fig1.set$activity <- as.factor(fig1.set$activity)
fig1.set$volant <- as.factor(fig1.set$volant)

fig1.set$gpt.pred <- predict(object = gpt.mod,
                             newdata = fig1.set,
                             n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                             type = "response") %>% exp()

fig1.set$move.pred <- predict(object = move.mod,
                              newdata = fig1.set,
                              n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                              type = "response") %>% exp() * 1000 / 60

fig1.set$corr.pred <- predict(object = corr.mod,
                              newdata = fig1.set,
                              n.trees = gbm.perf(corr.mod, method = "cv", plot.it = F),
                              type = "response") %>% exp()

germ.added.cols <- c("whole fruit", "greenhouse/nursery soil", "yes") %>% t() %>% as.data.frame()
colnames(germ.added.cols) <- germ.method.cols
fig1.set <- cbind(fig1.set, germ.added.cols)
fig1.set$compare.against <- factor(fig1.set$compare.against, levels = levels(germ.set$compare.against)) # "mechanically cleaned" "whole fruit" 
fig1.set$planting.location <- factor(fig1.set$planting.location, levels = levels(germ.set$planting.location)) # "field soil" "greenhouse/nursery soil" "other" "petri dish (or similar)"
fig1.set$feedtrial <- factor(fig1.set$feedtrial, levels = levels(germ.set$feedtrial)) # "both" "no" "yes"
fig1.set$germ.pred <- predict(object = germ.mod,
                              newdata = fig1.set,
                              n.trees = gbm.perf(germ.mod, method = "cv", plot.it = F),
                              type = "response")
fig1.set$dist <- with(fig1.set, gpt.pred * move.pred * corr.pred)



# We will use photos from wikimedia to show examples

# These get pulled in above

# Artibeus phaeotis https://en.m.wikipedia.org/wiki/Pygmy_fruit-eating_bat#/media/File%3AArtibeus_phaeotis.jpg #https://upload.wikimedia.org/wikipedia/commons/1/14/Artibeus_phaeotis.jpg
# Aceros corrugatus https://commons.wikimedia.org/wiki/File:Wrinkled_Hornbill.png # https://upload.wikimedia.org/wikipedia/commons/9/93/Wrinkled_Hornbill.png
# Loxodonta africana https://simple.wikipedia.org/wiki/African_elephant#/media/File:African_elephant.jpg # https://upload.wikimedia.org/wikipedia/commons/e/ec/African_elephant.jpg

# Tangara heinei https://upload.wikimedia.org/wikipedia/commons/0/09/Black-capped_Tanager_%28f%29_JCB.jpg
# Sus scrofa https://upload.wikimedia.org/wikipedia/commons/d/d1/20160208054949%21Wildschein%2C_N%C3%A4he_Pulverstampftor_%28cropped%29.jpg
# Varecia variegata https://upload.wikimedia.org/wikipedia/commons/0/03/Black_and_white_ruffed_lemur.jpg
# Vulpes vulpes https://sc.wikipedia.org/wiki/File:R%C3%B8d_r%C3%A6v_(Vulpes_vulpes).jpg # https://upload.wikimedia.org/wikipedia/commons/e/ef/R%C3%B8d_r%C3%A6v_%28Vulpes_vulpes%29.jpg



dev.off()
op <- par()
pdf("./figures/Functional outcomes.pdf", width = 7.25, height = 3.5)
par(mfrow=c(1,2), xpd = T)
par(mar = c(3, 4, 2.5, 0.5))

pt.cex <- 0.3
cex.axis <- 0.8
cex.axis.lab <- 0.9
jitterx <- 50
set.seed(4)
plot(log(fig1.set$dist) ~ jitter(fig1.set$body_mass_median, jitterx), 
     #log = "xy",
     xpd = T,
     las = 1,
     ylim = log(c(5,4000)),
     frame = F,
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     pch = 16,
     col = "grey40",
     cex = pt.cex)
mtext("Body mass (kg)", side = 1, line = 2, cex = cex.axis.lab)
mtext("Predicted dispersal distance (m)", side = 2, line = 2.75, cex = cex.axis.lab)
axis(1, at = log(c(5,500,50000,5000000)), labels = c("0.005","0.5","50","5000"),
     padj = -0.5,
     cex.axis = cex.axis)
axis(2, at = log(c(10,200,4000)), 
     labels = c("10","200","4000"), 
     las = 1,
     cex.axis = cex.axis)
put.fig.letter("a", cex = 1, font = 2, pos = 2)

pic.width <- 3.25

pic.bat.x <- log(2000)
pic.bat.y <- log(15)
pic.hornbill.x <- log(50)
pic.hornbill.y <- log(3000)
pic.elephant.x <- log(1000000)
pic.elephant.y <- log(150)

pic.radius <- pic.width / 2
pointer.radius <- 0.15
pic.lwd <- 3

pointerfun("Artibeus phaeotis", pic.bat.x, pic.bat.y)
pointerfun("Aceros corrugatus", pic.hornbill.x, pic.hornbill.y)
pointerfun("Loxodonta africana", pic.elephant.x, pic.elephant.y)

addImg(bat, x = pic.bat.x, y = pic.bat.y, width = pic.width)
addImg(hornbill, x = pic.hornbill.x, y = pic.hornbill.y, width = pic.width)
addImg(elephant, x = pic.elephant.x, y = pic.elephant.y, width = pic.width)

spp.lab.text <- 0.75

text(pic.hornbill.x + 1.5, pic.hornbill.y + 0.25, "Aceros\ncorrugatus",
     font = 3, cex = spp.lab.text,
     pos = 4)

text(pic.bat.x + 1.5, pic.bat.y - 0.5, "Artibeus\nphaeotis",
     font = 3, cex = spp.lab.text,
     pos = 4)

text(pic.elephant.x - 1.5, pic.elephant.y - 0.5, "Loxodonta\nafricana",
     font = 3, cex = spp.lab.text,
     pos = 2)


plot(fig1.set$germ.pred ~ jitter(fig1.set$body_mass_median, jitterx), 
     #log = "x",
     xpd = T,
     las = 1,
     ylim = c(-5,9),
     frame = F,
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     pch = 16,
     col = "grey40",
     cex = pt.cex)
mtext("Body mass (kg)", side = 1, line = 2, cex = cex.axis.lab)
mtext("Predicted germination effect (log odds)", side = 2, line = 2.25, cex = cex.axis.lab)
axis(1, at = log(c(5,500,50000,5000000)), labels = c("0.005","0.5","50","5000"),
     padj = -0.5,
     cex.axis = cex.axis)
abline(h = 0, lty = 2, lwd = 0.7, xpd = F)
axis(2, at = c(-4,0,4,8), las = 1,
     cex.axis = cex.axis)
put.fig.letter("b", cex = 1, font = 2, pos = 2)

# Pig, Lemur, Fox, Tanager 
pic.tanager.x <- log(400)
pic.tanager.y <- 8.5
pic.fox.x <- log(25)
pic.fox.y <- -3.2
pic.lemur.x <- log(700000)
pic.lemur.y <- 7.2
pic.pig.x <- log(1000000)
pic.pig.y <- -3.1


pointerfun("Tangara heinei", pic.tanager.x, pic.tanager.y, response.var = "germ.pred", log.y = F, single.plant = "Cecropia telealba")
pointerfun("Vulpes vulpes", pic.fox.x, pic.fox.y, response.var = "germ.pred", log.y = F, single.plant = "Ficus carica")
pointerfun("Varecia variegata", pic.lemur.x, pic.lemur.y, response.var = "germ.pred", log.y = F, single.plant = "Cryptocarya crassifolia")
pointerfun("Sus scrofa", pic.pig.x, pic.pig.y, response.var = "germ.pred", log.y = F, single.plant = "Dipteryx alata")


addImg(tanager, x = pic.tanager.x, y = pic.tanager.y, width = pic.width)
addImg(fox, x = pic.fox.x, y = pic.fox.y, width = pic.width)
addImg(lemur, x = pic.lemur.x, y = pic.lemur.y, width = pic.width)
addImg(pig, x = pic.pig.x, y = pic.pig.y, width = pic.width)


text(pic.tanager.x + 1.5, pic.tanager.y + 0.25, "Tangara\nheinei",
     font = 3, cex = spp.lab.text,
     pos = 4)

text(pic.fox.x + 1.5, pic.fox.y - 0.5, "Vulpes\nvulpes",
     font = 3, cex = spp.lab.text,
     pos = 4)

text(pic.lemur.x - 1.5, pic.lemur.y - 0.5, "Varecia\nvariegata",
     font = 3, cex = spp.lab.text,
     pos = 2)

text(pic.pig.x - 1.5, pic.pig.y - 0.5, "Sus\nscrofa",
     font = 3, cex = spp.lab.text,
     pos = 2)


dev.off()
par(op)





# Figure 2 novel network assembly ----------------------------------------------

pdf("./figures/map examples.pdf", height = 6, width = 7.5)

sp.cex.val <- 0.7

par(mfrow = c(3, 4))

layout.mat <- c(9,2,3,9,
                1,2,3,4,
                1,2,3,4,
                1,2,3,4,
                1,9,9,4,
                5,9,9,8,
                5,9,9,8,
                5,9,9,8,
                5,6,7,8,
                9,6,7,9,
                9,6,7,9,
                9,6,7,9)

layout(matrix(layout.mat, ncol = 4, byrow = T))

par(mar = c(0,3,3,1), xpd = F)
nets.to.plot <- c("Vizentin-Bugoni 2019 WAI",
                  "M_SD_031",
                  "Ramaswami 2017",
                  "Koike 2008",
                  "Nogales 2017",
                  "Toledo 2018",
                  "Ganzhorn 1988",
                  "Garcia 2014 Puhi-Puhi River")

plot.focal.net(nets.to.plot[1], xlim.max = 14, ylim.max = 22, panel.lab = "a", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # Hawaii # A
plot.focal.net(nets.to.plot[2], xlim.max = 14, ylim.max = 22, panel.lab = "b", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # Azores # C
plot.focal.net(nets.to.plot[3], xlim.max = 14, ylim.max = 22, panel.lab = "c", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # India # E
plot.focal.net(nets.to.plot[4], xlim.max = 14, ylim.max = 22, panel.lab = "d", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # Japan # G

plot.focal.net(nets.to.plot[5], xlim.max = 14, ylim.max = 22, panel.lab = "e", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # Galapagos # B
plot.focal.net(nets.to.plot[6], xlim.max = 14, ylim.max = 22, panel.lab = "f", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # Brazil # D
plot.focal.net(nets.to.plot[7], xlim.max = 14, ylim.max = 22, panel.lab = "g", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # F
plot.focal.net(nets.to.plot[8], xlim.max = 14, ylim.max = 22, panel.lab = "h", 
               sp.cex = sp.cex.val, col.srt = 90, sp.abbr = T) # NZ # H

grob.map <- base2grob(function(){
  blank.map(col="grey50", cex = 0.8, add.box = F, map.col = "grey70")
  text(metanet[nets.to.plot, c("longitude","latitude")], 
       letters[1:8],#letters[c(1,3,5,7,2,4,6,8)], 
       cex = 0.9)
  
  prob.val <- signif(seq(from = 0, to = 1, 
                         length.out = 6), 1)#[-1]
  leg.x <- 40
  leg.y <- 105
  leg.x.spread <- 25
  points(x = leg.x+((0.5+(-3:2))*leg.x.spread), 
         y = rep(leg.y * 0.975, length(prob.val)),
         pch = 15,
         cex = 1.2,
         col = rgb(0,0,1, prob.val))
  points(x = leg.x+((0.5+(-3:2))*leg.x.spread), 
         y = rep(leg.y * 0.975, length(prob.val)),
         pch = 0,
         cex = 1.25,
         col = rgb(0,0,0,0.3))
  text(x = leg.x, y = leg.y * 1.05,
       "Predicted interaction probability", cex = 0.7)
  # text(x = c(leg.x - leg.x.spread * 2, leg.x + leg.x.spread * 2), 
  #      y = leg.y * 0.85,
  #      c("low", "high"), 
  #      cex = 0.7, 
  #      font = 3)
  text(x = leg.x+((0.5+(-3:2))*leg.x.spread), 
       y = rep(leg.y * 0.9, length(prob.val)),
       prob.val,
       cex = sp.cex.val)
  
  points(x = leg.x+((0.5+(3))*leg.x.spread), 
         y = leg.y * 0.975,
         pch = 15,
         cex = 0.6,
         col = "black")
  text(x = leg.x+((0.5+(3))*leg.x.spread), y = leg.y * 0.95,
       pos = 4,
       font = 3,
       "Interaction\nobserved", 
       cex = sp.cex.val)
  
  
} )
vp <- viewport(0.45,0.47, 
               width = 0.6*1, # width = 0.6*1.05, 
               height = 1) # height = 1.05)
pushViewport(vp)
grid.draw(grob.map)

dev.off()





# Figure S3 estimating the lognormal standard deviation ------------------------

my.sdlog <- 0.6

pdf(file = "./figures/Figure S3.pdf",
    useDingbats = F,
    height = 5,
    width = 5.5)

op <- par()

par(mar = c(5.1,5.1,4.1,2.1))

plot(log(max.displacement) ~ log(displacement), data = max.displacement.dat,
     pch = 1,
     frame = F,
     las = 1,
     xlab = "Mean seed dispersal distance (m)",
     ylab = "Maximum seed dispersal\ndistance (m)",
     xaxt = "n",
     yaxt = "n",
     xlim = log(c(5, 3125)),
     ylim = log(c(25, 15625)))

axis(1, at = log(c(5, 25, 125, 625, 3125)),
     labels = c(5, 25, 125, 625, 3125))

axis(2, at = log(c(25, 125, 625, 3125, 15625)),
     labels = c(25, 125, 625, 3125, 15625),
     las = 1)


lines((-2):10,
      log(qlnorm(0.99, 
                 meanlog = (-2):10,
                 sdlog = my.sdlog)),
      lwd = 3,
      lty = 2,
      col = "grey30")

text(x = log(15),
     y = log(300),
     "99th percentile\nof lognormal\ndistribution",
     #pos = 1,
     cex = 0.9,
     font = 3)

dev.off()




# Figure 3 impacts on long-distance dispersal ----------------------------------

# 1. Determine which species to include in analyses
# Skip marine mammals
marine.mamm.to.skip <- filter(all.mamm.traits, marine == 1)$name.m

not.terrestrial <- matrix(colSums(m.mamm.current[-which(rownames(m.mamm.current) %in% marine.mamm.to.skip),]),
                          nrow = 142,
                          ncol = 360, byrow = T) == 0

# Only focus on mammal orders including frugivores (note that this
# overlooks a role for rodents and other ecologically similar small mammals
# in dispersal of fleshy-fruited plants)
orders.bats <- "Chiroptera"
orders.primates <- "Primates"
orders.mamm.carns <- c("Carnivora", "Cingulata", "Dasyuromorphia", 
                       "Didelphimorphia", "Eulipotyphla", "Microbiotheria", 
                       "Monotremata", "Paucituberculata", "Pholidota", 
                       "Tubulidentata")
orders.mamm.herbs <- c("Artiodactyla", "Cetartiodactyla", "Dermoptera",  # Note that the traits dataframe and net.long use different order names...
                       "Diprotodontia", "Hyracoidea", "Litopterna", 
                       "Notoungulata", "Perissodactyla", "Pilosa",
                       "Proboscidea")

orders.mamm.to.skip <- filter(all.mamm.traits, !Order.1.2 %in% c(orders.bats, orders.primates,
                                                                 orders.mamm.carns, orders.mamm.herbs) )$name.m

spp.mamm.to.skip <- c(marine.mamm.to.skip, orders.mamm.to.skip, "Homo sapiens")

m.mamm.current.no.intro <- m.mamm.current.no.intro[-which(rownames(m.mamm.current.no.intro) %in% spp.mamm.to.skip), ]
m.mamm.current <- m.mamm.current[-which(rownames(m.mamm.current) %in% spp.mamm.to.skip), ]
m.mamm.pres.nat <- m.mamm.pres.nat[-which(rownames(m.mamm.pres.nat) %in% spp.mamm.to.skip), ]


# 2. Get median trait values for the representative plant species
median.plant.traits <- apply(int.set[,19:29], 2, median) %>% t() %>% as.data.frame()

# Also going to make a function to calculate germination probabilities
# based on the log odds ratios predicted by the germ.mod.
# We will assume a 'baseline' germination of 40% based on the average
# germination percentage in the gut passage germination studies. 

pgerm <- function(yi, prob = 0.4){
  (exp(yi) * prob) / (1 + exp(yi) * prob - prob)
}



# 3. Assign threshold for long-distance dispersal.

# We will use 1000 meters
ldd.dist <- 1000


# 4. Create dataframes to be used for prediction of dispersal processes
# for each bird and mammal species combination with the median plant

# These are called x.bird and x.mamm - they will give the values of each
# dispersal component based on the animal traits and median plant traits

# Want to include animal group in this
all.bird.traits$animal.group <- "birds"
all.mamm.traits$animal.group <- ifelse(all.mamm.traits$Order.1.2 %in% orders.bats, "bats",
                                       ifelse(all.mamm.traits$Order.1.2 %in% orders.mamm.carns, "mamm.carns",
                                              ifelse(all.mamm.traits$Order.1.2 %in% orders.mamm.herbs, "mamm.herbs",
                                                     ifelse(all.mamm.traits$Order.1.2 %in% orders.primates, "primates", "other"))))

x.bird <- tibble(name.m = rownames(m.bird.current))
x.bird <- left_join(x.bird, all.bird.traits)
x.bird[,colnames(median.plant.traits)] <- median.plant.traits

x.mamm <- tibble(name.m = rownames(m.mamm.current))
x.mamm <- left_join(x.mamm, all.mamm.traits[!duplicated(all.mamm.traits$phyla.bin),])
x.mamm[,colnames(median.plant.traits)] <- median.plant.traits

# Also need to make additional columns for a few effects
germ.added.cols <- c("whole fruit", "greenhouse/nursery soil", "yes") %>% t() %>% as.data.frame()
germ.method.cols <- c("compare.against", "planting.location", "feedtrial")
colnames(germ.added.cols) <- germ.method.cols
x.bird[,colnames(germ.added.cols)] <- germ.added.cols
x.mamm[,colnames(germ.added.cols)] <- germ.added.cols

x.bird$compare.against <- factor(x.bird$compare.against, levels = levels(germ.set$compare.against)) # "mechanically cleaned" "whole fruit" 
x.bird$planting.location <- factor(x.bird$planting.location, levels = levels(germ.set$planting.location)) # "field soil" "greenhouse/nursery soil" "other" "petri dish (or similar)"
x.bird$feedtrial <- factor(x.bird$feedtrial, levels = levels(germ.set$feedtrial)) # "both" "no" "yes"
x.bird$animal.group <- factor(x.bird$animal.group, levels = levels(int.set$animal.group))

x.mamm$compare.against <- factor(x.mamm$compare.against, levels = levels(germ.set$compare.against)) # "mechanically cleaned" "whole fruit" 
x.mamm$planting.location <- factor(x.mamm$planting.location, levels = levels(germ.set$planting.location)) # "field soil" "greenhouse/nursery soil" "other" "petri dish (or similar)"
x.mamm$feedtrial <- factor(x.mamm$feedtrial, levels = levels(germ.set$feedtrial)) # "both" "no" "yes"
x.mamm$animal.group <- factor(x.mamm$animal.group, levels = levels(int.set$animal.group))



# 5. Actually make the predictions for each model

particip.preds.bird <- predict(object = particip.mod,
                               newdata = x.bird,
                               n.trees = gbm.perf(particip.mod, method = "cv", plot.it = F),
                               type = "response")
particip.preds.mamm <- predict(object = particip.mod,
                               newdata = x.mamm,
                               n.trees = gbm.perf(particip.mod, method = "cv", plot.it = F),
                               type = "response")

nd.preds.bird <- predict(object = animal.nd.mod,
                         newdata = x.bird,
                         n.trees = gbm.perf(animal.nd.mod, method = "cv", plot.it = F),
                         type = "response")
nd.preds.mamm <- predict(object = animal.nd.mod,
                         newdata = x.mamm,
                         n.trees = gbm.perf(animal.nd.mod, method = "cv", plot.it = F),
                         type = "response")

quant.preds.bird <- predict(object = quant.mod,
                            newdata = x.bird,
                            n.trees = gbm.perf(quant.mod, method = "cv", plot.it = F),
                            type = "response")
quant.preds.mamm <- predict(object = quant.mod,
                            newdata = x.mamm,
                            n.trees = gbm.perf(quant.mod, method = "cv", plot.it = F),
                            type = "response")

# now need to add animal nd to the x.bird & x.mamm predictor database
x.bird$animal.nd <- nd.preds.bird
x.mamm$animal.nd <- nd.preds.mamm

x.bird$plant.nd <- 1 # This is a somewhat arbitrary number basically saying a generalist plant species
x.mamm$plant.nd <- 1


int.preds.bird <- predict(object = int.mod,
                          newdata = x.bird,
                          n.trees = gbm.perf(int.mod, method = "cv", plot.it = F),
                          type = "response")
int.preds.mamm <- predict(object = int.mod,
                          newdata = x.mamm,
                          n.trees = gbm.perf(int.mod, method = "cv", plot.it = F),
                          type = "response")

gpt.preds.bird <- predict(object = gpt.mod,
                          newdata = x.bird,
                          n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                          type = "response")
gpt.preds.mamm <- predict(object = gpt.mod,
                          newdata = x.mamm,
                          n.trees = gbm.perf(gpt.mod, method = "cv", plot.it = F),
                          type = "response")

move.preds.bird <- predict(object = move.mod,
                           newdata = x.bird,
                           n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                           type = "response")
move.preds.mamm <- predict(object = move.mod,
                           newdata = x.mamm,
                           n.trees = gbm.perf(move.mod, method = "cv", plot.it = F),
                           type = "response")

corr.preds.bird <- predict(object = corr.mod,
                           newdata = x.bird,
                           n.trees = gbm.perf(corr.mod, method = "cv", plot.it = F),
                           type = "response")
corr.preds.mamm <- predict(object = corr.mod,
                           newdata = x.mamm,
                           n.trees = gbm.perf(corr.mod, method = "cv", plot.it = F),
                           type = "response")

germ.preds.bird <- predict(object = germ.mod,
                           newdata = x.bird,
                           n.trees = gbm.perf(germ.mod, method = "cv", plot.it = F),
                           type = "response")
germ.preds.mamm <- predict(object = germ.mod,
                           newdata = x.mamm,
                           n.trees = gbm.perf(germ.mod, method = "cv", plot.it = F),
                           type = "response")




# 6. Put these components of the dispersal process together for LDD index

# Index of how many seeds dispersed
x.bird$mu <- particip.preds.bird * exp(quant.preds.bird) * int.preds.bird * pgerm(germ.preds.bird) 

# Portion of those seeds that fall greater than the long distance threshold
x.bird$p.ldd <- 1-plnorm(q = ldd.dist, meanlog = log(exp(gpt.preds.bird) * exp(move.preds.bird)  * 1000 / 60 * exp(corr.preds.bird)),
                         sdlog = my.sdlog)

# The product of these two is the long-distance dispersal index
x.bird$prod <- x.bird$p.ldd * x.bird$mu 


# The same for mammals...
x.mamm$mu <- particip.preds.mamm * exp(quant.preds.mamm) * int.preds.mamm * pgerm(germ.preds.mamm) # Index of how many seeds dispersed
x.mamm$p.ldd <- 1-plnorm(q = ldd.dist, meanlog = log(exp(gpt.preds.mamm) * exp(move.preds.mamm)  * 1000 / 60 * exp(corr.preds.mamm)),
                         sdlog = my.sdlog)
x.mamm$prod <- x.mamm$p.ldd * x.mamm$mu 





# 7. We will map ldd index onto a raster (saved as a matrix object)

# Will do this for each scenario, given the bird and mammal assemblages there

s.bird.current <- x.bird$prod * m.bird.current
s.bird.current <- colSums(s.bird.current)
s.bird.current <- matrix(s.bird.current,
                         nrow = 142,
                         ncol = 360, byrow = T)

s.bird.current.no.intro <- x.bird$prod * m.bird.current.no.intro
s.bird.current.no.intro <- colSums(s.bird.current.no.intro)
s.bird.current.no.intro <- matrix(s.bird.current.no.intro,
                                  nrow = 142,
                                  ncol = 360, byrow = T)

s.bird.pres.nat <- x.bird$prod * m.bird.pres.nat
s.bird.pres.nat <- colSums(s.bird.pres.nat)
s.bird.pres.nat <- matrix(s.bird.pres.nat,
                          nrow = 142,
                          ncol = 360, byrow = T)



s.mamm.current <- x.mamm$prod * m.mamm.current
s.mamm.current <- colSums(s.mamm.current)
s.mamm.current <- matrix(s.mamm.current,
                         nrow = 142,
                         ncol = 360, byrow = T)

s.mamm.pres.nat <- x.mamm$prod * m.mamm.pres.nat
s.mamm.pres.nat <- colSums(s.mamm.pres.nat)
s.mamm.pres.nat <- matrix(s.mamm.pres.nat,
                          nrow = 142,
                          ncol = 360, byrow = T)

s.mamm.current.no.intro <- x.mamm$prod * m.mamm.current.no.intro
s.mamm.current.no.intro <- colSums(s.mamm.current.no.intro)
s.mamm.current.no.intro <- matrix(s.mamm.current.no.intro,
                                  nrow = 142,
                                  ncol = 360, byrow = T)


# 7a. Prep for panels b-d

# i. We want to first get the scenario where we no longer have 
# vulnerable and endangered species.

te.categories <- c("CR", "EN", "EP", "EW", "EX", "VU", "CR (PE)")

# mammals
te.mamm.spp <- filter(all.mamm.traits, iucn.category %in% te.categories)$name.m
m.mamm.no.te <- m.mamm.current
s.mamm.no.te <- x.mamm$prod * m.mamm.no.te
s.mamm.no.te <- s.mamm.no.te[!rownames(s.mamm.no.te) %in% te.mamm.spp,]
s.mamm.no.te <- colSums(s.mamm.no.te)
s.mamm.no.te <- matrix(s.mamm.no.te,
                       nrow = 142,
                       ncol = 360, byrow = T)

# birds
te.bird.spp <- filter(all.bird.traits, iucn.category %in% te.categories)$name.m
m.bird.no.te <- m.bird.current
s.bird.no.te <- x.bird$prod * m.bird.no.te
s.bird.no.te <- s.bird.no.te[!rownames(s.bird.no.te) %in% te.bird.spp,]
s.bird.no.te <- colSums(s.bird.no.te)
s.bird.no.te <- matrix(s.bird.no.te,
                       nrow = 142,
                       ncol = 360, byrow = T)

# Calculate total associated with both
s.both.current.no.intro <- s.mamm.current.no.intro + s.bird.current.no.intro
s.both.current <- s.mamm.current + s.bird.current
s.both.pres.nat <- s.mamm.pres.nat + s.bird.pres.nat
s.both.no.te <- s.mamm.no.te + s.bird.no.te




# 7b. Make panel a-d

pdf(file = "./figures/both.curr.pdf",
    useDingbats = F,
    height = 5.75,
    width = 4)

colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))

s2 <- s.both.current
s2[not.terrestrial] <- NA

s.both.stack <- stack(raster(s2), raster(s2), raster(s2))

lo.val <- 0
hi.val <- 6#max(maxValue(s.both.stack))

gr.val <- seq(0, round(hi.val), len = 5)

rasterVis::levelplot(s.both.stack, 
                     margin=FALSE,                       # suppress marginal graphics
                     colorkey=list(
                       space='bottom',                   # plot legend at bottom
                       labels=list(at=gr.val)      # legend ticks and labels 
                     ),    
                     par.settings=list(
                       panel.background=list(col="grey40"),
                       strip.border=list(col='transparent'),
                       strip.background=list(col='transparent'),
                       axis.line=list(col='transparent') # suppress axes and legend outline
                     ),
                     scales=list(draw=FALSE),            # suppress axis labels
                     col.regions=colr,                   # colour ramp
                     at=seq(lo.val, hi.val, len=101),
                     names.attr=rep('', nlayers(s.both.stack))) %>% 
  update(., aspect = 142/360)
lab.cex <- 0.75
lab.x <- 0.095#0.072
grid.text("Current", x = unit(lab.x,"npc"), y = unit(0.67,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

grid.text("a)", x = unit(lab.x,"npc"), y = unit(0.845,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))



dev.off()





pdf(file = "./figures/change from current three panel.pdf",
    useDingbats = F,
    height = 5.75,
    width = 4)

s.curr <- log(s.both.pres.nat / s.both.current)
s.curr[not.terrestrial] <- NA

s.no.te <- log(s.both.no.te / s.both.current)
s.no.te[not.terrestrial] <- NA

s.no.int <- log(s.both.current.no.intro / s.both.current)
s.no.int[not.terrestrial] <- NA



lo.cut <- log(0.125) # 0.15 # used to be 0.25
hi.cut <- log(8) # 6 # used to be 4

s.curr[s.curr < lo.cut] <- lo.cut
s.curr[s.curr > hi.cut] <- hi.cut

s.no.te[s.no.te < lo.cut] <- lo.cut
s.no.te[s.no.te > hi.cut] <- hi.cut

s.no.int[s.no.int < lo.cut] <- lo.cut
s.no.int[s.no.int > hi.cut] <- hi.cut


ratio.stack <- stack(raster(s.curr), raster(s.no.te), raster(s.no.int))

lo.val <- lo.cut
hi.val <- hi.cut

colr <- colorRampPalette(colors = c("#67001F",
                                    "#B2182B",
                                    "#F4A582",
                                    "#F7F7F7",
                                    "#92C5DE",
                                    "#2166AC",
                                    "#053061"))

gr.vals <- c(0.15, 0.5, 1, 2, 6) # used to be 0.25 and 4
gr.labs <- c("-85%", "-50%", "0%", "100%", "500%") # Used to be -75 and 300

rasterVis::levelplot(ratio.stack, 
                     margin=FALSE,                       # suppress marginal graphics
                     colorkey=list(
                       space='bottom',                   # plot legend at bottom
                       labels=list(at=log(gr.vals),
                                   labels = gr.labs)      # legend ticks and labels 
                     ),    
                     par.settings=list(
                       panel.background=list(col="grey40"),
                       strip.border=list(col='transparent'),
                       strip.background=list(col='transparent'),
                       axis.line=list(col='transparent') # suppress axes and legend outline
                     ),
                     scales=list(draw=FALSE),            # suppress axis labels
                     col.regions=colr,                   # colour ramp
                     at=seq(lo.val, hi.val, len=101),
                     names.attr=rep('', nlayers(ratio.stack))) %>% 
  update(., aspect = 142/360)
lab.cex <- 0.75
#lab.x <- 0.072
grid.text("Natural", x = unit(lab.x,"npc"), y = unit(0.67,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

grid.text("Extinction of\nendangered spp", x = unit(lab.x,"npc"), y = unit(0.395,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

grid.text("Extirpation of\nintroduced spp", x = unit(lab.x,"npc"), y = unit(0.12,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))


grid.text("b)", x = unit(lab.x,"npc"), y = unit(0.845,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

grid.text("c)", x = unit(lab.x,"npc"), y = unit(0.565,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

grid.text("d)", x = unit(lab.x,"npc"), y = unit(0.285,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))

dev.off()






# 8. Prep for panels e and f

ecoregion.df <- tibble(ecoregion.name = lapply(ecoregion.stack, names) %>% unlist(),
                       richness.pres.nat = NA,
                       richness.current = NA,
                       function.pres.nat = NA,
                       function.current = NA,
                       pixels = NA,
                       
                       mamm.richness.pres.nat = NA,
                       mamm.richness.current = NA,
                       mamm.function.pres.nat = NA,
                       mamm.function.current = NA,
                       
                       bird.richness.pres.nat = NA,
                       bird.richness.current = NA,
                       bird.function.pres.nat = NA,
                       bird.function.current = NA
)



rich.both.current <- colSums(rbind(m.bird.current, m.mamm.current))
rich.both.current <- matrix(rich.both.current,
                            nrow = 142,
                            ncol = 360, byrow = T)

rich.both.pres.nat <- colSums(rbind(m.bird.pres.nat, m.mamm.pres.nat))
rich.both.pres.nat <- matrix(rich.both.pres.nat,
                             nrow = 142,
                             ncol = 360, byrow = T)

rich.mamm.current <- colSums(m.mamm.current)
rich.mamm.current <- matrix(rich.mamm.current,
                            nrow = 142,
                            ncol = 360, byrow = T)

rich.mamm.pres.nat <- colSums(m.mamm.pres.nat)
rich.mamm.pres.nat <- matrix(rich.mamm.pres.nat,
                             nrow = 142,
                             ncol = 360, byrow = T)


rich.bird.current <- colSums(m.bird.current)
rich.bird.current <- matrix(rich.bird.current,
                            nrow = 142,
                            ncol = 360, byrow = T)

rich.bird.pres.nat <- colSums(m.bird.pres.nat)
rich.bird.pres.nat <- matrix(rich.bird.pres.nat,
                             nrow = 142,
                             ncol = 360, byrow = T)

# This takes quite a while - there must be something better

for(i in 1:length(ecoregion.stack)){
  ecoregion.set <- ecoregion.stack[[i]] %>% as.matrix(.) %>% is.na(.) %>% !.
  
  ecoregion.df$pixels[i] <- sum(ecoregion.set)
  
  ecoregion.df$mamm.richness.pres.nat[i] <- rich.mamm.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$mamm.richness.current[i] <- rich.mamm.current[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$mamm.function.pres.nat[i] <- s.mamm.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$mamm.function.current[i] <- s.mamm.current[ecoregion.set] %>% mean(., na.rm=T)
  
  ecoregion.df$bird.richness.pres.nat[i] <- rich.bird.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$bird.richness.current[i] <- rich.bird.current[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$bird.function.pres.nat[i] <- s.bird.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$bird.function.current[i] <- s.bird.current[ecoregion.set] %>% mean(., na.rm=T)
  
  
  ecoregion.df$richness.pres.nat[i] <- rich.both.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$richness.current[i] <- rich.both.current[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$function.pres.nat[i] <- s.both.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$function.current[i] <- s.both.current[ecoregion.set] %>% mean(., na.rm=T)
  
  print(paste(i, ecoregion.df$ecoregion.name[i]))
  
}


ecoregion.df <- ecoregion.df %>% 
  mutate(mamm.richness.loss = (mamm.richness.current - mamm.richness.pres.nat) / mamm.richness.pres.nat,
         bird.richness.loss = (bird.richness.current - bird.richness.pres.nat) / bird.richness.pres.nat,
         richness.loss = (richness.current - richness.pres.nat) / richness.pres.nat,
         function.loss = (function.current - function.pres.nat) / function.pres.nat)


fd.ecoreg <- tibble(ecoregion.name = names(fd.output$FDis)[grep("current", names(fd.output$FDis))] %>% gsub("current.", "", ., fixed = T),
                    FDis.current = fd.output$FDis[grep("current", names(fd.output$FDis))],
                    FDis.pres.nat = fd.output$FDis[grep("pres.nat", names(fd.output$FDis))])

fd.ecoreg <- fd.ecoreg %>% mutate(FDis.loss = (FDis.current - FDis.pres.nat) / FDis.pres.nat)

ecoregion.df <- left_join(ecoregion.df, fd.ecoreg)

ecoregion.df1 <- filter(ecoregion.df, pixels > 10)

# Analysis

richness.sma <- smatr::sma(log(function.loss + 1) ~ log(richness.loss + 1), 
                           data = ecoregion.df1)

fd.sma <- smatr::sma(log(function.loss + 1) ~ log(FDis.loss + 1), 
                           data = ecoregion.df1)



pdf(file = "./figures/function decline two panel by ecoregion.pdf",
    useDingbats = F,
    height = 3.25,
    width = 4)

y.max.loss <- -1
y.min.loss <- 0.25

axis.cex <- 0.8
par(mfrow = c(1, 2),
    mar = c(4,4,1,0))
plot(function.loss ~ richness.loss, data = ecoregion.df1,
     #asp = 1,
     pch = 16,
     frame = F,
     las = 1,
     cex = 0.8,
     col = "grey40",
     #col = rgb(0,0,0,0.5),
     ylim = c(y.min.loss+0.05, y.max.loss),
     xlim = c(0+0.05, -0.3),
     xlab = "Change in\nspp. richness",
     ylab = "Change in long-distance dispersal",
     xaxt = "n",
     yaxt = "n",
     xpd = T)
curve(x * 1, add = T, lty = 3)
axis(1, at = c(seq(0,-0.3, len = 3)),
     labels = paste(seq(0,-0.3, len = 3)*100, "%", sep = ""),
     cex.axis = axis.cex)
axis(2, at = c(seq(y.min.loss,y.max.loss, len = 6)), las = 1,
     labels = paste(seq(y.min.loss,y.max.loss, len = 6)*100, "%", sep = ""),
     cex.axis = axis.cex)

abline(h = 0, lty = 2)

curve(exp(coef(richness.sma)[1] + coef(richness.sma)[2] * x) - 1, add = T, lwd = 2,
      from = 0, to = -0.13)


text(x = 0.08,
     y = -0.96,
     cex = 0.9,
     pos = 4,
     "e)")

text(x = -0.3,
     y = -0.35,
     cex = 0.75,
     pos = 2,
     font = 3,
     "1:1 line")

text(x = -0.33,
     y = -0.97,
     cex = 0.75,
     pos = 2,
     font = 3,
     "fitted line")



par(mar = c(4,3,1,1))

plot(function.loss ~ FDis.loss, data = ecoregion.df1,
     #asp = 1,
     pch = 16,
     frame = F,
     las = 1,
     cex = 0.8,
     col = "grey40",
     #col = rgb(0,0,0,0.5),
     ylim = c(y.min.loss+0.05, y.max.loss),
     xlim = c(0+0.05,-0.3),
     xlab = "Change in\nfunctional dispersion",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
curve(x * 1, add = T, lty = 3)
axis(1, at = c(seq(0,-0.3, len = 3)),
     labels = paste(seq(0,-0.3, len = 3)*100, "%", sep = ""),
     cex.axis = axis.cex)
axis(2, at = c(seq(y.min.loss,y.max.loss, len = 6)), las = 1,
     labels = rep("",6),
     cex.axis = axis.cex)

abline(h = 0, lty = 2)

curve(exp(coef(fd.sma)[1] + coef(fd.sma)[2] * x) - 1, add = T, lwd = 2,
      from = 0.018, to = -0.2)



text(x = 0.08,
     y = -0.96,
     cex = 0.9,
     pos = 4,
     "f)")

dev.off()


















# Table S2 animal role model results -------------------------------------------


# First get a sample size matrix for each of these models
table.s2.n <- matrix(NA, ncol = 3, nrow = 3) %>% as.data.frame()
rownames(table.s2.n) <- c("Participation", "Normalized degree", "Quantity")
colnames(table.s2.n) <- c("N", "N animal spp.", "N plant spp.")

n.fun <- function(set){
  c(set %>% dim() %>% .[1],
    set$animal.phylo.id %>% l.unique(),
    set$plant.phylo.id %>% l.unique())
}

table.s2.n[1:3, 1:3] <- sapply(list(animal.particip.set, animal.nd.set.no0, quant.set), n.fun) %>% t()



### Participation

# particip mod
particip.preds <- predict(object = particip.mod,
                          newdata = animal.particip.set,
                          n.trees = gbm.perf(particip.mod, method = "cv", plot.it = F),
                          type = "response")

particip.conf.mat <- confusionMatrix(factor(ifelse(particip.preds>0.5,1,0)), 
                                     factor(animal.particip.set[,particip.mod$response.name]))

particip.cv.conf.mat <- confusionMatrix(factor(ifelse(plogis(particip.mod$cv.fitted)>0.5,1,0)), 
                                        factor(animal.particip.set[,particip.mod$response.name]))

# particip based on mass
particip.mass.mod <- glm(animal.particip.set[,particip.mod$response.name] ~ animal.particip.set$body_mass_median,
                         family = "binomial")

particip.mass.preds <- predict(object = particip.mass.mod,
                               newdata = animal.particip.set,
                               type = "response")
particip.mass.conf.mat <- confusionMatrix(factor(ifelse(particip.mass.preds>0.5,1,0)), 
                                          factor(animal.particip.set[,particip.mod$response.name]))


# participeraction results for table
table.s2.particip <- matrix(NA, ncol = 3, nrow = 4) %>% as.data.frame()
rownames(table.s2.particip) <- c("AUC", "Accuracy", "Kappa", "TSS")
colnames(table.s2.particip) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s2.particip[,1] <- c(gbm.roc.area(animal.particip.set[,particip.mod$response.name], particip.preds),
                           particip.conf.mat$overall[c("Accuracy", "Kappa")],
                           BIOMOD::TSS.Stat(particip.conf.mat$table))

table.s2.particip[,2] <- c(gbm.roc.area(animal.particip.set[,particip.mod$response.name], plogis(particip.mod$cv.fitted)),
                           particip.cv.conf.mat$overall[c("Accuracy", "Kappa")],
                           BIOMOD::TSS.Stat(particip.cv.conf.mat$table))

table.s2.particip[,3] <- c(gbm.roc.area(animal.particip.set[,particip.mod$response.name], particip.mass.preds),
                           particip.mass.conf.mat$overall[c("Accuracy", "Kappa")],
                           BIOMOD::TSS.Stat(particip.mass.conf.mat$table))
table.s2.particip[1:4, 1:3] <- sapply(table.s2.particip, function(x) round(x, 3))






### ND

# animal.nd mod
animal.nd.preds <- predict(object = animal.nd.mod,
                           newdata = animal.nd.set.no0,
                           n.trees = gbm.perf(animal.nd.mod, method = "cv", plot.it = F),
                           type = "response")

# animal.nd based on mass
animal.nd.mass.mod <- lm(animal.nd.set.no0[,animal.nd.mod$response.name] ~ animal.nd.set.no0$body_mass_median)

animal.nd.mass.preds <- predict(object = animal.nd.mass.mod,
                                newdata = animal.nd.set.no0,
                                type = "response")


# animal.nd results for table
table.s2.animal.nd <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s2.animal.nd) <- c("R2", "NRMSE")
colnames(table.s2.animal.nd) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s2.animal.nd[,1] <- c(cor(animal.nd.set.no0[,animal.nd.mod$response.name], animal.nd.preds)^2,
                            hydroGOF::nrmse(animal.nd.preds, animal.nd.set.no0[,animal.nd.mod$response.name]))

table.s2.animal.nd[,2] <- c(cor(animal.nd.set.no0[,animal.nd.mod$response.name], animal.nd.mod$cv.fitted)^2,
                            hydroGOF::nrmse(animal.nd.mod$cv.fitted, animal.nd.set.no0[,animal.nd.mod$response.name]))

table.s2.animal.nd[,3] <- c(cor(animal.nd.set.no0[,animal.nd.mod$response.name], animal.nd.mass.preds)^2,
                            hydroGOF::nrmse(animal.nd.mass.preds, animal.nd.set.no0[,animal.nd.mod$response.name]))
table.s2.animal.nd[1, 1:3] <- round(table.s2.animal.nd[1,], 3)





### Seed disperal quantity

# quant mod
quant.preds <- predict(object = quant.mod,
                       newdata = quant.set,
                       n.trees = gbm.perf(quant.mod, method = "cv", plot.it = F),
                       type = "response")

# quant based on mass
quant.mass.mod <- lm(quant.set[,quant.mod$response.name] ~ quant.set$body_mass_median)

quant.mass.preds <- predict(object = quant.mass.mod,
                            newdata = quant.set,
                            type = "response")


# quant results for table
table.s2.quant <- matrix(NA, ncol = 3, nrow = 2) %>% as.data.frame()
rownames(table.s2.quant) <- c("R2", "NRMSE")
colnames(table.s2.quant) <- c("BRT value", "BRT cv", "Size matching GLM value")
table.s2.quant[,1] <- c(cor(quant.set[,quant.mod$response.name], quant.preds)^2,
                        hydroGOF::nrmse(quant.preds, quant.set[,quant.mod$response.name]))

table.s2.quant[,2] <- c(cor(quant.set[,quant.mod$response.name], quant.mod$cv.fitted)^2,
                        hydroGOF::nrmse(quant.mod$cv.fitted, quant.set[,quant.mod$response.name]))

table.s2.quant[,3] <- c(cor(quant.set[,quant.mod$response.name], quant.mass.preds)^2,
                        hydroGOF::nrmse(quant.mass.preds, quant.set[,quant.mod$response.name]))
table.s2.quant[1, 1:3] <- round(table.s2.quant[1,], 3)



### Now put together this whole table

table.s2 <- matrix(NA, nrow =(4*1 + 2*2), ncol = 7) %>% as.data.frame()
colnames(table.s2) <- c(colnames(table.s2.n), "Metric", colnames(table.s2.particip))


table.s2[c(1,5,7), 1:3] <- table.s2.n

table.s2[,5:7] <- rbind(table.s2.particip, 
                        table.s2.animal.nd, table.s2.quant)

table.s2[,4] <- list(table.s2.particip,
                     table.s2.animal.nd, table.s2.quant) %>% lapply(rownames) %>% unlist()
table.s2[is.na(table.s2)] <- ""

table.s2 %>% write.csv("./figures/table.s2.csv")





# Figure 4 climate outpacing dispersal -----------------------------------------

vclimmean <- calc(vclim, function(x){mean(x)})


# Will use the matrix representations here and will need one for the 
# scenario where threatened and endangered species are absent
m.bird.te0 <- m.bird.current
m.bird.te0[rownames(m.bird.te0) %in% te.bird.spp,] <- 0

m.mamm.te0 <- m.mamm.current
m.mamm.te0[rownames(m.mamm.te0) %in% te.mamm.spp,] <- 0




bird.mu <- particip.preds.bird * exp(quant.preds.bird) * int.preds.bird * pgerm(germ.preds.bird) # Index of how many seeds dispersed

bird.cv.current <- lapply(1:length(bird.mu), 
                          function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                                 meanlog = (log(exp(gpt.preds.bird) * exp(move.preds.bird)  * 1000 / 60 * exp(corr.preds.bird))[x]),
                                                 sdlog = my.sdlog)) * bird.mu[x] * matrix(m.bird.current[x,], nrow=142, byrow = T)})

bird.cv.current.sum <- Reduce('+', bird.cv.current)



bird.cv.pres.nat <- lapply(1:length(bird.mu), 
                           function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                                  meanlog = (log(exp(gpt.preds.bird) * exp(move.preds.bird)  * 1000 / 60 * exp(corr.preds.bird))[x]),
                                                  sdlog = my.sdlog)) * bird.mu[x] * matrix(m.bird.pres.nat[x,], nrow=142, byrow = T)})

bird.cv.pres.nat.sum <- Reduce('+', bird.cv.pres.nat)




bird.cv.te0 <- lapply(1:length(bird.mu), 
                      function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                             meanlog = (log(exp(gpt.preds.bird) * exp(move.preds.bird)  * 1000 / 60 * exp(corr.preds.bird))[x]),
                                             sdlog = my.sdlog)) * bird.mu[x] * matrix(m.bird.te0[x,], nrow=142, byrow = T)})

bird.cv.te0.sum <- Reduce('+', bird.cv.te0)







#

mamm.mu <- particip.preds.mamm * exp(quant.preds.mamm) * int.preds.mamm * pgerm(germ.preds.mamm) # Index of how many seeds dispersed

mamm.cv.current <- lapply(1:length(mamm.mu), 
                          function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                                 meanlog = (log(exp(gpt.preds.mamm) * exp(move.preds.mamm)  * 1000 / 60 * exp(corr.preds.mamm))[x]),
                                                 sdlog = my.sdlog)) * mamm.mu[x] * matrix(m.mamm.current[x,], nrow=142, byrow = T)})



mamm.cv.current.sum <- Reduce('+', mamm.cv.current)




mamm.cv.pres.nat <- lapply(1:length(mamm.mu), 
                           function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                                  meanlog = (log(exp(gpt.preds.mamm) * exp(move.preds.mamm)  * 1000 / 60 * exp(corr.preds.mamm))[x]),
                                                  sdlog = my.sdlog)) * mamm.mu[x] * matrix(m.mamm.pres.nat[x,], nrow=142, byrow = T)})

mamm.cv.pres.nat.sum <- Reduce('+', mamm.cv.pres.nat)



mamm.cv.te0 <- lapply(1:length(mamm.mu), 
                      function(x) {(1-plnorm(q = as.matrix(vclimmean), 
                                             meanlog = (log(exp(gpt.preds.mamm) * exp(move.preds.mamm)  * 1000 / 60 * exp(corr.preds.mamm))[x]),
                                             sdlog = my.sdlog)) * mamm.mu[x] * matrix(m.mamm.te0[x,], nrow=142, byrow = T)})

mamm.cv.te0.sum <- Reduce('+', mamm.cv.te0)




# Both

past.cv.loss <- (((bird.cv.current.sum + mamm.cv.current.sum) - (bird.cv.pres.nat.sum + mamm.cv.pres.nat.sum)) / (bird.cv.pres.nat.sum + mamm.cv.pres.nat.sum))
mean(past.cv.loss, na.rm = T)


te.cv.loss <- (((bird.cv.te0.sum + mamm.cv.te0.sum) - (bird.cv.current.sum + mamm.cv.current.sum)) / (bird.cv.current.sum + mamm.cv.current.sum))
mean(te.cv.loss, na.rm = T)




# For the purposes of the figures,

past.cv.loss.fig <- past.cv.loss
past.cv.loss.fig[past.cv.loss.fig > 0] <- -0.000001 # Won't plot gains...


pdf(file = paste("./figures/both cv",
                 "pdf", sep = "."),
    useDingbats = F,
    height = 4,# * (142/360),
    width = 6.1)


colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))

both.cv.raster <- raster(bird.cv.current.sum + mamm.cv.current.sum)# %>% log()
#both.cv.raster[both.cv.raster < -1] <- -1
both.cv.raster[both.cv.raster > 6] <- 6

lo.val <- 0#min(minValue(both.cv.raster))
hi.val <- 6#maxValue(both.cv.raster)

gr.val <- seq(0,round(hi.val), len = 5)

rasterVis::levelplot(both.cv.raster, 
                     margin=FALSE,                       # suppress marginal graphics
                     colorkey=list(
                       space='bottom',                   # plot legend at bottom
                       labels=list(at=gr.val)      # legend ticks and labels 
                     ),    
                     par.settings=list(
                       panel.background=list(col="grey40"),
                       strip.border=list(col='transparent'),
                       strip.background=list(col='transparent'),
                       axis.line=list(col='transparent') # suppress axes and legend outline
                     ),
                     scales=list(draw=FALSE),            # suppress axis labels
                     col.regions=colr,                   # colour ramp
                     at=seq(lo.val, hi.val, len=101),
                     names.attr=rep('', nlayers(both.cv.raster))) %>% 
  update(., aspect = 142/360)
lab.cex <- 0.75
#lab.x <- 0.072

grid.text("a)", x = unit(lab.x,"npc"), y = unit(0.845,"npc") + unit(1, "lines"), 
          just = c("left", "bottom"),
          gp = gpar(cex=lab.cex,
                    col = "white",
                    font = 1))



dev.off()






pdf(file = paste("./figures/bivariate cv",
                 #Sys.time(),
                 "pdf", sep = "."),
    useDingbats = F,
    height = 6 * (142/360),
    width = 6)


plotRGB(stack(list(raster((past.cv.loss.fig * -1)^(2/3)) * 255,
                   raster(past.cv.loss.fig * 0),
                   raster((te.cv.loss * -1)^(2/3)) * 255)), 1, 2, 3,
        #alpha = alpha.vals,
        asp = 142/360,
        #stretch = F,
        #interpolate = F,
        #maxpixels = 1000000,
        colNA = "grey40")

l.grid <- 20 # Originally 4

rgb.grid <- expand.grid(seq(0,255, len = l.grid),
                        seq(0,255, len = l.grid))

cols.rgb <- rgb((rgb.grid[,1]/255)^(2/3),
                0,
                (rgb.grid[,2]/255)^(2/3), maxColorValue = 1)


xx <- 1:l.grid * 2.2/360 + 0.1
yy <- 1:l.grid * 2.2/142 + 0.19
points(expand.grid(xx, yy),
       pch = 15,
       cex = 0.5,
       col = cols.rgb)


text("Current function\nendangered",
     srt = 90,
     pos = 3,
     col = "white",
     #font = 3,
     x = min(xx)*0.5, 
     y = mean(yy) * 0.85,
     cex = 0.7,
     0.15)

text("Past function\nlost",
     pos = 1,
     col = "white",
     #font = 3,
     x = mean(xx), 
     y = 0.14,
     cex = 0.7,
     0.15)

text(c("0%", "50%", "100%"),
     pos = 1,
     x = c(min(xx), mean(xx), max(xx)),
     y = min(yy),
     col = "white",
     cex = 0.5)

text(c("0%", "50%", "100%"),
     pos = 2,
     y = c(min(yy), mean(yy), max(yy)),
     x = min(xx)*1.05,
     col = "white",
     cex = 0.5)

dev.off()




# Table S3 output ecoregion results --------------------------------------------

# Want to add info about change in climate-outpacing seed dispersal

ecoregion.df$both.cv <- NA
ecoregion.df$cv.pixels <- NA


cv.both.pres.nat <- bird.cv.pres.nat.sum + mamm.cv.pres.nat.sum
cv.both.current <- bird.cv.current.sum + mamm.cv.current.sum

for(i in 1:length(ecoregion.stack)){
  ecoregion.set <- ecoregion.stack[[i]] %>% as.matrix(.) %>% is.na(.) %>% !.
  
  ecoregion.df$cv.pixels[i] <- sum(ecoregion.set & !is.na(cv.both.pres.nat))
  
  ecoregion.df$cv.pres.nat[i] <- cv.both.pres.nat[ecoregion.set] %>% mean(., na.rm=T)
  ecoregion.df$cv.current[i] <- cv.both.current[ecoregion.set] %>% mean(., na.rm=T)
  
  
  print(paste(i, 
              ecoregion.df$ecoregion.name[i]#, 
              #sum(ecoregion.set)
  ))
  
}


ecoregion.df <- ecoregion.df %>% 
  mutate(cv.loss = (cv.current - cv.pres.nat) / cv.pres.nat)
colnames(ecoregion.df)

table.s3 <- ecoregion.df[,c("pixels", "ecoregion.name", 
                            "function.current", "cv.current",
                            "richness.loss", "FDis.loss", "function.loss", "cv.loss")] %>% 
  filter(pixels > 10)

table.s3 <- table.s3 %>% mutate_if(is.numeric, round, digits = 2)

table.s3$richness.loss <- paste(table.s3$richness.loss * 100, "%", sep = "")
table.s3$FDis.loss <- paste(table.s3$FDis.loss * 100, "%", sep = "")
table.s3$function.loss <- paste(table.s3$function.loss * 100, "%", sep = "")
table.s3$cv.loss <- paste(table.s3$cv.loss * 100, "%", sep = "")

table.s3 %>% write.csv("./figures/table.s3.csv")



