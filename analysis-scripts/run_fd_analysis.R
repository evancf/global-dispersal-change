library("FD")

# Functional diversity analyses

# Get matrices of species presence by ecoregion (0 or 1)
presence.by.ecoregion.current <- cbind(t(m.mamm.current[,1:length(ecoregion.stack)]), 
                                       t(m.bird.current[,1:length(ecoregion.stack)]))
presence.by.ecoregion.current[1:length(ecoregion.stack),] <- NA

presence.by.ecoregion.pres.nat <- cbind(t(m.mamm.pres.nat[,1:length(ecoregion.stack)]), 
                                        t(m.bird.pres.nat[,1:length(ecoregion.stack)]))
presence.by.ecoregion.pres.nat[1:length(ecoregion.stack),] <- NA


mamm.inds <- 1:nrow(m.mamm.current)
bird.inds <- (length(mamm.inds) + 1): (nrow(m.mamm.current)+nrow(m.bird.current))

for(i in 1:length(ecoregion.stack)){
  ecoreg.inds <- which(as.vector(ecoregion.stack[[i]]) == 1)
  
  presence.by.ecoregion.current[i, mamm.inds] <- apply(m.mamm.current[,ecoreg.inds, drop = F], 1, function(x) any(x == 1))
  presence.by.ecoregion.pres.nat[i, mamm.inds] <- apply(m.mamm.pres.nat[,ecoreg.inds, drop = F], 1, function(x) any(x == 1))
  
  presence.by.ecoregion.current[i, bird.inds] <- apply(m.bird.current[,ecoreg.inds, drop = F], 1, function(x) any(x == 1))
  presence.by.ecoregion.pres.nat[i, bird.inds] <- apply(m.bird.pres.nat[,ecoreg.inds, drop = F], 1, function(x) any(x == 1))
  
  print(i)
  
}

ecoregion.names <- lapply(ecoregion.stack, names) %>% unlist()

presence.by.ecoregion.current <- presence.by.ecoregion.current[,order(colnames(presence.by.ecoregion.current))]
rownames(presence.by.ecoregion.current) <- paste("current", ecoregion.names, sep = ".")

presence.by.ecoregion.pres.nat <- presence.by.ecoregion.pres.nat[,order(colnames(presence.by.ecoregion.pres.nat))]
rownames(presence.by.ecoregion.pres.nat) <- paste("pres.nat", ecoregion.names, sep = ".")

presence.by.ecoregion <- rbind(presence.by.ecoregion.current,
                               presence.by.ecoregion.pres.nat)

presence.by.ecoregion <- presence.by.ecoregion[-which(rowSums(presence.by.ecoregion) == 0),]
presence.by.ecoregion <- presence.by.ecoregion[,-which(colSums(presence.by.ecoregion) == 0)]

# Get a trait dataframe
animal.trait.cols <- colnames(all.mamm.traits)[c(9:18, 21, 27:30)]
fd.trait <- bind_rows(all.mamm.traits[,c("name.m", animal.trait.cols)],
                      all.bird.traits[,c("name.m", animal.trait.cols)])
fd.trait <- fd.trait[order(fd.trait$name.m),]
dim(fd.trait)


fd.trait <- fd.trait %>% filter(name.m %in% colnames(presence.by.ecoregion))
fd.trait <- fd.trait[!duplicated(fd.trait$name.m),]

rownames(fd.trait) <- fd.trait$name.m
fd.trait <- dplyr::select(fd.trait, -name.m)


# Will manipulate these traits some

fd.trait$volant <- ordered(fd.trait$volant)
fd.trait$activity <- ordered(fd.trait$activity)

fd.trait$BodyMass.Value <- log(fd.trait$BodyMass.Value)
fd.trait$hab_breadth <- log(fd.trait$hab_breadth)

# Because there are the many diet columns, will weight these to equal 1 unit
fd.weights <- c(rep(0.1, 10), rep(1,5))

# Use dbFD to calculate FDis. 
# This takes quite a while - perhaps a couple days on a laptop.

time0 <- Sys.time()
fd.output <- dbFD(x = fd.trait,
                  presence.by.ecoregion,
                  corr = "lingoes",
                  w = fd.weights,
                  calc.FRic = F,
                  calc.CWM = F,
                  calc.FDiv = F,
                  calc.FGR = F)
time1 <- Sys.time()


# setwd("/Users/evanfricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/analysis")
# save(fd.output,
#      file = "fd.output.Rdata")


