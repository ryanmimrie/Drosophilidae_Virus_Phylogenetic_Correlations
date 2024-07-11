# ------------------------------------------------------------------------------
# ----- Imrie et al., 2024 Phylogenetic MCMCglmms ------------------------------
# ------------------------------------------------------------------------------

# ----- 0. Initialisation ------------------------------------------------------

# ----- 0.1 Dependencies -------------------------------------------------------

library(tidyverse); library(ape); library(MCMCglmm); library(patchwork)


# ----- 0.2 Data ---------------------------------------------------------------

setwd("/path/to/this/script")

data <- read_csv("data/data_vials.csv")

data <- as.data.frame(data)

# ----- 0.3 Phylogeny ----------------------------------------------------------

tree <- read.tree("data/phylogeny_host.nwk")

tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% unique(data$animal)])

plot(tree)


# ----- 1. Univeriate models ---------------------------------------------------

# ----- 1.1 Priors -------------------------------------------------------------

prior.her.uni <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

prior.rep.uni <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

# ----- 1.2 Models -------------------------------------------------------------

uni.models <- list()

i <- 1 # set to 1 for quick look, 100 for study iterations

for(vir in unique(data$virus)){
  
  model.her <- MCMCglmm(foldchange ~ 1, random = ~animal + host,
                        rcov = ~units, pedigree = tree, prior = prior.her.uni, data = filter(data, virus == vir),
                        nitt = 130000*i, thin = 50*i, burnin = 30000*i, pr = TRUE)
  
  model.rep <- MCMCglmm(foldchange ~ 1, random = ~animal,
                        rcov = ~units, pedigree = tree, prior = prior.rep.uni, data = filter(data, virus == vir),
                        nitt = 130000*i, thin = 50*i, burnin = 30000*i, pr = TRUE)
  
  
  uni.models[[vir]] <- list(model.her = model.her, model.rep = model.rep)
  
}

#save(uni.models, file = "models/Univeriate_Models.RData")


# ----- 2. Bivariate models ----------------------------------------------------

# ----- 2.1 Priors -------------------------------------------------------------

prior.her.bi <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000)),
  R = list(V = diag(2), nu = 0.002))

prior.rep.bi <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000)),
  R = list(V = diag(2), nu = 0.002))


# ----- 2.2 Models -------------------------------------------------------------

bi.models <- list()

i <- 1 # set to 1 for quick look, 100 for study iterations

x <- 1

for(vir1 in unique(data$virus)){
  for(vir2 in unique(data$virus)){
  
    if(vir1 == vir2) next
    
    print(sprintf("Running virus combination %s of 110", x))
    x <- x + 1
    
    model.her <- MCMCglmm(foldchange ~ virus, random = ~us(virus) : animal + us(virus) : host,
                          rcov = ~idh(virus) : units, pedigree = tree, prior = prior.her.bi, data = filter(data, virus %in% c(vir1, vir2)),
                          nitt = 130000*i, thin = 50*i, burnin = 30000*i, pr = TRUE)
    
    model.rep <- MCMCglmm(foldchange ~ virus, random = ~us(virus) : animal,
                          rcov = ~idh(virus) : units, pedigree = tree, prior = prior.rep.bi, data = filter(data, virus %in% c(vir1, vir2)),
                          nitt = 130000*i, thin = 50*i, burnin = 30000*i, pr = TRUE)
    
    
    bi.models[[paste(vir1, vir2, sep = ":")]] <- list(model.her = model.her, model.rep = model.rep)
    
  }
  
}

#save(bi.models, file = "models/Bivariate_Models.RData")


# ----- 3. Univariate Analyses -------------------------------------------------

# ----- 3.1 Load Models --------------------------------------------------------

#load("models/Univeriate_Models.RData")

# ----- 3.2 Extract univariate posterior densities -----------------------------

estimates <- data.frame(virus = character(0),
                        across.species.mean = numeric(0),
                        phylo.var.rep = numeric(0),
                        residual.var.rep = numeric(0),
                        phylo.coefficient = numeric(0),
                        phylo.var.her = numeric(0),
                        species.var.her = numeric(0))

for (virus in unique(data$virus)){
  
  model <- uni.models[[virus]]$model.rep
  model2 <- uni.models[[virus]]$model.her
  
  current <- data.frame(virus = virus,
                        across.species.mean = model$Sol[,1],
                        phylo.var.rep = model$VCV[,1],
                        residual.var.rep = model$VCV[,2],
                        phylo.coefficient = ((model$VCV[,1]^0.5)/model$Sol[,1])*100,
                        phylo.var.her = model2$VCV[,1],
                        species.var.her = model2$VCV[,2])
  
  colnames(current) = c("virus", "across.species.mean", "phylo.var.rep",
                        "residual.var.rep", "phylo.coefficient",
                        "phylo.var.her", "species.var.her")
  
  estimates <- rbind(estimates, current)
  
}

# ----- 3.3 Wrangle estimates --------------------------------------------------

estimates <- na.omit(estimates)

estimates$virus <- factor(estimates$virus, levels = rev(c("CrPV-OG", "CrPV-GRA", "CrPV-VIC",
                                                          "DCV-C", "DCV-EB", "DCV-M",
                                                          "FHV", "DAV", "IIV6", "DmelNv", "BFV")))

estimates$her <- estimates$phylo.var.her / (estimates$phylo.var.her + estimates$species.var.her)
estimates$rep <- estimates$phylo.var.rep / (estimates$phylo.var.rep + estimates$residual.var.rep)

estimate.summaries <- estimates %>% group_by(virus) %>%
  summarise(across.species.mean.mean = mean(across.species.mean),
            across.species.mean.low = HPDinterval(as.mcmc(across.species.mean))[,1],
            across.species.mean.high = HPDinterval(as.mcmc(across.species.mean))[,2],
            phylo.var.rep.mean = mean(phylo.var.rep),
            phylo.var.rep.low = HPDinterval(as.mcmc(phylo.var.rep))[,1],
            phylo.var.rep.high = HPDinterval(as.mcmc(phylo.var.rep))[,2],
            residual.var.rep.mean = mean(residual.var.rep),
            residual.var.rep.low = HPDinterval(as.mcmc(residual.var.rep))[,1],
            residual.var.rep.high = HPDinterval(as.mcmc(residual.var.rep))[,2],
            phylo.coefficient.mean = mean(phylo.coefficient, trim = 0.001),
            phylo.coefficient.low = HPDinterval(as.mcmc(phylo.coefficient))[,1],
            phylo.coefficient.high = HPDinterval(as.mcmc(phylo.coefficient))[,2],
            her.mean = mean(her),
            her.low = HPDinterval(as.mcmc(her))[,1],
            her.high = HPDinterval(as.mcmc(her))[,2],
            rep.mean = mean(rep),
            rep.low = HPDinterval(as.mcmc(rep))[,1],
            rep.high = HPDinterval(as.mcmc(rep))[,2])

# ----- 3.4 Plot estimates -----------------------------------------------------

p1 <- ggplot(estimate.summaries) +
  geom_point(aes(x = across.species.mean.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = across.species.mean.low,
                    xmax = across.species.mean.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Mean Across-Species Fold-change in Viral Load",
                     breaks = seq(-4, 6, 2),
                     labels = c(expression(10^-4),
                                expression(10^-2),
                                expression(0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6))) +
  theme_bw() +
  theme(axis.title.y = element_blank())

p2 <- ggplot(estimate.summaries) +
  geom_point(aes(x = phylo.var.rep.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = phylo.var.rep.low,
                    xmax = phylo.var.rep.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = expression("Phylogenetic Variance (V"[p]~")"),
                     breaks = seq(0, 12, 4),
                     labels = c(expression(0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6))) +
  theme_bw() +
  theme(axis.title.y = element_blank())

p3 <- ggplot(estimate.summaries) +
  geom_point(aes(x = residual.var.rep.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = residual.var.rep.low,
                    xmax = residual.var.rep.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = expression("Residual Variance (V"[r]~")"),
                     breaks = seq(0, 2, 1),
                     labels = c(expression(0),
                                expression(10^1),
                                expression(10^2)),
                     limits = c(0, 2)) +
  theme_bw() +
  theme(axis.title.y = element_blank())

p4 <- ggplot(estimate.summaries) +
  geom_point(aes(x = phylo.coefficient.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = phylo.coefficient.low,
                    xmax = phylo.coefficient.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Coefficient of Phylogenetic Variation (%)",
                     breaks = seq(0, 300, 100)) +
  theme_bw() +
  theme(axis.title.y = element_blank())

(p1 | p2) / (p3 | p4)

p5 <- ggplot(estimate.summaries) +
  geom_point(aes(x = her.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = her.low, xmax = her.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Phylogenetic Heritability") +
  theme_bw() +
  theme(axis.title.y = element_blank())

p6 <- ggplot(estimate.summaries) +
  geom_point(aes(x = rep.mean, y = virus), size = 2) +
  geom_errorbar(aes(y = virus, xmin = rep.low, xmax = rep.high), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Repeatability") +
  theme_bw() +
  theme(axis.title.y = element_blank())

(p5 | p6)


# ----- 4. Bivariate analyses --------------------------------------------------

# ----- 4.1 Load models --------------------------------------------------------

#load("models/Bivariate_Models.RData")

# ----- 4.2 Extract bivariate posterior densities ------------------------------

estimates.bi <- data.frame(viruses = character(0),
                           cov = numeric(0),
                           var1 = numeric(0),
                           var2 = numeric(0))

for(viruses in names(bi.models)) {
  
  model <- bi.models[[viruses]]$model.rep
  
  current <- data.frame(viruses = viruses,
                        cov = model$VCV[,2],
                        var1 = model$VCV[,1],
                        var2 = model$VCV[,4])
  
  estimates.bi <- rbind(estimates.bi, current)

}

colnames(estimates.bi) <- c("viruses", "cov", "var1", "var2")

estimates.bi$correlation <- estimates.bi$cov/sqrt(estimates.bi$var1*estimates.bi$var2)
estimates.bi$slope <- estimates.bi$cov/estimates.bi$var1

estimate.summaries.bi <- estimates.bi %>% group_by(viruses) %>%
  summarise(correlation.mean = mean(as.mcmc(correlation)),
            correlation.low = HPDinterval(as.mcmc(correlation))[,1],
            correlation.high = HPDinterval(as.mcmc(correlation))[,2],
            slope.mean = mean(as.mcmc(slope)),
            slope.low = HPDinterval(as.mcmc(slope))[,1],
            slope.high = HPDinterval(as.mcmc(slope))[,2])

estimate.summaries.bi[, 2:7] <- round(estimate.summaries.bi[, 2:7], digits = 2)

estimate.summaries.bi <- estimate.summaries.bi %>%
  mutate(
    correlation = paste(correlation.mean, " (", correlation.low, ", ", correlation.high, ")", sep = ""),
    slope = paste(slope.mean, "(", slope.low, ", ", slope.high, ")", sep = "")
  ) %>%
  select(viruses, correlation, slope)

View(estimate.summaries.bi)
