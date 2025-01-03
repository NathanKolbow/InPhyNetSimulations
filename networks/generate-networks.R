args <- commandArgs(trailingOnly = TRUE)
ntaxa <- as.numeric(args[1])

library(SiPhyNetwork)
nets <- sim.bdh.taxa.ssa(n = ntaxa, numbsim = 1, lambda = 1, mu = 0.2, nu = 0.025,  # placeholder so that we can write to output
                        hybprops = c(0.5, 0.25, 0.25), hyb.inher.fxn = make.beta.draw(10, 10))
set.seed(42)

j <- 1
while(length(nets) < 100) {
    iter_net <- sim.bdh.taxa.ssa(n = ntaxa, numbsim = 1, lambda = 1, mu = 0.2, nu = 0.025,
                                 hybprops = c(0.5, 0.25, 0.25), hyb.inher.fxn = make.beta.draw(10, 10))[[1]]
    if(is.numeric(iter_net)) { next }
    if(length(iter_net$tip.label) != ntaxa) { next }
    if(getNetworkLevel(iter_net) > 1) { next }
    if(nrow(iter_net$reticulation) == 0) { next }
    
    nets[[j]] <- iter_net
    j <- j + 1
    cat(paste0("\rFound #", length(nets)+1, " - nretic = ", nrow(iter_net$reticulation)))
}

write.net(nets, paste0("n", ntaxa, ".netfile"))





