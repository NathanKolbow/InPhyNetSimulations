args <- commandArgs(trailingOnly = TRUE)

library(SiPhyNetwork)
n <- as.numeric(args[1])
ils <- args[2]
seed <- as.numeric(args[3])
output <- args[4]


minCycleSize <- function(network) {
    rt <- as.integer(length(network$tip.label)+1)
    edges <- rbind(network$edge, network$reticulation)
    mode(edges) <- 'integer'
    nNode <- length(network$tip.label)+network$Nnode
    blobs <- biconnectedComponents(edges, rt, nNode)
    min(unlist(lapply(blobs, length)))
}

nu <- 100.0
if(n == 50) {
    nu <- 0.025
} else if(n == 100) {
    nu <- 0.0035
} else if(n == 200) {
    nu <- 0.001
} else if(n == 500) {
    nu <- 0.00015
} else {
    cat("ERROR: n =", n, "not a valid parameter.")
}

iseed <- seed
set.seed(seed)
net <- sim.bdh.taxa.ssa(
    n = n, numbsim = 1, lambda = 1, mu = 0, nu = nu,
    hybprops = c(0.5, 0.25, 0.25), hyb.inher.fxn = make.beta.draw(10, 10)
)[[1]]
while(getNetworkLevel(net) > 1 || minCycleSize(net) <= 3) {
    seed <- seed + 1
    set.seed(seed)
    net <- sim.bdh.taxa.ssa(
        n = n, numbsim = 1, lambda = 1, mu = 0, nu = nu,
        hybprops = c(0.5, 0.25, 0.25), hyb.inher.fxn = make.beta.draw(10, 10)
    )[[1]]
    if(seed > iseed + 10000) {
        cat("ERROR: Looped 10,000 times without finding suitable network.\n")
        break
    }
}
cat("Generated network after", seed-iseed+1, "attempts.\n")
write.net(net, file=output)