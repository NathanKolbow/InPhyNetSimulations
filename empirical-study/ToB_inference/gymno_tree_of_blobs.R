library(ape)
library(MSCquartets)

setwd("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/empirical-study/ToB_inference")
gts <- read.tree(file="./gymno_gts.tre")

relevant_taxa <- c("Pinta_v2.0", "GTHK", "TOXE", "VDAO", "AQFM",
                   "AREG", "AWQB", "DZQM", "GAMH", "GGEA", "IIOL",
                   "IOVS", "JUWL", "MFTM", "NPRL", "VSRH", "WVWN",
                   "JBND", "ACWS", "AIGO", "AUDE", "BBDD", "BTTS",
                   "BUWV", "CDFR", "CGDN", "EFMS", "EGLZ", "ETCJ",
                   "FHST", "FMWZ", "FRPM", "GKCZ", "GMHZ", "HBGV",
                   "HILW", "HQOM", "IAJW", "IFLI", "IZGN", "JDQB",
                   "JRNA", "JZVE", "KLGF", "MHGD", "MIXZ", "NKIN",
                   "NRXL", "NVGZ", "OVIJ", "OWFC", "OXGJ", "QCGM",
                   "QFAE", "QHBI", "QNGJ", "QSNJ", "RMMV", "ROWR",
                   "RSCE", "SCEB", "UEVI", "UUJS", "VFYZ", "VGSX",
                   "WWSS", "XIRK", "XLGK", "XMGP", "XQSG", "XTZO",
                   "YFZK", "YLPM", "YYPE", "ZQVF", "ZQWM", "ZYAX",
                   "URDJ", "FZJL", "NWMY", "VZCI", "PZRT", "WTKZ",
                   "GNQG", "KAWQ", "WLIC", "XZUY", "SGTW", "RICC",
                   "MTGC", "NWWI", "JKAA", "Selmo_v1.0", "GAON")

qt <- quartetTableParallel(gts, numCores=6, taxonnames=relevant_taxa)
tt <- TINNIK(qt, alpha=0.01, beta=0.99, plot=FALSE)
write.tree(tt$ToB, file="gymnosperm_tob.tre")
save.image("gymnosperm_results.RData")
