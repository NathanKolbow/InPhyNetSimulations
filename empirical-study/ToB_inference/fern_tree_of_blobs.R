library(ape)
library(MSCquartets)

setwd("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/empirical-study/ToB_inference")
gts <- read.tree(file="./fern_gts.tre")

relevant_taxa <- c("DFHO", "NHCM", "UGNK", "BEGM", "DJSE", "EEAQ",
                   "ALVQ", "QVMR", "GANB", "PNZO", "EWXK", "UWOD",
                   "MEKP", "QIAD", "WGTU", "RFMZ", "UOMY", "VIBO",
                   "AFPO", "BMIF", "BMJR", "CJNT", "DCDT", "FCHS",
                   "FLTD", "FQGQ", "GSXD", "GYFU", "HEGQ", "HNDZ",
                   "HTFH", "IXLH", "JBLI", "KJZG", "LHLE", "MROH",
                   "MTGC", "NDUV", "NOKI", "NWWI", "OCZL", "OQWW",
                   "ORJE", "PIVW", "POPJ", "PSKY", "RFRB", "RICC",
                   "SKYV", "UFJN", "UJTT", "UJWU", "URCP", "VITX",
                   "VVRN", "WCLG", "WQML", "XDDT", "XXHP", "YCKE",
                   "YIXP", "YJJY", "YLJA", "YOWV", "YQEC", "ZQYU",
                   "ZXJO", "CVEG", "KIIX", "CQPW", "PBUU", "UTRE",
                   "PRIQ", "QWRA", "XDLL", "WRSL", "MFTM", "OWFC",
                   "QNGJ", "IZGN", "EGLZ", "TOXE", "QHBI", "GKCZ",
                   "AQFM", "NKIN", "PKOX", "GKAG", "GAON", "XNXF",
                   "YHZW", "Selmo_v1.0", "ABIJ", "UPMJ", "KJYC", "PYHZ")

qt <- quartetTableParallel(gts, numCores=8, taxonnames=relevant_taxa)
tt <- TINNIK(qt, alpha=0.01, beta=0.99, plot=FALSE)
write.tree(tt$ToB, file="fern_tob.tre")
save.image("fern_results.RData")
