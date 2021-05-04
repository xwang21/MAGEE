### R code from vignette source 'MAGEE.Rnw'

###################################################
### code chunk number 1: installation (eval = FALSE)
###################################################
## ## try http:// if https:// URLs are not supported
## ## remove "doMC" below if you are running Windows
## install.packages(c("devtools", "RcppArmadillo", "CompQuadForm", "doMC", 
##         "foreach", "Matrix", "GMMAT", "BiocManager", "testthat"), 
## 	repos = "https://cran.r-project.org/")
## BiocManager::install(c("SeqArray", "SeqVarTools"))
## devtools::install_github("https://github.com/xwang21/MAGEE")


###################################################
### code chunk number 2: convert2GDS (eval = FALSE)
###################################################
## SeqArray::seqVCF2GDS("VCF_file_name", "GDS_file_name")
## SeqArray::seqBED2GDS("BED_file_name", "FAM_file_name", "BIM_file_name", 
##         "GDS_file_name")


###################################################
### code chunk number 3: loading (eval = FALSE)
###################################################
## library(MAGEE)


###################################################
### code chunk number 4: help (eval = FALSE)
###################################################
## ?MAGEE


###################################################
### code chunk number 5: MAGEEglmmkingds (eval = FALSE)
###################################################
## GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "MAGEE")
## GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
## model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM,
##                   id = "id", family = binomial(link = "logit"))


###################################################
### code chunk number 6: MAGEEgeigds (eval = FALSE)
###################################################
## infile <- system.file("extdata", "geno.gds", package = "MAGEE")
## glmm.gei(model0, interaction='sex', geno.file = infile, 
##          outfile = "glmm.gei.gds.testoutfile.txt")


###################################################
### code chunk number 7: MAGEEgeibgen (eval = FALSE)
###################################################
## infile <- system.file("extdata", "geno.bgen", package = "MAGEE")
## samplefile <- system.file("extdata", "geno.sample", package = "MAGEE")
## glmm.gei(model0, interaction='sex', geno.file = infile, 
##          outfile = "glmm.gei.bgen.testoutfile.txt",
##          bgen.samplefile = samplefile)


###################################################
### code chunk number 8: MAGEEmageegds (eval = FALSE)
###################################################
## geno.file <- system.file("extdata", "geno.gds", package = "MAGEE")
## group.file <- system.file("extdata", "SetID.withweights.txt", 
##                           package = "MAGEE")
## out <- MAGEE(model0, interaction='sex', geno.file, group.file, 
##              group.file.sep = "\t", tests=c("JV", "JF", "JD"))


###################################################
### code chunk number 9: MKL (eval = FALSE)
###################################################
## Sys.setenv(MKL_NUM_THREADS = 1)


###################################################
### code chunk number 10: RhpcBLASctlL (eval = FALSE)
###################################################
## #install.packages("RhpcBLASctl")
## library(RhpcBLASctl)
## blas_set_num_threads(1)


