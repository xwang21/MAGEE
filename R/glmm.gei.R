glmm.gei <- function(null.obj, interaction, geno.file, outfile, bgen.samplefile=NULL, interaction.covariates=NULL, meta.output=F, center=T, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch=100, ncores = 1){
  if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
    warning("The package doMC is not available on Windows... Switching to single thread...")
    ncores <- 1
  }
  
  if(!grepl("\\.gds$|\\.bgen$", geno.file[1])) stop("Error: only .gds and .bgen format is supported in geno.file!")
  
  missing.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
  if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"omit\".")
  if(!class(interaction) %in% c("integer", "numeric", "character")) stop("Error: \"interaction\" should be an integer, numeric, or character vector.")
  if(any(duplicated(null.obj$id_include))) {
    J <- Matrix(sapply(unique(null.obj$id_include), function(x) 1*(null.obj$id_include==x)), sparse = TRUE)
    residuals <- as.numeric(as.matrix(crossprod(J, null.obj$scaled.residuals)))
    if(!is.null(null.obj$P)) null.obj$P <- as.matrix(crossprod(J, crossprod(null.obj$P, J)))
    else {
      null.obj$Sigma_iX <- crossprod(J, null.obj$Sigma_iX)
      null.obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(null.obj$Sigma_i,J)))
      null.obj$Sigma_i <- Matrix(null.obj$Sigma_i, sparse = TRUE)
    }
    rm(J)
  } else {
    residuals <- null.obj$scaled.residuals
  }
  
  qi <- length(interaction.covariates)
  ei <- length(interaction)
  if(class(interaction)=="character") {
    if (is.null(interaction.covariates)) {
      if (!all(interaction %in% colnames(null.obj$X))) {stop("there are interactions not in column name of covariate matrix.")}
      E <- as.matrix(null.obj$X[,interaction])
      E <- scale(E, scale = FALSE)
    } else {
      if (any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
      interaction <- c(interaction.covariates, interaction)
      if (!all(interaction %in% colnames(null.obj$X))) {stop("there are interaction and interaction.covariates not in column name of covariate matrix.")}
      E <- as.matrix(null.obj$X[,interaction])
      E <- scale(E, scale = FALSE)
    }
  } else {
    if (is.null(interaction.covariates)) {
      E <- as.matrix(null.obj$X[,interaction+1])
      E <- scale(E, scale = FALSE)
    } else {
      interaction <- c(interaction.covariates, interaction)
      E <- as.matrix(null.obj$X[,interaction+1])
      E <- scale(E, scale = FALSE)
    }
  }
  
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  
  if(grepl("\\.gds$", geno.file[1])) {
    if (class(geno.file)[1] != "SeqVarGDSClass") {
      gds <- SeqArray::seqOpen(geno.file) 
      sample.id <- SeqArray::seqGetData(gds, "sample.id")
      variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
    } else {
      sample.id <- SeqArray::seqGetData(geno.file, "sample.id")
      variant.idx.all <- SeqArray::seqGetData(geno.file, "variant.id")
    }
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    match.id <- match(sample.id, null.obj$id_include)
    residuals <- residuals[match.id]
    if(!is.null(null.obj$P)) {
      null.obj$P <- null.obj$P[match.id, match.id]
    } else {
      null.obj$Sigma_iX <- null.obj$Sigma_iX[match.id, , drop = FALSE]
      null.obj$Sigma_i <- null.obj$Sigma_i[match.id, match.id]
    }
    E <- as.matrix(E[match.id, , drop = FALSE])
    strata <- apply(E, 1, paste, collapse = ":")
    strata <- if(length(unique(strata))>length(strata)/100) NULL else as.numeric(as.factor(strata))
    if(!is.null(strata)) strata.list <- lapply(unique(strata), function(x) which(strata==x))
    if (class(geno.file)[1] != "SeqVarGDSClass") {
      SeqArray::seqClose(gds)
    }
    p.all <- length(variant.idx.all)
    
    if(ncores > 1) {
      doMC::registerDoMC(cores = ncores)
      p.percore <- (p.all-1) %/% ncores + 1
      n.p.percore_1 <- p.percore * ncores - p.all
      foreach(b=1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
        variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
        p <- length(variant.idx)
        if (class(geno.file)[1] != "SeqVarGDSClass") {
          gds <- SeqArray::seqOpen(geno.file)
        } else {
          gds <- geno.file
        }
        SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
        rm(sample.id)
        nbatch.flush <- (p-1) %/% 100000 + 1
        ii <- 0
        for(i in 1:nbatch.flush) {
          gc()
          tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
          SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
          MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
          AF <- 1 - SeqVarTools::alleleFrequency(gds)
          include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
          if(sum(include) == 0) {
            write.table("", paste0(outfile, "_tmp.", b), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=(ii > 1), na=".")
            next
          }
          ii <- ii + 1
          tmp.variant.idx <- tmp.variant.idx[include]
          tmp.p <- length(tmp.variant.idx)
          SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
          SNP <- SeqArray::seqGetData(gds, "annotation/id")
          SNP[SNP == ""] <- NA
          out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
          rm(SNP)
          alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
          out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
          out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
          out$MISSRATE <- MISSRATE[include]
          out$AF <- AF[include]
          rm(alleles.list, include)
          tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
            tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
            SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
            
            geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
            ng <- ncol(geno)
            
            N <- nrow(geno) - colSums(is.na(geno))
            AF.strata.min <- AF.strata.max <- rep(NA, ncol(geno))
            if(!is.null(strata)) { # E is not continuous
              #freq_strata <- apply(geno,2,function(x) range(tapply(x,strata,mean,na.rm=TRUE)/2))
              freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
              if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp))
              AF.strata.min <- freq_strata[1,]
              AF.strata.max <- freq_strata[2,]
            }
            
            if(center) geno <- scale(geno, scale = FALSE)
            
            miss.idx <- which(is.na(geno))
            if(length(miss.idx)>0) {
              geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
            }
            
            K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
            if(!is.null(interaction.covariates)) {
              geno <- cbind(geno, do.call(cbind, sapply(1:qi, function(xx) geno*E[,xx], simplify = FALSE), envir = environment()))
            }
            
            U <- as.vector(crossprod(geno, residuals))
            if(!is.null(null.obj$P)) {
              PG <- crossprod(null.obj$P, geno)
            } else {
              GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
              PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
            }
            
            
            GPG <- as.matrix(crossprod(geno, PG)) * (matrix(1, 1+qi, 1+qi) %x% diag(ng))
            GPG_i <- try(solve(GPG), silent = TRUE)
            if(class(GPG_i)[1] == "try-error") GPG_i <- MASS::ginv(GPG)
            if(is.null(interaction.covariates)) {
              V_i <- diag(GPG_i)
            } else {
              V <- diag(GPG)[1:ng]
              V_i <- ifelse(V < .Machine$double.eps, 0, 1/V)
            }
            
            V.MAIN.adj <- diag(GPG_i)[1:ng]
            BETA.MAIN.adj <- as.vector(crossprod(GPG_i, U))[1:ng]
            STAT.MAIN.adj <- ifelse(V.MAIN.adj > 0, BETA.MAIN.adj^2/V.MAIN.adj, NA)
            
            BETA.MAIN <- V_i * U[1:ng]
            SE.MAIN <- sqrt(V_i)
            STAT.MAIN <- BETA.MAIN * U[1:ng]
            PVAL.MAIN <- ifelse(V_i>0, pchisq(STAT.MAIN, df=1, lower.tail=FALSE), NA)
            
            if(!is.null(null.obj$P)) {
              KPK <- crossprod(K,crossprod(null.obj$P,K))
            } else {
              KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
              KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
            }
            KPK <- as.matrix(KPK) * (matrix(1, ei, ei) %x% diag(ng))
            KPG <- as.matrix(crossprod(K,PG)) * (matrix(1, ei, 1+qi) %x% diag(ng))
            
            if(!is.null(interaction.covariates)) {
              IV.U <- (rep(1, ei) %x% diag(ng)) * as.vector(crossprod(K,residuals)-tcrossprod(KPG, t(crossprod(GPG_i, U))))
              IV.V_i <- try(solve(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i))), silent = TRUE)
              if(class(IV.V_i)[1] == "try-error") IV.V_i <- MASS::ginv(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i)))
              BETA.INT <- crossprod(IV.V_i, IV.U)
              STAT.INT <- diag(crossprod(IV.U,  BETA.INT))
              PVAL.INT <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
              PVAL.JOINT <- ifelse(is.na(STAT.MAIN.adj), NA, pchisq(STAT.MAIN.adj + STAT.INT, df=1+ei, lower.tail=FALSE))
              
            } else {
              IV.U <- (rep(1, ei) %x% diag(ncol(geno))) * as.vector(crossprod(K,residuals)-tcrossprod(KPG, t(BETA.MAIN)))
              IV.V_i <- try(solve(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i))), silent = TRUE)
              if(class(IV.V_i)[1] == "try-error") IV.V_i <- MASS::ginv(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i)))
              BETA.INT <- crossprod(IV.V_i, IV.U)
              STAT.INT <- diag(crossprod(IV.U, BETA.INT))
              PVAL.INT <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
              PVAL.JOINT <- ifelse(is.na(PVAL.MAIN), NA, pchisq(STAT.MAIN + STAT.INT, df=1+ei, lower.tail=FALSE))
            }
            
            
            if (meta.output) {
              return(rbind(N, AF.strata.min, AF.strata.max, BETA.MAIN[1:ng], SE.MAIN, 
                           do.call(cbind, lapply(1:ng, function(x) {idx <- seq(x, ng*ei, ng); BETA.INT[idx,x]})),
                           do.call(cbind, lapply(1:ng, function(x) {idx <- seq(x, ng*ei, ng); IV.V_i[idx,idx][1:(ei*ei)]})),
                           PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
            } else {
              return(rbind(N, AF.strata.min, AF.strata.max, BETA.MAIN[1:ng], SE.MAIN,
                           PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
            }
          })
          if (meta.output) {
            cov.header = matrix(paste(rep(paste0("Cov_Gx", interaction[1:ei]), each = ei), interaction[1:ei], sep = "_Gx"), ei, ei)
            diag(cov.header) <- paste0("VAR.BETA.Gx", interaction[1:ei])
            ss.header = c(paste0("BETA.Gx", interaction[1:ei]), cov.header)
            tmp.out <- matrix(unlist(tmp.out), ncol = 9 + ei + ei*ei, byrow = TRUE, dimnames = list(NULL, c("N", "AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", ss.header, "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT")))
            if (ei != 1) {
              ss.header = c(ss.header[1:ei], diag(cov.header), cov.header[lower.tri(cov.header)])
            }
            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,c("MISSRATE","AF")], tmp.out[,c("AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", ss.header, "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT"), drop = F])
          } else {
            tmp.out <- matrix(unlist(tmp.out), ncol = 9, byrow = TRUE, dimnames = list(NULL, c("N", "AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT")))
            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,c("MISSRATE","AF")], tmp.out[,c("AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT"), drop = F])
          }
          
          names(out)[6] <- "N"
          rm(tmp.out)
          if(b == 1) {
            write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
          } else {
            write.table(out, paste0(outfile, "_tmp.", b), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=(ii > 1), na=".")
          }
          rm(out)
        }
        SeqArray::seqClose(gds)
      }
      for(b in 2:ncores) {
        system(paste0("cat ", outfile, "_tmp.", b, " >> ", outfile))
        unlink(paste0(outfile, "_tmp.", b))
      }
    } else { # use a single core
      variant.idx <- variant.idx.all
      rm(variant.idx.all)
      p <- length(variant.idx)
      if (class(geno.file)[1] != "SeqVarGDSClass") {
        gds <- SeqArray::seqOpen(geno.file)
      } else {
        gds <- geno.file
      }
      SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
      rm(sample.id)
      nbatch.flush <- (p-1) %/% 100000 + 1
      ii <- 0
      for(i in 1:nbatch.flush) {
        gc()
        tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
        AF <- 1 - SeqVarTools::alleleFrequency(gds)
        include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        if(sum(include) == 0) {
          next
        }
        ii <- ii + 1
        tmp.variant.idx <- tmp.variant.idx[include]
        tmp.p <- length(tmp.variant.idx)
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        SNP <- SeqArray::seqGetData(gds, "annotation/id")
        SNP[SNP == ""] <- NA
        out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
        rm(SNP)
        alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
        out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
        out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
        out$MISSRATE <- MISSRATE[include]
        out$AF <- AF[include]
        rm(alleles.list, include)
        tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
          tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
          SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
          
          geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
          ng <- ncol(geno)
          
          N <- nrow(geno) - colSums(is.na(geno))
          AF.strata.min <- AF.strata.max <- rep(NA, ncol(geno))
          if(!is.null(strata)) { # E is not continuous
            #freq_strata <- apply(geno,2,function(x) range(tapply(x,strata,mean,na.rm=TRUE)/2))
            freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
            if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp))
            AF.strata.min <- freq_strata[1,]
            AF.strata.max <- freq_strata[2,]
          }
          
          if(center) geno <- scale(geno, scale = FALSE)
          miss.idx <- which(is.na(geno))
          if(length(miss.idx)>0) {
            geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
          }
          
          K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
          if(!is.null(interaction.covariates)) {
            geno <- cbind(geno, do.call(cbind, sapply(1:qi, function(xx) geno*E[,xx], simplify = FALSE), envir = environment()))
          }
          
          U <- as.vector(crossprod(geno, residuals))
          if(!is.null(null.obj$P)) {
            PG <- crossprod(null.obj$P, geno)
          } else {
            GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
            PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
          }
          
          GPG <- as.matrix(crossprod(geno, PG)) * (matrix(1, 1+qi, 1+qi) %x% diag(ng))
          GPG_i <- try(solve(GPG), silent = TRUE)
          if(class(GPG_i)[1] == "try-error") GPG_i <- MASS::ginv(GPG)
          if(is.null(interaction.covariates)) {
            V_i <- diag(GPG_i)
          } else {
            V <- diag(GPG)[1:ng]
            V_i <- ifelse(V < .Machine$double.eps, 0, 1/V)
          }
          
          V.MAIN.adj <- diag(GPG_i)[1:ng]
          BETA.MAIN.adj <- as.vector(crossprod(GPG_i, U))[1:ng]
          STAT.MAIN.adj <- ifelse(V.MAIN.adj > 0, BETA.MAIN.adj^2/V.MAIN.adj, NA)
          
          BETA.MAIN <- V_i * U[1:ng]
          SE.MAIN <- sqrt(V_i)
          STAT.MAIN <- BETA.MAIN * U[1:ng]
          PVAL.MAIN <- ifelse(V_i>0, pchisq(STAT.MAIN, df=1, lower.tail=FALSE), NA)
          
          if(!is.null(null.obj$P)) {
            KPK <- crossprod(K,crossprod(null.obj$P,K))
          } else {
            KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
            KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
          }
          KPK <- as.matrix(KPK) * (matrix(1, ei, ei) %x% diag(ng))
          KPG <- as.matrix(crossprod(K,PG)) * (matrix(1, ei, 1+qi) %x% diag(ng))
          
          if(!is.null(interaction.covariates)) {
            IV.U <- (rep(1, ei) %x% diag(ng)) * as.vector(crossprod(K,residuals)-tcrossprod(KPG, t(crossprod(GPG_i, U))))
            IV.V_i <- try(solve(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i))), silent = TRUE)
            if(class(IV.V_i)[1] == "try-error") IV.V_i <- MASS::ginv(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i)))
            BETA.INT <- crossprod(IV.V_i, IV.U)
            STAT.INT <- diag(crossprod(IV.U,  BETA.INT))
            PVAL.INT <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
            PVAL.JOINT <- ifelse(is.na(STAT.MAIN.adj), NA, pchisq(STAT.MAIN.adj + STAT.INT, df=1+ei, lower.tail=FALSE))
            
          } else {
            IV.U <- (rep(1, ei) %x% diag(ncol(geno))) * as.vector(crossprod(K,residuals)-tcrossprod(KPG, t(BETA.MAIN)))
            IV.V_i <- try(solve(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i))), silent = TRUE)
            if(class(IV.V_i)[1] == "try-error") IV.V_i <- MASS::ginv(KPK - tcrossprod(KPG,tcrossprod(KPG, GPG_i)))
            BETA.INT <- crossprod(IV.V_i, IV.U)
            STAT.INT <- diag(crossprod(IV.U, BETA.INT))
            PVAL.INT <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
            PVAL.JOINT <- ifelse(is.na(PVAL.MAIN), NA, pchisq(STAT.MAIN + STAT.INT, df=1+ei, lower.tail=FALSE))
          }
          
          
          if (meta.output) {
            return(rbind(N, AF.strata.min, AF.strata.max, BETA.MAIN[1:ng], SE.MAIN, 
                         do.call(cbind, lapply(1:ng, function(x) {idx <- seq(x, ng*ei, ng); BETA.INT[idx,x]})),
                         do.call(cbind, lapply(1:ng, function(x) {idx <- seq(x, ng*ei, ng); IV.V_i[idx,idx][1:(ei*ei)]})),
                         PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
          } else {
            return(rbind(N, AF.strata.min, AF.strata.max, BETA.MAIN[1:ng], SE.MAIN,
                         PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
          }
          
        })
        
        SeqArray::seqClose(gds)
        if (meta.output) {
          cov.header = matrix(paste(rep(paste0("Cov_Gx", interaction[1:ei]), each = ei), interaction[1:ei], sep = "_Gx"), ei, ei)
          diag(cov.header) <- paste0("VAR.BETA.Gx", interaction[1:ei])
          ss.header = c(paste0("BETA.Gx", interaction[1:ei]), cov.header)
          
          tmp.out <- matrix(unlist(tmp.out), ncol = 9 + ei + ei*ei, byrow = TRUE, dimnames = list(NULL, c("N", "AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", ss.header, "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT")))
          if (ei != 1) {
            ss.header = c(ss.header[1:ei], diag(cov.header), cov.header[lower.tri(cov.header)])
          }
          out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,c("MISSRATE","AF")], tmp.out[,c("AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", ss.header, "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT"), drop = F])
        } else {
          tmp.out <- matrix(unlist(tmp.out), ncol = 9, byrow = TRUE, dimnames = list(NULL, c("N", "AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT")))
          out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,c("MISSRATE","AF")], tmp.out[,c("AF.strata.min", "AF.strata.max", "BETA.MAIN", "SE.MAIN", "PVAL.MAIN", "STAT.INT", "PVAL.INT", "PVAL.JOINT"), drop = F])
        }
        
        names(out)[6] <- "N"
        rm(tmp.out)
        write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
        rm(out)
      }
    }
    return(invisible(NULL))
  } else if (grepl("\\.bgen$", geno.file)) {
    
    bgenInfo <- .Call('bgenHeader', geno.file)
    
    if (bgenInfo$SampleIdFlag == 0) {
      if (is.null(bgen.samplefile)) {
        stop("Error: bgen file does not contain sample identifiers. A .sample file (bgen.samplefile) is needed.")
      }
      sample.id <- read.table(bgen.samplefile, header = TRUE, sep = " ")
      if ((nrow(sample.id)-1) != bgenInfo$N){
        stop(paste0("Error: Number of sample identifiers in BGEN sample file (", nrow(sample.id)-1, ") does not match number of samples in BGEN file (", bgenInfo$N,")."))
      }
      sample.id <- sample.id[-1, 1]
    } else {
      sample.id <- bgenInfo$SampleIds
    }
    
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of bgen sample file!")
    select <- match(sample.id, unique(null.obj$id_include))
    select[is.na(select)] <- 0
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    
    match.id <- match(sample.id, null.obj$id_include)
    residuals <- residuals[match.id]
    if(!is.null(null.obj$P)) {
      null.obj$P <- null.obj$P[match.id, match.id]
    } else {
      null.obj$Sigma_iX <- Matrix(null.obj$Sigma_iX[match.id, , drop = FALSE], sparse = TRUE)
      null.obj$Sigma_i <- Matrix(null.obj$Sigma_i[match.id, match.id], sparse = TRUE)
      null.obj$cov <- Matrix(null.obj$cov, sparse = TRUE)
    }
    E <- as.matrix(E[match.id, , drop = FALSE])
    strata <- apply(E, 1, paste, collapse = ":")
    strata <- if(length(unique(strata))>length(strata)/100) NULL else as.numeric(as.factor(strata))
    if(!is.null(strata)) {
      strata.list <- lapply(unique(strata), function(x) which(strata==x))
    } else {
      strata.list <- NULL
    }
    if (is.null(interaction.covariates)) {
      null.obj$E <- E
    } else {
      null.obj$E <- E[,(qi+1):(qi+ei), drop = F]
      null.obj$EC <- E[,1:qi, drop = F]
    }
    rm(sample.id, E)
    
    variant.idx.all <- 1:bgenInfo$M
    p.all <- length(variant.idx.all)
    
    if (ncores > bgenInfo$M) {
      ncores <- bgenInfo$M
      print(paste0("Warning: number of cores (", ncores,") is greater than number of variants in BGEN files (", bgenInfo$M,"). Using ", ncores, " instead."))
    }
    
    threadInfo <- .Call("getVariantPos", geno.file, bgenInfo$offset, bgenInfo$M, bgenInfo$N, bgenInfo$CompressionFlag, bgenInfo$LayoutFlag, ncores)
    center2 <- ifelse(center, 'c', 'n')
    missing.method <- substr(missing.method, 1, 1)
    
    if (ncores > 1) {
      doMC::registerDoMC(cores = ncores)
      foreach(i=1:ncores) %dopar% {
        if (bgenInfo$LayoutFlag == 2) {
          .Call("glmm_gei_bgen13", as.numeric(residuals), null.obj, geno.file, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, missing.method, nperbatch, ei, qi, is.null(null.obj$P), is.null(interaction.covariates), strata.list, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, meta.output)
        } else {
          .Call("glmm_gei_bgen11", as.numeric(residuals), null.obj, geno.file, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, missing.method, nperbatch, ei, qi, is.null(null.obj$P), is.null(interaction.covariates), strata.list, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, meta.output)
        }
      }
      
      outTmp <- file(outfile, "w")
      if (meta.output) {
        if (ei == 1) {
          cov.header = NULL
          ss.header  = c(paste0("BETA.Gx", interaction[1:ei]), paste0("VAR.BETA.Gx", interaction[1:ei]))
          
        } else {
          cov.header = matrix(paste(rep(paste0("Cov_Gx", interaction[1:ei]), each = ei), interaction[1:ei], sep = "_Gx"), ei, ei)
          cov.header = cov.header[lower.tri(cov.header)]
          ss.header  = c(paste0("BETA.Gx", interaction[1:ei]), paste0("VAR.BETA.Gx", interaction[1:ei]), cov.header)
        }
        
        writeLines(paste0("SNP\tRSID\tCHR\tPOS\tREF\tALT\tN\tMISSRATE\tAF\tAF.strata.min\tAF.strata.max\tBETA.MAIN\tSE.MAIN\t",paste0(ss.header, collapse = "\t"),"\tPVAL.MAIN\tSTAT.INT\tPVAL.INT\tPVAL.JOINT\n"), outTmp)
      } else {
        writeLines(paste0("SNP\tRSID\tCHR\tPOS\tREF\tALT\tN\tMISSRATE\tAF\tAF.strata.min\tAF.strata.max\tBETA.MAIN\tSE.MAIN\tPVAL.MAIN\tSTAT.INT\tPVAL.INT\tPVAL.JOINT\n"), outTmp)
      }
      
      for(i in 1:ncores){
        inTmp <- readLines(paste0(outfile, "_tmp.", i))
        writeLines(inTmp, outTmp)
        unlink(paste0(outfile, "_tmp.", i))
      }
      close(outTmp)
    } else {
      outTmp <- file(outfile, "w")
      if (meta.output) {
        if (ei == 1) {
          cov.header = NULL
          ss.header  = c(paste0("BETA.Gx", interaction[1:ei]), paste0("VAR.BETA.Gx", interaction[1:ei]))
          
        } else {
          cov.header = matrix(paste(rep(paste0("Cov_Gx", interaction[1:ei]), each = ei), interaction[1:ei], sep = "_Gx"), ei, ei)
          cov.header = cov.header[lower.tri(cov.header)]
          ss.header  = c(paste0("BETA.Gx", interaction[1:ei]), paste0("VAR.BETA.Gx", interaction[1:ei]), cov.header)
        }
        
        writeLines(paste0("SNP\tRSID\tCHR\tPOS\tREF\tALT\tN\tMISSRATE\tAF\tAF.strata.min\tAF.strata.max\tBETA.MAIN\tSE.MAIN\t",paste0(ss.header, collapse = "\t"),"\tPVAL.MAIN\tSTAT.INT\tPVAL.INT\tPVAL.JOINT\n"), outTmp)
      } else {
        writeLines(paste0("SNP\tRSID\tCHR\tPOS\tREF\tALT\tN\tMISSRATE\tAF\tAF.strata.min\tAF.strata.max\tBETA.MAIN\tSE.MAIN\tPVAL.MAIN\tSTAT.INT\tPVAL.INT\tPVAL.JOINT\n"), outTmp)
      }
      close(outTmp)
      if (bgenInfo$LayoutFlag == 2) {
        .Call("glmm_gei_bgen13", as.numeric(residuals), null.obj, geno.file, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, missing.method, nperbatch, ei, qi, is.null(null.obj$P), is.null(interaction.covariates), strata.list, select, threadInfo$begin[1], threadInfo$end[1], threadInfo$pos[1], bgenInfo$N, bgenInfo$CompressionFlag, meta.output)
      } else {
        .Call("glmm_gei_bgen11", as.numeric(residuals), null.obj, geno.file, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, missing.method, nperbatch, ei, qi, is.null(null.obj$P), is.null(interaction.covariates), strata.list, select, threadInfo$begin[1], threadInfo$end[1], threadInfo$pos[1], bgenInfo$N, bgenInfo$CompressionFlag, meta.output)
      }
      
    }
    return(invisible(NULL))
  }
}