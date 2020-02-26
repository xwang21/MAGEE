MAGEE <- function(null.obj, interaction, geno.file, group.file, group.file.sep = "\t", MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, missing.method = "impute2mean", method = "davies", tests = "JF", use.minor.allele = FALSE, auto.flip = FALSE, Garbage.Collection = FALSE, ncores = 1){
  if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
    warning("The package doMC is not available on Windows... Switching to single thread...")
    ncores <- 1
  }
  missing.method <- try(match.arg(missing.method, c("impute2mean", "impute2zero")))
  if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".")
  if(any(!tests %in% c("MV", "MF", "IV", "IF", "JV", "JF", "JD"))) stop("Error: \"tests\" should only include \"MV\" for the main effect variance component test, \"MF\" for the main effect combined test of the burden and variance component tests using Fisher\'s method, \"IV\" for the interaction variance component test, \"IF\" for the interaction combined test of the burden and variance component tests using Fisher\'s method, \"JV\" for the joint variance component test for main effect and interaction, \"JF\" for the joint combined test of the burden and variance component tests for main effect and interaction using Fisher\'s method, or \"JD\" for the joint combined test of the burden and variance component tests for main effect and interaction using double Fisher\'s method.")
  MV <- "MV" %in% tests
  MF <- "MF" %in% tests
  IV <- "IV" %in% tests
  IF <- "IF" %in% tests
  JV <- "JV" %in% tests
  JF <- "JF" %in% tests
  JD <- "JD" %in% tests
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
  } else residuals <- null.obj$scaled.residuals
  if(class(interaction)=="character") E <- as.matrix(null.obj$X[,which(colnames(null.obj$X) %in% interaction)])
  else E <- as.matrix(null.obj$X[,interaction+1])
  E <- scale(E, scale = FALSE)
  if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
  gds <- SeqArray::seqOpen(geno.file)
  sample.id <- SeqArray::seqGetData(gds, "sample.id")
  if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
  sample.id <- sample.id[sample.id %in% null.obj$id_include]
  if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
  match.id <- match(sample.id, null.obj$id_include)
  residuals <- residuals[match.id]
  if(!is.null(null.obj$P)) null.obj$P <- null.obj$P[match.id, match.id]
  else {
    null.obj$Sigma_iX <- null.obj$Sigma_iX[match.id, , drop = FALSE]
    null.obj$Sigma_i <- null.obj$Sigma_i[match.id, match.id]
  }
  E <- as.matrix(E[match.id, , drop = FALSE])
  strata <- apply(E, 1, paste, collapse = ":")
  strata <- if(length(unique(strata))>length(strata)/100) NULL else as.numeric(as.factor(strata))
  if(!is.null(strata)) strata.list <- lapply(unique(strata), function(x) which(strata==x))
  variant.idx <- SeqArray::seqGetData(gds, "variant.id")
  chr <- SeqArray::seqGetData(gds, "chromosome")
  pos <- SeqArray::seqGetData(gds, "position")
  alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
  ref <- unlist(lapply(alleles.list, function(x) x[1]))
  alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
  rm(alleles.list); gc()
  SeqArray::seqClose(gds)
  variant.id <- paste(chr, pos, ref, alt, sep = ":")
  rm(chr, pos, ref, alt); gc()
  group.info <- try(read.table(group.file, header = FALSE, col.names = c("group", "chr", "pos", "ref", "alt", "weight"), colClasses = c("character","character","integer","character","character","numeric"), sep = group.file.sep), silent = TRUE)
  if (class(group.info) == "try-error") {
    stop("Error: cannot read group.file!")
  }
  group.info <- group.info[!duplicated(paste(group.info$group, group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")), ]
  variant.id1 <- paste(group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")
  variant.idx1 <- variant.idx[match(variant.id1, variant.id)]
  group.info$variant.idx <- variant.idx1
  group.info$flip <- 0
  if(auto.flip) {
    cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
    variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    variant.idx2 <- variant.idx[match(variant.id2, variant.id)]
    if(any(!is.na(variant.idx1) & !is.na(variant.idx2))) {
      tmp.dups <- which(!is.na(variant.idx1) & !is.na(variant.idx2))
      cat("The following ambiguous variants were found:\n")
      cat("chr:", chr[tmp.dups], "\n")
      cat("pos:", pos[tmp.dups], "\n")
      cat("ref:", ref[tmp.dups], "\n")
      cat("alt:", alt[tmp.dups], "\n")
      cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
      variant.idx2[tmp.dups] <- NA
      rm(tmp.dups)
    }
    group.info$flip <- 1*(!is.na(variant.idx2))
    group.info$variant.idx[!is.na(variant.idx2)] <- variant.idx2[!is.na(variant.idx2)]
    rm(variant.id2, variant.idx2)
  }
  rm(variant.id, variant.id1, variant.idx1); gc()
  group.info <- subset(group.info, !is.na(variant.idx))
  groups <- unique(group.info$group)
  n.groups.all <- length(groups)
  group.info$group.idx <- as.numeric(factor(group.info$group))
  group.info <- group.info[order(group.info$group.idx, group.info$variant.idx), ]
  group.idx.end <- findInterval(1:n.groups.all, group.info$group.idx)
  group.idx.start <- c(1, group.idx.end[-n.groups.all] + 1)
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  if(ncores > 1) {
    doMC::registerDoMC(cores = ncores)
    n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
    b <- NULL
    out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
      idx <- if(b <= n.groups.percore_1) ((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1)) else (n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)
      n.groups <- length(idx)
      gds <- SeqArray::seqOpen(geno.file)
      SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
      n.variants <- rep(0,n.groups)
      miss.min <- rep(NA,n.groups)
      miss.mean <- rep(NA, n.groups)
      miss.max <- rep(NA, n.groups)
      freq.min <- rep(NA, n.groups)
      freq.mean <- rep(NA, n.groups)
      freq.max <- rep(NA, n.groups)
      freq.strata.min <- rep(NA, n.groups)
      freq.strata.max <- rep(NA, n.groups)
      if(MV | JV) MV.pval <- rep(NA, n.groups)
      if(IV | JV) IV.pval <- rep(NA, n.groups)
      if(JV) JV.pval <- rep(NA, n.groups)
      if(MF | JF | JD) MF.pval <- rep(NA, n.groups)
      if(IF | JF | JD) IF.pval <- rep(NA, n.groups)
      if(JF) JF.pval <- rep(NA, n.groups)
      if(JD) JD.pval <- rep(NA, n.groups)
      for(i in 1:n.groups) {
        tmp.idx <- group.idx.start[idx[i]]:group.idx.end[idx[i]]
        tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
        SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
        geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
        miss <- colMeans(is.na(geno))
        freq <- colMeans(geno, na.rm = TRUE)/2
        include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
        if(!is.null(strata)) { # E is not continuous
          #freq_strata <- apply(geno,2,function(x) range(tapply(x,strata,mean,na.rm=TRUE)/2))
          freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
          if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp))
          include <- include & !is.na(freq_strata[1,]) & !is.na(freq_strata[2,]) & freq_strata[1,] >= MAF.range[1] & freq_strata[2,] <= 1-MAF.range[1]
          rm(freq.tmp)
        }
        n.p <- sum(include)
        if(n.p == 0) next
        tmp.group.info <- tmp.group.info[include, , drop = FALSE]
        miss <- miss[include]
        freq <- freq[include]
        geno <- geno[, include, drop = FALSE]
        if(!is.null(strata)) freq_strata <- freq_strata[, include, drop = FALSE]
        N <- nrow(geno) - colSums(is.na(geno))
        if(sum(tmp.group.info$flip) > 0) {
          freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
          geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
          if(!is.null(strata)) freq_strata[, tmp.group.info$flip==1] <- 1 - freq_strata[, tmp.group.info$flip==1]
        }
        if(max(miss)>0) {
          miss.idx <- which(is.na(geno))
          geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
        }
        U <- as.vector(crossprod(geno, residuals))
        if(!is.null(null.obj$P)) PG <- crossprod(null.obj$P, geno)
        else {
          GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
          PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
        }
        V <- as.matrix(crossprod(geno, PG))
        if(IV | IF | JV | JF | JD) {
          K <- do.call(cbind, sapply(1:ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
          if(!is.null(null.obj$P)) KPK <- crossprod(K,crossprod(null.obj$P,K))
          else {
            KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
            KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
          }
          V_i <- MASS::ginv(V)
          KPG <- crossprod(K,PG)
          IV.U <- crossprod(K,residuals)-tcrossprod(tcrossprod(KPG,V_i),t(U))
          IV.V <- KPK - tcrossprod(tcrossprod(KPG,V_i),KPG)
        }
        if(use.minor.allele) {
          tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
          if(!is.null(strata)) freq_strata[, freq > 0.5] <- 1 - freq_strata[, freq > 0.5]
          freq[freq > 0.5] <- 1 - freq[freq > 0.5]
        }
        tmp.group.info$weight <- tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2])
        n.variants[i] <- n.p
        miss.min[i] <- min(miss)
        miss.mean[i] <- mean(miss)
        miss.max[i] <- max(miss)
        freq.min[i] <- min(freq)
        freq.mean[i] <- mean(freq)
        freq.max[i] <- max(freq)
        if(!is.null(strata)) {
          freq.strata.min[i] <- min(freq_strata)
          freq.strata.max[i] <- max(freq_strata)
        }
        U <- U*tmp.group.info$weight
        V <- t(V*tmp.group.info$weight)*tmp.group.info$weight
        if(IV | IF | JV | JF | JD) {
          IV.U <- IV.U*rep(tmp.group.info$weight, ncol(E))
          IV.V <- t(IV.V*rep(tmp.group.info$weight, ncol(E)))*rep(tmp.group.info$weight, ncol(E))
        }
        if(MV | JV) MV.pval[i] <- tryCatch(.quad_pval(U = U, V = V, method = method), error = function(e) { NA })
        if(IV | JV) IV.pval[i] <- tryCatch(.quad_pval(U = IV.U, V = IV.V, method = method), error = function(e) { NA })
        if(JV) JV.pval[i] <- tryCatch(fisher_pval(c(MV.pval[i], IV.pval[i])), error = function(e) { MV.pval[i] })
        if(MF | JF | JD) {
          MF.BU <- sum(U)
          MF.BV <- sum(V)
          MF.Bp <- pchisq(MF.BU^2/MF.BV,df=1,lower.tail=FALSE)
          V.rowSums <- rowSums(V)
          MF.U <- U - V.rowSums * MF.BU / MF.BV
          MF.V <- V - tcrossprod(V.rowSums) / MF.BV
          if(mean(abs(MF.V)) < sqrt(.Machine$double.eps)) MF.p <- NA
          else MF.p <- tryCatch(.quad_pval(U = MF.U, V = MF.V, method = method), error = function(e) { NA })
          MF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p)), error = function(e) { MF.Bp })
        }
        if(IF | JF | JD) {
          IF.BU <- sum(IV.U)
          IF.BV <- sum(IV.V)
          IF.Bp <- pchisq(IF.BU^2/IF.BV,df=1,lower.tail=FALSE)
          IV.V.rowSums <- rowSums(IV.V)
          IF.U <- IV.U - IV.V.rowSums * IF.BU / IF.BV
          IF.V <- IV.V - tcrossprod(IV.V.rowSums) / IF.BV
          if(mean(abs(IF.V)) < sqrt(.Machine$double.eps)) IF.p <- NA
          else IF.p <- tryCatch(.quad_pval(U = IF.U, V = IF.V, method = method), error = function(e) { NA })
          IF.pval[i] <- tryCatch(fisher_pval(c(IF.Bp, IF.p)), error = function(e) { IF.Bp })
        }
        if(JF) JF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p, IF.Bp, IF.p)), error = function(e) { MF.Bp })
        if(JD) JD.pval[i] <- tryCatch(fisher_pval(c(MF.pval[i], IF.pval[i])), error = function(e) { MF.pval[i] })
        rm(geno)
        if(Garbage.Collection) gc()
      }
      SeqArray::seqClose(gds)
      tmp.out <- data.frame(group=unique(group.info$group)[idx], n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max,freq.strata.min=freq.strata.min, freq.strata.max=freq.strata.max)
      if(MV | JV) tmp.out$MV.pval <- MV.pval
      if(MF | JF | JD) tmp.out$MF.pval <- MF.pval
      if(IV | JV) tmp.out$IV.pval <- IV.pval
      if(IF | JF | JD) tmp.out$IF.pval <- IF.pval
      if(JV) tmp.out$JV.pval <- JV.pval
      if(JF) tmp.out$JF.pval <- JF.pval
      if(JD) tmp.out$JD.pval <- JD.pval
      tmp.out
    }
  } else { # use a single core
    n.groups <- n.groups.all
    gds <- SeqArray::seqOpen(geno.file)
    SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
    n.variants <- rep(0,n.groups)
    miss.min <- rep(NA,n.groups)
    miss.mean <- rep(NA, n.groups)
    miss.max <- rep(NA, n.groups)
    freq.min <- rep(NA, n.groups)
    freq.mean <- rep(NA, n.groups)
    freq.max <- rep(NA, n.groups)
    freq.strata.min <- rep(NA, n.groups)
    freq.strata.max <- rep(NA, n.groups)
    if(MV | JV) MV.pval <- rep(NA, n.groups)
    if(IV | JV) IV.pval <- rep(NA, n.groups)
    if(JV) JV.pval <- rep(NA, n.groups)
    if(MF | JF | JD) MF.pval <- rep(NA, n.groups)
    if(IF | JF | JD) IF.pval <- rep(NA, n.groups)
    if(JF) JF.pval <- rep(NA, n.groups)
    if(JD) JD.pval <- rep(NA, n.groups)
    for(i in 1:n.groups) {
      tmp.idx <- group.idx.start[i]:group.idx.end[i]
      tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
      SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
      geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
      miss <- colMeans(is.na(geno))
      freq <- colMeans(geno, na.rm = TRUE)/2
      include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
      if(!is.null(strata)) { # E is not continuous
        #freq_strata <- apply(geno,2,function(x) range(tapply(x,strata,mean,na.rm=TRUE)/2))
        freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
        if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp))
        include <- include & !is.na(freq_strata[1,]) & !is.na(freq_strata[2,]) & freq_strata[1,] >= MAF.range[1] & freq_strata[2,] <= 1-MAF.range[1]
        rm(freq.tmp)
      }
      n.p <- sum(include)
      if(n.p == 0) next
      tmp.group.info <- tmp.group.info[include, , drop = FALSE]
      miss <- miss[include]
      freq <- freq[include]
      geno <- geno[, include, drop = FALSE]
      if(!is.null(strata)) freq_strata <- freq_strata[, include, drop = FALSE]
      N <- nrow(geno) - colSums(is.na(geno))
      if(sum(tmp.group.info$flip) > 0) {
        freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
        geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
        if(!is.null(strata)) freq_strata[, tmp.group.info$flip==1] <- 1 - freq_strata[, tmp.group.info$flip==1]
      }
      if(max(miss)>0) {
        miss.idx <- which(is.na(geno))
        geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
      }
      U <- as.vector(crossprod(geno, residuals))
      if(!is.null(null.obj$P)) PG <- crossprod(null.obj$P, geno)
      else {
        GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
        PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
      }
      V <- as.matrix(crossprod(geno, PG))
      if(IV | IF | JV | JF | JD) {
        K <- do.call(cbind, sapply(1:ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
        if(!is.null(null.obj$P)) KPK <- crossprod(K,crossprod(null.obj$P,K))
        else {
          KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
          KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
        }
        V_i <- MASS::ginv(V)
        KPG <- crossprod(K,PG)
        IV.U <- crossprod(K,residuals)-tcrossprod(tcrossprod(KPG,V_i),t(U))
        IV.V <- KPK - tcrossprod(tcrossprod(KPG,V_i),KPG)
      }
      if(use.minor.allele) {
        tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
        if(!is.null(strata)) freq_strata[, freq > 0.5] <- 1 - freq_strata[, freq > 0.5]
        freq[freq > 0.5] <- 1 - freq[freq > 0.5]
      }
      tmp.group.info$weight <- tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2])
      n.variants[i] <- n.p
      miss.min[i] <- min(miss)
      miss.mean[i] <- mean(miss)
      miss.max[i] <- max(miss)
      freq.min[i] <- min(freq)
      freq.mean[i] <- mean(freq)
      freq.max[i] <- max(freq)
      if(!is.null(strata)) {
        freq.strata.min[i] <- min(freq_strata)
        freq.strata.max[i] <- max(freq_strata)
      }
      U <- U*tmp.group.info$weight
      V <- t(V*tmp.group.info$weight)*tmp.group.info$weight
      if(IV | IF | JV | JF | JD) {
        IV.U <- IV.U*rep(tmp.group.info$weight, ncol(E))
        IV.V <- t(IV.V*rep(tmp.group.info$weight, ncol(E)))*rep(tmp.group.info$weight, ncol(E))
      }
      if(MV | JV) MV.pval[i] <- tryCatch(.quad_pval(U = U, V = V, method = method), error = function(e) { NA })
      if(IV | JV) IV.pval[i] <- tryCatch(.quad_pval(U = IV.U, V = IV.V, method = method), error = function(e) { NA })
      if(JV) JV.pval[i] <- tryCatch(fisher_pval(c(MV.pval[i], IV.pval[i])), error = function(e) { MV.pval[i] })
      if(MF | JF | JD) {
        MF.BU <- sum(U)
        MF.BV <- sum(V)
        MF.Bp <- pchisq(MF.BU^2/MF.BV,df=1,lower.tail=FALSE)
        V.rowSums <- rowSums(V)
        MF.U <- U - V.rowSums * MF.BU / MF.BV
        MF.V <- V - tcrossprod(V.rowSums) / MF.BV
        if(mean(abs(MF.V)) < sqrt(.Machine$double.eps)) MF.p <- NA
        else MF.p <- tryCatch(.quad_pval(U = MF.U, V = MF.V, method = method), error = function(e) { NA })
        MF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p)), error = function(e) { MF.Bp })
      }
      if(IF | JF | JD) {
        IF.BU <- sum(IV.U)
        IF.BV <- sum(IV.V)
        IF.Bp <- pchisq(IF.BU^2/IF.BV,df=1,lower.tail=FALSE)
        IV.V.rowSums <- rowSums(IV.V)
        IF.U <- IV.U - IV.V.rowSums * IF.BU / IF.BV
        IF.V <- IV.V - tcrossprod(IV.V.rowSums) / IF.BV
        if(mean(abs(IF.V)) < sqrt(.Machine$double.eps)) IF.p <- NA
        else IF.p <- tryCatch(.quad_pval(U = IF.U, V = IF.V, method = method), error = function(e) { NA })
        IF.pval[i] <- tryCatch(fisher_pval(c(IF.Bp, IF.p)), error = function(e) { IF.Bp })
      }
      if(JF) JF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p, IF.Bp, IF.p)), error = function(e) { MF.Bp })
      if(JD) JD.pval[i] <- tryCatch(fisher_pval(c(MF.pval[i], IF.pval[i])), error = function(e) { MF.pval[i] })
      rm(geno)
      if(Garbage.Collection) gc()
    }
    SeqArray::seqClose(gds)
    out <- data.frame(group=unique(group.info$group), n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max,freq.strata.min=freq.strata.min, freq.strata.max=freq.strata.max)
    if(MV | JV) out$MV.pval <- MV.pval
    if(MF | JF | JD) out$MF.pval <- MF.pval
    if(IV | JV) out$IV.pval <- IV.pval
    if(IF | JF | JD) out$IF.pval <- IF.pval
    if(JV) out$JV.pval <- JV.pval
    if(JF) out$JF.pval <- JF.pval
    if(JD) out$JD.pval <- JD.pval
  }
  return(out[match(groups, out$group),])
}

.Q_pval <- function(Q, lambda, method = "davies") {
  if(method == "davies") {
    tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
    pval <- tmp$Qq
    if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
  }
  if(method == "kuonen") {
    pval <- .pKuonen(x = Q, lambda = lambda)
    if(is.na(pval)) method <- "liu"
  }
  if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
  return(pval)
}

.quad_pval <- function(U, V, method = "davies") {
  Q <- sum(U^2)
  lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
  lambda <- lambda[lambda > 0]
  pval <- .Q_pval(Q, lambda, method = method)
  return(pval)
}

.pKuonen<-function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
  delta <- delta[lambda != 0]
  df <- df[lambda != 0]
  lambda <- lambda[lambda != 0]
  if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
  if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
  if (length(lambda) == 1) {
    pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE)
  }
  d <- max(lambda)
  lambda <- lambda/d
  x <- x/d
  k0 <- function(zeta) {
    -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
                                                     zeta)/(1 - 2 * zeta * lambda))
  }
  kprime0 <- function(zeta) {
    sapply(zeta, function(zz) {
      sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                                                                 zz * lambda^2)/(1 - 2 * zz * lambda)^2)
    })
  }
  kpprime0 <- function(zeta) {
    sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
          delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
  }
  if (any(lambda < 0)) {
    lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
  }
  else if (x > sum((df+delta)*lambda)) {
    lmin <- -0.01
  }
  else {
    lmin <- -length(lambda)*max(df+delta)/(2 * x)
  }
  lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
                     upper = lmax, tol = 1e-08)$root
  w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v <- hatzeta * sqrt(kpprime0(hatzeta))
  if (abs(hatzeta) < 1e-04)
    NA
  else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

MAF.weights.beta.fun <- function(freq, beta1, beta2) {
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}

fisher_pval <- function(p) {
  is.valid.p <- !is.na(p) & p > 0 & p <= 1
  if(sum(is.valid.p) == 0) return(NA)
  p <- p[is.valid.p]
  pchisq(-2*sum(log(p)), df = 2*length(p), lower.tail = FALSE)
}

