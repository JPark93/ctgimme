ctsgimme = function(varnames = NULL, dataframe = NULL,
                    id = NULL, time = NULL,
                    cores = 1, directory = NULL, 
                    sig.thrsh = 0.55, sub.sig.thrsh = 1.00,
                    Galpha = 0.05, ben.hoch = TRUE, S.Galpha = 0.05, Ialpha = 0.01,
                    ME.var = 1e-8, PE.var = diag(1, nvar),
                    subgroup.model = FALSE,
                    time.intervals = c(1)){
  JPmx = function (model, matrices = NA, full = TRUE){
    OpenMx:::warnModelCreatedByOldVersion(model)
    if (OpenMx:::single.na(matrices)) {
      matrices = names(model$matrices)
      if (is(model$expectation, "MxExpectationRAM")) {
        matrices = setdiff(matrices, model$expectation$F)
      }
    }
    if (imxHasWLS(model)) {
      stop("modification indices not implemented for WLS fitfunction")
    }
    param = omxGetParameters(model)
    param.names = names(param)
    gmodel = omxSetParameters(model, free = FALSE, labels = param.names)
    mi.r = NULL
    mi.f = NULL
    EPCs = NULL
    a.names = NULL
    new.models = list()
    for (amat in matrices) {
      matObj = model[[amat]]
      freemat = matObj$free
      sym.sel = upper.tri(freemat, diag = TRUE)
      notSymDiag = !(is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                                   "SymmMatrix"))
      for (i in 1:length(freemat)) {
        if (freemat[i] == FALSE && (notSymDiag || sym.sel[i] == 
                                    TRUE)) {
          tmpLab = gmodel[[amat]]$labels[i]
          plusOneParamModel = model
          if (length(tmpLab) > 0 && !is.na(tmpLab)) {
            gmodel = omxSetParameters(gmodel, labels = tmpLab, 
                                       free = TRUE)
            plusOneParamModel = omxSetParameters(plusOneParamModel, 
                                                  labels = tmpLab, free = TRUE)
          }
          else {
            gmodel[[amat]]$free[i] = TRUE
            plusOneParamModel[[amat]]$free[i] = TRUE
          }
          if (is(gmodel[[amat]])[1] %in% c("ZeroMatrix")) {
            cop = gmodel[[amat]]
            newSingleParamMat = mxMatrix("Full", nrow = nrow(cop), 
                                          ncol = ncol(cop), values = cop$values, free = cop$free, 
                                          labels = cop$labels, name = cop$name, lbound = cop$lbound, 
                                          ubound = cop$ubound, dimnames = dimnames(cop))
            bop = plusOneParamModel[[amat]]
            newPlusOneParamMat = mxMatrix("Full", nrow = nrow(bop), 
                                           ncol = ncol(bop), values = bop$values, free = bop$free, 
                                           labels = bop$labels, name = bop$name, lbound = bop$lbound, 
                                           ubound = bop$ubound, dimnames = dimnames(bop))
          }
          else if (is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                                "SymmMatrix")) {
            cop = gmodel[[amat]]
            newSingleParamMat = mxMatrix("Symm", nrow = nrow(cop), 
                                          ncol = ncol(cop), values = cop$values, free = (cop$free | 
                                                                                           t(cop$free)), labels = cop$labels, name = cop$name, 
                                          lbound = cop$lbound, ubound = cop$ubound, 
                                          dimnames = dimnames(cop))
            bop = plusOneParamModel[[amat]]
            newPlusOneParamMat = mxMatrix("Symm", nrow = nrow(bop), 
                                           ncol = ncol(bop), values = bop$values, free = (bop$free | 
                                                                                            t(bop$free)), labels = bop$labels, name = bop$name, 
                                           lbound = bop$lbound, ubound = bop$ubound, 
                                           dimnames = dimnames(bop))
          }
          else {
            newSingleParamMat = gmodel[[amat]]
            newPlusOneParamMat = plusOneParamModel[[amat]]
          }
          gmodel[[amat]] = newSingleParamMat
          plusOneParamModel[[amat]] = newPlusOneParamMat
          custom.compute = mxComputeSequence(list(mxComputeNumericDeriv(checkGradient = FALSE), 
                                                   mxComputeReportDeriv()))
          gmodel = mxModel(gmodel, custom.compute)
          grun = try(mxRun(gmodel, silent = TRUE, suppressWarnings = FALSE, 
                            unsafe = TRUE))
          nings = TRUE
          if (is(grun, "try-error")) {
            gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                       free = FALSE)
            next
          }
          grad = grun$output$gradient
          hess = grun$output$hessian
          modind = 0.5 * grad^2/hess
          if (full == TRUE) {
            custom.compute.smart = mxComputeSequence(list(mxComputeNumericDeriv(knownHessian = model$output$hessian, 
                                                                                 checkGradient = FALSE), mxComputeReportDeriv()))
            plusOneParamRun = mxRun(mxModel(plusOneParamModel, 
                                             custom.compute.smart), silent = TRUE, suppressWarnings = FALSE, 
                                     unsafe = TRUE)
            grad.full = plusOneParamRun$output$gradient
            grad.full[is.na(grad.full)] = 0
            hess.full = plusOneParamRun$output$hessian
            modind.full = 0.5 * t(matrix(grad.full)) %*%
              solve(hess.full) %*% matrix(grad.full)
            EPC = modind.full/grad.full[grad.full!=0]
          }
          else {
            modind.full = NULL
          }
          n.names = names(omxGetParameters(grun))
          if (length(modind) > 0) {
            a.names = c(a.names, n.names)
            mi.r = c(mi.r, modind)
            mi.f = c(mi.f, modind.full)
            EPCs = c(EPCs, EPC)
            new.models = c(new.models, plusOneParamModel)
          }
          gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                     free = FALSE)
        }
      }
      names(mi.r) = a.names
      if (full == TRUE) {
        names(mi.f) = a.names
        names(EPCs) = a.names
      }
      names(new.models) = a.names
    }
    if (length(model$submodels) > 0) {
      for (asubmodel in names(model$submodels)) {
        ret = c(ret, JPmx(asubmodel))
      }
    }
    # message(paste0("The N is ", N))
    return(list(MI = mi.r, MI.Full = mi.f, plusOneParamModels = new.models, EPC = -EPCs))
  }
  safe_read_vector = function(file, element, param_names) {
    out = rep(NA_real_, length(param_names))
    names(out) = param_names
    
    vec = tryCatch({
      raw = readRDS(file)[[element]]
      vec = c(raw)
      if (!is.null(names(vec)) && length(vec) != length(names(vec))) {
        stop("Length of vector and names do not match")
      }
      vec
    }, error = function(e) {
      message("Failed to read ", element, " from ", file, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(vec)) {
      intersect_names = intersect(names(vec), param_names)
      out[intersect_names] = vec[intersect_names]
    }
    
    return(out)
  }
  
  #  Creating Directory---###---###---###
  ###---###---###---###---###---###---###
  dir.create(directory, showWarnings = FALSE)
  dir.create(paste0(directory, "/MIs/"), showWarnings = FALSE)
  dir.create(paste0(directory, "/Models/"), showWarnings = FALSE)
  dir.create(paste0(directory, "/Models/Individuals/"), showWarnings = FALSE)  
  # Specifying Cores
  message("Scaling Variables for Analysis")
  ids = unique(dataframe[,id])
  scales = NULL
  for(ID in ids){
    temp = subset(dataframe, id == ID)
    temp[,varnames] = scale(temp[,varnames])
    scales = rbind(scales, temp[,varnames])
  }
  dataframe[,varnames] = scales
  nvar = length(varnames)
  if(length(ids) < cores){
    cores = length(ids)
    message("Cores adjusted to your sample-size!")
  }else{cores = cores}
  ###---###---###---###---###---###---###
  # Setting Up Multiple Model Estimation#
  ###---###---###---###---###---###---###
  library(parallel)
  library(dynr)
  library(OpenMx)
  library(igraph)
  library(qgraph)
  cl = makeCluster(cores, type = "PSOCK")
  clusterExport(cl, c("dataframe", "JPmx",
                      "varnames", "id", "time", "directory",
                      "nvar", "ME.var", "PE.var"),
                envir = environment())
  clusterEvalQ(cl, {
    packs = list('dynr', 'OpenMx', 'qgraph')
    invisible(lapply(packs, require, character.only = T))
  })
  ###---###---###---###---###---###---###
  # Running Step 1 in Parallel---###---##
  ###---###---###---###---###---###---###
  parLapply(cl, ids, function(i) {
    DRIFT = diag(paste0("A_", 1:nvar, ",", 1:nvar), nvar)
    subset_dat = subset(dataframe, id == i)
    # subset_dat[,varnames] = scale(subset_dat[,varnames])
    amat = mxMatrix("Full", nvar, nvar,
                     free   = DRIFT != "0",
                     name   = "A", 
                    ubound = 30, lbound = -30)
    bmat = mxMatrix('Zero', nvar, nvar, name='B')
    cdim = list(varnames, paste0('F', 1:nvar))
    cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
    dmat = mxMatrix('Zero', nvar, nvar, name='D')
    qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
    rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
        xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
    pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
    umat = mxMatrix('Zero', nvar, 1, name='u')
    tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
    osc = mxModel("OUMod", 
                  amat, bmat, cmat, dmat, qmat, 
                  rmat, xmat, pmat, umat, tmat,
                  mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                    'R', 'x0', 'P0', 'u', 'time'),
                  mxFitFunctionML(),
                  mxData(subset_dat, 'raw'))  
    analysis_result = tryCatch({
      fit = mxTryHard(osc)
    }, error = function(e) {
      message("Error for subject ", i, ": ", e$message)
    })
    if (!is.null(analysis_result)) {
      MIs = JPmx(analysis_result, matrices = "A")
      saveRDS(
        object = MIs,
        file   = paste0(directory, "/MIs/MI_", i, ".RDS")
      )
    }
  })
  stopCluster(cl)
  ###---###---###---###---###---###---###
  # Iteration Steps
  ###---###---###---###---###---###---###
  iterate = 0; count = 1
  m = (nvar^2)-nvar
  ks = matrix(NA, m, 1)
  for(k in 1:m){
    ks[k,] = (k/m) * Galpha
  }
  ks = matrix(sort(ks, TRUE), m, 1)
  if(ben.hoch == FALSE){
    ks = matrix(Galpha, nrow(ks), 1)
  }
  DRIFT = diag(paste0("A_", 1:nvar, 1:nvar), nvar)
  while(iterate < 1){
    # rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
    # files = NULL
    # EPCs = NULL
    # for (file in rdss) {
    #   file_id = gsub("MI_|\\.RDS", "", basename(file))
    #   files = cbind(files, tryCatch({
    #     c(readRDS(file)$"MI.Full")
    #   }, error = function(e) {
    #     message("Failed to read ", file, ": ", e$message)
    #     NULL
    #   }))
    #   message(length(c(readRDS(file)$"MI.Full")))
    #   EPCs = cbind(EPCs, tryCatch({
    #     c(readRDS(file)$"EPC")
    #   }, error = function(e) {
    #     message("Failed to read ", file, ": ", e$message)
    #     NULL
    #   }))
    # }
    # sig1 = rowSums(pchisq(files, 1, lower.tail = FALSE) <= ks[count,])/ncol(files)
    # sig2 = rowMeans(abs(EPCs))
    # sigs = matrix(cbind(sig1, sig2), nrow(files))
    param_names = character(0)
    for (j in 1:nvar) {
      for (i in 1:nvar) {
        # if (i != j) {
          param_names = c(param_names, sprintf("OUMod.A[%d,%d]", i, j))
        # }
      }
    }
    # param_names = sort(param_names)
    rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
    files = NULL
    EPCs = NULL
    for (file in rdss) {
      file_id = gsub("MI_|\\.RDS", "", basename(file))
      mi_full = abs(safe_read_vector(file, "MI.Full", param_names))
      epc = safe_read_vector(file, "EPC", param_names)
      files = cbind(files, mi_full)
      EPCs = cbind(EPCs, epc)
    }
    sig1 = rowSums(pchisq(files, 1, lower.tail = FALSE) <= ks[count, ], na.rm = TRUE) / ncol(files)
    sig2 = rowMeans(abs(EPCs), na.rm = TRUE)
    sigs = matrix(cbind(sig1, sig2), nrow(files))
    rownames(sigs) = rownames(files)
    SigThresh = sig1[which.max(sig1)] >= sig.thrsh
    if(SigThresh){
      param.to.add = which(names(which.max(sigs[which(sigs[,1] == max(sigs[,1])),2])) == rownames(sigs))
      cells = as.numeric(unlist(regmatches(rownames(files)[param.to.add], 
                                           gregexpr("\\d+", rownames(files)[param.to.add]))))
      DRIFT[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
      message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
      message(paste0("Completed Step ", count))
      unlink(paste0(directory, "/MIs/", "*"), recursive = TRUE, force = TRUE)
      count = count + 1
    }else{
      prune = 0
      message("IT'S PRUNING TIME!")
      while(prune < 1){
        models = list()
        rdss1 = list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE)
        for(prn in 1:length(rdss1)){
          file_id1 = gsub("Model_|\\.RDS", "", basename(rdss1[prn]))
          temp1 = tryCatch({readRDS(rdss1[prn])}, error = function(e) {
            message("Failed to read ", prn, ": ", e$message)
            NULL})
          drifts = subset(summary(temp1)$parameters, matrix == 'A')
          cells = matrix(
            as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
            ncol = 2, byrow = TRUE
          )
          row_indices = cells[, 1]
          col_indices = cells[, 2]
          temp.mat1 = matrix(NA, nvar, nvar)
          cells_str = as.character(cells)
          for(i in 1:nrow(cells)){
            temp.mat1[row_indices[i], col_indices[i]] = ifelse(abs(drifts[i,"Estimate"])/
                                                                 (drifts[i,"Std.Error"]*qnorm(0.95)) > 
                                                                 qnorm(0.95), TRUE, FALSE)
          }
          models[[prn]] = cbind(temp.mat1)
        }
        arr = simplify2array(models)
        true_count = apply(arr, c(1, 2), function(x) sum(x == TRUE, na.rm = TRUE))/length(rdss1)
        diag(true_count) = 1.00
        true_count = ifelse(true_count == 0, NA, true_count)
        true_count[DRIFT == "0"] = NA
        if(!any(true_count <= 0.70*sig.thrsh, na.rm = TRUE)){
          prune = 1
        }
        if(true_count[which.min(true_count)] <= 0.70*sig.thrsh){
          DRIFT[which.min(true_count)] = "0"
        }
        cl = makeCluster(cores, type = "PSOCK")
        clusterExport(cl, c("dataframe", "JPmx", "ME.var", "PE.var",
                            "varnames", "id", "time", "directory",
                            "nvar", "DRIFT", "count"),
                      envir = environment())
        clusterEvalQ(cl, {
          packs = list('dynr', 'OpenMx', 'qgraph')
          invisible(lapply(packs, require, character.only = T))
        })
        parLapply(cl, ids, function(i) {
          subset_dat = subset(dataframe, id == i)
          subset_dat[,varnames] = scale(subset_dat[,varnames])
          amat = mxMatrix("Full", nvar, nvar,
                          free   = DRIFT != "0",
                          name   = "A")
          bmat = mxMatrix('Zero', nvar, nvar, name='B')
          cdim = list(varnames, paste0('F', 1:nvar))
          cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
          dmat = mxMatrix('Zero', nvar, nvar, name='D')
          qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
          rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
          xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
          pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
          umat = mxMatrix('Zero', nvar, 1, name='u')
          tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
          osc = mxModel("OUMod", 
                        amat, bmat, cmat, dmat, qmat, 
                        rmat, xmat, pmat, umat, tmat,
                        mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                          'R', 'x0', 'P0', 'u', 'time'),
                        mxFitFunctionML(),
                        mxData(subset_dat, 'raw'))  
          analysis_result = tryCatch({
            fit = mxTryHard(osc)
          }, error = function(e) {
            message("Error for subject ", i, ": ", e$message)
          })
          if (!is.null(analysis_result)) {
            MIs = JPmx(analysis_result, matrices = "A")
            saveRDS(
              object = MIs,
              file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
            )
            saveRDS(
              object = analysis_result,
              file   = paste0(directory, "/Models/", "/Model_", i, ".RDS")
            )
          }
        })
        stopCluster(cl)
      }
      G.DRIFT = DRIFT
      output_path = file.path(directory, "Group Paths.png")
      png(filename = output_path, width = 800, height = 800)
      # Plot the community structure
      qgraph(t(G.DRIFT != "0") * 1, layout = "circle", labels = varnames, 
             edge.width = 5, diag = TRUE, edge.labels = "GROUP")
      dev.off()
      ks = ks[c(count:nrow(ks)),]
      message("Group Search Complete.")
      iterate = 1
    }
    cl = makeCluster(cores, type = "PSOCK")
    clusterExport(cl, c("dataframe", "JPmx", "ME.var", "PE.var",
                        "varnames", "id", "time", "directory",
                        "nvar", "DRIFT", "count"),
                  envir = environment())
    clusterEvalQ(cl, {
      packs = list('dynr', 'OpenMx', 'qgraph')
      invisible(lapply(packs, require, character.only = T))
    })
    parLapply(cl, ids, function(i) {
      subset_dat = subset(dataframe, id == i)
      subset_dat[,varnames] = scale(subset_dat[,varnames])
      amat = mxMatrix("Full", nvar, nvar,
                       free   = DRIFT != "0",
                       name   = "A")
      bmat = mxMatrix('Zero', nvar, nvar, name='B')
      cdim = list(varnames, paste0('F', 1:nvar))
      cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
      dmat = mxMatrix('Zero', nvar, nvar, name='D')
      qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
      rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
          xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
      pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
      umat = mxMatrix('Zero', nvar, 1, name='u')
      tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
      osc = mxModel("OUMod", 
                    amat, bmat, cmat, dmat, qmat, 
                    rmat, xmat, pmat, umat, tmat,
                    mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                      'R', 'x0', 'P0', 'u', 'time'),
                    mxFitFunctionML(),
                    mxData(subset_dat, 'raw'))  
      analysis_result = tryCatch({
        fit = mxTryHard(osc)
      }, error = function(e) {
        message("Error for subject ", i, ": ", e$message)
      })
      if (!is.null(analysis_result)) {
        MIs = JPmx(analysis_result, matrices = "A")
        saveRDS(
          object = MIs,
          file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
        )
        saveRDS(
          object = analysis_result,
          file   = paste0(directory, "/Models/", "/Model_", i, ".RDS")
        )
      }
    })
    stopCluster(cl)
  }
  ###---###---###---###---###---###---###
  # Subgrouping Stage
  ###---###---###---###---###---###---###
  if(sub.sig.thrsh == 1.00){
    memb = rep(1, length(ids))
    message("Subgrouping Disabled for Testing.")
  }else{
    # message("Beginning Subgrouping Stage")
    # rdss1 = list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE)
    # rdss2 = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
    # rdss = c(rdss1, rdss2)
    # models = list()
    # for (file in 1:(length(rdss)/2)) {
    #   file_id1 = gsub("Model_|\\.RDS", "", basename(rdss1[file]))
    #   file_id2 = gsub("MI_|\\.RDS", "", basename(rdss2[file]))
    #   temp1 = tryCatch({readRDS(rdss1[file])}, error = function(e) {
    #     message("Failed to read ", file, ": ", e$message)
    #     NULL})
    #   temp2 = tryCatch({readRDS(rdss2[file])}, error = function(e) {
    #     message("Failed to read ", file, ": ", e$message)
    #     NULL})$MI.Full
    #   temp3 = tryCatch({readRDS(rdss2[file])}, error = function(e) {
    #     message("Failed to read ", file, ": ", e$message)
    #     NULL})$EPC
    #   drifts = subset(summary(temp1)$parameters, matrix == 'A')
    #   cells = matrix(
    #     as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
    #     ncol = 2, byrow = TRUE
    #   )
    #   row_indices = cells[, 1]
    #   col_indices = cells[, 2]
    #   MI.cells = matrix(as.numeric(unlist(regmatches(names(c(temp2)), 
    #                                                  gregexpr("\\d+", names(c(temp2)))))), 
    #                     ncol = 2, byrow = TRUE)
    #   MI.cells = cbind(MI.cells, temp2, temp3)
    #   MI.cells[,3] = ifelse(MI.cells[,3] > qchisq(0.995, 1), MI.cells[,4], 0)
    #   temp.mat1 = temp.mat2 = matrix(NA, nvar, nvar)
    #   cells_str = as.character(cells)
    #   for(i in 1:nrow(cells)){
    #     temp.mat1[row_indices[i], col_indices[i]] = ifelse(abs(drifts[i,"Estimate"])/
    #                                                          (drifts[i,"Std.Error"]*qnorm(0.99)) > 
    #                                                          qnorm(0.99), drifts[i,"Estimate"], 0)
    #   }
    #   for(i in 1:nrow(MI.cells)){
    #     temp.mat2[MI.cells[i,1], MI.cells[i,2]] = MI.cells[i,4]
    #   }
    #   models[[file]] = cbind(temp.mat1, temp.mat2)
    # }
    # adj.mat = matrix(NA, length(models), length(models))
    # for(i in 1:length(models)){
    #   for(j in 1:length(models)){
    #     # Determine based on sign:
    #     adj.mat[i,j] = sum(sign(c(models[[i]])) == sign(c(models[[j]])), na.rm = TRUE)
    #     # Determine based on Distances:
    #     # adj.mat[i,j] = 1/(1 + norm(models[[i]] - models[[j]], type = "F"))
    #     # Determine based on magnitude:
    #     # adj.mat[i,j] = 
    #   }
    # }
    message("Beginning Subgrouping Stage")
    param_names = character(0)
    for (j in 1:nvar) {
      for (i in 1:nvar) {
          param_names = c(param_names, sprintf("OUMod.A[%d,%d]", i, j))
      }
    }
    # param_names = sort(param_names)
    safe_read_vector = function(file, element, param_names) {
      out = rep(NA_real_, length(param_names))
      names(out) = param_names
      vec = tryCatch({ c(readRDS(file)[[element]]) }, error = function(e) {
        message("Failed to read ", element, " from ", file, ": ", e$message)
        return(NULL)
      })
      if (!is.null(vec)) {
        intersect_names = intersect(names(vec), param_names)
        out[intersect_names] = vec[intersect_names]
      }
      return(out)
    }
    rdss1 = list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE)
    rdss2 = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
    models = list()
    for (file in 1:length(rdss1)) {
      file_id = gsub("Model_|\\.RDS", "", basename(rdss1[file]))
      model_obj = tryCatch({ readRDS(rdss1[file]) }, error = function(e) {
        message("Failed to read ", rdss1[file], ": ", e$message)
        return(NULL)
      })
      if (is.null(model_obj)) next
      mi_vec = abs(safe_read_vector(rdss2[file], "MI.Full", param_names))
      epc_vec = safe_read_vector(rdss2[file], "EPC", param_names)
      drift_table = subset(summary(model_obj)$parameters, matrix == "A")
      drift_cells = matrix(as.numeric(unlist(regmatches(drift_table$name, gregexpr("\\d+", drift_table$name)))), 
                           ncol = 2, byrow = TRUE)
      temp.mat1 = matrix(NA, nvar, nvar)
      for (i in 1:nrow(drift_cells)) {
        r = drift_cells[i, 1]
        c = drift_cells[i, 2]
        est = drift_table$Estimate[i]
        se = drift_table$Std.Error[i]
        zthr = qnorm(0.975)
        temp.mat1[r, c] = ifelse(abs(est / (se * zthr)) > zthr, est, 0)
      }
      
      temp.mat2 = matrix(NA, nvar, nvar)
      for (name in names(mi_vec)) {
        indices = as.numeric(unlist(regmatches(name, gregexpr("\\d+", name))))
        if (length(indices) == 2) {
          r = indices[1]
          c = indices[2]
          mi_val = mi_vec[name]
          epc_val = epc_vec[name]
          temp.mat2[r, c] = ifelse(!is.na(mi_val) && mi_val > qchisq(0.99, 1), epc_val, 0)
        }
      }
      
      models[[file_id]] = cbind(temp.mat1, temp.mat2)
    }
    
    adj.mat = matrix(NA, length(models), length(models))
    
    for (i in 1:length(models)) {
      for (j in 1:length(models)) {
        adj.mat[i, j] = sum(sign(c(models[[i]])) == sign(c(models[[j]])), na.rm = TRUE)
      }
    }
    
    adj.mat = adj.mat - min(adj.mat)
    diag(adj.mat) = 0
    g = graph_from_adjacency_matrix(adj.mat, mode = "undirected", weighted = TRUE, diag = FALSE)
    walktrap_comm = cluster_walktrap(g)
    output_path = file.path(directory, "walktrap_community_plot.png")
    png(filename = output_path, width = 800, height = 800)
    # Plot the community structure
    plot(walktrap_comm, g, vertex.size = 15, 
         vertex.label = 1:length(unique(dataframe[,id])), 
         main = "Walktrap Community Detection")
    dev.off()
    memb = membership(walktrap_comm)
  }
  print(memb)
  ###---###---###---###---###---###---###
  # Iteration Steps - Subgroup
  ###---###---###---###---###---###---###
  for(subgroup in sort(unique(memb))){
    DRIFT = G.DRIFT
    dir.create(paste0(directory, "/Models/Subgroup ", subgroup, "/"), showWarnings = FALSE)
    iterate = 0
    count = 1
    m = (nvar^2)-nvar
    sg.ks = matrix(NA, m, 1)
    for(k in 1:m){
      sg.ks[k,] = (k/m) * S.Galpha
    }
    sg.ks = matrix(sort(sg.ks, TRUE), m, 1)
    if(ben.hoch == FALSE){
      sg.ks = matrix(S.Galpha, nrow(sg.ks), 1)
    }
    memb.id = cbind(unique(dataframe$id), memb)
    while(iterate < 1){
      new.data = subset(dataframe, id %in% subset(memb.id[,1], memb == subgroup))
      # param_names = sort(param_names)
      valid_ids = unique(new.data$id)
      all_rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
      rdss = all_rdss[gsub("MI_|\\.RDS", "", basename(all_rdss)) %in% valid_ids]
      # rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
      files = NULL
      EPCs = NULL
      for (file in rdss) {
        file_id = gsub("MI_|\\.RDS", "", basename(file))
        mi_full = abs(safe_read_vector(file, "MI.Full", param_names))
        epc = safe_read_vector(file, "EPC", param_names)
        files = cbind(files, mi_full)
        EPCs = cbind(EPCs, epc)
      }
      sig1 = rowSums(pchisq(files, 1, lower.tail = FALSE) <= sg.ks[count, ], na.rm = TRUE) / ncol(files)
      sig2 = rowMeans(abs(EPCs), na.rm = TRUE)
      sigs = matrix(cbind(sig1, sig2), nrow(files))
      rownames(sigs) = rownames(files)
      SigThresh = sig1[which.max(sig1)] >= sub.sig.thrsh
      if(SigThresh){
        param.to.add = which(names(which.max(sigs[which(sigs[,1] == max(sigs[,1])),2])) == rownames(sigs))
        cells = as.numeric(unlist(regmatches(rownames(files)[param.to.add], 
                                             gregexpr("\\d+", rownames(files)[param.to.add]))))
        DRIFT[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
        message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
        message(paste0("Completed Step ", count))
        count = count + 1
      }else{
      if(count > 1 & !sub.sig.thrsh == 1.00){
        prune = 0
        message("IT'S PRUNING TIME!")
        while(prune < 1){
          models = list()
          rdss1 = list.files(paste0(directory, "/Models/Subgroup ", subgroup, "/"), pattern = "\\.RDS$", full.names = TRUE)
          for(prn in 1:length(rdss1)){
            file_id1 = gsub("Model_|\\.RDS", "", basename(rdss1[prn]))
            temp1 = tryCatch({readRDS(rdss1[prn])}, error = function(e) {
              message("Failed to read ", prn, ": ", e$message)
              NULL})
            drifts = subset(summary(temp1)$parameters, matrix == 'A')
            cells = matrix(
              as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
              ncol = 2, byrow = TRUE
            )
            row_indices = cells[, 1]
            col_indices = cells[, 2]
            temp.mat1 = matrix(NA, nvar, nvar)
            cells_str = as.character(cells)
            for(pp in 1:nrow(cells)){
              temp.mat1[row_indices[pp], col_indices[pp]] = ifelse(abs(drifts[pp,"Estimate"])/
                                                                   (drifts[pp,"Std.Error"]*qnorm(0.975)) > 
                                                                   qnorm(0.975), TRUE, FALSE)
            }
            models[[prn]] = cbind(temp.mat1)
          }
          arr = simplify2array(models)
          true_count = apply(arr, c(1, 2), function(x) sum(x == TRUE, na.rm = TRUE))/length(rdss1)
          true_count[which(DRIFT != "0")] = 1.00
          true_count = ifelse(true_count == 0, NA, true_count)
          true_count[DRIFT == "0"] = NA
          if(!any(true_count <= 0.70*sub.sig.thrsh, na.rm = TRUE)){
            prune = 1
          }
          if(true_count[which.min(true_count)] <= 0.70*sub.sig.thrsh){
            DRIFT[which.min(true_count)] = "0"
          }
          cl = makeCluster(cores, type = "PSOCK")
          clusterExport(cl, c("new.data", "JPmx", "ME.var", "PE.var",
                              "varnames", "id", "time", "directory",
                              "nvar", "DRIFT", "count",
                              "subgroup"),
                        envir = environment())
          clusterEvalQ(cl, {
            packs = list('dynr', 'OpenMx', 'qgraph')
            invisible(lapply(packs, require, character.only = T))
          })
          parLapply(cl, valid_ids, function(i) {
            subset_dat = subset(new.data, id == i)
            amat = mxMatrix("Full", nvar, nvar,
                            free   = DRIFT != "0",
                            name   = "A")
            bmat = mxMatrix('Zero', nvar, nvar, name='B')
            cdim = list(varnames, paste0('F', 1:nvar))
            cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
            dmat = mxMatrix('Zero', nvar, nvar, name='D')
            qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
            rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
            xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
            pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
            umat = mxMatrix('Zero', nvar, 1, name='u')
            tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
            osc = mxModel("OUMod", 
                          amat, bmat, cmat, dmat, qmat, 
                          rmat, xmat, pmat, umat, tmat,
                          mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                            'R', 'x0', 'P0', 'u', 'time'),
                          mxFitFunctionML(),
                          mxData(subset_dat, 'raw'))  
            analysis_result = tryCatch({
              fit = mxTryHard(osc)
            }, error = function(e) {
              message("Error for subject ", i, ": ", e$message)
            })
            if (!is.null(analysis_result)) {
              MIs = JPmx(analysis_result, matrices = "A")
              saveRDS(
                object = MIs,
                file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
              )
              saveRDS(
                object = analysis_result,
                file   = paste0(directory, "/Models/Subgroup ", subgroup, "/Model_", i, ".RDS")
              )
            }
          })
          stopCluster(cl)
        }
      }
        message(paste0("Subgroup Search ", subgroup," of ", max(memb)," Complete."))
        if(subgroup.model == TRUE){
          subgroup.data = NULL
          for(i in unique(new.data[,id])){
            temp = subset(new.data, id == i)
            temp = rbind(temp, rep(NA, ncol(temp)))
            subgroup.data = rbind(subgroup.data, temp)
          }
          ###---###---###---###---###---###---###---###---
          ## ONLY WORKS FOR OMID DATA###---###---###---###
          ###---###---###---###---###---###---###---###---
          # subgroup.data$time = 0:(nrow(subgroup.data)-1)
          nsubjs = nrow(subgroup.data)
          days = floor((nsubjs - 1)/5)
          offsets = seq(0, by = 1/8, length.out = 5)
          subgroup.data$Time = as.vector(sapply(0:days, function(d) d + offsets))[1:nsubjs]
          ###---###---###---###---###---###---###---###---
          ## ONLY WORKS FOR OMID DATA###---###---###---###
          ###---###---###---###---###---###---###---###---          
          message(paste0("Fitting Parameterized Model of Subgroup ", subgroup))
          amat = mxMatrix('Full', nvar, nvar, DRIFT != "0", 
                          name = 'A')
          bmat = mxMatrix('Zero', nvar, nvar, name='B')
          cdim = list(varnames, paste0('F', 1:nvar))
          cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
          dmat = mxMatrix('Zero', nvar, nvar, name='D')
          qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
          rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
              xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
          pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
          umat = mxMatrix('Zero', nvar, 1, name='u')
          tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
          osc = mxModel("OUMod", 
                        amat, bmat, cmat, dmat, qmat, 
                        rmat, xmat, pmat, umat, tmat,
                        mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                          'R', 'x0', 'P0', 'u', 'time'),
                        mxFitFunctionML(),
                        mxData(subgroup.data, 'raw'))  
          fit = mxTryHard(osc)
          sum.fit = summary(fit)
          effects = matrix(0, nvar, nvar)
          for(i in 1:nrow(sum.fit$parameters)){
            effects[sum.fit$parameters$col[i],sum.fit$parameters$row[i]] = sum.fit$parameters$Estimate[i]
          }
          # Determine significance stars
          sig = ifelse(abs(sum.fit$parameters$Estimate) > qnorm(0.975) * sum.fit$parameters$Std.Error, "*", "ns")
          
          # Format edge labels
          vals = cbind(round(sum.fit$parameters$Estimate, 2), sig)
          edge_labs = paste0(vals[, 2], " (", vals[, 1], ")")
          
          # Logical masks
          shared = G.DRIFT != "0" & DRIFT != "0"
          group_only = G.DRIFT == "0" & DRIFT != "0"
          colors = character(length = length(DRIFT))
          colors[shared] = "gray"
          colors[group_only & effects > 0] = "blue"
          colors[group_only & effects < 0] = "red"
          colors[group_only & effects == 0] = "black"
          colors[DRIFT == "0"] = "black"
          colors = c(colors)
          
          # Output path
          output_path = file.path(directory, paste0("Models/Subgroup ", subgroup, "/Subgroup ", subgroup, " Params.png"))
          
          # Plot
          png(filename = output_path, width = 800, height = 800)
          qgraph(effects, layout = "circle", labels = varnames, 
                 edge.width = 1, diag = TRUE, edge.labels = edge_labs,
                 theme = "colorblind", edge.color = colors, fade = TRUE)
          dev.off()
          
          for(ints in time.intervals){
            delt = (round(expm(effects * ints), 3))
            output_path = file.path(paste0(directory, "/Models/Subgroup ", subgroup, "/Subgroup ", subgroup, " Delta_t = ", ints,".png"))
            png(filename = output_path, width = 800, height = 800)
            # Plot the subgroup structure
            qgraph(delt, layout = "circle", labels = varnames, fade = TRUE,
                   edge.width = 1, diag = TRUE, edge.labels = delt, maximum = 1.00,
                   theme = "colorblind", title = paste0("Subgroup ", subgroup, "; Delta_t = ", ints))
            dev.off()
          }  
          if (!is.null(fit)) {
            saveRDS(
              object = fit,
              file   = paste0(directory, "/Models/Subgroup ", subgroup,"/Subgroup_", subgroup, "Model.RDS")
            )
          }
        }
        message(paste0("Beginning Individual Model Fitting for Subgroup Members"))
        m = (nvar^2)-sum(DRIFT!="0")
        nks = matrix(NA, m, 1)
        for(k in 1:m){
          nks[k,] = (k/m) * Ialpha
        }
        nks = matrix(sort(nks, TRUE), m, 1)
        SG.DRIFT = DRIFT
        output_path = file.path(paste0(directory, "/Models/Subgroup ", subgroup, "/Subgroup ", subgroup, " Paths.png"))
        png(filename = output_path, width = 800, height = 800)
        # Plot the subgroup structure
        qgraph(t(abs(((G.DRIFT != "0") * 1) - ((DRIFT != "0") * 1))), layout = "circle", labels = varnames, 
               edge.width = 5, diag = TRUE, edge.labels = paste0("SG-", subgroup))
        dev.off()
        cl = makeCluster(cores, type = "PSOCK")
        clusterExport(cl, c("new.data", "JPmx", "ME.var", "PE.var",
                            "varnames", "id", "time", "directory",
                            "nvar", "SG.DRIFT", "count",
                            "subgroup", "nks"),
                      envir = environment())
        clusterEvalQ(cl, {
          packs = list('dynr', 'OpenMx', 'qgraph')
          invisible(lapply(packs, require, character.only = T))
        })
        parLapply(cl, valid_ids, function(i) {
          DRIFT = SG.DRIFT
          subset_dat = subset(new.data, id == i)
          amat = mxMatrix("Full", nvar, nvar,
                           free   = DRIFT != "0",
                           name   = "A")
          bmat = mxMatrix('Zero', nvar, nvar, name='B')
          cdim = list(varnames, paste0('F', 1:nvar))
          cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
          dmat = mxMatrix('Zero', nvar, nvar, name='D')
          qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
          rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
              xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
          pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
          umat = mxMatrix('Zero', nvar, 1, name='u')
          tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
          osc = mxModel("OUMod", 
                        amat, bmat, cmat, dmat, qmat, 
                        rmat, xmat, pmat, umat, tmat,
                        mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                          'R', 'x0', 'P0', 'u', 'time'),
                        mxFitFunctionML(),
                        mxData(subset_dat, 'raw'))  
          fit = mxTryHard(osc)
          fit2 = fit
          optimization = 0
          count = 1
          while(optimization < 1){
            if(count > 1){
              fit2 = fit
              fit = mxTryHard(osc)
              if(any(is.na(summary(fit)$parameters$Std.Error))){
                fit = fit2
                optimization = 1
                break
              }
            }
            MIs = JPmx(fit, matrices = "A")
          #   if(is.null(MIs) | is.null(MIs$MI.Full)){optimization = 1; fit = fit2}
          #   if(abs(MIs$MI.Full)[which.max(abs(MIs$MI.Full))] >= qchisq(1-nks[count,], df = 1)){
          #     cells = as.numeric(unlist(regmatches(names(which.max(MIs$MI.Full)),
          #                                          gregexpr("\\d+", names(which.max(MIs$MI.Full))))))
          #     
          #     osc$A$free[cells[1], cells[2]] = TRUE
          #     osc$A$labels[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
          #     message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
          #     MIs = NULL
          #     count = count + 1
          #     if(sum(osc$A$free) == nvar^2){
          #       optimization = 1
          #     }
          #   }else{
          #     optimization = 1
          #   }
          # }
            
            
            if (is.null(MIs) || is.null(MIs$MI.Full)) {
              optimization = 1
              fit = fit2
            } else if (length(MIs$MI.Full) == length(names(MIs$MI.Full))) {
              max_idx = which.max(abs(MIs$MI.Full))
              max_val = abs(MIs$MI.Full[max_idx])
              
              if (max_val >= qchisq(1 - nks[count, ], df = 1)) {
                max_name = names(MIs$MI.Full)[max_idx]
                cells = as.numeric(unlist(regmatches(max_name, gregexpr("\\d+", max_name))))
                
                osc$A$free[cells[1], cells[2]] = TRUE
                osc$A$labels[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
                
                message(paste0("Adding drift parameter A[", cells[1], ",", cells[2], "]"))
                
                MIs = NULL
                count = count + 1
                
                if (sum(osc$A$free) == nvar^2) {
                  optimization = 1
                }
              } else {
                optimization = 1
              }
            } else {
              message(paste0("Malformed MI.Full — names and values mismatch. Skipping subject ", i))
              optimization = 1
              fit = fit2
            }
            
          }
          message("Pruning Stage.")
          prune = 0
          count = 1
          stat.sig1 = summary(fit)$parameters
          while (prune < 1) {
            stat.sig = summary(fit)$parameters
            ests = matrix(0, nvar, nvar)
            for(jj in 1:nrow(stat.sig)){
              ests[stat.sig[,3][jj], stat.sig[,4][jj]] = stat.sig[jj,5]
            }
            prunable = stat.sig[stat.sig[,5] %in% ests[setdiff(which(ests != 0), 
                                                             which(ests != 0 & DRIFT != "0"))],]
            if (nrow(prunable) == 0 | is.null(nrow(prunable))) {
              message("No unprotected parameters left to prune.")
              break
            }
            se = prunable[, 6] * qnorm(0.975)
            z.scores = abs(prunable[, 5]) / se
            min_z_index = which.min(z.scores)
            if(length(z.scores[min_z_index]) == 0){
              message("All Prunable Paths Removed.")
              break
              prune = 1
            }
            if (z.scores[min_z_index] < qnorm(0.975)) {
              this = prunable[min_z_index,]
              cells = matrix(c(this$row, this$col), 1, 2)
              osc$A$free[cells[1,1], cells[1,2]] = FALSE
              fit = mxTryHard(osc)
              message(paste0("NOTE: Pruned drift parameter: A[", cells[1,1], ",", cells[1,2], "]!"))
            } else {
              message("No Pruning Conducted.")
              prune = 1
            }
          }
          if (!is.null(fit)) {
            saveRDS(
              object = fit,
              file   = paste0(directory, "/Models/Individuals/FinalModel_", i, ".RDS")
            )
            stat.sig = summary(fit)$parameters
            ests = matrix(0, nvar, nvar)
            for(jj in 1:nrow(stat.sig)){
              ests[stat.sig[,3][jj], stat.sig[,4][jj]] = stat.sig[jj,5]
            }
            ests = t(ests)
            png(filename = paste0(directory, "/Models/Individuals/FinalModel_", i, ".PNG"), 
                width = 800, height = 800)
            qgraph(ests, layout = "circle", labels = varnames, 
                   edge.width = 1, diag = TRUE, edge.labels = round(c(ests), 2),
                   theme = "colorblind", fade = FALSE)
            dev.off()
            unlink(list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE))
            unlink(list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE))
            unlink(list.files(paste0(directory, "/Models/Subgroup ",subgroup,"/"), pattern = "\\.RDS$", full.names = TRUE))
          }
        })
        stopCluster(cl)
        iterate = 1
      }
      cl = makeCluster(cores, type = "PSOCK")
      clusterExport(cl, c("new.data", "JPmx", "ME.var", "PE.var",
                          "varnames", "id", "time", "directory",
                          "nvar", "DRIFT", "count",
                          "subgroup"),
                    envir = environment())
      clusterEvalQ(cl, {
        packs = list('ctsem', 'ctsemOMX', 'dynr', 'OpenMx', 'qgraph')
        invisible(lapply(packs, require, character.only = T))
      })
      # qmat = mxMatrix('Diag', nvar, nvar, FALSE, PE.var, name='Q')
      # rmat = mxMatrix('Diag', nvar, nvar, FALSE, ME.var, name='R')
      parLapply(cl, valid_ids, function(i) {
        subset_dat = subset(new.data, id == i)
        amat = mxMatrix("Full", nvar, nvar,
                         free   = DRIFT != "0",
                         name   = "A")
        bmat = mxMatrix('Zero', nvar, nvar, name='B')
        cdim = list(varnames, paste0('F', 1:nvar))
        cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
        dmat = mxMatrix('Zero', nvar, nvar, name='D')
        qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
        rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
            xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
        pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
        umat = mxMatrix('Zero', nvar, 1, name='u')
        tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
        osc = mxModel("OUMod", 
                      amat, bmat, cmat, dmat, qmat, 
                      rmat, xmat, pmat, umat, tmat,
                      mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                        'R', 'x0', 'P0', 'u', 'time'),
                      mxFitFunctionML(),
                      mxData(subset_dat, 'raw'))  
        analysis_result = tryCatch({
          fit = mxTryHard(osc)
        }, error = function(e) {
          message("Error for subject ", i, ": ", e$message)
        })
        if (!is.null(analysis_result)) {
          MIs = JPmx(analysis_result, matrices = "A")
          saveRDS(
            object = MIs,
            file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
          )
          saveRDS(
            object = analysis_result,
            file   = paste0(directory, "/Models/Subgroup ", subgroup, "/Model_", i, ".RDS")
          )
        }
      })
      stopCluster(cl)
    }
  }
  message(paste0("Subgrouping with Continuous-Time GIMME Complete. Find networks in ", directory, "."))
  unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
  if(sub.sig.thrsh == 1){
    return(message("Continuous-Time S-GIMME Complete."))
  }else{return(walktrap_comm)}
}

