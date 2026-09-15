# tmb_name <- "nsddm_TMBExports"
# tmb_flags <- commandArgs(trailingOnly = TRUE)

# if(file.exists(paste0(tmb_name, ".cpp"))) {
#   if(length(tmb_flags) == 0) tmb_flags <- ""
#   TMB::compile(file = paste0(tmb_name, ".cpp"),
#                PKG_CXXFLAGS = tmb_flags,
#                safebounds = FALSE, safeunload = FALSE)
#   file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
#             to = "..", overwrite = TRUE)
# }

tmb_name  <- "nsddm_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
    
    sysname <- Sys.info()["sysname"]
    
    if( sysname == "Darwin" ) {
        openmp_flags <- "-Xclang -fopenmp"
        openmp_libs  <- "-lomp"
    } else {
        # Windows & Linux:
        openmp_flags <- "-fopenmp"
        openmp_libs  <- "-lgomp"
    }
    
    if ( sysname == "Windows" ) {
        openmp_flags <- paste(openmp_flags, "-Wa,-mbig-obj")
    }
    
    if(length(tmb_flags) == 0) tmb_flags <- openmp_flags
    
    TMB::compile(file = paste0(tmb_name, ".cpp"),
                 PKG_CXXFLAGS = tmb_flags,
                 PKG_LIBS     = openmp_libs,
                 safebounds = FALSE, safeunload = FALSE)
    
    file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
              to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
