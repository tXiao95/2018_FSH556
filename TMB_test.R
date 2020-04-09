######################
# Simulate data for a linear mixed model with random intercepts
######################


#' Steps to install TMB: https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Steps-to-install-TMB
#' Rtools download: https://cran.r-project.org/bin/windows/Rtools/
#' R handholding Rtools guide: https://thecoatlessprofessor.com/programming/cpp/installing-rtools-for-compiled-code-via-rcpp/
#' RStudio issues: https://github.com/rstudio/rstudio/issues/3563

set.seed(1)
Factor = rep( 1:10, each=10)
Z = rnorm( length(unique(Factor)), mean=0, sd=1)

X0 = 0
Y = Z[Factor] + X0 + rnorm( length(Factor), mean=0, sd=1)

######################
# Run in TMB
######################

#install.packages("TMB")
library(TMB)
Version = "linear_mixed_model"

# Download CPP file
setwd( tempdir() )
download.file( url="https://raw.githubusercontent.com/James-Thorson/mixed-effects/master/linear_mixed_model/linear_mixed_model.cpp", destfile="linear_mixed_model.cpp", method="auto")


# TH: This is the line that ends up causing problems. It basically finds the C++ compiler to compile 
#' this .cpp files, and if you're on Windows and get an error, then Rtools was not installed correctly. 
#' When I originally tried to run this, I kept getting an error because my Rtools installation was in 
#' C:\RBuildTools\3.5\mingw_64\bin, which I think is the default filepath of the Rtools installation when 
#' done through the pkgbuild package. However, my Makeconf file for R (which I cannot write to even as 
#' an admin because the file permissions are 555) has the BINPREF environemnt variable, which is where 
#' R looks for the compiler, at 'C:\Rtools\mingw_64\bin' so I ended up just re-installing it and downloading
#' from source rather than pkgbuild. Link above.
TMB::compile( paste0(Version,".cpp") )

# Generate inputs for TMB
Data = list( "n_data"=length(Y), "n_factors"=length(unique(Factor)), "Factor"=Factor-1, "Y"=Y)
Parameters = list( "X0"=-10, "log_SD0"=2, "log_SDZ"=2, "Z"=rep(0,Data$n_factor) )
Random = c("Z")

# Build TMB object
dyn.load( dynlib(Version) )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Check that TMB is working - value should be 313.4137
Obj$fn( Obj$par )
