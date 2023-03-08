\name{PT41F.t15.tmb.MP-package}
\alias{PT41F.t15.tmb.MP-package}
\alias{PT41F.t15.tmb.MP}
\docType{package}
\title{
\packageTitle{PT41F.t15.tmb.MP}
}
\description{
\packageDescription{PT41F.t15.tmb.MP}
}
\details{
}

\author{
  \packageAuthor{PT41F.t15.tmb.MP}

  Maintainer: \packageMaintainer{PT41F.t15.tmb.MP}
}

\seealso{
  \code{\link{exampleBET_jointCPUE_DataFilePath}}
  \code{\link{assessmentCPUE}}
  \code{\link{assessMP}}
  \code{\link{runShinyMP}}
}

\examples{
# Create assessment CPUE series from region based CPUE series
exportFileName  <- "BETOMcpue2022.csv"
MPDatFile       <- "BETMPcpue2022.csv"
IsSeasonal      <- FALSE

assessmentCPUE(exampleBET_jointCPUE_DataFilePath(), exportFileName, MPDatFile, IsSeasonal)

# Get path to example dat file
mp_example_data_file <- system.file("extdata", "mp_example_data.csv", package="PT41F.t15.tmb.MP")

print(mp_example_data_file)

# Assess the MP recomendation using the example data
results <- assessMP(mp_example_data_file)

# list the results variables
ls(results)

# print the MP recommendation results
print(results$TAC)
print(results$B)
print(results$q)
print(results$Depletion)

# display the TMB model data
print(results$ModelData)

# plot the MP recommendation evaluation plots
print(results$plots)

\dontrun{
# run the shiny app incorporating the MP
runShinyMP()
}
}

\keyword{package}
