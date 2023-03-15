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
library(PT41F.t15.tmb.MP)

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
# read catch data from IOTC (assumes readxl package is installed) and create a csv file
# containing necessary data to run the MP, and then run the MP assessment
library(readxl)

download.file('https://iotc.org/sites/default/files/documents/2022/09/IOTC-2022-WPTT24AS-DATA03-NC.zip',"catch_data.zip")
unzip('catch_data.zip')
catch_data <- read_xlsx("IOTC-2022-WPTT24(AS)-DATA03-NC.xlsx", sheet="Summary_Resume", skip = 9)

# Only take the data we want
catch <- subset(catch_data, select = c("Year/An","BET"))

# need to exclude the "Grand Total" row
rows <- which(catch[["Year/An"]] != "Grand Total")
catch <- catch[rows,]

# only data after 1978
catch <- catch[catch[["Year/An"]] > 1978,]

# coerce to the correct data type
catch[["Year/An"]] <- as.numeric(catch[["Year/An"]])
catch$BET <- as.numeric(catch$BET)

# cpue data
cpue <- read.csv("BETMPcpue2022.csv", header=T)

# combine data to make input table for MP
MP.Input <- data.frame(Year=1979:2021, Catch=NA, CPUE=NA, TAC=NA)
MP.Input[catch[["Year/An"]] - 1978,]$Catch <- catch$BET
MP.Input[cpue$yr - 1978,]$CPUE <- cpue$cpue

# Add a final year TAC value so that the MP can determine the +/-15% change limit
LastTAC <- MP.Input[MP.Input$Year == max(MP.Input$Year), ]$Catch
MP.Input[MP.Input$Year == max(MP.Input$Year), ]$TAC <- LastTAC

# Save the MP data file
write.csv(MP.Input, "mp_data_file.csv", row.names=F)

# run MP assessment and print some results
results2 <- assessMP("mp_data_file.csv")

print(results2$TAC)
print(results2$plots)

# run the shiny app incorporating the MP
runShinyMP()
}
}

\keyword{package}
