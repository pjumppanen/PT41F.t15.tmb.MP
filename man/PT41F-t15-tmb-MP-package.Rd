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

# plot the MP cpue fit plot
print(results$plots$cpue_plot)

\dontrun{
# read catch data from IOTC (assumes readxl package is installed) and create a csv file
# containing necessary data to run the MP, and then run the MP assessment
library(readxl)
library(ggplot2)
library(reshape2)

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

# Add a guestimate of final year TAC value so that the MP can determine the +/-15% change limit
TAC.years <- (max(MP.Input$Year)-3):max(MP.Input$Year)
LastTAC <- mean(MP.Input[MP.Input$Year \%in\% TAC.years, ]$Catch)
MP.Input[MP.Input$Year \%in\% TAC.years, ]$TAC <- LastTAC

# Save the MP data file
write.csv(MP.Input, "mp_data_file.csv", row.names=F)

# run MP assessment and print some results
results2 <- assessMP("mp_data_file.csv")

print(results2$TAC)

data <- MP.Input
RecommendedTAC <- rep(NA, nrow(data))
RecommendedTAC[length(RecommendedTAC)] <- results2$TAC
data <- cbind(data, "Recommended TAC"=RecommendedTAC)

data_melt <- melt(data[, c("Year", "Catch", "TAC", "Recommended TAC")], id.vars='Year', value.name='Catch')
colors <- c("Catch"="#00345D", "TAC"="#00A9CE", "Recommended TAC"="#000080")
shapes <- c("Catch"=NA, "TAC"=NA, "Recommended TAC"=1)
types <- c("Catch"=1, "TAC"=1, "Recommended TAC"=0)

print(ggplot(data_melt, aes(x=Year, y=Catch / 1000, linetype=variable, color=variable, shape=variable)) +
             geom_line(size=2) +
             geom_point(size=6) + 
             scale_linetype_manual(values=types) + 
             scale_shape_manual(values=shapes) + 
             scale_color_manual(values=colors) + 
             scale_y_continuous(limits=c(0, 1.2 * max(data_melt$Catch / 1000)), expand=c(0, 0)) +
             xlab('') + 
             ylab('Catch (x 1000)\n ') + 
             theme_bw() + 
             theme(legend.title=element_blank()))

print(results2$plots$cpue_plot)

# run the shiny app incorporating the MP
runShinyMP()
}
}

\keyword{package}
