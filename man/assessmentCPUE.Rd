\name{assessmentCPUE}
\alias{assessmentCPUE}
\title{
  Construct the assessment OM and MP CPUE csv files. 
}
\description{
  Calling \code{assessmentCPUE} constructs assessment OM and MP CPUE csv files.
}
\usage{
  assessmentCPUE(CPUEDataFile, exportFileName, MPDatFile, IsSeasonal=FALSE)
}
\arguments{
  \item{CPUEDataFile}{
    \code{CPUEDataFile} is a csv data file containing raw CPUE data for the Indian Ocean BET fishery.
  }
  \item{exportFileName}{
    \code{exportFileName} is the output file and path where the OM CPUE csv data is written to.
  }
  \item{MPDatFile}{
    \code{MPDatFile} is the output file and path where the MP CPUE csv data is written to.
  }
  \item{IsSeasonal}{
    \code{IsSeasonal} is a boolean that when \code{TRUE} results in seasonal CPUE output series.
  }
}
\details{
  Calling \code{assessmentCPUE} constructs assessment OM and MP CPUE csv files from the raw Indian Ocean BET fishery CPUE data. 
  It expects the contents of \code{CPUEDataFile} to contain CPUE data for four regions of the BET fishery with columns \code{YrQtr},
  \code{R1N}, \code{R1S}, \code{R2} and \code{R3}. It returns a single unified CPUE series that is constructed by a weighted 
  combination of those regional CPUE series with a 1\% catchability trend applied. \code{exportFileName} is intended as the CPUE 
  series input file to the \emph{Stock Synthesis} operating models used to construct the MSE evaluation model set and \code{MPDatFile}
  is the CPUE series input file to be used for the \emph{niMSE-IO-BET-YFT} MSE MP evaluation. See the function source code for a 
  complete understanding of how the series are combined.
}
\value{
  This function does not return a result.
}
\author{
  \packageAuthor{PT41F.t15.tmb.MP}

  Maintainer: \packageMaintainer{PT41F.t15.tmb.MP}
}

\seealso{
  \code{\link{exampleBET_jointCPUE_DataFilePath}}
  \code{\link{PT41F.t15.tmb.MP-package}}
}

\note{
  see \code{\link{PT41F.t15.tmb.MP-package}} for examples of use.
}
