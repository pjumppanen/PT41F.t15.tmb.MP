\name{exampleBET_jointCPUE_DataFilePath}
\alias{exampleBET_jointCPUE_DataFilePath}
\title{
  Return the fully qualified path to an example BET joint CPUE csv file
}
\description{
  Calling \code{exampleBET_jointCPUE_DataFilePath} returns the fully qualified path to an example BET joint CPUE csv file.
}
\usage{
  exampleBET_jointCPUE_DataFilePath()
}
\arguments{
}
\details{
  returns the fully qualified path to the \code{IOTC_BET_JointCPUE_submitted.csv} file for 2022 bundled with this package. The returned file path
  can then be used as the \code{CPUEDataFile} argument to the \code{\link{assessmentCPUE}} function.
}
\value{
  Fully qualified path to IOTC_BET_JointCPUE_submitted.csv.
}
\author{
  \packageAuthor{PT41F.t15.tmb.MP}

  Maintainer: \packageMaintainer{PT41F.t15.tmb.MP}
}

\seealso{
  \code{\link{PT41F.t15.tmb.MP-package}}
}

\note{
  see \code{\link{PT41F.t15.tmb.MP-package}} for examples of use.
}
