\name{BuildHMDP1}
\alias{BuildHMDP1}
\title{Build the HMDP (3 levels with shared linking) using the C++ binary writer.}
\usage{
BuildHMDP1(filePrefix, param, dlms, dglmParam)
}
\arguments{
  \item{filePrefix}{Prefix used by the binary files storing
  the MDP model.}

  \item{param}{Model parameters a list created using
  \code{\link{setParameters}}.}

  \item{dlms}{A list a long as the number of rations where
  each entry is a DLM created using \code{\link{iniDLM}}
  and \code{\link{buildDLM}}.}

  \item{dglmParam}{Model parameters a list created using
  \code{\link{iniDGLM}}.}
}
\value{
Build log (character).
}
\description{
The MDP specified at level 3 is generated for each ration
for a dummy state at the first stage at level 2.
}
\author{
Lars Relund \email{lars@relund.dk}
}

