\name{BD_Obs}
\alias{BD_Obs}
\docType{data}
\title{Observed Data from Ground Stations}
\description{
  A dataset containing daily precipitation observations from multiple ground stations.
}

\format{
  A \code{data.table} or \code{data.frame} with \emph{n} rows (dates) and \emph{m+1} columns (stations + Date):
  \itemize{
    \item \strong{Date}: A \code{Date} object representing the observation date.
    \item \strong{STxxx}: Numeric values representing observed precipitation measurements for each station, where \code{STxxx} is the unique station identifier.
  }
}

\usage{data(BD_Obs)}

\examples{
  # Load dataset
  data(BD_Obs)

  # Check structure
  str(BD_Obs)

  # Display first rows
  head(BD_Obs)
}

\details{
  The dataset follows the structure:
  \preformatted{
    > BD_Obs
    # A data.table or data.frame with n rows (dates) and m+1 columns (stations + Date)
       Date        ST001  ST002  ST003  ST004  ...
       <date>      <dbl>  <dbl>  <dbl>  <dbl>  ...
    1  2015-01-01    0      0      0      0    ...
    2  2015-01-02    0      0      0     0.2   ...
    3  2015-01-03   0.1     0      0     0.1   ...
  }

  - Each station column contains numeric values representing observed measurements.
  - The column names (station identifiers) must be unique and match those in \code{BD_Coord$Cod} to ensure proper spatial referencing.
}

\seealso{\code{BD_Coord} for station metadata.}
