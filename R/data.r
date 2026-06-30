#' Entrecampos Air Quality Dataset
#'
#' This dataset contains interval data of air pollutants' concentrations, including min-max values and microdata.
#' This air quality dataset was obtained from a monitoring station in Entrecampos, Lisbon.
#' It is composed of 9 pollutants' concentration measures in µg/m3 during the years 2019, 2020, and 2021: sulphur dioxide (SO2), particles < 10µm, ozone (O3), nitrogen dioxide (NO2), carbon monoxide (CO), benzene (C6H6), particles < 2.5µm, nitrogen oxides (NOx), and nitrogen monoxide (NO).
#' For the `microdata_transformed`, `min_max`, and `intData`, the pollutant "benzene" was removed due to a high number of missing values.
#' The aggregation of the microdata was done by day.
#' 
#' @docType data
#' @name entrecampos_air_quality
#' @usage data(entrecampos_air_quality)
#' 
#' @references This data was retrieved from the Portuguese Environment Agency database available at \url{https://qualar.apambiente.pt/}.
#' 
#' @format A list with the following components:
#' \describe{
#'   \item{\code{microdata_raw}}{A data frame with `26304` rows and `11` columns. It contains the raw microdata, with individual measurements of each variable for all observations.}
#'   \item{\code{microdata_transformed}}{A data frame with `26304` rows and `10` columns. It contains the microdata, with individual measurements of each variable for all observations. Logarithmic transformations were applied to all variables and interpolation to deal with missing values.}
#'   \item{\code{min_max}}{A data frame with `1096` rows and `17` columns. Each row corresponds to a different observation, and each column gives the minimum and maximum values for each variable. The first column corresponds to the day, the next 8 to the minimum and the last 8 to the maximum.}
#'   \item{\code{intData}}{An \code{\linkS4class{intData}} object, constructed using KDE for estimating the parameters of the latent distributions.}
#' }
#' 
#' @examples
#' data(entrecampos_air_quality)
#' head(entrecampos_air_quality$microdata_raw)
#' head(entrecampos_air_quality$microdata_transformed)
#' head(entrecampos_air_quality$min_max)
#' head(entrecampos_air_quality$intData)
#'
#' @keywords datasets
"entrecampos_air_quality"

#' Credit Card Dataset
#'
#' This dataset contains interval data of credit card expenses, including min-max values, centers and ranges, microdata, and an \code{\linkS4class{intData}} object.
#' It is composed of 5 variables: Food, Social, Travel, Gas, and Clothes. It was aggregated by person-month.
#' 
#' @docType data
#' @name creditcard
#' @usage data(creditcard)
#' 
#' @references This data was retrieved from Billard, L. and Diday, E. (2006). Symbolic Data Analysis: Conceptual Statistics and Data Mining. 
#' John Wiley & Sons. \doi{10.1002/9780470090183}.
#' 
#' @format A list with the following components:
#' \describe{
#'   \item{\code{microdata}}{A data frame with `1000` rows and `7` columns. It contains the microdata, with individual measurements of each variable for all observations.}
#'   \item{\code{min_max}}{A data frame with `36` rows and `10` columns. Each row corresponds to a different observation, and each column gives the minimum and maximum values for each variable.}
#'   \item{\code{centers_ranges}}{A data frame with `36` rows and `10` columns. Each row corresponds to the centers and ranges of the interval data.}
#'   \item{\code{intData}}{An \code{\linkS4class{intData}} object with `36` interval-valued observations and `5` variables, constructed assuming the microdata follow symmetric triangular distributions.}
#' }
#' 
#' @examples
#' data(creditcard)
#' head(creditcard$min_max)
#' head(creditcard$microdata)
#' head(creditcard$intData)
#'
#' @keywords datasets
"creditcard"

#' Spotify Tracks Dataset
#'
#' This dataset contains interval data of Spotify tracks' audio features, including min-max values and trimmed intervals, as well as the microdata.
#' It is composed of 11 audio features: duration, danceability, energy, loudness, speechiness, acousticness, instrumentalness, liveness, valence, tempo, and popularity.
#' The aggregation of the microdata was done by track genre.
#' 
#' @docType data
#' @name spotify_tracks
#' @usage data(spotify_tracks)
#' 
#' @references This data was retrieved from Kaggle (DOI:10.34740/KAGGLE/DSV/4372070; Spotify Tracks Dataset by Maharshi Pandya).
#' 
#' @format A list with the following components:
#' \describe{
#'      \item{\code{microdata}}{A data frame with `81033` rows and `20` columns. It contains the microdata, with individual measurements of each variable for all observations.}
#'      \item{\code{microdata_transformed}}{A data frame with `81033` rows and `20` columns. It contains the transformed microdata, with individual measurements of each variable for all observations. Logarithmic transformations were applied to "loudness" and "tempo". "duration_ms" in milliseconds was converted to "duration" in minutes. "popularity" was scaled to the range `[0,1]`.}
#'      \item{\code{intData_minmax}}{An \code{\linkS4class{intData}} object with `111` interval-valued observations and `11` variables, constructed using min-max aggregation based on the transformed microdata.}
#'      \item{\code{intData_trimmed}}{An \code{\linkS4class{intData}} object with `111` interval-valued observations and `11` variables, constructed using trimmed aggregation (`1\%` trimming) based on the transformed microdata.}
#' }
#' 
#' @examples
#' data(spotify_tracks)
#' head(spotify_tracks$intData_minmax)
#' head(spotify_tracks$intData_trimmed)
#' head(spotify_tracks$microdata)
#' head(spotify_tracks$microdata_transformed)
#' 
#' @keywords datasets
"spotify_tracks"

#' Cars Dataset
#' 
#' This dataset contains interval data of car specifications, including min-max values.
#' It is composed of 5 variables: Engine Capacity, Top Speed, Acceleration, Price and Class.
#' The aggregation of the microdata was done by car model.
#' 
#' @docType data
#' @name intCars
#' @usage data(intCars)
#' 
#' @references This data was retrieved from the \code{MAINT.Data} package, available at \url{https://cran.r-project.org/package=MAINT.Data}.
#' 
#' @format A list with the following components:
#' \describe{
#'  \item{\code{microdata}}{A data frame with `27` rows and `9` columns. It contains the lower and upper bounds for each variable.}
#'  \item{\code{intData}}{An \code{\linkS4class{intData}} object with `27` interval-valued observations and `4` variables. The variable "Price" was log-transformed into "lnPrice". The microdata are not available, thus the default parameters of the latent distributions were used assuming a uniform distribution.}
#' }
#' 
#' @examples
#' data(intCars)
#' head(intCars$microdata)
#' head(intCars$intData)
#' 
#' @keywords datasets
"intCars"