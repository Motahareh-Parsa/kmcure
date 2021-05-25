#' Data of Heart Failure Patients
#'
#' A dataset containing the survival time, event/censoring status, and 11 explanatory variables of heart failure patients who were admitted to the Institute of Cardiology and Allied hospital Faisalabad-Pakistan during April-December 2015.
#'
#' This dataset is adapted from [an open access article (DOI: 10.1371/journal.pone.0181001)](https://doi.org/10.1371/journal.pone.0181001) with Creative Commons Attribution License [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/) which was published in the PLOS ONE journal in 2017.
#'
#' Dataset citation: Ahmad T, Munir A, Bhatti SH, Aftab M, Raza MA (2017) Survival analysis of heart failure patients: A case study. PLoS ONE 12(7): e0181001. https://doi.org/10.1371/journal.pone.0181001
#'
#' @format A data frame with 299 rows and 13 variables:
#' \describe{
#'   \item{TIME}{Time to death duo to Cardiovascular Heart Disease (CHD) by day}
#'   \item{Event}{Event status which is coded as 1 for event and 0 for censoring}
#'   \item{Gender}{Gender which is coded as 0 for women and 1 for men}
#'   \item{Smoking}{Smoking status}
#'   \item{Diabetes}{Diabetes status}
#'   \item{BP}{Blood Pressure (BP) status}
#'   \item{Anaemia}{Anaemia status}
#'   \item{Age}{Age by years}
#'   \item{Ejection.Fraction}{Ejection Fraction (EF)}
#'   \item{Sodium}{Sodium levels}
#'   \item{Creatinine}{Creatinine}
#'   \item{Pletelets}{Pletelets}
#'   \item{CPK}{Creatinine Phosphokinase}
#' }
#' @source \url{https://doi.org/10.1371/journal.pone.0181001}
"hfp"
