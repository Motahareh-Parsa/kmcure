#' Data of Heart Failure Patients
#'
#' A dataset containing the survival time, event/censoring status, and 11 novel explanatory variables of heart failure patients who were admitted to the Institute of Cardiology and Allied hospital Faisalabad-Pakistan during April-December 2015.
#'
#' This dataset is adapted from [an open access article (DOI: 10.1371/journal.pone.0181001)](https://doi.org/10.1371/journal.pone.0181001) with Creative Commons Attribution License [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/) which was published in the PLOS ONE journal in 2017.
#'
#' Dataset citation: Ahmad T, Munir A, Bhatti SH, Aftab M, Raza MA (2017) Survival analysis of heart failure patients: A case study. PLoS ONE 12(7): e0181001. https://doi.org/10.1371/journal.pone.0181001
#'
#' @format A data frame with 299 rows and 20 variables:
#' \describe{
#'   \item{Time}{Time to death duo to Cardiovascular Heart Disease (CHD) per day}
#'   \item{Event}{Event status which is coded as 1 for event and 0 for censoring}
#'   \item{Gender}{Gender variable which is coded 1 for men and 0 for women}
#'   \item{Smoking}{Smoking status}
#'   \item{Diabetes}{Diabetes status}
#'   \item{BP}{Blood Pressure (BP) status}
#'   \item{Anaemia}{Anaemia status}
#'   \item{Age}{Age by years}
#'   \item{EF}{Ejection Fraction (EF)}
#'   \item{Sodium}{serum Sodium}
#'   \item{Creatinine}{serum Creatinine}
#'   \item{Platelets}{Platelets}
#'   \item{CPK}{Creatinine Phosphokinase}
#'   \item{EFlow}{Ejection Fraction less than or equal 30}
#'   \item{EFmed}{Ejection Fraction greater than 30 but less than or equal 45}
#'   \item{EFhigh}{Ejection Fraction greater than 45}
#'   \item{PLTlow}{Platelets less than or equal the first quartile}
#'   \item{PLTmed}{Platelets greater than the first quartile but less than or equal the third quartile}
#'   \item{PLThigh}{Platelets greater than the third quartile}
#'   \item{logCPK}{Natural logarithms of Creatinine Phosphokinase}
#' }
#' @source \url{https://doi.org/10.1371/journal.pone.0181001}
"hfp"
