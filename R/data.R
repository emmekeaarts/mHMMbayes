#' Nonverbal communication of patients and therapist
#'
#' A dataset containing the nonverbal communication of 10 patient-therapist couples,
#' recorded for 15 minutes at a frequency of 1 observation per second (= 900 observations per couple).
#'
#' @format A matrix with 10 * 900 rows and 5 variables:
#' \describe{
#' \item{id}{id variable of patient - therapist couple to distinguish which
#' observation belongs to which couple}
#' \item{p_verbalizing}{verbalizing behavior of the patient, consisting of 1 =
#' not verbalizing, 2 = verbalizing, 3 = back channeling}
#' \item{p_looking}{looking behavior of the patient, consisting of 1 = not
#' looking at therapist, 2 = looking at therapist}
#' \item{t_verbalizing}{verbalizing behavior of the therapist, consisting of 1
#' = not verbalizing, 2 = verbalizing, 3 = back channeling}
#' \item{t_looking}{looking behavior of the therapist, consisting of 1 = not
#' looking at patient, 2 = looking at patient}
#' }
"nonverbal"

#' Predictors of nonverbal communication
#'
#' A dataset containing predictors of nonverbal communication of 10 patient-therapist couples.
#'
#' @format A matrix with 10 rows and 3 variables:
#' \describe{
#' \item{diagnosis}{Diagnosis of the patient, consisting of 0 = depression, 1 = anxiety}
#' \item{std_CDI_change}{Change in measure for depression (CDI) before and after therapy, standardized scale}
#' \item{std_SCA_change}{Change in measure for anxiety (SCARED) before and after therapy, standardized scale}
#' }
"nonverbal_cov"
