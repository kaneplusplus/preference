#' @title Design and Analysis of Two-stage Randomized Clinical Trials
#' @name preference
#' @docType package
#' @description The \pkg{preference} package is used for the design and analysis of
#' two-stage randomized trials with a continuous outcome measure. In this study,
#' patients are first randomized to either a random or choice arm. Patients 
#' initially randomized to the choice arm are allowed to select their preferred
#' treatment from the available treatment options; patients initially randomized
#' to the random arm undergo a second randomization procedure to one of the 
#' available treatment options. The design has also been extended to include
#' important stratification variables; the functions provided in this package
#' can accomodate both the unstratified and stratified designs.
#' 
#' In this study, there are three effects that may be of interest. The treatment
#' effect captures the difference in outcome between patients randomized to 
#' treatment A and treatment B (similar to a traditional RCT). The selection 
#' effect captures the difference in outcome between patients that prefer 
#' treatment A and patients that prefer treatment B, regardless of the 
#' treatment that is actually received. Finally, the preference effect compares
#' the outcomes of patients who receive their preferred treatment (either 
#' treatment A or treatment B) and patients who do not receive their preferred
#' treatment. 
#' 
#' To aid in the design of these two-stage randomized studies, sample size 
#' functions are provided to determine the necessary sample size to 
#' detect a particular selection, preference, and/or treatment effect. If the 
#' sample size is fixed prior to the start of the study, functions are provided
#' to calculate the study power to detect each effect. Finally, the 
#' \code{theta_optim} function can be used to determine the optimal proportion
#' of patients randomized to the choice arm in the initial randomization. 
#' 
#' To analyze the data from the two-stage randomized trial, two analysis
#' functions are provided. The function \code{analysis_raw} computes the test
#' statistic and p-value for each effect given provided raw study data. The
#' function \code{analysis_summary} uses provided summary data (mean, variance,
#' and sample size) of each study group to compute the test statistic and 
#' p-value of each effect. 
#' 
#' Sample Size Function calls
#' \itemize{
#' \item{n_sel: }{required sample size to detect a given selection effect}
#' \item{n_pref: }{required sample size to detect a given preference effect}
#' \item{n_trt: }{required sample size to detect a given treatment effect}
#' \item{n_overall: }{required sample size to detect a given set of selection, 
#' preference, and treatment effects}
#' }
#' 
#' Power Function Calls
#' \itemize{
#' \item{sel_pwr: }{study power to detect a given selection effect}
#' \item{pref_pwr: }{study power to detect a given preference effect}
#' \item{trt_pwr: }{study power to detect a given treatment effect}
#' \item{pwr_overall: }{study power to detect a given set of selection, 
#' preference, and treatment effects}
#' }
#' 
#' Analysis Function Calls
#' \itemize{
#' \item{analysis_raw: }{computes test statistic and p-value for observed
#' selection, preference, and treatment effects using provided raw data}
#' \item{analysis_summary: }{computes test statistic and p-value for observed
#' selection, preference, and treatment effects using provided summary data 
#' (mean, variance, sample size)}
#' }
#' 
#' Other Function Calls
#' \itemize{
#' \item{trt_effect: }{computes the treatment effect that can be detected given a
#' specified sample size and power}
#' \item{theta_optim: }{computes the optimal proportion randomized to choice arm
#' (defined for unstratified design only)}
#' \item{calc_effects: }{computes the treatment, selection, and preference effect
#' sizes provided the study means in each treatment arm}
#' }
#' 
#' Data Sets
#' \itemize{
#' \item{imap: }{summary SF36 outcome data for the two-stage randomized IMAP 
#' study}
#' \item{imap_strat: }{summary SF36 outcome data for the two-stage randomized 
#' IMAP study stratified by high vs. low STAI score}
#' }
#' @references Rucker G (1989). "A two-stage trial design for testing treatment,
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @references McCaffery et al. (2010) "Psychosocial outcomes of three triage 
#' methods for the management of borderline abnormal cervical smears: an open 
#' randomised trial." \emph{BMJ}, \strong{340}:b4491.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2827716/}{PubMed})
#' @references Walter et. al. (2011). "Optimal allocation of participants for
#' the estimation of selection, preference and treatment effects in the 
#' two-stage randomised trial design." \emph{Stat Med}, 
#' \strong{31}(13):1307-1322.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22362374}{PubMed})
#' @references McCaffery et al. (2011) "Determining the Impact of Informed 
#' Choice: Separating Treatment Effects from the Effects of Choice and Selection
#' in Randomized Trials." \emph{Med Decis Making}, \strong{31}(2):229-236.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/21041538}{PubMed})
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
NULL


#' Summary data for the IMAP study
#'
#' @name imap
#' @docType data
#' @description Data set is from the Improving Management of Abnormal Pap Smears
#' study, which used a doubly randomized preference trial design to evaluate 
#' psychosocial outcomes in women found to have atypical cells in a Pap Smear. 
#' Two systems for managing the atypical cells were tested (repeated Pap smears 
#' or HCV triage) and a doubly randomized design was used to evaluate the role 
#' of patient preference. The data set provides mean, standard deviation and 
#' sample sizes of the SF36 outcome for each treatment in both the choice and 
#' random arms.
#' @references McCaffery et al. (2010) "Psychosocial outcomes of three triage 
#' methods for the management of borderline abnormal cervical smears: an open 
#' randomised trial." \emph{BMJ}, \strong{340}:b4491.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2827716/}{PubMed})
#' @references McCaffery et al. (2011) "Determining the Impact of Informed 
#' Choice: Separating Treatment Effects from the Effects of Choice and Selection
#' in Randomized Trials." \emph{Med Decis Making}, \strong{31}(2):229-236.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/21041538}{PubMed})
#' @keywords data
NULL


#' Stratified summary data for the IMAP study
#'
#' @name imap_strat
#' @docType data
#' @description Data set is from the Improving Management of Abnormal Pap Smears
#' study, which used a doubly randomized preference trial design to evaluate 
#' psychosocial outcomes in women found to have atypical cells in a Pap Smear. 
#' Two systems for managing the atypical cells were tested (repeated Pap smears
#' or HCV triage) and a doubly randomized design was used to evaluate the role 
#' of patient preference. 
#' 
#' In this data set, patients are stratified according to their baseline score
#' on the six-item abbreviated Spielberger State-Trait Anxiety Inventory (STAI).
#' This assessment is used to assess the level of anxiety experienced by 
#' patients. Stratum 1 includes patients with low STAI scores (<=10 on STAI 
#' averaged across whole study), which indicates low levels of anxiety. 
#' Stratum 2 includes patients with high STAI  scores (>10), indicating higher
#' levels of anxiety. The data set provides mean, standard deviation and 
#' sample sizes of the SF36 outcome for each treatment in both the choice and 
#' random arms for both strata.
#' 
#' @references McCaffery et al. (2010) "Psychosocial outcomes of three triage 
#' methods for the management of borderline abnormal cervical smears: an open 
#' randomised trial." \emph{BMJ}, \strong{340}:b4491.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2827716/}{PubMed})
#' @references McCaffery et al. (2011) "Determining the Impact of Informed 
#' Choice: Separating Treatment Effects from the Effects of Choice and Selection
#' in Randomized Trials." \emph{Med Decis Making}, \strong{31}(2):229-236.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/21041538}{PubMed})
#' @references Marteau TM, Bekker H. (1992) "The development of a six-item 
#' short-form of the state scale of the Spielberger State-Trait Anxiety 
#' Inventory (STAI)." \emph{Br J Clin Psychol}, \strong{31}:301-306.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/1393159#}{PubMed})
#' @keywords data
NULL