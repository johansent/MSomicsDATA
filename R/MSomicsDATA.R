#' MSomicsDATA
#'
#' The MSomicsDATA package contains four sets of example mass spectrometry data: proteomics (both peptide- and protein-level data), metabolomics, and lipidomics.
#'
#'
#' @section Data types:
#'
#' This package contains accurate mass and time (AMT) tag mass spectrometry (MS) data for an experiment investigating the human Calu-3 cell response to an infectious clone of Middle Eastern Respiratory Syndrome coronavirus. Data at the peptide (LS-MS/MS), metabolite (GC-MS), and lipid (LC-MS/MS) level are available. Peptide-level data has also been quantified to the protein level using Bayesian proteoform modeling (Webb-Robertson, 2014). Each dataset is from the same experiment, where samples were either infected ("Infection" group) or mock infected ("Mock" group). All samples were collected at the same time point after infection (either the true infection or the mock infection).
#'
#'
#' @section Data formats:
#'
#' Each type of data is provided in two formats. The first format is an S3 object class used by the R pacakge \code{MSomicsQC}. Data available as S3 objects 'pepData', 'proData', 'metabData', and 'lipidData' are created by \code{\link[MSomicsQC]{as.pepData}}, \code{\link[MSomicsQC]{as.proData}}, \code{\link[MSomicsQC]{as.metabData}}, or \code{\link[MSomicsQC]{as.lipidData}}, respectively. The second format corresponds to the individual components of the S3 object classes, \code{e_data} (required), \code{f_data} (required), and \code{e_meta} (optional). See \code{MSomicsQC} for more details.
#'
#'
#' @section Support:
#'
#' Data generation was supported by the National Institutes of Health (NIH)/National Institute of Allergy and Infectious Diseases (NIAID) through grant U19 A1-106772 and computational work was supported by the NIH/National Cancer Institute through grant U01–1CA184783.  Data were collected and processed in the Environmental Molecular Sciences Laboratory (EMSL). EMSL is a national scientific user facility supported by the Department of Energy. All work was performed at PNNL, which is a multiprogram national laboratory operated by Battelle for the U.S. Department of Energy under contract DE-AC06-76RL01830.
#'
#' @references Webb-Robertson BJ, Matzke MM, Datta S, Payne SH, Kang J, Bramer LM, Nicora CD, Shukla AK, Metz TO, Rodland KD, Smith RD, Tardiff MF, McDermott JE, Pounds JG, Waters KM (2014), \emph{Bayesian proteoform modeling improves protein quantification of global proteomic measurements}. Molecular & Cellular Proteomics. doi: 10.1074/mcp.M113.030932.
#'
#' @docType package
#' @name MSomicsDATA
#' @rdname MSomicsDATA
#'
#' @seealso \code{\link[MSomicsQC]{as.proData}}
#' @seealso \code{\link[MSomicsQC]{as.lipidData}}
#' @seealso \code{\link[MSomicsQC]{as.metabData}}
#' @seealso \code{\link[MSomicsQC]{as.pepData}}
NULL
