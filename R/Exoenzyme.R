
#' Structure of the S4 class "Exoenzyme"
#'
#' @aliases Exoenzyme
#'
#' @exportClass Exoenzyme
#'
#' @slot id Character for exoenzyme ID
#' @slot name Characer of exoenzyme name
#' @slot D Numeric of the enzymes' diffusion coefficients. Unit: µm^2 per
#' sec
#' @slot lambda Numeric for enzyme's decay constants. Unit: per hr
#' @slot Kcat Numeric for turnover rate (kcat) in MM-kinetics. Unit: 1/s
#' @slot Km Numeric for Km in MM-kinetics. Unit: mM
#' @slot mets Character vector of metabolite IDs for compounds in the reaction.
#' Main substrate should be the first element.
#' @slot stoich Numeric vector with the stoichiometries of `mets` in exoenzyme's
#' reaction.
setClass("Exoenzyme",

         slots = c(
           id = "character",
           name = "character",
           D = "numeric",
           lambda = "numeric",
           Kcat = "numeric",
           Km = "numeric",
           mets = "character",
           stoich = "numeric"

         )
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#  Adding an Exoenzyme          #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Add an exoenzyme to organism and simulation
#'
#' @description Adds a new exoenzyme to the simulation and links it to a specific
#' organism.
#'
#' @param object S4-object of type \link{growthSimulation}
#' @param organism Character indicating the name of the organism, to which the
#' exoenzyme production is linked
#' @param id Character defining the ID used for the specific exoenzyme
#' @param mets Character vector with the compound IDs that participate in reaction
#' catalyzed by the exoenzyme.
#' @param stoich Numeric vector of the same length and order as `mets` specifying
#' the stoichiometries of compounds. First entry (the substrate's coefficient)
#' should be -1.
#' @param production.rate Numeric indicating the production rate of the enzyme
#' by the organism. Unit: nmol enzyme catalytic centers per gDW cells per hr.
#' Default: 0.01
#' @param name Character with an optional name for the enzyme.
#' @param D Diffusion coefficient of the enzyme. Unit µm^2/s. Default: 10
#' @param lambda Numeric indicating the decay rate of the enzyme. Unit: per hour.
#' Default: 0.4 . The enzyme's half life can be calculated by ln(2)/lambda
#' @param Kcat Numeric for enzyme's turnover rate. Unit: 1/s . Default: 10000
#' @param Km Numeric for Menten-Michaelis Km value. Unit: mM . Default: 100
#' @param init.conc Numeric indicating the initial concentration of the enzyme
#' in the growth environment. Unit: nM . Default: 0
#'
#' @details
#' Exoenzymes can differ markedly in their kinetic parameters. The defaults
#' provided here do not represent a typical enzymes, but are chosen based
#' on data from inveratases (EC 3.2.1.26) from Zymomonas mobilis. Kcat (~10000 s^-1)
#' and Km (100 mM) were obtained from https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.26&Suchword=&reference=&UniProtAcc=&organism%5B%5D=Zymomonas+mobilis&show_tm=0.
#'
#' @importFrom methods new
#'
#' @export
add_exoenzyme <- function(object, organism, id, mets, stoich,
                          production.rate = 0.01,
                          name = NULL,
                          D = 10,
                          lambda = 0.4,
                          Kcat = 10000,
                          Km = 100,
                          init.conc = 0) {

  ind <- length(object@environ@exoenzymes) + 1

  # Sanity checks
  if(id %in% names(object@environ@exoenzymes)) {
    stop(paste0("Exoenzmye with ID '",id,"' is already present."))
  }

  if(stoich[1] != -1)
    stop("First enty of 'stoich' represents the stoichiometric coefficient of the substrate and should be -1.")

  if(production.rate < 0)
    stop("Production rate should be postive.")

  if(is.null(name))
    name <- id

  # If mets are not yet part of the environment, add them here
  new_mets <- mets[!(mets %in% object@environ@compounds)]
  if(length(new_mets) > 0)
    object <- add_compounds(object, new_mets,
                            concentrations = rep(0, length(new_mets)))

  # create Exoenzyme object
  exec.obj <- new("Exoenzyme",
                  id = id,
                  name = name,
                  D = D,
                  lambda = lambda,
                  Kcat = Kcat,
                  Km = Km,
                  mets = mets,
                  stoich = stoich)
  object@environ@exoenzymes[[ind]] <- exec.obj
  names(object@environ@exoenzymes)[ind] <- id

  # Add new column to exoenzyme concentration matrix
  object@environ@exoenzymes.conc <- cbind(object@environ@exoenzymes.conc,
                                          rep(init.conc, object@environ@nfields))

  # Add Exoenzyme info to Organism object
  object@models[[organism]]@exoenzymes <- c(object@models[[organism]]@exoenzymes,
                                            id)
  object@models[[organism]]@exoenzymes.prod <- c(object@models[[organism]]@exoenzymes.prod,
                                                 production.rate)

  return(object)
}

