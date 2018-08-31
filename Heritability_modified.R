##---Heritability function modified from Hilary's

#' Estimate Broad Sense Heritability of a Trait
#'
#' @param responses A character vector of responses for which to calculate heritability
#' @param random.predictors A character vector of random predictor variables. Genotype should be included. Genet (or other genetic variable) should be listed first to recieve correct heritability estimates.
#' @param fixed.predictors A character vector of fixed covariates to include.
#' @param data The data object containing columns corresponding to responses and predictors.
#' @param boxcox.lambda A single value by which to transform the response variable by (default=1).
#'
#' @return The broad sense heritability estimates of each response. genet variance/total variance (extracted from linear mixed effects model lme4::lmer()).
#' @export
#'
#'
#' @examples
#' heritability(responses = c("Phyllocolpa","Harmandia","Phloem_Feeding","Endophagous"),
#'     random.predictors = c("Genet","Block"),
#'     fixed.predictors = c("survey.event","Observer.initials"),
#'     boxcox.lambda = c(1,.5),
#'     data=Master.CPUE.Insect.Data
#' )
#'
#'

heritability <- function(
  ## Variables
  responses, #character vector of response variables
  random.predictors, #character vector of predictor variables with random effects
  fixed.predictors, #character vector of predictor variables with fixed effects
  data=Master.CPUE.Insect.Data, #data set to apply model to
  boxcox.lambda=1 #boxcox transormation to apply to response (response^boxcox.lambda)
  ){
  require(lme4) #needed for mixed effects modeling
  ## H2 as df
  # H2 <- as.data.frame(matrix(nrow = length(responses),ncol = 1)) #empty df for results
  # rownames(H2) <- responses #name the values according to responses
  # names(H2) <- c("H2") #name the column according to output

  ## H2 as vector
  H2 <- vector(mode = "numeric",length = length(responses)) #assign an empty vector of appropriate length
  names(H2) <- responses #name the values according to the responses
  for(response in responses){ #loop through for each response variable
    ## formula objects

    # if(length(boxcox.lambda)>1){
    #   i <- which(responses==response)
    # } else {i <- 1}

    random.form<- paste(paste("(1|",random.predictors,")",
                              sep=""),collapse = "+") #random effects (1|X)+(1|Y)
    fixed.form <- paste(fixed.predictors,collapse = "+") #fixed effects A+B
    predictors.mixed <- paste(random.form,fixed.form,sep = "+") #all predictors (1|X)+A
    form.mixed <- formula(paste(response,"^",boxcox.lambda,"~",
                                predictors.mixed,sep = "")) #full formula
    ## mixed effects model
    m <- lmer(form.mixed,data = data) #mixed effects model using formla
    ## calculate heritability
    VarianceComponents <- as.data.frame(VarCorr(m),comp = "Variance")#extract the variance components from the lmer model
    H2[response] <- VarianceComponents[1,4] / sum(VarianceComponents$vcov) #genet variance (or covariance) divided by total variance (or covariacne)

    ## H2 as df
    # H2[response,"H2"] <- VarianceComponents[1,4] / sum(VarianceComponents$vcov)#calculate the broad-sense heritability (genotypic variance / phenotypic variance)
  }
  return(H2)
}
