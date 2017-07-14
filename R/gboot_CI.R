#' @name gboot_CI
#' @aliases gboot_CI
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stats quantile
#' @importFrom utils capture.output
#'
#' @title Bootstrap Confidence Interval
#'
#' @description Provides a two-sided bootstrap confidence interval.
#' @usage gboot_CI(x,alpha=0.05,digits=3)
#'
#' @author Diogo Francisco Rossoni \email{dfrossoni@uem.br}
#' @author Vinicius Basseto Felix \email{felix_prot@hotmail.com}
#'
#' @param x object generate by functions \code{\link[geotoolsR]{gboot_block}},
#'  \code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
#'   \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}
#' @param alpha significance level (Default=0.05).
#' @param digits number of decimal places.
#'
#' @details Examples of this function can be found in \code{\link[geotoolsR]{gboot_block}},
#'  \code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
#'   \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}
#'
#' @keywords Bootstrap CI
#' @export

gboot_CI<-function(x,alpha=0.05,digits=3){
  
  Parameter=Value=Bound=NULL

  #Auxiliary functions
  quiet<-function(x){
    invisible(capture.output(x))}

  #Confidence interval
  CI_pars<-data.frame(Parameter=c("Nugget","Sill","Contribution",
                                  "Range","Practical Range"),
                      Lower=numeric(5),
                      Estimate=c(x[[4]][1,1],x[[4]][1,2],
                                 x[[4]][1,3],x[[4]][1,4],x[[4]][1,5]),
                      Upper=numeric(5))

  CI_pars[,2]<-apply(x[[3]],2,function(x)round(quantile(x,alpha/2),digits))

  CI_pars[,4]<-apply(x[[3]],2,function(x)round(quantile(x,1-alpha/2),digits))

  CI_pars[,3]<-round(CI_pars[,3],digits)

  cat("Parameters confidence interval (",
      100*(1-alpha),"%): \n",sep="")
  CI_pars

  #Plot
  CI<-gather(x[[3]],Parameter,Value)%>%
    mutate(Parameter=factor(Parameter,
                            c("Nugget",
                              "Sill",
                              "Contribution",
                              "Range",
                              "Practical Range"))) %>%
    group_by(Parameter)%>%
    summarise(Lower=quantile(Value,probs = alpha/2),
              Upper=quantile(Value,probs = 1-alpha/2))%>%
    gather(Bound,Value,-Parameter)

    quiet(p1<-gather(x[[3]],Parameter,Value) %>%
          mutate(Parameter=factor(Parameter,
                                  c("Nugget",
                                    "Sill",
                                    "Contribution",
                                    "Range",
                                    "Practical Range"))) %>%
    ggplot(aes(Value))+
    geom_histogram(col="black",bins=30)+
    theme_minimal()+
    geom_vline(data=CI,aes(xintercept = Value),
               linetype="dashed",col="red")+
    facet_wrap(~Parameter,scales = "free",ncol=3)+
    labs(x="Value",y="Frequency"))

    #Output

    quiet(print(p1))

    return(CI_pars)

}

