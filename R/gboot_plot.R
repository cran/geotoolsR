#' @name gboot_plot
#' @aliases gboot_plot
#' 
#' @import ggplot2
#'
#' @title Bootstrap plot
#'
#' @description A graphic with the original variogram and each of the B bootstrap variograms.
#' @usage gboot_plot(x)
#'
#' @author Diogo Francisco Rossoni \email{dfrossoni@uem.br}
#' @author Vinicius Basseto Felix \email{felix_prot@hotmail.com}
#'
#' @param x object generate by functions \code{\link[geotoolsR]{gboot_block}},
#'  \code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
#'   \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}
#'
#'
#' @details Examples of this function can be found in \code{\link[geotoolsR]{gboot_block}},
#'  \code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
#'   \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}
#' 
#' @export


# gboot_plot --------------------------------------------------------------

gboot_plot<-function(x){
    
    Distance=Semivariance=NULL

    p1<-ggplot(x[[1]],aes(x=Distance,y=Semivariance))+
    geom_point(alpha=.5,col="gray")+
    theme_minimal()+
    labs(x="Distance",y="Semivariance")+
    scale_x_discrete(breaks=x[[2]]$Distance,
                     labels=round(x[[2]]$Length))+
    geom_point(data=x[[2]],aes(x=Distance,y=Semivariance))

  p1

}
