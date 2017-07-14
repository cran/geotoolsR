#' @name gboot_cross
#' @aliases gboot_cross
#'
#' @import geoR
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @importFrom utils capture.output
#'
#'
#' @title Cross-validation bootstrap
#'
#' @description Performs a boostrap based on error from the cross-validation
#' @usage gboot_cross(data,var,model,B)
#'
#' @author Diogo Francisco Rossoni \email{dfrossoni@uem.br}
#' @author Vinicius Basseto Felix \email{felix_prot@hotmail.com}
#'
#' @param data object of the class geodata.
#' @param var object of the class variogram.
#' @param model object of the class variomodel.
#' @param B number of the bootstrap that will be performed (default B=1000).
#'
#' @return \bold{variogram_boot} gives the variogram of each bootstrap.
#' @return \bold{variogram_or} gives the original variogram.
#' @return \bold{pars_boot} gives the estimatives of the nugget, sill, contribution, range and practical range for each bootstrap.
#' @return \bold{pars_or} gives the original estimatives of the nugget, sill, contribution, range and practical range.
#' @return Invalid arguments will return an error message.
#'
#' @details We can define the error of prediction by \eqn{\epsilon({s_i})=Z({s_i})-\hat Z({s_i})},
#' where \eqn{\hat Z({s_i})} are obtained from cross-validation. The steps of the algorithm are:
#' @details
#' \enumerate{
#' \item Set \eqn{{s_i}^*={s_i}};
#' \item Obtain \eqn{\hat Z({s_i})} from \eqn{\hat Z({s_i})=\sum\limits_{j \ne i}^{n - 1}{{\lambda _j}Z({s_j})}};
#' \item Calculate \eqn{\epsilon({s_i})=Z({s_i})-\hat Z({s_i})}
#' \item Sample with replacement \eqn{\epsilon^*(s_i)} from \eqn{\epsilon (s_i) - \bar \epsilon (s_i)};
#' \item The new data will be \eqn{Z^*({s_i})=\hat Z({s_i})+ \epsilon^*(s_i)};
#' \item Calculate the new variogram;
#' \item Calculate and save the statistics of interest;
#' \item Return to step 4 and repeat the process at least 1000 times.
#' }
#'
#'
#'
#' @keywords Spatial Bootstrap Cross-validation
#' @examples
#' 
#' \dontrun{
#' # Example 1
#'
#' ## transforming the data.frame in an object of class geodata
#' data<- as.geodata(soilmoisture)
#'
#' points(data) ## data visualization
#'
#' var<- variog(data, max.dist = 140) ## Obtaining the variogram
#' plot(var)
#'
#' ## Fitting the model
#' mod<- variofit(var,ini.cov.pars = c(2,80),nugget = 2,cov.model = "sph")
#' lines(mod, col=2, lwd=2) ##fitted model
#'
#' ## Bootstrap procedure
#'
#' boot<- gboot_cross(data,var,mod,B=10)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#'
#' # Example 2
#'
#' ## transforming the data.frame in an object of class geodata
#' data<- as.geodata(NVDI)
#'
#' points(data) ## data visualization
#'
#' var<- variog(data, max.dist = 18) ## Obtaining the variogram
#' plot(var)
#'
#' ## Fitting the model
#' mod<- variofit(var,ini.cov.pars = c(0.003,6),nugget = 0.003,cov.model = "gaus")
#' lines(mod, col=2, lwd=2) ##fitted model
#'
#' ## Bootstrap procedure
#'
#' boot<- gboot_cross(data,var,mod,B=10)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#' }
#'
#' @export






# gboot_cross ----------------------------------------------------------------

gboot_cross<-function(data,var,model,B=1000){

  Distance=Semivariance=NULL
  
  #Testing
  if(is.geodata(data) == T){
  }else{
    stop("Object data is not of the class geodata")
  }
  if(isTRUE(class(var) == "variogram")){
  }else{
    stop("Object var is not of the class variogram")
  }
  if(isTRUE(class(model)[1] == "variomodel") & isTRUE(class(model)[2] == "variofit")){
  }else{
    stop("Object model is not of the class variomodel/variofit")
  }
  if(B >0 ){
  }else{
    stop("Object B must be positive")
  }

  #Auxiliary functions
  quiet<-function(x){
    invisible(capture.output(x))}

  #Settings
  model_name<-substr(model$cov.model,1,3)

  max_dist<-var$max.dist

  x<-var$u

  c0<-model$nugget

  c1<-model$cov.pars[1]

  a<-model$cov.pars[2]

  pars_or<-data.frame(C0=c0,
                      Sill=c0+c1,
                      C1=c1,
                      a=a,
                      PR=model$practicalRange)

  quiet(v_cruzada<- xvalid(data,
                           model=model))

  bin<-length(var$u)

  var_df<-matrix(0,nrow=B,ncol=bin)

  pars<-data.frame(C0=rep(0,B),
                   C1=rep(0,B),
                   a=rep(0,B),
                   Sill=rep(0,B),
                   `Pratical Range`=rep(0,B))

  #Bootstrap
  for(i in 1:B ){
    error_new<-sample((v_cruzada$error-mean(v_cruzada$error)),
                     length(v_cruzada$error),replace = T)

    df_new<-as.geodata(data.frame(data$coords[,1],
                                  data$coords[,2],
                                  v_cruzada$predicted+error_new))

    quiet(var_new<-variog(df_new,
                          max.dist=max_dist))

    var_df[i,]<-var_new$v

    quiet(mod_new<-variofit(var_new,
                            ini.cov.pars=c(c1,a),
                            nugget=c0,
                            cov.model=model_name))

    pars[i,]<-c(as.numeric(summary(mod_new)$estimated.pars[1]),
                sum(as.numeric(c(summary(mod_new)$estimated.pars)[1:2])),
                as.numeric(c(summary(mod_new)$estimated.pars[2:3])),
                mod_new$practicalRange)
    }

  var_df<-as.data.frame(var_df)

  names(var_df)<-paste("Class",letters[1:bin])

  var_df<-gather(var_df,Distance,Semivariance)

  var_df$B<-rep(1:B,bin)

  var_aux<-data.frame(Distance=paste("Class",letters[1:bin]),Semivariance=var$v)

  var_aux$Length<-var$u

  names(pars)<-c("Nugget","Sill","Contribution","Range","Practical Range")

  names(pars_or)<-c("Nugget","Sill","Contribution","Range","Practical Range")

  return(list(variogram_boot=var_df,
              variogram_or=var_aux,
              pars_boot=pars,
              pars_or=pars_or))
}







