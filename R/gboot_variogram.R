#' @name gboot_variogram
#' @aliases gboot_variogram
#'
#' @import geoR
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @importFrom utils capture.output
#'
#'
#' @title Variogram bootstrap
#'
#' @description Perform a boostrap based on error from the fitted model of the variogram.
#' @usage gboot_variogram(data,var,model,B)
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
#' @details The algorithm for the bootstrap variogram is the same presented for
#' Davison and Hinkley (1997) for the non linear regression. We can write the
#' variogram as \eqn{\hat \gamma(h) = \gamma_{mod}(h)+\epsilon}, where \eqn{\gamma_{mod}(h)}
#' is the fitted model. The steps of the algorithm are:
#' @details
#' \enumerate{
#' \item Set \eqn{h^*=h};
#' \item Sample with replacement \eqn{\epsilon^*} from \eqn{\epsilon - \bar \epsilon};
#' \item The new variogram will be \eqn{\gamma^*(h^*) = \gamma_{mod}(h)+\epsilon^*};
#' \item Calculate and save the statistics of interest;
#' \item Return to step 2 and repeat the process at least 1000 times.
#' }
#'
#'
#'
#' @references DAVISON, A.C.; HINKLEY, D. V. Bootstrap Methods and their Application. [s.l.] Cambridge University Press, 1997. p. 582
#'
#' @keywords Spatial Bootstrap Variogram
#' @examples
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
#' boot<- gboot_variogram(data,var,mod,B=10)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#'
#' \donttest{
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
#' boot<- gboot_variogram(data,var,mod,B=10)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#' }
#' @export




# gboot_variogram ----------------------------------------------------------------

gboot_variogram<-function(data,var,model,B=1000){

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

  #Testing model

  model_name<-substr(model$cov.model,1,3)

  models<-c("sph","exp","gau")

  if(model_name %in% models){
  }else{
    stop("This functions only accepts the following models: gaussian, spherical and exponential")
  }

  #Auxiliary functions

  quiet<-function(x){
    invisible(capture.output(x))}

  #Predict
  if(model_name=="gau" ){
    gboot_predict<-function(u,c0,c1,a){
      ifelse(u==0,0,c0+c1*(1-exp(-3*(u/a)^2)))
    }
  }
  if(model_name=="exp"){
    gboot_predict<-function(u,c0,c1,a){
      ifelse(u==0,0,c0+c1*(1-exp(-3*u/a)))
    }
  }
  if(model_name=="sph"){
    gboot_predict<-function(u,c0,c1,a){
      ifelse(u==0,0,
             ifelse(u>a,c0+c1,c0+c1*(1.5*u/a-.5*(u/a)^3)))
    }
  }

  #Settings
  max_dist<-var$max.dist

  x<-var$u

  c0<-model$nugget

  c1<-model$cov.pars[1]

  a<-model$cov.pars[2]

  y<-var$v

  pars_or<-data.frame(C0=c0,
                      Sill=c0+c1,
                      C1=c1,
                      a=a,
                      PR=model$practicalRange)

  bin<-length(var$u)

  var_df<-matrix(0,
                 nrow=B,
                 ncol=bin)

  pars<-data.frame(C0=rep(0,B),
                   C1=rep(0,B),
                   a=rep(0,B),
                   Sill=rep(0,B),
                   `Pratical Range`=rep(0,B))

  error_or<-y-gboot_predict(x,c0,c1,a)

  var_new<-var

  #Bootstrap
  for(i in 1:B ){
    error_new<-sample((error_or-mean(error_or)),length(error_or),replace = T)

    var_new$v<-gboot_predict(x,c0,c1,a)+error_new

    var_df[i,]<-var_new$v

    quiet(mod_new<-variofit(var_new,ini.cov.pars=c(c1,a),
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







