#' @name gboot_cloud
#' @aliases gboot_cloud
#'
#' @import geoR
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @importFrom utils capture.output
#'
#'
#' @title Bootstrap of the variogram cloud
#'
#' @description Performs a boostrap based on the variogram cloud
#' @usage gboot_cloud(data,var,model,B)
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
#' @details The variogram cloud is computed by the function \code{\link[geoR]{variog}}.
#' It provides all the possible pairs that will generate the classical variogram.
#' The algorithm performs a classical bootstrap in each lag of the variogram.
#' The steps are:
#'
#' @details
#' \enumerate{
#' \item Calculate the variogram cloud;
#' \item Obtain the number of lags (See details in \code{\link[geoR]{variog}}: defining the bins);
#' \item Sample with replacement in each lag;
#' \item Create a new variogram using the average of all pairs in each lag;
#' \item Calculate and save the statistics of interest;
#' \item Return to step 3 and repeat the process at least 1000 times.
#' }
#'
#' @keywords Spatial Bootstrap Variogram Cloud
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
#' boot<- gboot_cloud(data,var,mod,B=10)
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
#' boot<- gboot_cloud(data,var,mod,B=10)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#'}
#' @export





# gboot_cloud ----------------------------------------------------------------

gboot_cloud<-function(data,var,model,B=1000){

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

  sample_geo_boot<-function(x){
    sample(x,length(x),replace=T)
  }

  #Settings
  model_name<-substr(model$cov.model,1,3)

  max_dist<-var$max.dist

  x<-var$u

  c0<-model$nugget

  c1<-model$cov.pars[1]

  a<-model$cov.pars[2]

  y<-var$v

  pars_or<-data.frame(C0=c0,
                      Sill=c0+c1,
                      C1=c1,a=a,
                      PR=model$practicalRange)

  quiet(var<-variog(data,
                    max.dist=max_dist,
                    bin.cloud=T))

  bin<-length(var$u)

  var_df<-matrix(0,
                 nrow=B,
                 ncol=bin)

  pars<-data.frame(C0=rep(0,B),
                   C1=rep(0,B),
                   a=rep(0,B),
                   Sill=rep(0,B),
                   `Pratical Range`=rep(0,B))

  quiet(var_new<-variog(data,
                        max.dist = max_dist))


  #Bootstrap
  for(i in 1:B ){
    var_df[i,]<-sapply(lapply(var$bin.cloud,sample_geo_boot),mean)

    var_new$v<-var_df[i,]

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







