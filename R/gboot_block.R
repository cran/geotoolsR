#' @name gboot_block
#' @aliases gboot_block
#'
#' @import geoR
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @importFrom utils capture.output
#'
#' @title Block bootstrap
#'
#' @description Performs a bootstrap based on subdivision of data in blocks
#' @usage gboot_block(data,var,model,B,L1,L2)
#' @author Diogo Francisco Rossoni \email{dfrossoni@uem.br}
#' @author Vinicius Basseto Felix \email{felix_prot@hotmail.com}
#'
#' @param data object of the class geodata.
#' @param var object of the class variogram.
#' @param model object of the class variomodel.
#' @param B number of the bootstrap that will be performed (default B=1000).
#' @param L1 number of cuts in the vertical (L1xL2 blocks).
#' @param L2 number of cuts in the horizontal (L1xL2 blocks).
#'
#' @return \bold{variogram_boot} gives the variogram of each bootstrap.
#' @return \bold{variogram_or} gives the original variogram.
#' @return \bold{pars_boot} gives the estimatives of the nugget, sill, contribution, range and practical range for each bootstrap.
#' @return \bold{pars_or} gives the original estimatives of the nugget, sill, contribution, range and practical range.
#' @return Invalid arguments will return an error message.
#'
#' @details The algorithm for the block bootstrap is an adaptation of the time series bootstrap.
#' Consider that your data presents the second order stationarity, so, we can subdivide them into small blocks.
#' The steps of the algorithm are:
#' @details
#' \enumerate{
#' \item Subdivide the data into L1xL2 blocks;
#' \item Realocate each block with probability\eqn{\frac{1}{L1L2}} ;
#' \item Calculate the new variogram from the new data;
#' \item Calculate and save the statistics of interest;
#' \item Return to step 2 and repeat the process at least 1000 times.
#' }
#'
#'
#'
#' @references DAVISON, A.C.; HINKLEY, D. V. Bootstrap Methods and their Application. [s.l.] Cambridge University Press, 1997. p. 582
#'
#' @keywords Spatial Bootstrap Block
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
#' boot<- gboot_block(data,var,mod,B=10, L1=2, L2=2)
#' ## For better Confidence Interval, try B=1000
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
#' boot<- boot<- gboot_block(data,var,mod,B=10, L1=2, L2=2)
#' ## For better Confidence interval, try B=1000
#'
#' gboot_CI(boot,digits = 4) ## Bootstrap Confidence Interval
#'
#' gboot_plot(boot) ## Bootstrap Variogram plot
#' }
#'
#' @export




# gboot_block ----------------------------------------------------------------

gboot_block<-function(data,var,model,B=1000,L1=2,L2=2){

  lat=long=block=value=Distance=Semivariance=NULL
  
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

  bin<-length(var$u)

  var_df<-matrix(0,
                 nrow=B,
                 ncol=bin)

  pars<-data.frame(C0=rep(0,B),
                   C1=rep(0,B),
                   a=rep(0,B),
                   Sill=rep(0,B),
                   `Pratical Range`=rep(0,B))

  x<-var$u

  c0<-model$nugget

  c1<-model$cov.pars[1]

  a<-model$cov.pars[2]

  pars_or<-data.frame(C0=c0,
                      Sill=c0+c1,
                      C1=c1,
                      a=a,
                      PR=model$practicalRange)

  df_new <- data.frame(lat=data[[1]][,1],
                       long=data[[1]][,2],
                       value=data[[2]])

  df_new <- within(df_new, {
    x = cut(lat , L1, labels = FALSE)
    y = cut(long, L2, labels = FALSE)
  })

  df_new$block<-paste0(df_new$y,df_new$x)

  # p1<-df_new %>%
  #   ggplot(aes(lat,long,col=block))+
  #   geom_point()+
  #   labs(x="Latitude",y="Longitude",col="Block:")+
  #   coord_equal()+
  #   theme_bw()

  df_new %>%
    group_by(x) %>%
    mutate(lat=max(lat)) %>%
    group_by(y) %>%
    mutate(long=max(long)) ->teste

  p1<-df_new %>%
    ggplot(aes(lat,long,col=block))+
    geom_point()+
    labs(x="Latitude",y="Longitude",col="Block:")+
    coord_equal()+
    theme_minimal()+
    geom_vline(data=teste,aes(xintercept=lat),linetype="dashed")+
    geom_hline(data=teste,aes(yintercept=long),linetype="dashed")+
    ylim(min(df_new$long),max(df_new$long))+
    xlim(min(df_new$lat),max(df_new$lat))

  sample<-data

  #Bootstrap
  for(i in 1:B ){
    sample$data<- df_new %>%
      group_by(block) %>%
      nest() %>%
      sample_frac() %>%
      unnest() %>%
      select(value)

    quiet(var_new<-variog(sample,
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

  print(p1)

  return(list(variogram_boot=var_df,
              variogram_or=var_aux,
              pars_boot=pars,
              pars_or=pars_or))
}







