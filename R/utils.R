#' Create_Column
#'
#' Worker Function
#' @param knot is a data frame of the data locations
#' @param resolution is a spam object of values of the same length as there are rows in the dataframe
#' @param locations are the initial knots
#' @param max_dist are the initial knots
#' @keywords multiresolution
#' @export
#' @examples
#' create_column()

create_column<-function(knot,resolution,locations,max_dist,nu)
{
  temp=rdist(locations,knot)
  K_kernel_temp=(1-(temp/max_dist[resolution])^2)^nu
  K_kernel_temp[is.nan(K_kernel_temp)]<-0
  K_kernel_temp[which(temp>max_dist[resolution])]<-0
  K_kernel_temp=as.spam(K_kernel_temp)
  return(K_kernel_temp)
}
#' create_knots_2d
#'
#' Worker Function
#' @param knots is a data frame of the data locations
#' @param resolutions is a spam object of values of the same length as there are rows in the dataframe
#' @keywords multiresolution
#' @export
#' @examples
#' create_knots_2d()
create_knots_2d<-function(knots,resolutions,p,resolution_dist)
{
  repknot=knots[rep(1:nrow(knots),each=p),]
  temp=.5*resolution_dist[c(resolutions+1)]
  tempx=rep(temp,each=p)*rep(c(-1,-1,1,1),nrow(knots))
  tempy=rep(temp,each=p)*rep(c(-1,1,-1,1),nrow(knots))
  return(cbind(repknot+cbind(tempx,tempy),rep(resolutions+1,each=p)))
}

#' create_knots_1d
#'
#' Worker Function
#' @param knots is a data frame of the data locations
#' @param resolutions is a spam object of values of the same length as there are rows in the dataframe
#' @keywords multiresolution
#' @export
#' @examples
#' create_knots_2d()
create_knots_1d<-function(knots,resolutions,p,resolution_dist)
{
  repknot=knots[rep(1:length(knots),each=p)]
  temp=.5*resolution_dist[c(resolutions+1)]
  tempx=rep(temp,each=p)*rep(c(-1,1),length(knots))
  return(cbind(repknot+tempx,rep(resolutions+1,each=p)))
}

#' calculate_prior_fixed_pi
#'
#' Worker Function
#' @param params is a list with first element a data frame of the knots with resolutions, second element prior inclusion prob, and third element dimension
#' @keywords multiresolution
#' @export
#' @examples
#' calculate_prior_fixed_pi()
calculate_prior<-function(params,pi_method)
{
  if(pi_method=="beta binomial pi"){return(calculate_prior_beta_pi(params))}
  if(pi_method=="fixed pi"){return(calculate_prior_fixed_pi(params))}
  if(pi_method=="giantdata"){return(calculate_prior_tiny_pi(params))}
}


#' calculate_prior_fixed_pi
#'
#' Worker Function
#' @param params is a list with first element a data frame of the knots with resolutions, second element prior inclusion prob, and third element dimension
#' @keywords multiresolution
#' @export
#' @examples
#' calculate_prior_fixed_pi()
calculate_prior_fixed_pi<-function(params)
{
  knot_resolutions=params[[1]]
  pii=params[[2]]
  d=params[[3]]
  #number of knots at different resolutions
  p_r=as.numeric(table(knot_resolutions))
  p_r=c(p_r,0)
  log_prior_prob=0
  for(r in 2:length(p_r))
  {
    log_prior_prob=log_prior_prob+dbinom(p_r[r],p_r[r-1]*2^d,pii,log=T)
  }
  return(log_prior_prob)
}
#' calculate_prior_tiny_pi
#'
#' Worker Function
#' @param params is a list with first element a data frame of the knots with resolutions, second element prior inclusion prob, and third element dimension
#' @keywords multiresolution
#' @export
#' @examples
#' calculate_prior_tiny_pi()
calculate_prior_tiny_pi<-function(params)
{
  knot_resolutions=params[[1]]
  d=params[[3]]
  a_pi=params[[2]]
  #number of knots at different resolutions
  p_r=as.numeric(table(knot_resolutions))
  p_r=c(p_r,0)
  log_prior_prob=0
  for(r in 2:length(p_r))
  {
    log_prior_prob=log_prior_prob+lchoose(p_r[r-1]*2^d,p_r[r])-p_r[r]*a_pi
  }
  return(log_prior_prob)
}

#' calculate_prior_beta_pi
#'
#' Worker Function
#' @param params is a list with first element a data frame of the knots with resolutions, second and third a and b parameters of the beta prior, and fourth element dimension
#' @keywords multiresolution
#' @export
#' @examples
#' calculate_prior_beta_pi()
calculate_prior_beta_pi<-function(params)
{
  knot_resolutions=params[[1]]
  a=params[[2]]
  b=params[[3]]
  d=params[[4]]
  #number of knots at different resolutions
  p_r=as.numeric(table(knot_resolutions))
  p_r=c(p_r,0)
  a_new=a+sum(p_r[-1])
  b_new=b
  for(r in 2:length(p_r))
  {
    b_new=b_new+p_r[r-1]*2^d-p_r[r]
  }
  log_prior_prob=lbeta(a_new,b_new)-lbeta(a,b)
  return(log_prior_prob)
}