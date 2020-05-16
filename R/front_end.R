multires.r1_create<-function(locs,spatial_dimension,n_x=10,buffer=0)
{
  if(spatial_dimension==1)
  {
    knots_r1=seq(floor(min(locs)-buffer),ceiling(max(locs)+buffer),length=(n_x))
    return(knots_r1)
  }
  else if(spatial_dimension==2)
  {
    knots_r1_x=seq(floor(min(locs[,1])-buffer),ceiling(max(locs[,1]+buffer)),length=(n_x))
    domain_width_x=knots_r1_x[length(knots_r1_x)]-knots_r1_x[1]
    partition_width_x=domain_width_x/(length(knots_r1_x)-1)
    knots_r1_y=seq(floor(min(locs[,2])-buffer),
                ceiling(max(locs[,2])+buffer),by=partition_width_x)
    return(expand.grid(knots_r1_x,knots_r1_y))
  }
  else{
    warning("spatial dimension must be 1 or 2")
  }
}