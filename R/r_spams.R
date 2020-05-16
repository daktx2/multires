#################
.verif_enum <- function(arg,ename,msg) {
  defName = paste(".__E___", ename, sep = "")
  l = eval.parent(parse(text = sprintf("l <- %s",defName)))
  if (! (arg %in% names(l)))
    stop("ERROR : bad enum value ",msg,"\n")
}
###########  linalg ##############

multires.conjGrad <- function(A,b,x0 = NULL,tol = 1e-10,itermax = NULL) {
  n = ncol(A)
  if(is.null(x0)) {
    x = as.vector(matrix(c(0),ncol = n))
  } else {
    x = c(x0)
  }
  if(is.null(itermax)) {
    itermax = n
  }
  conjugateGradient(A,b,x,tol,itermax)
  return(x)
}

multires.fistaTree <- function(Y,X,W0,tree,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
              fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
              a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
              compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
              resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
              log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
              inner_weights=c(0.),size_group=1,sqrt_step=TRUE,transpose=FALSE,linesearch_mode=0,choler=Matrix(0,2,2),eigenvecinv=matrix(0,2,2),eigenvalinv=matrix(0,2,2)) {

  
  eta_g = tree[['eta_g']]
  groups = tree[['groups']]
  own_variables = tree[['own_variables']]
  N_own_variables = tree[['N_own_variables']]
  m = nrow(W0)
  n = ncol(W0)
#  W = matrix(rep(0,m * n),nrow = m,ncol = n)
  if (length(tree) != 4) {
    stop("fistaTree : tree should be a list of 4 elements")
  }
  W = matrix(c(0),nrow = m,ncol = n)
  optim_info = fistaTree(Y,X,W0,W,eta_g,groups,own_variables,N_own_variables,numThreads ,max_it ,L0,fixed_step,gamma,lambda1,delta,lambda2,lambda3,a,b,c,tol,it0,max_iter_backtracking,compute_gram,lin_admm,admm,intercept,resetflow,regul,loss,verbose,pos,clever,log,ista,subgrad,logName,is_inner_weights,inner_weights,size_group,sqrt_step,transpose,linesearch_mode,choler,eigenvecinv,eigenvalinv)
  if(return_optim_info == TRUE)
    return(list(W,optim_info))
  else
    return (W)
}