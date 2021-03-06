\encoding{latin1}
\name{spams.fistaGraph}
\alias{spams.fistaGraph}
\title{
spams.fistaGraph
}
\description{
    spams.fistaGraph solves sparse regularized problems.
}
\usage{
spams.fistaGraph(Y,X,W0,graph,return_optim_info = FALSE,numThreads =-1,max_it =1000,L0=1.0,
                 fixed_step=FALSE,gamma=1.5,lambda1=1.0,delta=1.0,lambda2=0.,lambda3=0.,
                 a=1.0,b=0.,c=1.0,tol=0.000001,it0=100,max_iter_backtracking=1000,
                 compute_gram=FALSE,lin_admm=FALSE,admm=FALSE,intercept=FALSE,
                 resetflow=FALSE,regul="",loss="",verbose=FALSE,pos=FALSE,clever=FALSE,
                 log=FALSE,ista=FALSE,subgrad=FALSE,logName="",is_inner_weights=FALSE,
                 inner_weights=c(0.),size_group=1,sqrt_step=TRUE,transpose=FALSE,
                 linesearch_mode=0)
}
\arguments{
\item{Y}{%
double dense m x n matrix
\preformatted{%
}}
\item{X}{%
double dense or sparse m x p matrix   
\preformatted{%
}}
\item{W0}{%
double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
\preformatted{%
initial guess
}}
\item{graph}{%
struct (see documentation of spams.proximalGraph)
\preformatted{%
}}
\item{return_optim_info}{%
\preformatted{%
if true the function will return 2 matrices in a list.
}}
\item{loss}{%
(choice of loss, see above)
\preformatted{%
}}
\item{regul}{%
(choice of regularization, see function spams.proximalFlat)
\preformatted{%
}}
\item{lambda1}{%
(regularization parameter)
\preformatted{%
}}
\item{lambda2}{%
(optional, regularization parameter, 0 by default)
\preformatted{%
}}
\item{lambda3}{%
(optional, regularization parameter, 0 by default)
\preformatted{%
}}
\item{verbose}{%
(optional, verbosity level, false by default)
\preformatted{%
}}
\item{pos}{%
(optional, adds positivity constraints on the
\preformatted{%
coefficients, false by default)
}}
\item{numThreads}{%
(optional, number of threads for exploiting
\preformatted{%
multi-core / multi-cpus. By default, it takes the value -1,
which automatically selects all the available CPUs/cores).
}}
\item{max_it}{%
(optional, maximum number of iterations, 100 by default)
\preformatted{%
}}
\item{it0}{%
(optional, frequency for computing duality gap, every 10 iterations by default)
\preformatted{%
}}
\item{tol}{%
(optional, tolerance for stopping criteration, which is a relative duality gap
\preformatted{%
if it is available, or a relative change of parameters).
}}
\item{gamma}{%
(optional, multiplier for increasing the parameter L in fista, 1.5 by default)
\preformatted{%
}}
\item{L0}{%
(optional, initial parameter L in fista, 0.1 by default, should be small enough)
\preformatted{%
}}
\item{fixed_step}{%
(deactive the line search for L in fista and use L0 instead)
\preformatted{%
}}
\item{compute_gram}{%
(optional, pre-compute X^TX, false by default).
\preformatted{%
}}
\item{intercept}{%
(optional, do not regularize last row of W, false by default).
\preformatted{%
}}
\item{ista}{%
(optional, use ista instead of fista, false by default).
\preformatted{%
}}
\item{subgrad}{%
(optional, if not ista, use subradient descent instead of fista, false by default).
\preformatted{%
}}
\item{a}{%
\preformatted{%
}}
\item{b}{%
(optional, if subgrad, the gradient step is a/(t+b)
\preformatted{%
also similar options as spams.proximalTree

the function also implements the ADMM algorithm via an option admm=true. It is not documented
and you need to look at the source code to use it.
}}
\item{delta}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{c}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{max_iter_backtracking}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{lin_admm}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{admm}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{resetflow}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{clever}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{log}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{logName}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{is_inner_weights}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{inner_weights}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{sqrt_step}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{size_group}{%
undocumented; modify at your own risks!
\preformatted{%
}}
\item{transpose}{%
undocumented; modify at your own risks!
\preformatted{%
}}
}
\details{
\preformatted{%

        X is a design matrix of size m x p
        X=[x^1,...,x^n]', where the x_i's are the rows of X
        Y=[y^1,...,y^n] is a matrix of size m x n
        It implements the algorithms FISTA, ISTA and subgradient descent.
        
        It implements the algorithms FISTA, ISTA and subgradient descent for solving
        
          min_W  loss(W) + lambda1 psi(W)
          
        The function psi are those used by spams.proximalGraph (see documentation)
        for the loss functions, see the documentation of spams.fistaFlat
        
        This function can also handle intercepts (last row of W is not regularized),
        and/or non-negativity constraints on W.
}
}
\value{
\item{W}{%
double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
\preformatted{%
}}
\item{optim}{%
optional, double dense 4 x n matrix.
\preformatted{%
first row: values of the objective functions.
third row: values of the relative duality gap (if available)
fourth row: number of iterations
}}
\item{optim_info}{%
vector of size 4, containing information of the optimization.
\preformatted{%
W <- spams.fistaGraph(Y,X,W0,graph,return_optim_info = FALSE,...)
v <- spams.fistaGraph(Y,X,W0,graph,return_optim_info = TRUE,...)
W <- v[[1]]
optim_info <- v[[2]]
}}
}
\author{
Julien MAIRAL, 2010 (spams, matlab interface and documentation)
Jean-Paul CHIEZE 2011-2012 (R interface)
}
\note{
    Valid values for the regularization parameter (regul) are:
      "l0", "l1", "l2", "linf", "l2-not-squared", "elastic-net", "fused-lasso",
      "group-lasso-l2", "group-lasso-linf", "sparse-group-lasso-l2",
      "sparse-group-lasso-linf", "l1l2", "l1linf", "l1l2+l1", "l1linf+l1",
      "tree-l0", "tree-l2", "tree-linf", "graph", "graph-ridge", "graph-l2",
      "multi-task-tree", "multi-task-graph", "l1linf-row-column", "trace-norm",
      "trace-norm-vec", "rank", "rank-vec", "none"
}
