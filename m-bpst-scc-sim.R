library(Triangulation)
library(BPST)
library(ImageSCC)
library(tidyverse)
library(MASS)
library(pracma)
library(VGAM)
library(argparser)
library(Rcpp)
library(RcppArmadillo)

##### Parameters for Slurm job #####
args <- arg_parser('Robust FDA with BPST simulation') %>%
  add_argument('--results-dir', type = 'character', help = 'Result path', default = file.path('m-bpst-scc-sim')) %>%
  add_argument('--job', type = 'integer', help = 'job ID.') %>%
  add_argument('--task', type = 'integer', help = 'Task number.') %>%
  parse_args()

# let the program sleep such that the C++ cache is not accessed concurrently 
Sys.sleep(args$task)

source("robust_bpst_fda_functions.R")
# sourceCpp("cpp_func_cls.cpp")
sourceCpp("cpp_func_cls.cpp", cacheDir = "/scratch/ylong5/rcpp-cache/robust-fda")

# parameters from Slurm file to create subtasks
nsim <- 20
con_rate <- c(0, 0.1)
n <- c(100, 200, 400)
start <- nsim * c(0:(100/nsim-1))

func_type <- c('quadratic', 'exponential', 'cubic')
cont_type <- c('stripe', 'square2')
# tri <- c(2, 3)
# d <- c(4, 5, 6)
settings <- expand.grid(n=n, start=start, func_type=func_type, con_rate=con_rate, cont_type=cont_type)

# for sensitivity analysis
# settings <- expand.grid(n=n, start=start, con_rate=con_rate, cont_type=cont_type, tri=tri, d=d) %>%
#   filter((tri == 2 & d == 4) | (tri == 2 & d == 6) | (tri == 3 & d == 4) | (tri == 3 & d == 5))

cont_rate <- settings$con_rate[args$task]  # contamination rate
n <- settings$n[args$task]
cont_type <- settings$cont_type[args$task]
true_func_type <- settings$func_type[args$task]
start <- settings$start[args$task]
# d.est <- settings$d[args$task]
# triangulation <- settings$tri[args$task]

# # Triangulation information
# if (triangulation == 2) {
#   V.est=as.matrix(Brain.V2)
#   Tr.est=as.matrix(Brain.Tr2)
#   V.band=as.matrix(Brain.V2)
#   Tr.band=as.matrix(Brain.Tr2)
# } else {
#   V.est=as.matrix(Brain.V3)
#   Tr.est=as.matrix(Brain.Tr3)
#   V.band=as.matrix(Brain.V3)
#   Tr.band=as.matrix(Brain.Tr3)
# }

V.est=as.matrix(Brain.V2)
Tr.est=as.matrix(Brain.Tr2)
V.band=as.matrix(Brain.V2)
Tr.band=as.matrix(Brain.Tr2)


# Sample location information
n1=40; n2=40;
# n1=79; n2=79;
npix=n1*n2
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
uu.mtx=matrix(uu,n2,n1)
vv.mtx=matrix(vv,n2,n1)
Z=as.matrix(cbind(uu,vv))
ind.inside=inVT(V.est,Tr.est,Z[,1],Z[,2])$ind.inside

# parameter for simulation
d.est <- 5
d.band <- 2
r <- 1
# n <- 200 # of images
lam1 <- 0.5
lam2 <- 0.2
lambda <- 10^{c(-5:2)} # penalty parameter
alpha <- 0.05 # alpha for SCB

# contamination settings
a_0 <- c(0.3, 0.6) # location of stripe contamination
a_1 <- 0.3 # location of square contamination (min)
a_2 <- 0.5 # location of square contamination (max)
a_square <- matrix(c(0.3, 0.4, 0.6, 0.7), nrow=2, byrow=TRUE)
clean_sd <- .2
cont_dist_min <- 10
cont_dist_max <- 20

# model parameter
huber_const <- .1

n
nsim
cont_rate
cont_type
true_func_type
huber_const
clean_sd

n_b <- 1000
bw_ls <- 0
bw_m_bs <- 0
coverage_m_bs_pw <- c()
coverage_ls_pw <- c()
coverage_ls <- 0
# seed_all <- c()

# initialize the matrices to reduce repetitive computation
Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
B.est <- Bfull.est$B
ind.inside.est <- Bfull.est$Ind.inside
Q2.est <- Bfull.est$Q2
K.est <- Bfull.est$K
B_til <- B.est%*%Q2.est
D <- t(Q2.est)%*%K.est%*%Q2.est

for (i in (start+1):(start+nsim)) {

  set.seed(i)
  
  cat('Iteration No.',i,'\n')
  ptm0=Sys.time()
  gen_sample <- generate_sample(
    Z, ind.inside, n=n, type=true_func_type, 
    error_type='bpst', lam1=lam1, lam2=lam2, 
    con_rate=cont_rate, con_type=cont_type, clean_sd=clean_sd,
    a_0=a_0,
    a_1=a_1, a_2=a_2, 
    a_square=a_square,
    rand_cont_ct=192,
    cont_dist_min=cont_dist_min,
    cont_dist_max=cont_dist_max)
  true_func <- gen_sample$true_func
  Y <- gen_sample$y
  Y_inside <- Y[,ind.inside]
  
  # bootstrap method
  rscc_fit <- robust_imagescc(
    B.est, ind.inside, Q2.est, K.est, Y, lambda=lambda, 
    huber_const=huber_const,
    maxit = 1000, tol = 1e-6
  )
  
  cc.l <- rscc_fit$cc.l
  cc.u <- rscc_fit$cc.u
  
  coverage_m_bs_pw <- rbind(
    coverage_m_bs_pw,
    (cc.l<true_func[ind.inside]) & (cc.u>true_func[ind.inside]))

  bw_m_bs <- bw_m_bs + rscc_fit$bw
  
  out_ls=scc.image(
    Ya=Y,Z=Z,V.est.a=V.est,Tr.est.a=Tr.est,V.band.a=V.band,Tr.band.a=Tr.band,
    d.est=d.est,d.band=d.band,r=r,
    penalty=TRUE,lambda=lambda,alpha.grid=alpha,adjust.sigma=TRUE)
  
  coverage_ls_pw <- rbind(
    coverage_ls_pw,
    (out_ls$scc[,1,1]<true_func[ind.inside]) & (out_ls$scc[,2,1]>true_func[ind.inside]))
  
  coverage_ls <- coverage_ls + 
    apply(out_ls$scc,3,FUN=function(scc){
      (sum((scc[,1]<true_func[ind.inside]) & (scc[,2]>true_func[ind.inside]))/length(ind.inside)) == 1
    })
  bw_ls <- bw_ls + out_ls$bw
  
  ptm1=Sys.time()
  cat('Time: ',ptm1-ptm0,'\n')
}

input <- list(
  # seed_all = seed_all,
  sample_size = n,
  # func_type = true_func_type,
  cont_rate = cont_rate,
  n_sim = nsim,
  huber_const = huber_const,
  # parameter for sample generation
  cont_type = cont_type, # location of stripe contamination
  a_0 = a_0,
  a_1 = a_1, # location of square contamination (min)
  a_2 = a_2, # location of square contamination (max)
  a_square = matrix(c(0.3, 0.4, 
                      0.6, 0.7), nrow=2, byrow=TRUE),
  clean_sd = clean_sd,
  cont_dist_min = cont_dist_min,
  cont_dist_max = cont_dist_max,
  triangulation = triangulation,
  d.est = d.est,
  n1=n1,
  n2=n2
)


output <- list(
  coverage_ls = coverage_ls,
  coverage_ls_pw = coverage_ls_pw,
  coverage_m_bs_pw = coverage_m_bs_pw,
  bw_ls = bw_ls/nsim,
  bw_m_bs = bw_m_bs/nsim
)

all_results <- list(input = input, output = output)

dir.create(args$results_dir, showWarnings = F)
job_id <- sprintf('job=%01d-task=%02d', args$job, args$task)
job_file <- file.path(args$results_dir, paste(job_id, 'rds', sep = '.'))

saveRDS(all_results, file = job_file)
