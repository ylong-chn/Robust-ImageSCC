# Load the latest version of R in Git Bash
# 0. check the latest version of GNU module: module avail 
# 1. load the latest version of GNU module: module load gnu10/10.3.0-ya
# 2. load the latest version of R: module load r/4.1.2-dx

# Install required packages in Hopper
# tidyverse, VGAM, fields, ImageSCC (requires BPST and Triangulation), robustbase

library(Triangulation)
library(BPST)
library(mgcv)
library(ImageSCC)
library(tidyverse)
library(keras3)
library(MASS)
library(pracma)
library(VGAM)
library(tensorflow)
library(argparser)
library(Rcpp)
library(RcppArmadillo)

##### Parameters for Slurm job #####
args <- arg_parser('Robust FDA with BPST simulation') %>%
  add_argument('--results-dir', type = 'character', help = 'Result path', default = file.path('m-bpst-est-sim-r2-hs')) %>%
  add_argument('--job', type = 'integer', help = 'job ID.') %>%
  add_argument('--task', type = 'integer', help = 'Task number.') %>%
  # add_argument('--nsim', type = 'integer', help = 'Number of simulation runs.') %>%
  parse_args()

# let the program sleep such that the C++ cache is not accessed concurrently 
Sys.sleep(args$task)

source("robust_bpst_fda_functions.R")
# sourceCpp("cpp_func_cls.cpp")
sourceCpp("cpp_func_cls.cpp", cacheDir = "/scratch/ylong5/rcpp-cache/robust-fda")
# source("rdnn.R")

# seed <- args$job + args$task
# set.seed(seed)

cat('task:', args$task, '...\n')
# cat('seed:', seed, '...\n')

##### Parameters #####
# parameters from Slurm file to create sub-tasks
con_rate <- c(0.05, 0.1, 0.15, 0)
# con_rate <- c(0.15, 0)
# con_rate <- 0
n <- c(100, 200, 400)
# func_type <- c('quadratic', 'exponential', 'cubic', 'sine', 'horseshoe')
# func_type <- c('quadratic', 'exponential', 'cubic', 'sine')
# func_type <- c('sine', 'horseshoe')
# func_type <- c('quadratic', 'exponential', 'cubic')
# func_type <- 'exponential'
func_type <- 'horseshoe'
cont_type <- c('stripe', 'square2', 'cauchy', 'slash')
# cont_type <- c('stripe', 'square2', 'cauchy')
# cont_type <- 'stripe'
nsim <- 100
start <- nsim * c(0:(100/nsim-1))

# for sensitivity analysis, add more variables to the grid expansion
# tri <- c(2, 3)
# d <- c(4, 5, 6)
tri <- 2
d <- 5
# huber_const <- c(0.01, 0.1, 1)
huber_const <- 0.1
# cont_type <- c('square2')
settings <- expand.grid(n=n, start=start, tri=tri, d=d, huber_const=huber_const, cont_type=cont_type, func_type=func_type, con_rate=con_rate) %>% 
  filter((tri == 2 & d == 4) | (tri == 2 & d == 5) | (tri == 2 & d == 6) | (tri == 3 & d == 4) | (tri == 3 & d == 5)) %>% 
  filter(con_rate %in% c(0.05, 0.1, 0.15) | (con_rate==0 & cont_type=='stripe'))
# settings <- expand.grid(n=n, start=start, mix_weight=mix_weight, cont_type=cont_type, tri=tri, d=d, huber_const=huber_const) %>%
  # filter((tri == 2 & d == 4) | (tri == 2 & d == 5) | (tri == 2 & d == 6) | (tri == 3 & d == 4) | (tri == 3 & d == 5))
# settings <- expand.grid(n=n, start=start, mix_weight=mix_weight, cont_type=cont_type, tri=tri, d=d, huber_const=huber_const) %>%
#   filter((tri == 2 & d == 4) | (tri == 2 & d == 6) | (tri == 3 & d == 4) | (tri == 3 & d == 5))
# settings <- expand.grid(n=n, start=start, func_type=func_type, cont_type=cont_type, con_rate=con_rate)

cont_rate <- settings$con_rate[args$task]  
n <- settings$n[args$task]
cont_type <- settings$cont_type[args$task]
true_func_type <- settings$func_type[args$task]
start <- settings$start[args$task]
d.est <- settings$d[args$task]
triangulation <- as.numeric(settings$tri[args$task])
fixed_k <- settings$huber_const[args$task]

# for sensitivity analysis
# true_func_type <- 'exponential'
# 
# # Triangulation information

# triangulation 
if (true_func_type == 'horseshoe') {
  data('horseshoe')
  # n_tri <- 5 # fineness of the triangulation 
  n_tri <- triangulation + 3
  VT <- TriMesh(horseshoe, n_tri)
  V.est <- VT$V
  Tr.est <- VT$Tr
} else {
  # # triangulation 
  # data(Brain.V2)
  # data(Brain.Tr2)
  # V.est <- Brain.V2
  # Tr.est <- Brain.Tr2
  if (triangulation == 1) {
    V.est=as.matrix(Brain.V1)
    Tr.est=as.matrix(Brain.Tr1)
    V.band=as.matrix(Brain.V1)
    Tr.band=as.matrix(Brain.Tr1)
  } else if (triangulation == 2) {
    V.est=as.matrix(Brain.V2)
    Tr.est=as.matrix(Brain.Tr2)
    V.band=as.matrix(Brain.V2)
    Tr.band=as.matrix(Brain.Tr2)
  } else {
    V.est=as.matrix(Brain.V3)
    Tr.est=as.matrix(Brain.Tr3)
    V.band=as.matrix(Brain.V3)
    Tr.band=as.matrix(Brain.Tr3)
  }
}

# image resolution
# n1 <- 80
# n2 <- 80
n1 <- 40
n2 <- 40

# parameters for bivariate spline over triangulation
# d.est <- 5
r <- 1
# r <- 0 

# parameter for simulation
lam1 <- 0.5
lam2 <- 0.2
# n_sim <- 100

if (true_func_type == 'horseshoe') {
  xm <- seq(-1, 3.5, length=n1)
  yn <- seq(-1, 1, length=n2)
  gridpoints <- expand.grid(xm, yn) %>% arrange_all() %>% as.matrix()
  
  z <- gridpoints
} else {
  z <- generate_grid(n1, n2)
}

# Basis function and preliminary computation
basis_start <- Sys.time()
Bfull <- basis(V=V.est, Tr=Tr.est, d=d.est, r=r, Z=z)
B <- Bfull$B
ind_inside <- Bfull$Ind.inside
Q2 <- Bfull$Q2
K <- Bfull$K
B_til <- B%*%Q2
D <- t(Q2)%*%K%*%Q2
basis_end <- Sys.time()
basis_time <- as.numeric(difftime(basis_end, basis_start, units = 'mins'))

# huber const
# fixed_k <- 0.1

# parameter for sample generation
lambda <- 10^{c(-6:3)} # penalty parameter

# contamination settings
a_0 <- c(0.3, 0.6) # location of stripe contamination
# segment <- FALSE # entire stripe or only segments on the stripe is contaminated 
a_1 <- 0.3 # location of square contamination (min)
a_2 <- 0.5 # location of square contamination (max)
a_square <- matrix(c(0.3, 0.4, 0.6, 0.7), nrow=2, byrow=TRUE)
clean_sd <- .2
cont_dist_min <- 10
cont_dist_max <- 20

if (cont_type %in% c('cauchy', 'slash')) {
  mix_weight <- cont_rate
  cont_rate <- 1
} else mix_weight <- 0

cat('n =', n, '...\n')
cat('cont_type =', as.character(cont_type), '...\n')
cat('cont_rate =', cont_rate, '...\n')
cat('mix_weight =', mix_weight, '...\n')
cat('start =', start, '...\n')
cat('number of simulation =', nsim, '...\n')
cat('d =', d.est, '...\n')
cat('triangulation =', triangulation, '...\n')
cat('N =', n1*n2, '...\n')
cat('Huber constant =', fixed_k, '...\n')

# other parameters for rdnn functions
Y_dim <- 2 # dimension of the problem
# Grid <- list(xm, yn)
N <- c(n1, n2)
L <- 1 # number of hidden layers
p <- 500 # width within each layer
s <- 0.5 # dropout rate
epoch <- 100
batch <- 256

# boundary for soap-film smoothing
add_boundary_column <- function(coordinates) {
  # Convert to a data frame for easier manipulation with dplyr
  coordinates_df <- as.data.frame(coordinates)
  colnames(coordinates_df) <- c("x", "y")
  
  # Find boundary points for each unique x
  boundary_points_x <- coordinates_df %>%
    group_by(x) %>%
    filter(y == min(y) | y == max(y)) %>%
    ungroup()
  
  # Find boundary points for each unique y
  boundary_points_y <- coordinates_df %>%
    group_by(y) %>%
    filter(x == min(x) | x == max(x)) %>%
    ungroup()
  
  # Combine boundary points and remove duplicates
  boundary_points <- bind_rows(boundary_points_x, boundary_points_y) %>%
    distinct() %>%
    mutate(boundary = 1)  # Mark boundary points with 1
  
  # Left join with original coordinates and replace NA with 0 for non-boundary points
  coordinates_with_boundary <- coordinates_df %>%
    left_join(boundary_points, by = c("x", "y")) %>%
    mutate(boundary = ifelse(is.na(boundary), 0, boundary))
  
  return(coordinates_with_boundary)
}

z_inside <- data.frame(z[ind_inside,])
# boundary <- fs.boundary()
boundary <- fs.boundary(r0=0.05, l = 3.05)
x <- z_inside[,1]
y <- z_inside[,2]
inside <- inSide(boundary, x, y)
knots <- z_inside[inside,]
# y_max <- max(knots[,2])
# y_min <- min(knots[,2])
# knots <- knots[!(knots[,2] %in% c(y_max, y_min)),]

z_boundary <- add_boundary_column(knots)
ind.internal <- z_boundary$boundary == 0
knot_place <- knots[ind.internal,]

# inflate_boundary <- function(boundary, scale=1.01) {
#   cx <- mean(boundary$x)
#   cy <- mean(boundary$y)
#   data.frame(
#     x = cx + scale*(boundary$x - cx),
#     y = cy + scale*(boundary$y - cy)
#   )
# }
# 
# boundary_inf <- inflate_boundary(boundary, scale=1.1)

# boundary <- fs.boundary(r0=0.05, l = 3.05)
# plot(boundary, type='l')
# points(df$x, df$y, col='black')
# boundary_r0 <- fs.boundary(r0=0.05, l = 3.05)
# lines(boundary_r0, col='blue')
# plot(boundary_inf$x, boundary_inf$y, type='l')
# points(df$x, df$y, col='black')
# lines(boundary$x, boundary$y, col='blue')
# points(knot_place, col='red')

##### Simulation #####
sim_result <- c()
# chosen_ks <- c()

# for (i in c(1:52, 54:101)) {
for(i in 1:100) {
  
  set.seed(i)
  
  cat('Iteration No.',i,'\t')
  ptm0 <- Sys.time()

  # generate sample 
  gen_sample <- generate_sample(
    z, ind_inside, n=n, type=true_func_type, 
    error_type='bpst', lam1=lam1, lam2=lam2, 
    con_rate=cont_rate, con_type=cont_type, clean_sd=clean_sd,
    a_0=a_0,
    a_1=a_1, a_2=a_2, 
    a_square=a_square,
    rand_cont_ct=192,
    cont_dist_min=cont_dist_min,
    cont_dist_max=cont_dist_max,
    mix_weight=mix_weight)
  
  true_func <- gen_sample$true_func
  Y <- gen_sample$y
  Y_inside <- Y[,ind_inside]
  
  
  
  # LS estimator
  mfit1 <- fit.image(matrix(apply(Y,2,mean),nrow=1),Z=z,V.est,Tr.est,d.est,r,lambda=lambda,Bfull=Bfull)
  
  # M estimator
  huber_sol <- mean_func_m_est(
    B, ind_inside, Q2, K, Y, lambda=lambda,
    huber_const=fixed_k,
    maxit = 1000, tol = 1e-6
  )
  
  # thin-plate
  mfit.tp <- gam(drop(apply(Y,2,mean)) ~ s(z[,1], z[,2], bs = "tp", k = 800))
  
  # soap-film
  # mfit.sf <- gam(drop(apply(Y,2,mean)) ~ s(z[,1], z[,2], bs = "so", xt = list(bnd = list(boundary))), 
  #                knots = list(x = knot_place$x, y = knot_place$y))
  df <- data.frame(
    x    = z_inside[,1],
    y    = z_inside[,2],
    Ybar = drop(apply(Y_inside, 2, mean))
  )
  
  mfit.sf <- gam(Ybar ~ s(x, y, bs = "so", xt = list(bnd = list(boundary))), 
                 data = df, 
                 knots = list(x = knot_place[,1], y = knot_place[,2]))

  # tensor-product
  mfit.bts <- gam(drop(apply(Y,2,mean)) ~ te(z[,1], z[,2], bs = "bs", m = 2, k = 13))

  # NN estimator
  rdnn_sol <- rFDADNN(Data=NULL, Y_ins=Y[,ind_inside], z_ins=z[ind_inside,],
                      Y_dim, Grid=NULL, N, n, L, p, s, epoch, batch, huber_const=fixed_k, quantile=NULL)
  
  
  ptm1 <- Sys.time()
  cat('Time: ',ptm1-ptm0,'\n')
  
  # append result
  sim_result <- rbind(
    sim_result,
    data.frame(
      func_type=true_func_type,
      cont_type=cont_type,
      cont_rate=cont_rate,
      mix_weight=mix_weight,
      huber_const=fixed_k,
      triangulation=triangulation,
      d.est=d.est,
      r=r,
      sample_size=n,
      sim_ct=i,
      N=n1,
      l2_ls=mean((mfit1$Yhat-true_func[ind_inside])^2),
      l2_m=mean((drop(huber_sol$est_func_robust)-true_func[ind_inside])^2),
      l2_tp=mean((drop(mfit.tp$fitted.values)-true_func[ind_inside])^2),
      ls_sf=mean((drop(mfit.sf$fitted.values)-true_func[ind_inside])^2),
      l2_bts=mean((drop(mfit.bts$fitted.values)-true_func[ind_inside])^2),
      l2_nn=mean((as.numeric(rdnn_sol$estimation)-true_func[ind_inside])^2)
    )
  )
}

input <- list(
  sample_size = n,
  # parameter for sample generation
  a_0 = a_0, # location of stripe contamination
  a_1 = a_1, # location of square contamination (min)
  a_2 = a_2, # location of square contamination (max)
  a_square = a_square,
  # mix_weight = mix_weight, # weight of fat-tail distribution
  clean_sd = clean_sd,
  cont_dist_min = cont_dist_min,
  cont_dist_max = cont_dist_max,
  n1=n1,
  n2=n2
)

output <- sim_result
all_results <- list(input = input, output = output)

##### Save Results #####
dir.create(args$results_dir, showWarnings = F)
job_id <- sprintf('job=%01d-task=%02d', args$job, args$task)
job_file <- file.path(args$results_dir, paste(job_id, 'rds', sep = '.'))
saveRDS(all_results, file = job_file)
