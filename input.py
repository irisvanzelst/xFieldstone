# input for xFieldstone 

# ========== constants ========== #
year        = 365.25 * 24. * 60. * 60.  
cm          = 0.01
kelvin      = 273.15 
g           = 9.81 
Rgas        = 8.3145
eps         = 1e-9 # a little epsilon equivalent to zero

# ========== model specific constants (for the entire model run) ========== #
Ttop        = 0.      #C 
Tmax        = 1300.   #C
alpha       = 0.005

# ========== temperature dependent variables ========== # 

option_k    = 0 
# 0 = constant according to Parsons & Slater 1977; 
# 1 = T-dependent with k from (simplified) Hofmeister 1999 as posed by McKenzie et al., 2005
# 2 = T-dependent with k from Xu et al., 2004 

option_C_p  = 0 
# 0 = constant according to Parsons & Slater 1977
# 1 = 100% forsterite Berman 1988
# 2 = 100% fayalite Berman 1988
# 3 = 89% forsterite & 11% fayalite Berman 1988
# 4 = 100% forsterite Berman & Aranovich 1996
# 5 = 100% fayalite Berman & Aranovich 1996
# 6 = 89% forsterite & 11% fayalite Berman & Aranovich 1996

option_rho  = 0 
# 0 = constant according to Parsons & Slater 1977
# 1 = T-dependent based on McKenzie et al., 2005

# ======== 1D thermal structure of oceanic lithosphere ======== # 

recalculate = 0
# do you want to recalculate the lookup table for the LHS of the 2D problem?
# 0 = no 
# 1 = yes 

# model setup 
dz          = 1000.                # m 
dt          = 1000. * year        # year
end_time    = 200. * 1e6 * year    # million years 
nz          = 107                  # size of the spatial array (in km; max depth, assuming dz = 1km) 
depth_melt1 = 0e3 

# how to solve the model 
option_1D_solve = 2
# 1 = brute force
# 2 = tridiagonal banded matrix solve

# visualisation 
visualisation_1D = 0

# ======== 2D kinematic subduction with dynamic wedge ======== #

# resolution 
nelx = 66 #264
nely = 60 #240

# model setup
slab_vel = 5. * cm / year   # velocity of the subducting slab 
slab_age = 50.e6 * year   # slab age

left_boundary_condition = 0
# 0 = Van Keken half-space cooling model 
# 1 = 1D FD model of oceanic lithosphere cooling 

# ----- rheology ----- #
case = '2a' 
# 1a = analytical cornerflow model 
# 1b = dynamical flow in isoviscous wedge I 
# 1c = dynamical flow in isoviscous wedge II 
# 2a = dynamical flow with diffusion creep 
# 2b = dynamical flow with dislocation creep 
# 2c = dynamical flow with diffusion and dislocation creep (need to be made!)

# diffusion creep 
Q_diff  =   335e3
A_diff  =   1.32043e9

# dislocation creep
Q_disl  =   540e3
A_disl  =   28968.6
n_disl  =   3.5

# ----- numerical considerations ----- # 
eta_max =   1.e26   # maximum viscosity in the model 
eta_def =   1.e21   # default viscosity 
relax   =   0.8 
niter   =   50       # number of iterations 
tol     =   1e-5    # tolerance 












