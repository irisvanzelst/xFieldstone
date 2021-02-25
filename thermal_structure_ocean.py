# calculate the thermal structure of the ocean according to Richards et al., 2018 JGR: Solid Earth 
# we use a time- and space-centered Crank-Nicholson finite-difference scheme with a predictor-corrector step (Press et al., 1992)

# import all the constants and defined model setup parameters 
import input

# modules
import numpy as np
from numpy import load
import matplotlib.pyplot as plt
import time as timing
import scipy.linalg as la 
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
import thermal_parameters

start       = timing.time()
zeros       = np.zeros
power       = np.power
exp         = np.exp
transpose   = np.transpose 
array       = np.array
full        = np.full
all         = np.all
diagonal    = np.diagonal
solve_banded= la.solve_banded
        
# we will solve the system Mx = D with
# x = the temperature vector (our unknowns)
# M = our matrix of coefficients in front of the temperature (known)
# D = the right hand side vector (known)

# first, we need to build M and D
# M and D are different in the predictor and corrector step 

# ========== time ========== 

nt = int(input.end_time / input.dt + 1)

t = zeros((1,nt))

# ========== initial conditions ========== 

z = np.arange(0,input.nz*input.dz,input.dz)

for i in range(0,input.nz):
    if z[i] == input.depth_melt1:
        index_depth_melt1 = i

Told = zeros((1,input.nz))
for i in range(0,input.nz):
    slope = (input.alpha * input.g) / thermal_parameters.heat_capacity(input.option_C_p,Told[0,i-1]+input.kelvin,0,0)
    if z[i] == 0:
        Told[0,i] = 0. + input.kelvin 
    elif z[i] < input.depth_melt1:  
        Tdepth_melt1 = slope * input.depth_melt1 + input.Tmax  - slope * z[-1] + input.kelvin 
        Told[0,i] = Told[0,i-1] + Tdepth_melt1/(index_depth_melt1+1)
    else:
        Told[0,i] = slope * z[i] + input.Tmax  - slope * z[-1] + input.kelvin


temperature = zeros((nt,input.nz))
temperature[0,:] = Told
print(Told)

dummy = zeros((1,input.nz))

options = zeros((3,input.nz))
options[0,:] = input.option_k
options[1,:] = input.option_C_p
options[2,:] = input.option_rho

for time in range(1,nt):
    print(time)
    t[0,time] = t[0,time] + time*input.dt
    print(t[0,time]/1e6/input.year)

    M = zeros((input.nz, input.nz))     # pre-allocate M array

    A               = zeros((input.nz, 1))
    D               = zeros((input.nz, 1))     # pre-allocate D column vector
    Q               = zeros((input.nz, 1)) 
    B               = zeros((input.nz, 1)) 

    for step in range(0,2):
        # 0 = predictor step 
        # 1 = corrector step 

        # first make the 1D thermal parameters, so we can use them to calculate M later         
        if step == 0:
            # predictor step
            # m = n
            k_m                 = transpose(array(list(map(thermal_parameters.heat_conductivity,options[0,:],Told,dummy,dummy))))
            density_m           = transpose(array(list(map(thermal_parameters.density,options[2,:],Told,dummy,dummy))))
            heat_capacity_m     = transpose(array(list(map(thermal_parameters.heat_capacity,options[1,:],Told,dummy,dummy))))
            dz_m                = full((input.nz,1), input.dz) # current assumption: incompressible
        elif step == 1: 
            # corrector step 
            # m = n+1/2
            k_m                 = (transpose(array(list(map(thermal_parameters.heat_conductivity,options[0,:],Told,dummy,dummy)))) + transpose(array(list(map(thermal_parameters.heat_conductivity,options[0,:],Tdummy,dummy,dummy)))) ) / 2.
            density_m           = (transpose(array(list(map(thermal_parameters.density,options[2,:],Told,dummy,dummy)))) + transpose(array(list(map(thermal_parameters.density,options[2,:],Tdummy,dummy,dummy)))) ) / 2.
            heat_capacity_m     = (transpose(array(list(map(thermal_parameters.heat_capacity,options[1,:],Told,dummy,dummy)))) + transpose(array(list(map(thermal_parameters.heat_capacity,options[1,:],Tdummy,dummy,dummy)))) ) / 2.
            dz_m                = (full((input.nz,1), input.dz) + full((input.nz,1), input.dz) ) / 2.  # current assumption: incompressible
            
        for i in range(0,input.nz):
            if i == 0: 
                start_loop = i
                end_loop   = i+1
            elif i == input.nz-1:
                start_loop = i-1
                end_loop   = i  
            else: 
                start_loop = i-1
                end_loop   = i+1  
                
            for j in range(start_loop,end_loop+1): 
        
                if j == 0:
                    # calculate A  
                    A[j] = input.dt / ( density_m[j] * heat_capacity_m[j] * ( dz_m[j] + dz_m[j] ) )
                elif j > 0:
                    # calculate A 
                    A[j] = input.dt / ( density_m[j] * heat_capacity_m[j] * ( dz_m[j] + dz_m[j-1] ) )   
                
                # ========== BUILD M ========== - coefficient matrix 
                # ========== boundary conditions ==========
                if (i == 0 and j == 0): 
                    # boundary condition at the top 
                    M[i,j] = 1. 
                    D[j]   = input.Ttop + input.kelvin 
                elif (i == input.nz-1 and j == input.nz-1):
                    # boundary condition at the bottom 
                    M[i,j] = 1. 
                    D[j]   = input.Tmax + input.kelvin  
                else:
                    # ========== BUILD M ========== - coefficient matrix 
                    if i - j == 1 and j < input.nz - 2:
                        # T_j+1 
                        if j > 0:
                            M[i,j] = -A[j] * ( (k_m[j] + k_m[j+1] ) / 2. ) / dz_m[j]
                        elif j == 0:
                            M[i,j] = -A[j] * ( (k_m[j] + k_m[j] ) / 2. ) / dz_m[j]
                    elif i == j:
                        # diagonal: T_j
                        M[i,j] = 1. + A[j] * ( ( (k_m[j] + k_m[j+1] ) / 2. ) / dz_m[j] + ( (k_m[j] + k_m[j-1] ) / 2. ) / dz_m[j-1])
                    elif i - j == -1 and j > 1 and j < input.nz:
                        # T_j-1
                        if j < input.nz - 1:
                            M[i,j] = -A[j] * ( (k_m[j] + k_m[j-1] ) / 2. ) / dz_m[j-1]
                        elif j == input.nz - 1:
                            M[i,j] = -A[j] * ( (k_m[j] + k_m[j] ) / 2. ) / dz_m[j]
            
                    # ========== BUILD D ========== - right hand side vector 
                    # D consists of multiple components
                    # we say D = T + A * Q + B
                    if j > 0 and j < input.nz-1:
                        Q[j] = ( ( (k_m[j] + k_m[j+1] ) / 2. ) / dz_m[j] )*Told[0,j+1] - (( ( (k_m[j] + k_m[j+1] ) / 2. ) / dz_m[j] ) + ( ( (k_m[j] + k_m[j-1] ) / 2. ) / dz_m[j-1] )) * Told[0,j] + ( ( (k_m[j] + k_m[j-1] ) / 2. ) / dz_m[j-1] )*Told[0,j-1]
            
                        # B - correction that represents the second term on the right-hand side on the equation 
                        if step == 0:
                            # B - predictor step 
                            B[j] = -Told[0,j] * ( thermal_parameters.density(input.option_rho,Told[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Told[0,j],0.,0.) - thermal_parameters.density(input.option_rho,temperature[time-1,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,temperature[time-1,j],0.,0.)) / (thermal_parameters.density(input.option_rho,Told[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Told[0,j],0.,0.))
                        elif step == 1: 
                            # B - corrector step 
                            B[j] = - ( (Tdummy[0,j] + Told[0,j]) * ( thermal_parameters.density(input.option_rho,Tdummy[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Tdummy[0,j],0.,0.) - thermal_parameters.density(input.option_rho,Told[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Told[0,j],0.,0.)) ) / ( thermal_parameters.density(input.option_rho,Tdummy[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Tdummy[0,j],0.,0.) + thermal_parameters.density(input.option_rho,Told[0,j],0.,0.) * thermal_parameters.heat_capacity(input.option_C_p,Told[0,j],0.,0.))

                        D[j] = Told[0,j] + A[j] * Q[j] + B[j]
      

        # ========== Solve system of equations using numpy.linalg.solve ==========
        
        if input.option_1D_solve == 1:
            # brute force solve
            Tnew = transpose(np.linalg.solve(M, D))    # solve Mx = D for x
        elif input.option_1D_solve == 2:    
            # make ordered banded diagonal matrix MD from M 
            upper = 1
            lower = 1
            n = M.shape[1]
            assert(all(M.shape ==(n,n)))
    
            ab = zeros((2*n-1, n))
    
            for i in range(n):
                ab[i,(n-1)-i:] = diagonal(M,(n-1)-i)
        
            for i in range(n-1): 
                ab[(2*n-2)-i,:i+1] = diagonal(M,i-(n-1))

            mid_row_inx = int(ab.shape[0]/2)
            upper_rows = [mid_row_inx - i for i in range(1, upper+1)]
            upper_rows.reverse()
            upper_rows.append(mid_row_inx)
            lower_rows = [mid_row_inx + i for i in range(1, lower+1)]
            keep_rows = upper_rows+lower_rows
            ab = ab[keep_rows,:]
        
            Tnew = transpose(solve_banded((1,1),ab, D))
        
        if step == 0:
            Tdummy = Tnew 
            #temperature[time,:] = Tnew
            #print(Tnew)
            #Told = Tnew
        elif step == 1:
            #print(time) 
            temperature[time,:] = Tnew
            Told = Tnew
            Tdummy = []
            print(Tnew)
        
        # end for-loop predictor-corrector step 

print("time: %.3f s" % (timing.time() - start))
print(Tnew)

# ========== save thermal structure ==========

np.savez('databases/slab_k{0}_Cp{1}_rho{2}.npz'.format(input.option_k,input.option_C_p,input.option_rho),temperature=temperature,t=t,z=z)
#data = np.load('thermal_structure_ocean_k{0}_Cp{1}_rho{2}.npz'.format(input.option_k,input.option_C_p,input.option_rho))

# ========== figure - visualisation ==========

if input.visualisation_1D == 1:

    X, Y  = np.meshgrid(t/1e6/input.year, z/1e3)     # make a mesh of X and Y for easy plotting in the right units (i.e., t in Myr and z in km) 
    Zaxis = transpose(temperature-input.kelvin)                   # convert temperature back into C for easy plotting

    from IPython.display import set_matplotlib_formats
    #matplotlib inline
    set_matplotlib_formats('png')

    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    # Initialize plot objects
    rcParams['figure.figsize'] = 10, 5 # sets plot size
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # define contours
    levels = array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300.])

    import matplotlib.cm as cm # matplotlib's color map library
    cpf = ax.contourf(X,Y,Zaxis, len(levels), cmap='RdBu_r')

    line_colors = ['black' for l in cpf.levels]

    # Generate a contour plot
    cp = ax.contour(X, Y, Zaxis, levels=levels, colors=line_colors)
    #plt.colorbar(orientation='vertical') 
    ax.clabel(cp, fontsize=10, colors=line_colors)
    ax.set_xlabel('Age (Ma)')
    ax.set_ylabel('Depth below seafloor (km)')
    plt.gca().invert_yaxis()
    plt.show()

