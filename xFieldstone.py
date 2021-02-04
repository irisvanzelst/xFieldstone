# overarching main xFieldstone script that depends on and calls all other .py scripts in this folder

# import modules specific to xFieldstone 
import input
import thermal_parameters

# standard modules
import os

# ======== 1D thermal structure of oceanic lithosphere ======== # 
# the 1D model provides us with a lookup table necessary to define the LHS of the 2D problem

# check if lookup table exists 
lookup_table = os.path.isfile('slab_k{0}_Cp{1}_rho{2}.npz'.format(input.option_k,input.option_C_p,input.option_rho))

# run the 1D model to create the lookup table if (and only if) the required lookup table with the desired T-dependence of variables does not exist yet, or the user wants to recalculate the lookup table (for example with a different spatial or temporal resolution)
if input.left_boundary_condition == 1 and (lookup_table == 0 or input.recalculate == 1):
    import thermal_structure_ocean

# ======== 2D kinematic subduction with dynamic wedge ======== # 
import main_solver_Stokes_T_grid
