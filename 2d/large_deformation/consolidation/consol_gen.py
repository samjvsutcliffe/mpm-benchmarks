import pycbg.preprocessing as utl
from pycbg.mesh import Mesh
import numpy as np
import math
dt = 1e-3

# The usual start of a PyCBG script:
sim = utl.Simulation(title="consol")

length = 50
resolution = 0.5

resolutions = [resolution,resolution]
# Creating the mesh:

#shelf_length = resolution
particle_dims   = (resolution,length)
domain_dims     = (resolution,length + (2*resolution))

#length = 1
#particle_dims = (resolution,resolution*2)
#domain_dims =   (resolution,resolution*3)


#node_type = "ED2Q4"
#node_type = "ED2Q16G"
node_type = "ED2GIMP"
#node_type = "ED3H64G"
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)],cell_type = node_type)
#gen_gimp_mesh(domain_dims,[x//r for x,r in zip(domain_dims,resolutions)],"consol",sim.mesh)
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,0), ncells=[x//r for x,r in zip(particle_dims,resolutions)],cell_type=node_type)

mps_per_cell = 2
sim.particles = utl.Particles(pmesh,mps_per_cell,directory=sim.directory,particle_type="FS")

#mps_array = [1,1]
#eff_res = [r/(mps_per_cell*mps) for r,mps in zip(resolutions,mps_array)]
#origin = [r/(2*mps_per_cell*mps) for r,mps in  zip(resolutions,mps_array)]
#
#sim.particles.positions = np.mgrid[
#        origin[0]:particle_dims[0]:eff_res[0],
#        origin[1]:particle_dims[1]:eff_res[1]
#        #0.5:1:0.5
#        #origin[2]:particle_dims[2]:eff_res[2],
#        ].reshape(2,-1).T 

sim.init_entity_sets()

sim.particles._filename = "particles.txt"
sim.particles.write_file()



# Creating entity sets (the 2 materials), using lambda functions:
maxwell_particles = sim.entity_sets.create_set(lambda x,y: True, typ="particle")

E = 1e6
nu = 0.0
density = 800
density_water = 999

# The materials properties:
sim.materials.create_LinearElastic(pset_id=maxwell_particles,density=density,youngs_modulus=E,poisson_ratio=nu)

walls = []

walls.append([sim.entity_sets.create_set(lambda x,y: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
#walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [0, sim.mesh.l2]])
for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]

# Other simulation parameters (gravity, number of iterations, time step, ..):
sim.set_gravity([0,-10])
time = 100
nsteps = time//dt
sim.set_analysis_parameters(dt=dt,type="MPMExplicit2D", nsteps=nsteps, 
        output_step_interval=nsteps/100,
        damping=0.00)
        

sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": E*1e-3}
sim.analysis_params["velocity_update"] = True
#sim.analysis_params["damping"] = {"type": "Cundall", "damping_factor": 0.05}
sim.post_processing["vtk"] = ["stresses","volume"]

# Save user defined parameters to be reused at the postprocessing stage:
sim.add_custom_parameters({"maxwell_particles": maxwell_particles,
							"walls": walls})
# Final generation of input files:
sim.write_input_file()

l = E * nu / ((1+nu)*(1-2*nu))
mu = E / (2*(1+nu))

bulk_stiffness = ((3 * l)+(2*mu))/3
min_x = min(resolutions[0:1])
speed_of_sound = math.sqrt(bulk_stiffness/density)

courant = speed_of_sound * dt/min_x

print("Ice properties")
print("bulk modulus {}".format(bulk_stiffness))
print("speed of sound {}".format(speed_of_sound))
print("Courant number of {}".format(courant))

print(
"""
Simulation properties:
Domain size {}
Particle size {}
dt {}
resolutions {}
E {}
nu {}
density {}
""".format(domain_dims,particle_dims,dt,resolutions,E,nu,density))

