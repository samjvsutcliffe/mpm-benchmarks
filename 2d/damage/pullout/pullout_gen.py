import pycbg.preprocessing as utl
from pycbg.mesh import Mesh
import numpy as np
import math
dt = 1e-3

def create_Glen2D(self, pset_id=0,density=900,
        youngs_modulus=1e6,
        poisson_ratio=0.3,
        viscosity=1e8,
        viscous_power=1):

	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "Glen2D",
	"density": density,
	"youngs_modulus": youngs_modulus,
	"poisson_ratio": poisson_ratio,
	"viscosity": viscosity,
	"viscous_power":viscous_power 
        })
	
# The usual start of a PyCBG script:
sim = utl.Simulation(title="pullout")

resolution = 50
resolutions = [resolution,resolution ]

# Creating the mesh:


shelf_length = 500
particle_dims = (shelf_length,100.)
domain_dims = (1000.,400.)


node_type = "ED2Q4"
#node_type = "ED2Q16G"
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)],cell_type = node_type)
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,0), ncells=[x//r for x,r in zip(particle_dims,resolutions)],cell_type =node_type)

#Pseudo-2K
#sim.create_mesh(dimensions=domain_dims, ncells=(domain_dims[0]//resolution,domain_dims[1]//resolution,1))
#particle_dims = (500,100,1)
#domain_dims = (600.,200.,1.)

#pmesh = utl.Mesh(dimensions=particle_dims, ncells=(particle_dims[0]//resolution,particle_dims[1]//resolution,1))
# Creating Material Points, could have been done by filling an array manually:
#sim.create_particles(npart_perdim_percell=1)
mps_per_cell = 2

sim.particles = utl.Particles(pmesh,mps_per_cell,directory=sim.directory,particle_type = "")

mps_array = [1,1]
eff_res = [r/(mps_per_cell*mps) for r,mps in zip(resolutions,mps_array)]
origin = [r/(2*mps_per_cell*mps) for r,mps in  zip(resolutions,mps_array)]

sim.particles.positions = np.mgrid[
        origin[0]:particle_dims[0]:eff_res[0],
        origin[1]:particle_dims[1]:eff_res[1]
        #0.5:1:0.5
        #origin[2]:particle_dims[2]:eff_res[2],
        ].reshape(2,-1).T 

sim.init_entity_sets()

sim.particles._filename = "particles.txt"
sim.particles.write_file()



# Creating entity sets (the 2 materials), using lambda functions:
maxwell_particles = sim.entity_sets.create_set(lambda x,y: True, typ="particle")

E = 1e8
nu = 0.325
density = 900
density_water = 999

# The materials properties:
#sim.materials.create_MohrCoulomb3D(pset_id=lower_particles)
sim.materials.create_LinearElastic(pset_id=maxwell_particles,density=density,youngs_modulus=E,poisson_ratio=nu)
#sim.materials.create_Newtonian3D(pset_id=maxwell_particles)

#create_Glen2D(sim.materials,pset_id=maxwell_particles,
#        density=900,
#        youngs_modulus=E,
#        poisson_ratio=nu,
#        #viscosity=2.24e-24,
#        viscosity=visc,
#        viscous_power=3
#        )

# Boundary conditions on nodes entity sets (blocked displacements):
walls = []
walls.append([sim.entity_sets.create_set(lambda x,y: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
#walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [100]])
#walls.append([sim.entity_sets.create_set(lambda x,y: z==lim, typ="node") for lim in [0, sim.mesh.l2]])
for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]

#sim.add_velocity_condition(1,0.0,sim.entity_sets.create_set(lambda x,y: x<=shelf_length and y==100, typ="node"))
#sim.add_friction_condition(1,-1,0.9,sim.entity_sets.create_set(lambda x,y: y==0 and x <=250, typ="node"))
#
#sim.add_force(1,9.81 * (density_water) * (resolution**2),sim.entity_sets.create_set(lambda x,y: x>=shelf_length and y<100, typ="node"))
#sim.add_velocity_condition(0,0.0,sim.entity_sets.create_set(lambda x,y: x==500 and y<100, typ="node"))

# Other simulation parameters (gravity, number of iterations, time step, ..):
sim.set_gravity([0,-9.81])
time = 10
nsteps = time//dt
sim.set_analysis_parameters(dt=dt,type="MPMExplicit2D", nsteps=nsteps, 
        output_step_interval=nsteps/100,
        damping=0.00)
        

#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": 1e8}
sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": E*1e-3}
#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": K*1e-2}
#sim.analysis_params["damping"] = {"type": "Cundall", "damping_factor": 0.05}
sim.post_processing["vtk"] = ["stresses"]

# Save user defined parameters to be reused at the postprocessing stage:
sim.add_custom_parameters({"maxwell_particles": maxwell_particles,
							"walls": walls})
# Final generation of input files:
sim.write_input_file()

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

