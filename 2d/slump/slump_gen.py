import pycbg.preprocessing as utl
import numpy as np
import math
dt = 1e-3
def create_Maxwell3D(self, pset_id=0,density=900,elasticity=1e6,viscosity=1e8):
	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "Maxwell3D",
	"density": density,
	"elasticity": elasticity,
	"viscosity": viscosity})
	
def create_Newtonian3D(self, pset_id=0, density=1.225, 
									bulk_modulus=1.42e5, 
									dynamic_viscosity=1.81e-5):
	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
		"type": "Newtonian3D",
		"densimy": density,
		"bulk_modulus": bulk_modulus * dt,
		"dynamic_viscosity": dynamic_viscosity})
# The usual start of a PyCBG script:
sim = utl.Simulation(title="slump")

resolution = 10
resolutions = [20,20,1]

# Creating the mesh:

particle_dims = (500.,100.,1)
domain_dims = (600.,200.,1)
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)])
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,0), ncells=[x//r for x,r in zip(particle_dims,resolutions)])

#Pseudo-2d
#sim.create_mesh(dimensions=domain_dims, ncells=(domain_dims[0]//resolution,domain_dims[1]//resolution,1))
#particle_dims = (500,100,1)
#domain_dims = (600.,200.,1.)
#pmesh = utl.Mesh(dimensions=particle_dims, ncells=(particle_dims[0]//resolution,particle_dims[1]//resolution,1))
# Creating Material Points, could have been done by filling an array manually:
#sim.create_particles(npart_perdim_percell=1)
sim.particles = utl.Particles(pmesh,4,directory=sim.directory)
mps_per_cell = 1
eff_res = [r/mps_per_cell for r in resolutions]
origin = [r/(2*mps_per_cell) for r in resolutions]

sim.particles.particles = np.mgrid[
        origin[0]:particle_dims[0]:eff_res[0],
        origin[1]:particle_dims[1]:eff_res[1],
        0.5:1:0.5].reshape(3,-1).T 

sim.particles._filename = "particles.txt"
sim.particles.write_file(filename=sim._file_prefix+"particles")



# Creating entity sets (the 2 materials), using lambda functions:
sim.init_entity_sets()
maxwell_particles = sim.entity_sets.create_set(lambda x,y,z: True, typ="particle")

# The materials properties:
#sim.materials.create_MohrCoulomb3D(pset_id=lower_particles)
sim.materials.create_LinearElastic3D(pset_id=maxwell_particles,density=900,youngs_modulus=9.5e9,poisson_ratio=.3)
#sim.materials.create_Newtonian3D(pset_id=maxwell_particles)
#create_Newtonian3D(sim.materials,pset_id=maxwell_particles)
#create_Maxwell3D(sim.materials,pset_id=maxwell_particles)

# Boundary conditions on nodes entity sets (blocked displacements):
walls = []
walls.append([sim.entity_sets.create_set(lambda x,y,z: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
walls.append([sim.entity_sets.create_set(lambda x,y,z: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
walls.append([sim.entity_sets.create_set(lambda x,y,z: z==lim, typ="node") for lim in [0, sim.mesh.l2]])
for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]

# Other simulation parameters (gravity, number of iterations, time step, ..):
sim.set_gravity([0,-9.81,0])
time = 10
nsteps = time//dt
sim.set_analysis_parameters(dt=dt,type="MPMExplicit3D", nsteps=nsteps, output_step_interval=nsteps/100)
sim.post_processing["vtk"] = ["stresses"]

# Save user defined parameters to be reused at the postprocessing stage:
sim.add_custom_parameters({"maxwell_particles": maxwell_particles,
							"walls": walls})
# Final generation of input files:
sim.write_input_file()

E = 9.5e9
nu = 0.3
density = 900
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

