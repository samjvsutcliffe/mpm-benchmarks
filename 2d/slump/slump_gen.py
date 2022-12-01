import pycbg.preprocessing as utl
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
resolutions = [20,20,20]

# Creating the mesh:

particle_dims = (500.,100.,100.)
domain_dims = (600.,200.,300)
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)])
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,100), ncells=[x//r for x,r in zip(particle_dims,resolutions)])
.
#Pseudo-2d
#sim.create_mesh(dimensions=domain_dims, ncells=(domain_dims[0]//resolution,domain_dims[1]//resolution,1))
#particle_dims = (500,100,1)
#domain_dims = (600.,200.,1.)
#pmesh = utl.Mesh(dimensions=particle_dims, ncells=(particle_dims[0]//resolution,particle_dims[1]//resolution,1))
# Creating Material Points, could have been done by filling an array manually:
#sim.create_particles(npart_perdim_percell=1)
sim.particles = utl.Particles(pmesh,1,directory=sim.directory)
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
nsteps = 1.5e4
sim.set_analysis_parameters(dt=dt,type="MPMExplicit3D", nsteps=nsteps, output_step_interval=nsteps/100)
sim.post_processing["vtk"] = ["stresses"]

# Save user defined parameters to be reused at the postprocessing stage:
sim.add_custom_parameters({"maxwell_particles": maxwell_particles,
							"walls": walls})
# Final generation of input files:
sim.write_input_file()
