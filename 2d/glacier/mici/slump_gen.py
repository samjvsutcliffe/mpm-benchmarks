import pycbg.preprocessing as utl
from pycbg.mesh import Mesh
import numpy as np
import math
dt = 1e-2

def create_LinearElasticDamage(self, pset_id=0,density=900,
        youngs_modulus=1e6,
        poisson_ratio=0.3,
        critical_stress=1e5,
        critical_damage=1,
        damage_rate=1
        ):
	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "LinearElasticDamage2D",
	"density": density,
	"youngs_modulus": youngs_modulus,
	"poisson_ratio": poisson_ratio,
	"critical_stress": critical_stress,
	"critical_damage": critical_damage,
        "local_length": 50,
        "local_length_damaged": 0.1,
        "degredation_function": "deviatoric",
        "stress_measure": "VM",
	"damage_rate": damage_rate
        })
def create_GlenDamage(self, pset_id=0,density=900,
        youngs_modulus=1e6,
        poisson_ratio=0.3,
        viscosity=1e6,
        viscous_power=1.0,
        critical_stress=1e5,
        critical_damage=1,
        damage_rate=1
        ):
	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "GlenDamage2D",
	"density": density,
	"youngs_modulus": youngs_modulus,
	"poisson_ratio": poisson_ratio,
	"viscosity": viscosity,
	"viscous_power":viscous_power,
	"critical_stress": critical_stress,
	"critical_damage": critical_damage,
        "local_length": 50,
        "local_length_damaged": 50,
	"damage_rate": damage_rate
        })

def create_Glen2D(self, pset_id=0,density=900,
        youngs_modulus=1e9,
        poisson_ratio=0.3,
        viscosity=1e6,
        viscous_power=1.0):

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
sim = utl.Simulation(title="slump")

resolution = 50
resolutions = [resolution,resolution ]

# Creating the mesh:


shelf_ratio = 4

shelf_length = 500
#particle_dims = (shelf_length,500.)
shelf_height = 500
datum = shelf_height * 0.5
shelf_length = shelf_height * shelf_ratio
particle_dims = (shelf_length,shelf_height)
domain_dims = (shelf_length+shelf_height,shelf_height + 100.)


#node_type = "ED2Q4"
node_type = "ED2GIMP"
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

sim.particles = utl.Particles(pmesh,mps_per_cell,directory=sim.directory,particle_type = "_DAMAGE")

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

density = 900
density_water = 1000
g = -9.8
height = sim.particles.positions[:,1]
p = density_water * np.maximum(datum-height,0) * g 
s_patch = np.tile(np.array([1,1,1,0,0,0]),(len(height),1))
print(np.transpose(np.tile(height,(6,1))))
print(np.multiply(s_patch,np.transpose(np.tile(p,(6,1)))))
sim.set_initial_particles_stresses(np.multiply(s_patch,np.transpose(np.tile(p,(6,1)))))
#sim.set_initial_particles_volumes(np.ones(len(sim.particles.positions)) * (np.power(resolution / mps_per_cell,2)))
sim.init_entity_sets()
sim.particles._filename = "particles.txt"

#sim.particles.type = "inject"
#sim.particles.nparticles_per_dir = mps_per_cell
#inflow_rate = 1
#sim.particles._particle_velocity = [inflow_rate,0]
#sim.particles._particle_duration = [0,2000]
#points = sim.mesh.nodes
#ids = []
#for i, p in enumerate(points):
#    if (lambda x,y: y>100 and y<200 and x<50)(*p): 
#        ids.append(i)
#
#sim.particles.cset_id = sim.entity_sets.create_set(lambda *args:any(arg in ids for arg in args), typ="cell")

sim.particles.write_file()



# Creating entity sets (the 2 materials), using lambda functions:
maxwell_particles = sim.entity_sets.create_set(lambda x,y: True, typ="particle")

E = 1e9
nu = 0.325
visc = 111e6
crit = 0.33e6
damage_rate = 1e-08
critical_damage = 0.60

# The materials properties:
#sim.materials.create_MohrCoulomb3D(pset_id=lower_particles)
#sim.materials.create_Newtonian3D(pset_id=maxwell_particles)
#create_Newtonian3D(sim.materials,pset_id=maxwell_particles)
#create_Maxwell3D(sim.materials,pset_id=maxwell_particles)
#create_NortonHoff2D(sim.materials,pset_id=maxwell_particles,
#        density=900,
#        youngs_modulus=E,
#        poisson_ratio=0.45,
#        viscosity=1e-30,viscous_power=3)

#sim.materials.create_Newtonian(pset_id=maxwell_particles)
#create_Glen2D(sim.materials,pset_id=maxwell_particles,
#       density=900,
#       youngs_modulus=E,
#       poisson_ratio=nu,
#       viscosity=visc,
#       viscous_power=3.0
#       )
#sim.materials.create_LinearElastic(pset_id=maxwell_particles,density=900,youngs_modulus=1e7,poisson_ratio=nu)

create_LinearElasticDamage(sim.materials,pset_id=maxwell_particles,density=density,youngs_modulus=E,poisson_ratio=nu,critical_stress=crit,damage_rate=damage_rate,critical_damage=critical_damage)
#create_GlenDamage(sim.materials,pset_id=maxwell_particles,
#        density=density,
#        youngs_modulus=E,
#        poisson_ratio=nu,
#        viscosity=visc,
#        viscous_power=3.0,
#        critical_stress=crit,
#        damage_rate=damage_rate,
#        critical_damage=critical_damage)
#create_Glen2D(sim.materials,pset_id=maxwell_particles,
#        density=900,
#        youngs_modulus=E,
#        poisson_ratio=nu,
#        viscosity=1e-35,
#        viscous_power=3)

# Boundary conditions on nodes entity sets (blocked displacements):
walls = []
walls.append([sim.entity_sets.create_set(lambda x,y: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
#walls.append([sim.entity_sets.create_set(lambda x,y: y==lim, typ="node") for lim in [100]])
#walls.append([sim.entity_sets.create_set(lambda x,y: z==lim, typ="node") for lim in [0, sim.mesh.l2]])
for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]
#sim.add_velocity_condition(0, 0., walls[1][0])
#sim.add_friction_condition(1,-1,0.5,walls[1][0])
sim.add_buoyancy_condition(datum,density_water,[0,8000,0,1000])

#sim.add_velocity_condition(1,0.0,sim.entity_sets.create_set(lambda x,y: x<=shelf_length and y==100, typ="node"))
#sim.add_friction_condition(1,-1,0.9,sim.entity_sets.create_set(lambda x,y: y==0 and x <=250, typ="node"))
#
#sim.add_force(1,9.81 * (density_water) * (resolution**2),sim.entity_sets.create_set(lambda x,y: x>=shelf_length and y<100, typ="node"))
#sim.add_velocity_condition(0,0.0,sim.entity_sets.create_set(lambda x,y: x==500 and y<100, typ="node"))

l = E * nu / ((1+nu)*(1-2*nu))
mu = E / (2*(1+nu))

bulk_stiffness = ((3 * l)+(2*mu))/3
min_x = min(resolutions[0:1])
speed_of_sound = math.sqrt(bulk_stiffness/density)

mass_scale = 1e4

pmod = (E / (1+nu)*(1-nu))
dt = min_x * np.sqrt(density/pmod) * np.sqrt(mass_scale)
dt *= 0.1
#dt = 1e-2
print("Pmod dt est {}".format(dt))
# Other simulation parameters (gravity, number of iterations, time step, ..):
sim.set_gravity([0,-9.81])
time = 10000
nsteps = time//dt
sim.set_analysis_parameters(dt=dt,type="MPMExplicit2D", nsteps=nsteps, 
        output_step_interval=nsteps/100,
        damping=0.00)
        

#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": 1e8}
#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": E*1e-4}
sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": 0.01}
#sim.analysis_params["damping"] = {"type": "Cundall", "damping_factor": 0.05}
sim.post_processing["vtk"] = ["stresses"]
sim.analysis_params["damage_enable"] = True
sim.analysis_params["damage_removal"] = False
sim.analysis_params["damage_nonlocal"] = True
sim.analysis_params["mass_scale"] = mass_scale
sim.post_processing["vtk"] = ["stresses","damage","ybar","displacements"]

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

