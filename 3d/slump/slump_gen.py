import pycbg.preprocessing as utl
from pycbg.mesh import Mesh
import numpy as np
import math
dt = 1e-2

def sdf_circle(xi,circle_pos,radius):
    xi = np.array(xi)
    circle_pos = np.array(circle_pos)
    diff = xi - circle_pos
    dist = np.power(diff,2).sum(axis=-1)
    return np.sqrt(dist) - radius

def length(x):
    return np.sqrt(np.sum(x*x,axis=-1))

def sdf_rect(xi,pos,size):
    xi = np.array(xi)
    pos = np.array(pos)
    size = np.array(size)
    diff = xi - pos
    q = np.abs(diff) - size
    return length(np.maximum(q,0)) + min(0,np.max(q))

def remove_sdf(sim,sdf):
    ids = sdf(sim.particles.positions) > 0
    print(ids)
    sim.particles.positions = np.array(sim.particles.positions[ids])

def create_Glen3D(self, pset_id=0,density=900,
        youngs_modulus=1e9,
        poisson_ratio=0.3,
        viscosity=1e6,
        viscous_power=1.0,
        critical_stress=1e5,
        damage_rate=1e-9
        ):

	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "Glen3D",
	"density": density,
	"youngs_modulus": youngs_modulus,
	"poisson_ratio": poisson_ratio,
	"viscosity": viscosity,
	"viscous_power":viscous_power,
	"critical_stress": critical_stress,
	"damage_rate": damage_rate,
	"local_length": 50
        })

def create_LinearElasticDamage(self, pset_id=0,density=900,
        youngs_modulus=1e6,
        poisson_ratio=0.3,
        critical_stress=1e5,
        damage_rate=1
        ):
	self.pset_ids.append(pset_id)
	self.materials.append({"id": len(self.materials),
	"type": "LinearElasticDamage3D",
	"density": density,
	"youngs_modulus": youngs_modulus,
	"poisson_ratio": poisson_ratio,
	"critical_stress": critical_stress,
	"damage_rate": damage_rate,
	"local_length": 50
        })
	
# The usual start of a PyCBG script:
sim = utl.Simulation(title="slump")

resolution = 20
resolutions = [resolution,resolution,resolution ]

# Creating the mesh:


shelf_length = 500
particle_dims = (shelf_length,100.,200)
domain_dims = (2000.,200.,200)


#node_type = "ED2GIMP"
node_type = "ED3H8"
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)],cell_type = node_type)
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,0), ncells=[x//r for x,r in zip(particle_dims,resolutions)],cell_type =node_type)

mps_per_cell = 2

sim.particles = utl.Particles(pmesh,mps_per_cell,directory=sim.directory,particle_type = "_DAMAGE")

#print(type(sim.particles.positions))
remove_sdf(sim,lambda xi: sdf_rect(xi,np.array([250,100,100]),np.array([20,20,100])))
#remove_sdf(sim,lambda xi: sdf_circle(xi,np.array([250,100,100]),10))
#print(sim.particles.positions)

sim.init_entity_sets()

sim.particles._filename = "particles.txt"
sim.particles.write_file()

maxwell_particles = sim.entity_sets.create_set(lambda x,y,z: True, typ="particle")


E = 1e8
nu = 0.325
visc = 11e6
density = 900
density_water = 999
crit = 0.33e6
damage_rate = 1e-7

#create_LinearElasticDamage(sim.materials,
#        pset_id=maxwell_particles,
#        density=density,
#        youngs_modulus=E,
#        poisson_ratio=nu,
#        critical_stress=crit,
#        damage_rate=damage_rate)

create_Glen3D(sim.materials,pset_id=maxwell_particles,
       density=900,
       youngs_modulus=E,
       poisson_ratio=nu,
       viscosity=visc,
       viscous_power=3.0,
       critical_stress=crit,
       damage_rate=damage_rate
       )

# Boundary conditions on nodes entity sets (blocked displacements):
walls = []
walls.append([sim.entity_sets.create_set(lambda x,y,z: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
walls.append([sim.entity_sets.create_set(lambda x,y,z: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
walls.append([sim.entity_sets.create_set(lambda x,y,z: z==lim, typ="node") for lim in [0, sim.mesh.l2]])

for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]

sim.set_gravity([0,-9.81,0])
time = 100
nsteps = time//dt
sim.set_analysis_parameters(dt=dt,type="MPMExplicit3D", nsteps=nsteps, 
        output_step_interval=nsteps/100,
        damping=0.00)
        

#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": 1e8}
#sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": E*1e-4}
sim.analysis_params["damping"] = {"type": "Viscous", "damping_factor": E*1e-3}
#sim.analysis_params["damping"] = {"type": "Cundall", "damping_factor": 0.05}
sim.post_processing["vtk"] = ["stresses","damage"]

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

