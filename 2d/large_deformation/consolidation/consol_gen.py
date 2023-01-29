import pycbg.preprocessing as utl
from pycbg.mesh import Mesh
import numpy as np
import math
dt = 1e-2
def in_bounds(x,y,cells):
    return (x>=0) and (y>=0) and (x<=cells[0]) and (y<=cells[1])
def id_from_pos(x,y,cells):
    return x + (y*(cells[0]+1))
def id_from_pos_cell(x,y,cells):
    return x + (y*(cells[0]))
gimp_order = [
        [0, 0],
        [1 ,0] ,
        [1 ,1]  ,
        [0, 1] ,
        [-1, 0],
        [2 ,0 ],
        [1 ,-1] ,
        [1 ,2 ] ,
        [2 ,1 ] ,
        [-1, 1] ,
        [0, 2 ],
        [0, -1],
        [-1, -1],
        [2 ,-1 ],
        [2 ,2  ],
        [-1, 2 ]
        ]
def gen_gimp_mesh(length,cells,outdir="consol",mesh=None):
    print(length)
    print(cells)
    nd = 2
    length = np.array(length)
    cells = np.array(cells)
    cellvs = cells+1
    resolution = length/(cells)
    verts = np.zeros((cellvs[0]*cellvs[1],2))
    for x in range(cellvs[0]):
        for y in range(cellvs[1]):
            verts[id_from_pos(x,y,cells),0] = x*resolution[0]
            verts[id_from_pos(x,y,cells),1] = y*resolution[1]
    elements = np.zeros(((cells[0])*(cells[1]),16))
    for e_x in range(cells[0]):
        for e_y in range(cells[1]):
            ni = 0
            ei = id_from_pos_cell(e_x,e_y,cells)
            #for dx in range(-1,3):
            #    for dy in range(-1,3):
            for ni in range(16):
                dx,dy = gimp_order[ni]
                #print("{} {}".format(e_x+dx,e_y+dy))
                if in_bounds(e_x+dx,e_y+dy,cells):
                    elements[ei,ni] = id_from_pos(e_x+dx,e_y+dy,cells)
                else:
                    elements[ei,ni] = -1
    print(verts)
    print(elements)
    with open("{}/mesh.txt".format(outdir),"w") as f:
        f.write("{}   {}\n".format(len(verts),len(elements)))
        for v in verts:
            f.write("{} {}\n".format(v[0],v[1]))
        for e in elements:
            for ni in e:
                f.write("{} ".format(int(ni)))
            f.write("\n")
    if mesh:
        mesh.cells = []
        for e in elements:
            mesh.cells.append(e.astype(int).tolist()[0:4])
        mesh.nodes = verts
    return elements




# The usual start of a PyCBG script:
sim = utl.Simulation(title="consol")

resolution = 10
resolutions = [resolution,resolution]
# Creating the mesh:

#shelf_length = resolution
length = 50
particle_dims   = (resolution,length)
domain_dims     = (resolution,length + (2*resolution))

#length = 1
#particle_dims = (resolution,resolution*2)
#domain_dims =   (resolution,resolution*3)


#node_type = "ED2Q4"
node_type = "ED2Q16G"
#node_type = "ED3H64G"
sim.create_mesh(dimensions=domain_dims, ncells=[x//r for x,r in zip(domain_dims,resolutions)],cell_type = node_type)
gen_gimp_mesh(domain_dims,[x//r for x,r in zip(domain_dims,resolutions)],"consol",sim.mesh)
pmesh = utl.Mesh(dimensions=particle_dims,origin=(0,0,0), ncells=[x//r for x,r in zip(particle_dims,resolutions)],cell_type=node_type)

mps_per_cell = 4
sim.particles = utl.Particles(pmesh,mps_per_cell,directory=sim.directory,particle_type="")

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

E = 1e7
nu = 0.0
density = 80
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

