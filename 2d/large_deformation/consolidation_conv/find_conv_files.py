import matplotlib.pyplot as plt
import pandas as pd
import os
import shutil
import re
import numpy as np
for mesh_res in [2**x for x in range(1,13)]:
    # The usual start of a PyCBG script:
    sim_name = "consol_conv_{}".format(mesh_res)
    data_dir = "./{}/results/{}/".format(sim_name,sim_name)
    files = os.listdir(data_dir)
    p = re.compile('.*\.h5') 
    h5_files = [f for f in files if p.match(f)]
    terminus = []
    time = []
    shutil.copyfile(data_dir+h5_files[-1], "conv_files/conv_{}.h5".format(mesh_res))
    #for f in h5_files:
    #    df = pd.read_hdf(data_dir+f)
    #df = pd.read_hdf(data_dir+h5_files[-1])
    #df.to_hdf("conv_files/conv_{}".format(mesh_res))
    #plt.title("Terminus evolution over time")
    #plt.xlabel("Time (h)")
    #plt.ylabel("u_x (m)")
    #plt.plot(time/(60**2),terminus_displacement)
    #plt.show()
