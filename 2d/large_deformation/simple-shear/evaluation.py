import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import numpy as np
dt = 1e-2
data_dir = "./simple_shear/results/simple_shear/"
files = os.listdir(data_dir)
p = re.compile('.*\.h5') 
h5_files = [f for f in files if p.match(f)]

numbers = re.compile("\d*")
#element_counts = [[int(y) for y in x if y][0] for x in map(numbers.findall,h5_files)]
times = sorted([[int(y) for y in x if y][0] for x in map(numbers.findall,h5_files)])
h5_files = [x for _,x in sorted(zip(times,h5_files))]
times = sorted(times)

shear_stress = []
for f in h5_files:
    df = pd.read_hdf(data_dir+f)
    shear_stress.append(df["tau_xy"].mean())
time = np.array(times)
shear_stress = np.array(shear_stress)
plt.title("Shear stress over time")
plt.xlabel("Time (s)")
plt.ylabel("u_x (m)")
plt.plot(time,shear_stress)
plt.show()
