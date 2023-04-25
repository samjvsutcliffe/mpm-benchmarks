import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import numpy as np
dt = 1e-2 * 1e-3
data_dir = "./consol/results/consol/"
files = os.listdir(data_dir)
p = re.compile('.*\.h5') 
h5_files = [f for f in files if p.match(f)]
terminus = []
time = []
for f in h5_files:
    print(f)
    df = pd.read_hdf(data_dir+f)
    terminus.append(df["coord_y"].max())
    time.append(dt * float(re.findall("\d+",f)[0]))
time = np.array(time)
terminus = np.array(terminus)
terminus_displacement = terminus - terminus[0]
df = pd.DataFrame({"Time (s)":time,"Height":terminus_displacement})
df.to_csv("height_small.csv")
plt.title("Terminus evolution over time")
plt.xlabel("Time (h)")
plt.ylabel("u_x (m)")
plt.plot(time/(60**2),terminus_displacement)
print()
plt.show()
