import pandas as pd
import matplotlib.pyplot as plt
for f in ["1e6","2e6","4e6"]:
    visc = float(f)
    E = 1e9
    tau = E/visc
    df = pd.read_csv("terminus_position_{}.csv".format(f))
    plt.plot((0.1*df["Time (s)"]),df["Terminus displacement"],label="nu={}".format(f))
#E = 1e9
#df = pd.read_csv("terminus_position.csv")
#tau = 2e6/E
#plt.plot(df["Time (s)"],df["Terminus displacement"],label="2e6")
#df = pd.read_csv("terminus_position_1e6.csv")
#tau = 1e6/E
#plt.plot(df["Time (s)"],df["Terminus displacement"],label="1e6")
plt.xlabel("Time (s)")
plt.ylabel("Terminus displacement (m)")
plt.legend()
plt.show()
