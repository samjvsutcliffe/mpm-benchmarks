import pandas as pd
import matplotlib.pyplot as plt
for f in ["fs","ss"]:
    df = pd.read_csv("terminus_position_{}.csv".format(f))
    plt.plot((df["Time (s)"]),df["Terminus displacement"],label=f)
#E = 1e9
#df = pd.read_csv("terminus_position.csv")
#tau = 2e6/E
#plt.plot(df["Time (s)"],df["Terminus displacement"],label="2e6")
#df = pd.read_csv("terminus_position_1e6.csv")
#tau = 1e6/E
#plt.plot(df["Time (s)"],df["Terminus displacement"],label="1e6")
plt.legend()
plt.show()
