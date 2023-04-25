import matplotlib.pyplot as plt
import pandas as pd
for f in ["1e3","1e6"]:
    df = pd.read_csv("height_{}.csv".format(f))
    plt.plot(df["Time (s)"],df["Height"],label=f)
plt.legend()
plt.show()

