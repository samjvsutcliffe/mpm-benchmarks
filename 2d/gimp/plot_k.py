import pandas as pd
import matplotlib.pyplot as plt
for f in ["0325","045","049"]:
    df = pd.read_csv("terminus_position_{}.csv".format(f))
    plt.plot((df["Time (s)"]),df["Terminus displacement"],label=f)
plt.legend()
plt.show()
