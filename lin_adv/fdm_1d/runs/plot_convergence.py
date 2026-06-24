import pandas as pd
import matplotlib.pyplot as plt

# Read CSV
df = pd.read_csv("convergence_results.csv")
import numpy as np

df.columns = df.columns.str.strip()

# Plot
plt.figure(figsize=(8,6))

plt.loglog(df["Np"], df["LELE_CD6"], "o-", label="LELE CD6")
plt.loglog(df["Np"], df["LELE_CD8"], "s-", label="LELE CD8")
plt.loglog(df["Np"], df["LELE_CD10"], "^-", label="LELE CD10")

Np = df["Np"].values

# Reference order lines
ref6  = df["LELE_CD6"].iloc[0]  * (Np[0]/Np)**6
ref8  = df["LELE_CD8"].iloc[0]  * (Np[0]/Np)**8
ref10 = df["LELE_CD10"].iloc[0] * (Np[0]/Np)**10

plt.loglog(Np, ref6,  'k--', linewidth=1.5, label='O(h^6)')
plt.loglog(Np, ref8,  'k-.', linewidth=1.5, label='O(h^8)')
plt.loglog(Np, ref10, 'k:',  linewidth=2,   label='O(h^10)')

plt.xlabel("Number of Grid Points (Np)")
plt.ylabel("L2 Error")
plt.ylim(1e-15, 1e-2)  #double precision limit
plt.title("Grid Convergence of Compat Schemes")
plt.grid(True, which="both")
plt.legend()

plt.savefig("convergence_plot.png", dpi=300)
plt.show()
