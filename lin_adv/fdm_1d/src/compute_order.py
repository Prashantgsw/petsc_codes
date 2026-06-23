import pandas as pd
import numpy as np

df = pd.read_csv("convergence_results.csv")
df.columns = df.columns.str.strip()

for scheme in ["LELE_CD6", "LELE_CD8", "LELE_CD10"]:

    err = df[scheme].values

    print(f"\n{scheme}")

    for i in range(len(err)-1):
        p = np.log(err[i]/err[i+1]) / np.log(2)

        print(
            f"{int(df['Np'][i]):4d} -> "
            f"{int(df['Np'][i+1]):4d} : "
            f"p = {p:.4f}"
        )