import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter
import pandas as pd

# -------- plot O2(a1) --------
df1 = pd.read_csv("output_H2O2HE.csv")
df1.columns = [col.strip() for col in df1.columns]
t1 = 1e3*df1['Time(s)'].values
ppm1 = 1e6*df1['O2(A1)'].values/df1['N_gas(#/cm^3)'].values
plt.plot(t1, ppm1, linewidth=1.5, label='$\mathrm{O_2({a^1\Delta_g})}$, ChemPlasKin', color='b')
plt.xlabel("Time (ms)", fontsize=14)
plt.ylabel("Mole Fraction (ppm)", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(0, 160)
plt.legend(fontsize=12, loc='upper left')
plt.tight_layout()
plt.show()