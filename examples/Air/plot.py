import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file into a DataFrame
df = pd.read_csv("outputAir.csv")
df.columns = [col.strip() for col in df.columns]

# List of species to plot
species = ["O", "N2(A3)", "N2(C3)"]

# Create the plot
fig, ax1 = plt.subplots()

# Plot each species
for specie in species:
    ax1.plot(1E9*df["Time(s)"], df[specie], label=specie)

# Add labels and title
ax1.set_xlabel("Time [ns]")
ax1.set_ylabel("Number Density [$\mathrm{\#/cm^3]}$")
ax1.set_yscale("log")  
# ax1.set_ylim(1E12, 1E19)
ax1.legend(loc='lower center')

ax2 = ax1.twinx()
ax2.plot(1E9*df["Time(s)"], df["T_gas(K)"], label="Temeprature", color='r')
ax2.set_ylabel("Temeprature [K]")
ax2.legend(loc='center right')

plt.grid(True)

# Show the plot
# plt.savefig("species_density.png", dpi=300, bbox_inches='tight')
plt.show()
