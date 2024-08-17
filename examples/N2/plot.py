import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file into a DataFrame
df = pd.read_csv("outputN2.csv")
df.columns = [col.strip() for col in df.columns]

# List of species to plot
# species = ["e", "N", "N2(A)", "N2(B)", "N2(a')", "N2(C)", "N^+", "N2^+", "N3^+", "N4^+"]
species = ["N", "N2(A)", "N2(C)", "N2^+", "N4^+"]

# Create the plot
plt.figure(figsize=(10, 6))

# Plot each species
for specie in species:
    plt.plot(1000*df["Time(s)"], df[specie], label=specie)

# Add labels and title
plt.xlabel("Time [ms]")
plt.ylabel("Number Density [$\mathrm{\#/cm^3]}$")
plt.title("Species Number Density vs Time")
plt.yscale("log")  
plt.ylim(1E4, 1E14)
plt.legend(loc='lower right')
plt.grid(True)

# Show the plot
plt.savefig("species_density.png", dpi=300, bbox_inches='tight')
plt.show()
