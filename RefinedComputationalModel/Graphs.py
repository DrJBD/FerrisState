import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Simulated example metabolite pool time series data ---
# Simulating 1 hour at 10 sec intervals
timesteps = np.arange(0, 3601, 10)  # seconds
n_steps = len(timesteps)

# Simulated metabolite concentrations (example values in µM)
formate_series = np.log1p(timesteps) * 5  # example curve
ser_series = np.linspace(2000, 1900, n_steps)  # decreasing serine
thf_series = np.linspace(21, 19, n_steps)      # slight drop
nadpp_series = np.linspace(20, 10, n_steps)    # depletion
nadph_series = np.linspace(200, 300, n_steps)  # increasing
gly_series = formate_series * 1.9              # arbitrary correlation

# --- Data table ---
metabolite_data = pd.DataFrame({
    "Time (s)": timesteps,
    "Formate (µM)": formate_series,
    "Serine (µM)": ser_series,
    "THF (µM)": thf_series,
    "NADP⁺ (µM)": nadpp_series,
    "NADPH (µM)": nadph_series,
    "Glycine (µM)": gly_series
})

# --- Graph 1: Formate vs Time ---
plt.figure(figsize=(8, 5))
plt.plot(timesteps / 60, formate_series, label='Formate', linewidth=2)
plt.xlabel("Time (min)")
plt.ylabel("Formate Concentration (µM)")
plt.title("Formate Accumulation Over Time")
plt.grid(True)
plt.tight_layout()
plt.savefig("/mnt/data/formate_vs_time.png")

# --- Graph 2: Serine Dose-Response on Final Formate ---
ser_doses = [500, 1000, 1500, 2000, 2500]
final_formate = [50, 125, 200, 275, 300]  # Simulated

plt.figure(figsize=(7, 4.5))
plt.plot(ser_doses, final_formate, marker='o', linestyle='--', color='purple')
plt.xlabel("Initial Serine Concentration (µM)")
plt.ylabel("Final Formate Produced (µM)")
plt.title("Effect of Serine Input on Formate Yield")
plt.grid(True)
plt.tight_layout()
plt.savefig("/mnt/data/serine_vs_formate.png")

# Display the data table for download
import caas_jupyter_tools as caas
caas.display_dataframe_to_user(name="Metabolite Concentration Over Time", dataframe=metabolite_data)

("/mnt/data/formate_vs_time.png", "/mnt/data/serine_vs_formate.png")
