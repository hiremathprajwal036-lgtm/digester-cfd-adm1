import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Digester geometry
V_total = 500        # m³ total volume
HRT_nominal = 20.0   # days

# Dead zone scenarios
scenarios = {
    'No Dead Zones (0%)':   0.00,
    'Low Dead Zones (15%)': 0.15,
    'Baseline (31.8%)':     0.318,
    'High Dead Zones (50%)':0.50,
    'Severe (65%)':         0.65,
}

# ADM1 parameters
mu_max = 0.4
Ks = 0.5
kd = 0.02
CH4_per_VS = 0.35
OLR = 2.5

def ch4_yield(HRT_eff):
    mu = mu_max * (1 / (1 + Ks / (OLR * HRT_eff))) - kd
    efficiency = max(0, min(1, mu * HRT_eff / (1 + mu * HRT_eff)))
    return OLR * CH4_per_VS * efficiency * HRT_eff

results = []
for name, dz_frac in scenarios.items():
    V_active = V_total * (1 - dz_frac)
    HRT_eff = HRT_nominal * (1 - dz_frac)
    CH4 = ch4_yield(HRT_eff)
    results.append({
        'Scenario': name,
        'Dead Zone (%)': dz_frac * 100,
        'Active Volume (m³)': V_active,
        'Effective HRT (days)': round(HRT_eff, 2),
        'CH4 Yield (m³/day)': round(CH4, 2)
    })

df = pd.DataFrame(results)
print(df.to_string(index=False))

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('CFD Dead Zone Impact on Digester Performance', fontsize=14)

colors = ['green', 'yellowgreen', 'orange', 'tomato', 'darkred']

axes[0].bar(df['Dead Zone (%)'], df['Active Volume (m³)'], color=colors, width=8)
axes[0].set_xlabel('Dead Zone Fraction (%)')
axes[0].set_ylabel('Active Volume (m³)')
axes[0].set_title('Active Volume vs Dead Zones')
axes[0].grid(True, axis='y')

axes[1].bar(df['Dead Zone (%)'], df['CH4 Yield (m³/day)'], color=colors, width=8)
axes[1].set_xlabel('Dead Zone Fraction (%)')
axes[1].set_ylabel('CH₄ Yield (m³/day)')
axes[1].set_title('CH₄ Yield vs Dead Zones')
axes[1].grid(True, axis='y')

plt.tight_layout()
plt.savefig('../results/cfd_dead_zones.png', dpi=150, bbox_inches='tight')
print("\nPlot saved to results/cfd_dead_zones.png")

df.to_csv('../results/cfd_dead_zones.csv', index=False)
print("CSV saved to results/cfd_dead_zones.csv")