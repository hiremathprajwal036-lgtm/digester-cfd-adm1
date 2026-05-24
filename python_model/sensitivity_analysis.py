import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Base parameters
HRT_base = 20.0
T_base = 37.0
OLR_base = 2.5
mu_max = 0.4
Ks = 0.5
kd = 0.02
CH4_per_VS = 0.35

def ch4_yield(HRT, T, OLR):
    T_factor = 1 + 0.04 * (T - 37)
    mu = mu_max * T_factor * (1 / (1 + Ks / (OLR * HRT))) - kd
    efficiency = max(0, min(1, mu * HRT / (1 + mu * HRT)))
    return OLR * CH4_per_VS * efficiency * HRT

# Sensitivity ranges
HRT_range = np.linspace(10, 30, 50)
T_range   = np.linspace(30, 42, 50)
OLR_range = np.linspace(1, 5, 50)

ch4_HRT = [ch4_yield(h, T_base, OLR_base) for h in HRT_range]
ch4_T   = [ch4_yield(HRT_base, t, OLR_base) for t in T_range]
ch4_OLR = [ch4_yield(HRT_base, T_base, o) for o in OLR_range]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('ADM1 Sensitivity Analysis – CH₄ Yield', fontsize=14)

axes[0].plot(HRT_range, ch4_HRT, color='steelblue', linewidth=2)
axes[0].set_xlabel('HRT (days)'); axes[0].set_ylabel('CH₄ yield (m³/day)')
axes[0].set_title('Effect of HRT'); axes[0].grid(True)

axes[1].plot(T_range, ch4_T, color='tomato', linewidth=2)
axes[1].set_xlabel('Temperature (°C)'); axes[1].set_ylabel('CH₄ yield (m³/day)')
axes[1].set_title('Effect of Temperature'); axes[1].grid(True)

axes[2].plot(OLR_range, ch4_OLR, color='seagreen', linewidth=2)
axes[2].set_xlabel('OLR (kg VS/m³/day)'); axes[2].set_ylabel('CH₄ yield (m³/day)')
axes[2].set_title('Effect of OLR'); axes[2].grid(True)

plt.tight_layout()
plt.savefig('../results/sensitivity_analysis.png', dpi=150, bbox_inches='tight')
print("Plot saved to results/sensitivity_analysis.png")

df = pd.DataFrame({
    'HRT_days': HRT_range, 'CH4_HRT': ch4_HRT,
    'Temp_C': T_range, 'CH4_Temp': ch4_T,
    'OLR': OLR_range, 'CH4_OLR': ch4_OLR
})
df.to_csv('../results/sensitivity_results.csv', index=False)
print("CSV saved to results/sensitivity_results.csv")