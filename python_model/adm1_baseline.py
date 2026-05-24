import sys
sys.path.insert(0, '../PyADM1')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Reactor parameters ──────────────────────────────────────────
HRT_nominal = 20.0      # days - design HRT
f_dead      = 0.318     # dead zone fraction from CFD (OpenFOAM result)
HRT_eff     = HRT_nominal * (1 - f_dead)  # = 13.64 days

V_reactor   = 1000.0    # m³ - reactor volume
T           = 308.15    # K  - mesophilic (35°C)

# ── Substrate parameters (cattle manure, Batstone et al. 2002) ──
S_in = 10.0    # kg COD/m³ - influent substrate concentration
mu_max = 0.40  # 1/day - maximum specific growth rate
K_s    = 0.50  # kg COD/m³ - half-saturation constant
Y      = 0.040 # kg VSS/kg COD - yield coefficient
k_d    = 0.015 # 1/day - decay rate
k_m    = 0.35  # m³ CH4/kg COD - methane yield per COD removed

# ── Simple Monod ODE model ───────────────────────────────────────
def digester(t, y, HRT, S_in):
    S, X = y                          # substrate, biomass
    Q_in = V_reactor / HRT            # m³/day influent flow
    mu   = mu_max * S / (K_s + S)     # Monod growth rate
    dS   = (S_in - S) / HRT - mu * X / Y
    dX   = (mu - k_d) * X - X / HRT
    return [dS, dX]

# ── Initial conditions ───────────────────────────────────────────
y0   = [S_in * 0.5, 1.0]   # initial substrate and biomass
tspan = (0, 200)            # simulate 200 days to reach steady state
t_eval = np.linspace(0, 200, 1000)

# ── Run both scenarios ───────────────────────────────────────────
sol_nom = solve_ivp(digester, tspan, y0, t_eval=t_eval,
                    args=(HRT_nominal, S_in), method='RK45', rtol=1e-6)
sol_cfd = solve_ivp(digester, tspan, y0, t_eval=t_eval,
                    args=(HRT_eff, S_in),     method='RK45', rtol=1e-6)

# ── Methane production rate at steady state ──────────────────────
def ch4_rate(sol, HRT):
    S_ss = sol.y[0, -1]              # steady-state substrate
    removal = (S_in - S_ss) / S_in  # COD removal efficiency
    Q_in    = V_reactor / HRT
    return k_m * Q_in * (S_in - S_ss)  # m³ CH4/day

ch4_nom = ch4_rate(sol_nom, HRT_nominal)
ch4_cfd = ch4_rate(sol_cfd, HRT_eff)
ch4_loss = (ch4_nom - ch4_cfd) / ch4_nom * 100

print(f"Nominal HRT        : {HRT_nominal:.1f} days")
print(f"CFD-corrected HRT  : {HRT_eff:.2f} days")
print(f"Dead zone fraction : {f_dead*100:.1f}%")
print(f"CH4 nominal        : {ch4_nom:.1f} m³/day")
print(f"CH4 with dead zones: {ch4_cfd:.1f} m³/day")
print(f"Yield loss due to mixing: {ch4_loss:.1f}%")

# ── Plot ─────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(sol_nom.t, sol_nom.y[0], 'b-', label=f'Nominal HRT={HRT_nominal}d', linewidth=2)
ax1.plot(sol_cfd.t, sol_cfd.y[0], 'r--', label=f'CFD HRT={HRT_eff:.1f}d', linewidth=2)
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Substrate S (kg COD/m³)')
ax1.set_title('Substrate Degradation')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.bar(['Nominal HRT\n(20 days)', 'CFD-corrected\n(13.6 days)'],
        [ch4_nom, ch4_cfd], color=['steelblue', 'tomato'], width=0.5)
ax2.set_ylabel('CH₄ Production (m³/day)')
ax2.set_title('Impact of Dead Zones on Methane Yield')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('../results/adm1_baseline.png', dpi=150, bbox_inches='tight')
print("Plot saved to results/adm1_baseline.png")