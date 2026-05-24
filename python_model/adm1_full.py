import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_ivp

# ── ADM1 kinetic parameters ──────────────────────────────────────────────────
params = {
    'k_hyd':   0.10,   # hydrolysis rate constant        (1/day)
    'mu_ac':   0.40,   # max growth rate acetogens        (1/day)
    'mu_meth': 0.30,   # max growth rate methanogens      (1/day)
    'Ks_ac':   0.50,   # half-saturation acetogenesis     (g COD/L)
    'Ks_meth': 0.30,   # half-saturation methanogenesis   (g COD/L)
    'kd_ac':   0.02,   # decay rate acetogens              (1/day)
    'kd_meth': 0.015,  # decay rate methanogens            (1/day)
    'Y_ac':    0.04,   # yield acetogens
    'Y_meth':  0.03,   # yield methanogens
    'f_ch4':   0.35,   # CH4 fraction of COD degraded
}

# ── Initial conditions [S_part, S_ac, X_ac, X_meth, CH4_cum] ────────────────
y0 = [10.0, 0.5, 0.2, 0.1, 0.0]  # g COD/L

def adm1_ode(t, y, p):
    S_part, S_ac, X_ac, X_meth, CH4_cum = y

    # Hydrolysis
    r_hyd = p['k_hyd'] * S_part

    # Acetogenesis (VFA production)
    r_ac = p['mu_ac'] * (S_ac / (S_ac + p['Ks_ac'])) * X_ac

    # Methanogenesis
    r_meth = p['mu_meth'] * (S_ac / (S_ac + p['Ks_meth'])) * X_meth

    # ODEs
    dS_part = -r_hyd
    dS_ac   = r_hyd - r_ac / p['Y_ac']
    dX_ac   = (p['Y_ac'] * r_ac) - p['kd_ac'] * X_ac
    dX_meth = (p['Y_meth'] * r_meth) - p['kd_meth'] * X_meth
    dCH4    = p['f_ch4'] * r_meth

    return [dS_part, dS_ac, dX_ac, dX_meth, dCH4]

# ── Solve ODE ────────────────────────────────────────────────────────────────
t_span = (0, 60)
t_eval = np.linspace(0, 60, 300)
sol = solve_ivp(adm1_ode, t_span, y0, t_eval=t_eval,
                args=(params,), method='RK45', rtol=1e-6)

t = sol.t
S_part, S_ac, X_ac, X_meth, CH4_cum = sol.y

# ── Print summary ─────────────────────────────────────────────────────────────
print(f"Final particulate substrate : {S_part[-1]:.3f} g COD/L")
print(f"Final VFA (acetate)         : {S_ac[-1]:.3f} g COD/L")
print(f"Final acetogen biomass      : {X_ac[-1]:.3f} g COD/L")
print(f"Final methanogen biomass    : {X_meth[-1]:.3f} g COD/L")
print(f"Cumulative CH4 produced     : {CH4_cum[-1]:.3f} g COD/L")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle('Full ADM1 – Biological Process Dynamics (60 days)', fontsize=14)

axes[0,0].plot(t, S_part, color='steelblue', linewidth=2)
axes[0,0].set_title('Particulate Substrate'); axes[0,0].set_ylabel('g COD/L')
axes[0,0].set_xlabel('Time (days)'); axes[0,0].grid(True)

axes[0,1].plot(t, S_ac, color='tomato', linewidth=2)
axes[0,1].set_title('VFA / Acetate'); axes[0,1].set_ylabel('g COD/L')
axes[0,1].set_xlabel('Time (days)'); axes[0,1].grid(True)

axes[1,0].plot(t, X_ac, color='seagreen', label='Acetogens', linewidth=2)
axes[1,0].plot(t, X_meth, color='purple', label='Methanogens', linewidth=2)
axes[1,0].set_title('Microbial Biomass'); axes[1,0].set_ylabel('g COD/L')
axes[1,0].set_xlabel('Time (days)'); axes[1,0].legend(); axes[1,0].grid(True)

axes[1,1].plot(t, CH4_cum, color='darkorange', linewidth=2)
axes[1,1].set_title('Cumulative CH₄ Production'); axes[1,1].set_ylabel('g COD/L')
axes[1,1].set_xlabel('Time (days)'); axes[1,1].grid(True)

plt.tight_layout()
plt.savefig('../results/adm1_full.png', dpi=150, bbox_inches='tight')
print("Plot saved to results/adm1_full.png")

# ── Save CSV ──────────────────────────────────────────────────────────────────
df = pd.DataFrame({
    'Time_days': t, 'S_particulate': S_part, 'S_acetate': S_ac,
    'X_acetogens': X_ac, 'X_methanogens': X_meth, 'CH4_cumulative': CH4_cum
})
df.to_csv('../results/adm1_full.csv', index=False)
print("CSV saved to results/adm1_full.csv")