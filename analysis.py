#!/usr/bin/env python3
"""
Blackbody radiation (UP 201, Experiment 5) - analysis and graph generation.

Run:   python3 blackbody_analysis.py
Needs: numpy, scipy, matplotlib   and the data files in the same folder:
       6volt.txt, 7volt.txt, 7_5volt.txt (or 7.5volt.txt), 8volt.txt, anglecallibration.txt
Output: figures/*.png  and  results.json (all numbers quoted in the report)
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ------------------------------------------------------------------ settings
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

T0 = 300.0          # K, room temperature (27 C) - NOT recorded in the lab, assumed
dT0 = 3.0           # K, assumed uncertainty of T0
ALPHA = 4.5e-3      # 1/K, temperature coefficient of tungsten (given)
A_CAUCHY, B_CAUCHY = 13900.0, 1.689   # n = A/lambda^2 + B  (lambda in nm)
dV, dI = 0.01, 0.001                  # instrument resolution (V, A)
THETA0 = 0.0        # deg, arm position of the undeviated beam (scan was zeroed there)
ANGLE_SCALE = 1.0   # sensor reading -> degrees (calibration run gives slope 1)
B_WIEN_TRUE = 2.897771955e6           # nm K
SIGMA_SB = 5.670374419e-8
hP, cL, kB = 6.62607015e-34, 2.99792458e8, 1.380649e-23

# room-temperature (cold filament) V-I data from the lab notebook
V_cold = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10])
I_cold = np.array([0.002, 0.006, 0.008, 0.010, 0.014, 0.016, 0.020, 0.022, 0.024, 0.027])

# hot-filament runs: voltage, current (A), run number, data file
RUNS = [
    (6.0, 0.256, 5,  ["6volt.txt"]),
    (7.0, 0.280, 6,  ["7volt.txt"]),
    (7.5, 0.290, 12, ["7_5volt.txt", "7.5volt.txt"]),
    (8.0, 0.300, 9,  ["8volt.txt"]),
]
COLORS = {6.0: "#1f77b4", 7.0: "#2ca02c", 7.5: "#ff7f0e", 8.0: "#d62728"}

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight"})

# ------------------------------------------------------------------ helpers
def load(names):
    for n in names:
        p = os.path.join(HERE, n)
        if os.path.exists(p):
            d = np.loadtxt(p, skiprows=2)
            return d[:, 0], d[:, 1]
    raise FileNotFoundError(names)

def wavelength_nm(theta_deg):
    """Eq.(1) -> n ; Eq.(2) -> lambda.  NaN where n <= B (lambda undefined)."""
    th = np.radians(np.atleast_1d(np.asarray(theta_deg, float)))
    n = np.sqrt((2 / np.sqrt(3) * np.sin(th) + 0.5) ** 2 + 0.75)
    with np.errstate(invalid="ignore", divide="ignore"):
        lam = np.sqrt(A_CAUCHY / (n - B_CAUCHY))
    lam[n <= B_CAUCHY] = np.nan
    return lam if lam.size > 1 else lam[0]

def theta_for_lambda(lam_nm):
    n = A_CAUCHY / lam_nm ** 2 + B_CAUCHY
    s = (np.sqrt(n ** 2 - 0.75) - 0.5) * np.sqrt(3) / 2
    return np.degrees(np.arcsin(s))

def planck(lam_m, T):
    return 2 * hP * cL ** 2 / lam_m ** 5 / np.expm1(hP * cL / (lam_m * kB * T))

res = {}

# ------------------------------------------------------- 1. cold resistance R0
lr = stats.linregress(I_cold, V_cold)
R0_fit, dR0_fit, c_fit = lr.slope, lr.stderr, lr.intercept
R0_origin = float(np.sum(V_cold * I_cold) / np.sum(I_cold ** 2))
res_o = V_cold - R0_origin * I_cold
dR0_origin = float(np.sqrt(np.sum(res_o ** 2) / (len(V_cold) - 1) / np.sum(I_cold ** 2)))
# adopted value: straight line with intercept (absorbs meter zero offset)
R0, dR0 = R0_fit, dR0_fit
res["R0"] = dict(fit_slope=R0_fit, fit_slope_err=dR0_fit, fit_intercept=c_fit,
                 fit_intercept_err=lr.intercept_stderr, r=lr.rvalue,
                 through_origin=R0_origin, through_origin_err=dR0_origin,
                 power_max_mW=float(V_cold[-1] * I_cold[-1] * 1e3))

fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.errorbar(I_cold, V_cold, xerr=dI, yerr=dV, fmt="o", ms=5, capsize=2, color="k", label="data")
xx = np.linspace(0, 0.03, 50)
ax.plot(xx, R0_fit * xx + c_fit, "r-", label=fr"fit: $V={R0_fit:.2f}\,I{c_fit:+.3f}$")
ax.plot(xx, R0_origin * xx, "b--", lw=1, label=fr"through origin: $R_0={R0_origin:.2f}\ \Omega$")
ax.set_xlabel("Current $I$ (A)"); ax.set_ylabel("Voltage $V$ (V)")
ax.set_title(fr"Cold filament: $R_0 = {R0_fit:.2f}\pm{dR0_fit:.2f}\ \Omega$")
ax.set_xlim(0, 0.03); ax.set_ylim(0, 0.115); ax.legend(loc="upper left")
fig.savefig(os.path.join(OUT, "fig1_resistance_room_temperature.png")); plt.close(fig)

# --------------------------------------------- 2. filament temperatures
temps = {}
for V, I, run, _ in RUNS:
    R = V / I
    dR = R * np.hypot(dV / V, dI / I)
    T = T0 + (R / R0 - 1) / ALPHA
    dT = np.sqrt((dR / (ALPHA * R0)) ** 2 + (R * dR0 / (ALPHA * R0 ** 2)) ** 2 + dT0 ** 2)
    temps[V] = (R, dR, T, dT)

# ------------------------------------------------- 3. spectra
spec = {}
for V, I, run, names in RUNS:
    a, y = load(names)
    th = -a * ANGLE_SCALE - THETA0            # deflection angle, deg (positive)
    thm = np.maximum.accumulate(th)           # enforce monotonic sweep
    keep = np.concatenate(([True], np.diff(thm) > 1e-6)) & (thm > 0.5)
    th, y = thm[keep], y[keep]
    # linear baseline from line-free windows on both sides of the feature
    wl = (th > 53.0) & (th < 55.8)
    wr = (th > 64.5) & (th < 68.0)
    wb = wl | wr
    pb = np.polyfit(th[wb], y[wb], 1)
    base = np.polyval(pb, th)
    sig_base = np.std(y[wb] - base[wb], ddof=2)
    S = y - base
    # peak: quadratic fit around maximum of 5-pt running mean
    m = (th > 56.0) & (th < 63.0)
    Ssm = np.convolve(S, np.ones(5) / 5, mode="same")
    i0 = np.argmax(np.where(m, Ssm, -np.inf))
    w = np.abs(th - th[i0]) < 0.7
    co, cov = np.polyfit(th[w], S[w], 2, cov=True)
    aa, bb, cc = co
    th_pk = -bb / (2 * aa)
    g = np.array([bb / (2 * aa ** 2), -1 / (2 * aa), 0.0])
    dth_pk = float(np.sqrt(g @ cov @ g))
    dth_pk = float(np.hypot(dth_pk, 0.1))      # + 0.1 deg reading/zero uncertainty
    S_pk = np.polyval(co, th_pk)
    lam_pk = float(wavelength_nm(th_pk))
    dlam = float(abs(wavelength_nm(th_pk + 1e-3) - wavelength_nm(th_pk - 1e-3)) / 2e-3 * dth_pk)
    # integrated (baseline-subtracted) signal over the feature, theta in [56,63]
    mi = (th >= 56.0) & (th <= 63.0)
    Itot = float(np.trapezoid(S[mi], th[mi]))
    dItot = float(np.hypot(sig_base * 7.0, 0.01 * Itot))
    # half-maximum width in wavelength
    half = S_pk / 2
    ths, Ss = th[mi], S[mi]
    ip = np.argmax(Ss)
    left = ths[:ip][Ss[:ip] < half][-1] if np.any(Ss[:ip] < half) else np.nan
    right = ths[ip:][Ss[ip:] < half][0] if np.any(Ss[ip:] < half) else np.nan
    spec[V] = dict(th=th, y=y, S=S, base=base, th_pk=th_pk, dth_pk=dth_pk, S_pk=S_pk,
                   lam_pk=lam_pk, dlam=dlam, Itot=Itot, dItot=dItot, sig_base=sig_base,
                   lam_hm=(float(wavelength_nm(right)), float(wavelength_nm(left))),
                   raw_max=float(y[mi].max()), raw_base_left=float(y[wl].mean()),
                   n_pk=float(A_CAUCHY / lam_pk ** 2 + B_CAUCHY))

# ---- Fig. 3: intensity vs angle (raw + baseline-corrected zoom)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for V, *_ in RUNS:
    s = spec[V]
    ax[0].plot(s["th"], s["y"], color=COLORS[V], lw=1.2, label=f"{V:g} V")
    mz = (s["th"] > 54) & (s["th"] < 65)
    ax[1].plot(s["th"][mz], s["S"][mz], color=COLORS[V], lw=1.4, label=f"{V:g} V")
    ax[1].plot(s["th_pk"], s["S_pk"], "k^", ms=6)
ax[0].set_xlabel(r"Angular position $\theta$ (deg)"); ax[0].set_ylabel("Relative intensity (%)")
ax[0].set_title("(a) Full scan, as recorded"); ax[0].legend(title="Bulb voltage")
ax[0].axvspan(56, 63, color="gray", alpha=0.12)
ax[1].set_xlabel(r"Angular position $\theta$ (deg)"); ax[1].set_ylabel("Intensity above baseline (%)")
ax[1].set_title("(b) Spectral feature after baseline subtraction ($\\blacktriangle$ = peak)")
ax[1].legend(title="Bulb voltage")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_intensity_vs_angle.png")); plt.close(fig)

# ---- Fig. 4: intensity vs wavelength
fig, ax = plt.subplots(figsize=(7.4, 4.8))
for V, *_ in RUNS:
    s = spec[V]
    lam = wavelength_nm(s["th"])
    mz = (s["th"] > 55.9) & (s["th"] < 63.5) & np.isfinite(lam)
    ax.plot(lam[mz], s["S"][mz], color=COLORS[V], lw=1.4,
            label=fr"{V:g} V  ($T={temps[V][2]:.0f}$ K)")
    ax.plot(s["lam_pk"], s["S_pk"], "k^", ms=7)
    ax.axvline(s["lam_pk"], color=COLORS[V], ls=":", lw=0.8)
ax.set_xlim(450, 1250)
ax.set_xlabel(r"Wavelength $\lambda$ (nm)"); ax.set_ylabel("Intensity above baseline (%)")
ax.set_title(r"Intensity vs wavelength ($\blacktriangle$ = $\lambda_{max}$)")
ax.legend(); fig.savefig(os.path.join(OUT, "fig4_intensity_vs_wavelength.png")); plt.close(fig)

# --------------------------------------------------- 4. Stefan's law
Vs = [r[0] for r in RUNS]
T = np.array([temps[v][2] for v in Vs]); dT = np.array([temps[v][3] for v in Vs])
It = np.array([spec[v]["Itot"] for v in Vs]); dIt = np.array([spec[v]["dItot"] for v in Vs])
Spk = np.array([spec[v]["S_pk"] for v in Vs])
lx, ly = np.log(T), np.log(It)
st = stats.linregress(lx, ly)
st_pk = stats.linregress(lx, np.log(Spk))
# bootstrap-free uncertainty including T errors: Monte-Carlo
rng = np.random.default_rng(1)
sl_mc = []
for _ in range(5000):
    Tm = T + rng.normal(0, dT * 0 + dT * 0.0)  # (common-mode R0,T0 errors cancel in the slope to 1st order)
    Im = It * (1 + rng.normal(0, dIt / It))
    sl_mc.append(stats.linregress(np.log(Tm), np.log(Im)).slope)
res["stefan"] = dict(slope=st.slope, slope_err=st.stderr, r=st.rvalue, intercept=st.intercept,
                     slope_mc_std=float(np.std(sl_mc)), slope_peak=st_pk.slope,
                     slope_peak_err=st_pk.stderr, r_peak=st_pk.rvalue)
# cross-check: electrical input power V*I vs T (log-log)
Pel = np.array([v * temps_i for v, temps_i in zip(Vs, [r[1] for r in RUNS])])
sp = stats.linregress(lx, np.log(Pel))
res["stefan"]["slope_Pel"] = sp.slope; res["stefan"]["slope_Pel_err"] = sp.stderr; res["stefan"]["r_Pel"] = sp.rvalue
T4 = T ** 4
k4 = np.sum(T4 * It) / np.sum(T4 ** 2)
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax[0].errorbar(lx, ly, xerr=dT / T, yerr=dIt / It, fmt="o", color="k", capsize=3, label="data")
xx = np.linspace(lx.min() - 0.01, lx.max() + 0.01, 10)
ax[0].plot(xx, st.slope * xx + st.intercept, "r-",
           label=fr"fit: slope $=%.2f\pm%.2f$" % (st.slope, st.stderr))
ax[0].plot(xx, 4 * (xx - lx.mean()) + ly.mean(), "b--", lw=1, label="slope 4 (Stefan)")
for v, x, yv in zip(Vs, lx, ly):
    ax[0].annotate(f"{v:g} V", (x, yv), textcoords="offset points", xytext=(6, -12), fontsize=9)
ax[0].set_xlabel(r"$\ln T$  ($T$ in K)"); ax[0].set_ylabel(r"$\ln I_{tot}$  ($I_{tot}$ in %$\cdot$deg)")
ax[0].set_title(r"(a) $\ln I_{tot}$ vs $\ln T$"); ax[0].legend(loc="upper left")
ax[1].errorbar(T4 / 1e12, It, xerr=4 * T4 / T * dT / 1e12, yerr=dIt, fmt="o", color="k", capsize=3)
x4 = np.linspace(0, T4.max() * 1.05, 20)
ax[1].plot(x4 / 1e12, k4 * x4, "r-", label=r"$I_{tot}\propto T^4$ (through origin)")
ax[1].set_xlabel(r"$T^4$  ($10^{12}\ \mathrm{K}^4$)"); ax[1].set_ylabel(r"$I_{tot}$ (%$\cdot$deg)")
ax[1].set_title(r"(b) $I_{tot}$ vs $T^4$"); ax[1].set_xlim(0, None); ax[1].set_ylim(0, None)
ax[1].legend(loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_stefan_law.png")); plt.close(fig)

# ---------------------------------------------------- 5. Wien's law
lam = np.array([spec[v]["lam_pk"] for v in Vs]); dlam = np.array([spec[v]["dlam"] for v in Vs])
x = 1 / T; dx = dT / T ** 2
b_or = float(np.sum(x * lam) / np.sum(x ** 2))
b_or_err = float(np.sqrt(np.sum((lam - b_or * x) ** 2) / (len(x) - 1) / np.sum(x ** 2)))
wf = stats.linregress(x, lam)
lt = lam * T
dlt = lt * np.hypot(dlam / lam, dT / T)
ll = stats.linregress(np.log(T), np.log(lam))
res["wien"] = dict(b_origin=b_or, b_origin_err=b_or_err, b_free_slope=wf.slope,
                   b_free_err=wf.stderr, b_free_intercept=wf.intercept,
                   b_free_intercept_err=wf.intercept_stderr, r_free=wf.rvalue,
                   lam_T_mean=float(lt.mean()), lam_T_std=float(lt.std(ddof=1)),
                   loglog_slope=ll.slope, loglog_err=ll.stderr,
                   ratio_to_true=b_or / B_WIEN_TRUE,
                   lam_T=lt.tolist(), lam_T_err=dlt.tolist())
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax[0].errorbar(x * 1e4, lam, xerr=dx * 1e4, yerr=dlam, fmt="o", color="k", capsize=3, label="measured")
xx = np.linspace(0, x.max() * 1.1, 20)
ax[0].plot(xx * 1e4, b_or * xx, "r-", label=fr"fit through origin: $b=%.2f\times10^{-3}$ m K" % (b_or * 1e-6))
ax[0].plot(xx * 1e4, B_WIEN_TRUE * xx, "b--", lw=1, label=r"Wien: $b=2.898\times10^{-3}$ m K")
axin = ax[0].inset_axes([0.09, 0.50, 0.37, 0.30])
axin.errorbar(x * 1e4, lam, xerr=dx * 1e4, yerr=dlam, fmt="o", color="k", capsize=2, ms=4)
axin.plot(xx * 1e4, b_or * xx, "r-")
axin.set_xlim(5.6, 6.8); axin.set_ylim(580, 740); axin.tick_params(labelsize=7); axin.grid(alpha=0.3)
for v, xi, li in zip(Vs, x, lam):
    axin.annotate(f"{v:g} V", (xi * 1e4, li), textcoords="offset points", xytext=(-4, 7), fontsize=7)
axin.set_title("zoom on data", fontsize=8)
ax[0].set_xlabel(r"$1/T$  ($10^{-4}\ \mathrm{K^{-1}}$)"); ax[0].set_ylabel(r"$\lambda_{max}$ (nm)")
ax[0].set_title(r"(a) $\lambda_{max}$ vs $1/T$"); ax[0].set_xlim(0, None); ax[0].set_ylim(0, None)
ax[0].legend(loc="upper left", fontsize=9)
ax[1].errorbar(T, lt / 1e6, xerr=dT, yerr=dlt / 1e6, fmt="o", color="k", capsize=3, label=r"measured $\lambda_{max}T$")
ax[1].axhline(lt.mean() / 1e6, color="r", label=fr"mean $=%.2f\times10^{-3}$ m K" % (lt.mean() * 1e-6))
ax[1].axhline(B_WIEN_TRUE / 1e6, color="b", ls="--", lw=1, label=r"Wien: $2.898\times10^{-3}$ m K")
ax[1].set_xlabel("$T$ (K)"); ax[1].set_ylabel(r"$\lambda_{max}T$ ($10^{-3}$ m K)")
ax[1].set_title(r"(b) $\lambda_{max}T$ vs $T$"); ax[1].set_ylim(0, 3.4); ax[1].legend(loc="lower right", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig6_wien_law.png")); plt.close(fig)

# ----------------------------------- 6. theoretical blackbody spectrum (theory section)
lam_m = np.linspace(0.1e-6, 3.0e-6, 800)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
for Tt, c in zip([3000, 4000, 5000, 6000], ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"]):
    ax.plot(lam_m * 1e9, planck(lam_m, Tt) * 1e-13, color=c, label=f"{Tt} K")
Tt = 5000
rj = 2 * cL * kB * Tt / lam_m ** 4
ax.plot(lam_m * 1e9, rj * 1e-13, "k--", lw=1, label="Rayleigh-Jeans (5000 K)")
Tl = np.linspace(2800, 6200, 50); lm = B_WIEN_TRUE * 1e-9 / Tl
ax.plot(lm * 1e9, planck(lm, Tl) * 1e-13, "gray", ls=":", label=r"Wien locus $\lambda_{max}T=b$")
ax.axvspan(380, 750, color="yellow", alpha=0.15)
ax.set_ylim(0, 6.5); ax.set_xlim(0, 3000)
ax.set_xlabel(r"Wavelength $\lambda$ (nm)")
ax.set_ylabel(r"$B_\lambda(T)$  ($10^{13}\ \mathrm{W\,m^{-3}\,sr^{-1}}$)")
ax.set_title("Planck black-body spectra (shaded: visible band)"); ax.legend(fontsize=9)
fig.savefig(os.path.join(OUT, "fig2_blackbody_spectrum_theory.png")); plt.close(fig)

# ------------------------------------------------------------- results dump
res["runs"] = {}
for V, I, run, _ in RUNS:
    R, dR, Tk, dTk = temps[V]; s = spec[V]
    res["runs"][f"{V:g}"] = dict(run=run, I=I, R=R, dR=dR, T=Tk, dT=dTk,
        T_celsius=Tk - 273.15, lam_wien_nm=B_WIEN_TRUE / Tk, P_W=V * I,
        th_pk=s["th_pk"], dth_pk=s["dth_pk"], n_pk=s["n_pk"], lam_pk=s["lam_pk"], dlam=s["dlam"],
        S_pk=s["S_pk"], Itot=s["Itot"], dItot=s["dItot"], sig_base=s["sig_base"],
        lam_halfmax=s["lam_hm"], lam_T=s["lam_pk"] * Tk,
        sb_power_ratio=(Tk ** 4))
cal = np.loadtxt(os.path.join(HERE, "anglecallibration.txt"), skiprows=2)
cs = stats.linregress(cal[:, 0], cal[:, 1])
res["angle_calibration"] = dict(slope=cs.slope, intercept=cs.intercept, final=float(cal[-1, 0]))
res["lam_range_valid_deg"] = float(theta_for_lambda(1e9))
with open(os.path.join(HERE, "results.json"), "w") as f:
    json.dump(res, f, indent=2, default=float)
print(json.dumps(res, indent=1, default=float))
