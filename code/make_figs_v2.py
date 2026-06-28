#!/usr/bin/env python3
"""
make_figs_v2.py -- regenerate the two v2 figures for the erratum.

fig1_scalar_clock.pdf
  (a) master <r>-vs-CV curve for a ZERO-repulsion i.i.d. gamma renewal
      clock, with Poisson/GOE/GUE reference lines and the four provinces
      sitting ON the curve (their <r> is predicted by CV alone).
  (b) observed <r> vs SHUFFLE <r> per province (obs ~ shuf => marginal /
      scalar, not sequential), with the closed-data-null obs-shuf band.

fig2_age_error_superposition.pdf
  (a) age-error Monte Carlo: pooled <r> under age jitter per province,
      with Poisson/GOE/GUE lines and the v1 reported value marked.
  (b) the RETAINED real result: single (isolated) vs mixed (superposed)
      <r> -> superposition drives the local statistic to Poisson.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os, sys

RNG = np.random.default_rng(20260627)
POISSON_R, GOE_R, GUE_R = 0.386, 0.531, 0.603
DATADIR = sys.argv[1] if len(sys.argv) > 1 else '../data'
OUT = sys.argv[2] if len(sys.argv) > 2 else '../paper/figs'
os.makedirs(OUT, exist_ok=True)

C_OBS = '#b2182b'; C_SHUF = '#2166ac'; C_GAMMA = '#4d4d4d'
C_SINGLE = '#762a83'; C_MIXED = '#1b7837'

DATASETS = [
    ('chile_porphyry_expanded.csv', 'district', 'Andean\nporphyry Cu', 0.601),
    ('nanling_wsn.csv',             'deposit',  'Nanling\nW-Sn',        0.574),
    ('tethyan_porphyry.csv',        'district', 'Tethyan\nporphyry Cu', 0.712),
    ('orogenic_gold.csv',           'deposit',  'Orogenic gold\n(non-magmatic)', 0.678),
]

def r_stat(sp):
    sp = np.asarray(sp, float)
    if len(sp) < 3: return np.nan
    return float(np.mean(np.minimum(sp[:-1], sp[1:]) / np.maximum(sp[:-1], sp[1:])))

def shuffle_r(sp, n=3000):
    sp = np.asarray(sp, float)
    return float(np.mean([r_stat(RNG.permutation(sp)) for _ in range(n)]))

def v1_pool(df, dcol):
    pooled = []
    for _, g in df.groupby(dcol, sort=False):
        a = np.sort(g['age_ma'].values.astype(float)); iv = np.diff(a); iv = iv[iv > 0]
        if len(iv) >= 2: pooled.extend((iv / iv.mean()).tolist())
    return np.array(pooled)

def mixed_r(df):
    a = np.sort(df['age_ma'].values.astype(float)); iv = np.diff(a); iv = iv[iv > 0]
    return r_stat(iv / iv.mean())

def age_error_mc(df, dcol, n=4000):
    err = df['age_err'].values.astype(float) if 'age_err' in df.columns else np.full(len(df), 5.0)
    err = np.where(np.isfinite(err) & (err > 0), err, 5.0)
    ages = df['age_ma'].values.astype(float); grp = df[dcol].values
    rs = []
    for _ in range(n):
        d2 = pd.DataFrame({dcol: grp, 'age_ma': ages + RNG.normal(0, err)})
        p = v1_pool(d2, dcol)
        if len(p) >= 3: rs.append(r_stat(p))
    return np.array(rs)

# ---- gather data ----
labels, cvs, robs, rshuf, rmixed, mc, v1r = [], [], [], [], [], [], []
for fname, dcol, lab, v1 in DATASETS:
    df = pd.read_csv(os.path.join(DATADIR, fname), comment='#')
    p = v1_pool(df, dcol)
    labels.append(lab); v1r.append(v1)
    cvs.append(np.std(p)/np.mean(p)); robs.append(r_stat(p)); rshuf.append(shuffle_r(p))
    rmixed.append(mixed_r(df)); mc.append(age_error_mc(df, dcol))

# ---- master scalar-clock curve <r>(CV) for i.i.d. gamma, with finite-N band ----
cv_grid = np.linspace(0.18, 1.05, 32)
r_curve, r_lo, r_hi = [], [], []
Nrep = int(np.median([len(v1_pool(pd.read_csv(os.path.join(DATADIR, f), comment='#'), d))
                      for f, d, _, _ in DATASETS]))  # typical pooled N (~30)
for cv in cv_grid:
    k = 1.0/cv**2
    asymp = np.mean([r_stat(RNG.gamma(k, 1.0, size=40000)) for _ in range(3)])
    finite = np.array([r_stat(RNG.gamma(k, 1.0, size=Nrep)) for _ in range(2000)])
    r_curve.append(asymp); r_lo.append(np.percentile(finite, 2.5)); r_hi.append(np.percentile(finite, 97.5))
r_curve = np.array(r_curve); r_lo = np.array(r_lo); r_hi = np.array(r_hi)

# ============================ FIGURE 1 ============================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

# (a) master curve
axA.fill_between(cv_grid, r_lo, r_hi, color=C_GAMMA, alpha=0.16,
                 label=f'finite-N 95% band (N$\\approx${Nrep})')
axA.plot(cv_grid, r_curve, '-', color=C_GAMMA, lw=2.2,
         label='i.i.d. gamma renewal clock\n(zero repulsion)')
for ref, name, ls in [(POISSON_R, 'Poisson 0.386', ':'),
                      (GOE_R, 'GOE 0.531', '--'),
                      (GUE_R, 'GUE 0.603', '-.')]:
    axA.axhline(ref, color='grey', ls=ls, lw=1, alpha=0.8)
    axA.text(1.045, ref+0.004, name, ha='right', va='bottom', fontsize=8, color='grey')
for lab, cv, r in zip(labels, cvs, robs):
    axA.plot(cv, r, 'o', color=C_OBS, ms=9, mec='k', mew=0.6, zorder=5)
    axA.annotate(lab.replace('\n', ' '), (cv, r), textcoords='offset points',
                 xytext=(7, -2), fontsize=7.5)
axA.set_xlabel('CV of within-deposit spacings  (sigma/mu)')
axA.set_ylabel(r'spacing ratio  $\langle r\rangle$')
axA.set_title('(a)  $\\langle r\\rangle$ is set by CV alone', fontsize=11, loc='left')
axA.legend(loc='lower left', fontsize=8, framealpha=0.9)
axA.set_xlim(0.15, 1.07); axA.set_ylim(0.33, 0.78)

# (b) obs vs shuffle
x = np.arange(len(labels)); w = 0.36
axB.bar(x-w/2, robs, w, color=C_OBS, label=r'observed $\langle r\rangle$', edgecolor='k', lw=0.5)
axB.bar(x+w/2, rshuf, w, color=C_SHUF, label=r'shuffled $\langle r\rangle$', edgecolor='k', lw=0.5)
axB.axhline(POISSON_R, color='grey', ls=':', lw=1)
axB.text(len(labels)-0.5, POISSON_R+0.005, 'Poisson', fontsize=8, color='grey', ha='right')
for xi, ro, rsf in zip(x, robs, rshuf):
    axB.text(xi, max(ro, rsf)+0.012, f'{ro-rsf:+.3f}', ha='center', fontsize=7.5)
axB.set_xticks(x); axB.set_xticklabels([l.replace('\n',' ') for l in labels], fontsize=7.5, rotation=12)
axB.set_ylabel(r'$\langle r\rangle$')
axB.set_title('(b)  observed $\\approx$ shuffled  (marginal, not sequential)', fontsize=11, loc='left')
axB.legend(loc='lower right', fontsize=8); axB.set_ylim(0, 0.82)
axB.text(0.02, 0.96, r'$\langle r\rangle_{\rm obs}-\langle r\rangle_{\rm shuf}$ printed above bars',
         transform=axB.transAxes, fontsize=7.5, va='top', style='italic', color='#444')

fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig1_scalar_clock.pdf'), bbox_inches='tight')
print('wrote fig1_scalar_clock.pdf')

# ============================ FIGURE 2 ============================
fig2, (axC, axD) = plt.subplots(1, 2, figsize=(11, 4.4))

# (c) age-error MC
positions = np.arange(len(labels))
bp = axC.boxplot(mc, positions=positions, widths=0.55, patch_artist=True,
                 showfliers=False, medianprops=dict(color='k', lw=1.4),
                 whiskerprops=dict(color='#888'), capprops=dict(color='#888'))
for patch in bp['boxes']:
    patch.set_facecolor('#fddbc7'); patch.set_edgecolor('#b2182b'); patch.set_alpha(0.9)
for ref, name, ls in [(POISSON_R, 'Poisson', ':'), (GOE_R, 'GOE', '--'), (GUE_R, 'GUE', '-.')]:
    axC.axhline(ref, color='grey', ls=ls, lw=1)
    axC.text(len(labels)-0.45, ref+0.004, name, fontsize=8, color='grey', ha='right')
axC.plot(positions, v1r, '*', color='#1a1a1a', ms=15, zorder=6, label='v1 reported value')
axC.set_xticks(positions); axC.set_xticklabels([l for l in labels], fontsize=7.5)
axC.set_ylabel(r'pooled $\langle r\rangle$ under age jitter')
axC.set_title('(c)  age-error Monte Carlo: $\\langle r\\rangle$ unresolvable', fontsize=11, loc='left')
axC.legend(loc='lower left', fontsize=8); axC.set_ylim(0.30, 0.80)

# (d) superposition -> Poisson (retained real result)
w = 0.36
axD.bar(positions-w/2, robs, w, color=C_SINGLE, label='single source\n(isolated deposit)', edgecolor='k', lw=0.5)
axD.bar(positions+w/2, rmixed, w, color=C_MIXED, label='mixed\n(superposed)', edgecolor='k', lw=0.5)
axD.axhline(POISSON_R, color='grey', ls=':', lw=1.2)
axD.text(len(labels)-0.5, POISSON_R+0.006, 'Poisson 0.386', fontsize=8, color='grey', ha='right')
axD.set_xticks(positions); axD.set_xticklabels([l for l in labels], fontsize=7.5)
axD.set_ylabel(r'$\langle r\rangle$')
axD.set_title('(d)  superposition $\\to$ Poisson  (retained, real)', fontsize=11, loc='left')
axD.legend(loc='upper right', fontsize=7.5, ncol=1); axD.set_ylim(0, 0.82)

fig2.tight_layout()
fig2.savefig(os.path.join(OUT, 'fig2_age_error_superposition.pdf'), bbox_inches='tight')
print('wrote fig2_age_error_superposition.pdf')
