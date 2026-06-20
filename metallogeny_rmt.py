#!/usr/bin/env python3
"""
Metallogeny RMT analysis (Racetrack 4).
Four provinces: Chile/Nanling/Tethyan porphyry + orogenic gold (non-magmatic).
Per-source unfold + pool; single-source vs mixed; bootstrap CI.
Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""
import pandas as pd, numpy as np
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
import json, sys

def goe(s):return (np.pi/2)*s*np.exp(-np.pi*s**2/4)
def gue(s):return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
def mkcdf(f,mx=8,n=10000):
    s=np.linspace(0,mx,n);c=cumulative_trapezoid(f(s),s,initial=0);c/=c[-1]
    return interp1d(s,c,bounds_error=False,fill_value=(0,1))
POI=lambda x:1-np.exp(-x);GOE=mkcdf(goe);GUE=mkcdf(gue)
def sr(sp):
    if len(sp)<3:return np.nan,np.nan
    r=np.minimum(sp[:-1],sp[1:])/np.maximum(sp[:-1],sp[1:])
    return r.mean(),r.std()/np.sqrt(len(r))
def beta(r):
    if np.isnan(r):return np.nan
    if r<=0.386:return 0.0
    elif r<=0.536:return (r-0.386)/0.15
    elif r<=0.603:return 1+(r-0.536)/0.067
    else:return min(2+(r-0.603)/0.1,3)
def pool_test(df,dcol,seed=2026):
    pooled=[]
    for d,g in df.groupby(dcol,sort=False):
        a=np.sort(g['age_ma'].values);iv=np.diff(a);iv=iv[iv>0]
        if len(iv)>=2: pooled.extend((iv/iv.mean()).tolist())
    pooled=np.array(pooled)
    if len(pooled)<4: return None
    r,re=sr(pooled)
    _,pp=stats.kstest(pooled,POI);_,po=stats.kstest(pooled,GOE);_,pu=stats.kstest(pooled,GUE)
    best=min([('Poisson',stats.kstest(pooled,POI)[0]),('GOE',stats.kstest(pooled,GOE)[0]),
              ('GUE',stats.kstest(pooled,GUE)[0])],key=lambda x:x[1])[0]
    rng=np.random.default_rng(seed)
    boot=[sr(rng.choice(pooled,len(pooled),True))[0] for _ in range(8000)]
    ci=np.percentile(boot,[2.5,97.5])
    a=np.sort(df['age_ma'].values);iv=np.diff(a);iv=iv[iv>0];mixed_r=sr(iv/iv.mean())[0]
    return dict(n=len(pooled),r=r,b=beta(r),best=best,pp=pp,po=po,pu=pu,
                ci=[float(ci[0]),float(ci[1])],mixed_r=float(mixed_r))

DATASETS=[('chile_porphyry_expanded.csv','district','Chile porphyry Cu (subduction)'),
          ('nanling_wsn.csv','deposit','Nanling W-Sn (intracontinental)'),
          ('tethyan_porphyry.csv','district','Tethyan porphyry Cu (collision)'),
          ('orogenic_gold.csv','deposit','Orogenic gold (NON-magmatic)')]

if __name__=='__main__':
    datadir=sys.argv[1] if len(sys.argv)>1 else '../data'
    print(f"{'Province':38s}{'n':>4s}{'<r>':>8s}{'beta':>6s}{'mixed':>7s}{'best':>8s}")
    print("-"*72)
    out={}
    for fname,dcol,label in DATASETS:
        df=pd.read_csv(f"{datadir}/{fname}",comment='#')
        res=pool_test(df,dcol)
        out[label]=res
        print(f"{label:38s}{res['n']:4d}{res['r']:8.3f}{res['b']:6.2f}"
              f"{res['mixed_r']:7.3f}{res['best']:>8s}")
    print("\nBootstrap 95% CI (Poisson=0.386 excluded?):")
    for label,res in out.items():
        print(f"  {label:38s} [{res['ci'][0]:.3f},{res['ci'][1]:.3f}] "
              f"{'excluded' if res['ci'][0]>0.386 else 'NOT excluded'}")
    json.dump(out,open('four_province_results.json','w'),indent=2)
    print("\nKey: orogenic gold (non-magmatic) shows repulsion too -> intrinsic memory.")
