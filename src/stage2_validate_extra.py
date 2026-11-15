"""Extra LOYO variants (0.5 deg): B2 = B1 x dominant land-cover ratio; M1 without the FRP term."""
import sys; sys.path.insert(0, __import__('os').path.dirname(__file__))
import numpy as np, pandas as pd
import harmonize as H
from harmonize import *
from stage2_validate import cell_oblast_map
d=load_det(); outage,_=outage_weeks(d); res=0.5
cells=np.unique(cell_id(d.lat,d.lon,res))
_design=H.design
def design_nofrp(g):
    X=_design(g); return np.delete(X, 31, axis=1)  # drop log1p_modis_frp
rows=[]
for sats in (("Terra","Aqua"),("Terra",)):
    m,v=build_table(d,res,TRAIN_YEARS,modis_sats=sats)
    g=prepare(full_grid(cells,TRAIN_YEARS,m,v),lc_shares(res),cell_oblast_map(d,res),outage); g=g[~g.outage].copy()
    print(g.assign(l=dominant_lc(g)).groupby('l').m_n.sum())
    p=loyo(g,TRAIN_YEARS,{"B2":B2})
    for k,pp in p.items(): rows.append(dict(metrics(pp),model=k,modis="+".join(sats)))
    H.design=design_nofrp; M1.start=None
    a=M1().fit(g).alpha
    p=loyo(g,TRAIN_YEARS,{"M1_noFRP":lambda: M1(alpha=a)})
    for k,pp in p.items(): rows.append(dict(metrics(pp),model=k,modis="+".join(sats)))
    H.design=_design; M1.start=None
t=pd.DataFrame(rows); print(t[["modis","model","ukr_month_WAPE","ukr_month_bias","ukr_year_WAPE","oblast_month_WAPE","cell_month_WAPE"]].round(3).to_string(index=False))
t.assign(grid=0.5).to_csv('data/processed/validation_loyo_extra.csv',index=False)
