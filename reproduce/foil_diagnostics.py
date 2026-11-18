"""Run with the v0.3 NPZ, indexed-222 NPZ and upstream unweighted matrix as arguments."""
import hashlib
import json
import sys
from pathlib import Path

import contourpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm, Normalize
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import DisjointSet
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import lstsq
from scipy.ndimage import gaussian_filter
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

repo = Path(__file__).resolve().parents[1]
out = repo / "comparison"
config = json.loads((repo/"run25521/config/geometry.json").read_text())
source = json.loads((out/"atlas_sources.json").read_text())
index = pd.read_csv(repo/"indexing/cluster_index.csv")
versions = ["v0.3 clusters  /  upstream unweighted matrix", "Indexed-222  /  refitted unweighted matrix"]
colors = ["#a4522d", "#146679"]
foilnames = ["Foil 0  (+10 cm)", "Foil 1  (0 cm)", "Foil 2  (-10 cm)"]
delta_edges = [-25,-22,-20,-18,-16,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,22,25,30,35,40]
audit = {"matrix_checks": [], "cross_validation": [], "pages": []}
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def selected(f, eligible):
    n = min(np.count_nonzero(eligible & (f==i)) for i in range(3))
    return np.sort(np.concatenate([np.flatnonzero(eligible & (f==i))[:n] for i in range(3)]))

def compute(npz, matrix, side):
    ref = source[["left","right"][side]]
    assert digest(matrix) == ref["matrix_sha256"]
    if side == 0:
        assert digest(npz) == ref["source_sha256"]
    with np.load(npz, allow_pickle=False) as a:
        f,c,r,delta,bx,by,entry,obj = [a[k] for k in
            "foil_id hole_col hole_row delta_percent beam_x_cm beam_y_cm source_full_entry cluster_id".split()]
        parent = a["parent_before_fp_split_id"] if side else obj
        ok = ~a["index_review_required"] if side else np.ones(len(f),bool)
        fp = [a["fp_x_cm"]/100,a["fp_xp_rad"],a["fp_y_cm"]/100,a["fp_yp_rad"]]
    assert ok.sum() == [106058,105911][side] and len(entry)==106058 and (np.diff(entry)>0).all()
    z=(1-f)*10.
    theta=np.deg2rad(config["theta_deg"])
    xv=-by-config["mispoint_x_cm"]
    yv=-z*np.sin(theta)+bx*np.cos(theta)-config["mispoint_y_cm"]
    zv=z*np.cos(theta)+bx*np.sin(theta)
    length=config["sieve_z_cm"]
    hx=c*config["hole_pitch_cm"][0]; hy=r*config["hole_pitch_cm"][1]
    h1,h2=config["hb_delta_polynomial_cm"]; hb=(h1+h2*delta)*delta
    xp=(hx-xv)/(length-zv); yp=(hy-hb-yv)/(length-zv)
    xt=xv-xp*zv; yt=yv-yp*zv
    truth=np.column_stack([xp,yt/100,yp])
    rows=[t.split() for t in matrix.read_text().splitlines() if len(t.split())==5]
    coef=np.array([t[:3] for t in rows],float)
    exps=np.array([list(map(int,t[-1])) for t in rows])
    free=exps[:,-1]==0
    assert len(rows)==132 and free.sum()==126
    coords=fp+[xt/100]
    terms=np.column_stack([np.prod([v**k for v,k in zip(coords,e) if k],axis=0) if e.any()
                           else np.ones(len(f)) for e in exps])
    full=terms@coef
    inside=(delta>-15)&(delta<22)
    sel=selected(f,ok&inside)
    rmse=np.sqrt(np.mean(((full[sel]-truth[sel])*1000)**2,axis=0))
    assert len(sel)==ref["evaluated_events"] and np.allclose(rmse,ref["rmse_xptar_mrad_ytar_mm_yptar_mrad"],atol=1e-7)
    audit["matrix_checks"].append(dict(matrix_sha256=digest(matrix),data_sha256=digest(npz),selected=len(sel),
                                       indexed=len(f),used_for_residuals=int(ok.sum()),outside_fit_window=int((ok&~inside).sum()),rmse=rmse.tolist()))

    # Keep siblings and the same physical hole across all three foils in one fold.
    ds=DisjointSet()
    for col,row,par in np.unique(np.c_[c[ok],r[ok],parent[ok]],axis=0):
        h=(int(col),int(row)); p=int(par); ds.add(h); ds.add(p); ds.merge(h,p)
    labels={v:i for i,g in enumerate(ds.subsets()) for v in g}
    ii=np.flatnonzero(ok); groups=np.array([labels[key] for key in zip(c[ii],r[ii])])
    A=terms[:,free]; fixed=terms[:,~free]@coef[~free]
    held=np.full_like(full,np.nan)
    for fold,(train,test) in enumerate(GroupKFold(5).split(ii,groups=groups)):
        mask=np.zeros(len(f),bool); mask[ii[train]]=True
        tr=selected(f,mask&inside); te=ii[test]
        assert not set(groups[train]).intersection(groups[test]) and not np.intersect1d(tr,te).size
        design=A[tr]; target=(truth-fixed)[tr]
        # Direct unweighted SVD with the square-root equivalent of ROOT's normal-matrix threshold.
        beta,_,rank,_=lstsq(design,target,cond=np.sqrt(np.finfo(float).eps),lapack_driver='gelsd')
        held[te]=A[te]@beta+fixed[te]
        audit["cross_validation"].append(dict(side=side,fold=fold+1,training=len(tr),evaluated=len(te),rank=int(rank)))
    assert np.isfinite(held[ok]).all()
    def project(p):
        xp,yt,yp=p.T; yt=yt*100
        zz=(z+(yt+config["mispoint_y_cm"])*np.sin(theta))/(np.cos(theta)-yp*np.sin(theta))
        sx=xt+length*xp; sy=yt+length*yp+hb; fx=xt+xp*zz; fy=yt+yp*zz
        return dict(y=yt*10,dy=yt*10-truth[:,1]*1000,sx=sx,sy=sy,dx=(sx-hx)*10,ds=(sy-hy)*10,
                    fx=fx,fy=fy,dfx=(fx-xv)*10,dfy=(fy-yv)*10)
    ideal=project(truth)
    assert all(np.max(np.abs(ideal[k])) < 1e-8 for k in ('dy','dx','ds','dfx','dfy'))
    cv=project(held)
    return pd.DataFrame(dict(entry=entry,foil=f,col=c,row=r,delta=delta,obj=obj,parent=parent,ok=ok,inside=inside,
                 hx=hx,hy=hy,xv=xv,yv=yv,theory=truth[:,1]*1000,cvdy=cv['dy'],cvdx=cv['dx'],cvds=cv['ds'],**project(full)))

def empirical_boundary(points):
    # Smooth in locally whitened coordinates, then calibrate the level on actual member densities.
    med=np.median(points,axis=0)
    clamped=np.clip(points,*np.quantile(points,[.02,.98],axis=0))
    L=np.linalg.cholesky(np.cov(clamped.T)+np.eye(2)*.015**2)
    z=np.linalg.solve(L,(points-med).T).T
    bw=len(z)**(-1/6)
    lim=np.quantile(z,[.002,.998],axis=0)+np.array([[-1],[1]])*4*bw
    step=np.maximum((lim[1]-lim[0])/190,bw/5)
    edges=[np.arange(lim[0,j],lim[1,j]+step[j]*1.5,step[j]) for j in range(2)]
    h=np.histogram2d(*z.T,bins=edges)[0]
    sm=gaussian_filter(h,bw/step,mode="constant")
    centers=[(e[1:]+e[:-1])/2 for e in edges]
    density=RegularGridInterpolator(centers,sm,bounds_error=False,fill_value=0)(z)
    level=np.quantile(density,.1)
    assert level>0
    contours=contourpy.contour_generator(x=centers[0],y=centers[1],z=sm.T).lines(level)
    loops=[np.asarray(v)@L.T+med for v in contours if len(v)>3]
    inside=np.zeros(len(points),bool)
    for poly in loops:
        assert np.allclose(poly[0],poly[-1])
        inside ^= MplPath(poly).contains_points(points)
    assert loops and abs(inside.mean()-.9)<=max(.024,2/len(z))
    return loops,inside.sum()


plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'axes.labelcolor':'#344856','axes.titlecolor':'#183849','figure.facecolor':'white'})
def layout(title,sub,nr=3,nc=1,foot='',height=9):
    fig=plt.figure(figsize=(15.6,height))
    fig.text(.055,.955,title,size=19,weight='bold',color='#153647',va='top')
    fig.text(.055,.91,sub,size=10,color='#405869',va='top')
    gs=fig.add_gridspec(1,2,left=.066,right=.974,bottom=.092,top=.80,wspace=.17)
    panels=[]
    for i in range(2):
        fig.text([.066,.555][i],.845,versions[i],size=12,weight='bold',color=colors[i])
        panels.append(gs[i].subgridspec(nr,nc,hspace=.77 if nr>2 else .52,wspace=.38).subplots(squeeze=False))
    fig.text(.06,.022,foot or 'Full indexed populations: v0.3  106,058  |  Indexed-222  105,911 resolved events. All displayed statistics use unit event weights.',size=8.5,color='#405869')
    fig.text(.972,.022,f"Run 25521   |   {len(audit['pages'])+1}",ha='right',size=9)
    return fig,panels

def finish(fig,title):
    pdf.savefig(fig,dpi=150);plt.close(fig);audit['pages'].append(title)

def byfoil(side,foil):
    return usable[side][usable[side].foil==foil]

def edges(k,step):
    x=pd.concat([g[k] for g in usable])
    return np.arange(np.floor(x.min()/step)*step,(np.ceil(x.max()/step)+1)*step,step)
def sync(axs,axis):
    get,setter=('get_ylim','set_ylim') if axis=='y' else ('get_xlim','set_xlim')
    lo,hi=zip(*(getattr(ax,get)() for ax in axs))
    for ax in axs:getattr(ax,setter)(min(lo),max(hi))
def hist(ax,v,e,**kwargs):
    n=np.histogram(v,bins=e)[0]; assert n.sum()==len(v)
    ax.stairs(n,e,**kwargs)
def histpage(residual):
    title='Target Y residual by foil' if residual else 'Reconstructed and geometric target Y by foil'
    fig,pan=layout(title,'Event counts in 1 mm bins. Dashed vertical lines mark zero or the mean geometric target.')
    k='dy' if residual else 'y';e=edges(k,1)
    if not residual:
        et=edges('theory',1);e=np.arange(min(e[0],et[0]),max(e[-1],et[-1])+1)
    for f in range(3):
        for i in range(2):
            a=pan[i][f,0];g=byfoil(i,f)
            hist(a,g[k],e,color=colors[i],label='Reconstructed' if not residual else 'Residual')
            if not residual:hist(a,g.theory,e,color='#43545a',ls=':',label='Geometric target')
            a.axvline(0 if residual else g.theory.mean(),ls='--',color='.4',lw=.8)
            a.set(title=f'{foilnames[f]}   |   N = {len(g):,}',ylabel='Events / mm',
                  xlabel='Reconstructed minus geometric Ytar (mm)' if residual else 'Ytar (mm)')
            a.text(.98,.83,f'Mean {g[k].mean():+.2f} mm    SD {g[k].std(ddof=0):.2f} mm',transform=a.transAxes,ha='right',fontsize=8)
            if not residual and f==0:a.legend(frameon=False,loc='upper left',fontsize=8)
        sync([p[f,0] for p in pan],'y')
    finish(fig,title)

def maps(mode):
    cv=mode=='held';reco=mode=='reconstructed'
    title={'held':'Held-hole cross-validation of target Y','nominal':'Target Y bias and spread by assigned hole',
           'reconstructed':'Target Y bias and spread on the reconstructed sieve'}[mode]
    note=('Auxiliary five-fold unweighted refits of each indexing. Each physical hole and its sibling clusters stay in one fold.' if cv else
          'Each cell shows the mean residual or its population standard deviation. Blank cells have insufficient entries.')
    fig,pan=layout(title,note,3,2,foot=('Cross-validation uses auxiliary matrices; the other full-sample maps use the two named matrices. ' if cv else '')+
        ('0.5 cm bins; at least 10 events per cell. ' if reco else 'One cell per nominal hole. ')+'Shared color limits: bias -8 to +8 mm; SD 0 to 10 mm. Values beyond the limits saturate.')
    xe,ye=(np.arange(-16,16.01,.5),np.arange(-11,11.01,.5)) if reco else (np.arange(-5.5,6)*2.5,np.arange(-5.5,6)*1.64)
    for i in range(2):
        for f in range(3):
            g=byfoil(i,f); x,y=(g.sx,g.sy) if reco else (g.hx,g.hy);v=g.cvdy if cv else g.dy
            bins=(xe,ye);n=np.histogram2d(x,y,bins)[0];numer=np.histogram2d(x,y,bins,weights=v)[0]
            mean=numer/np.maximum(n,1);sd=np.sqrt(np.maximum(0,np.histogram2d(x,y,bins,weights=v*v)[0]/np.maximum(n,1)-mean**2))
            for j,(vals,label,cmap,l,h) in enumerate([(mean,'Bias','RdBu_r',-8,8),(sd,'SD','viridis',0,10)]):
                ax=pan[i][f,j];im=ax.pcolormesh(xe,ye,np.ma.masked_where(n<(10 if reco else 1),vals).T,cmap=cmap,vmin=l,vmax=h,rasterized=True)
                ax.set(title=f'{foilnames[f]}  |  {label}',xlabel='Sieve X (cm)',ylabel='Sieve Y (cm)')
                ax.set_facecolor('#eceff0');ax.set_aspect('equal')
                fig.colorbar(im,ax=ax,shrink=.80,label='mm',extend='both' if j==0 else 'max')
    finish(fig,title)

def trend(axis,k):
    delta=axis=='delta'
    unit={'dy':'target Y','dx':'sieve X','ds':'sieve Y'}[k];title=f'{unit.capitalize()} residual by '+({'row':'sieve row','col':'sieve column','delta':'momentum deviation'}[axis])
    fig,pan=layout(title,'Points show the mean. Bars span one population SD. The dashed line shows the mean for held-hole predictions.')
    for f in range(3):
        for i in range(2):
            g=byfoil(i,f);groups=g.groupby(pd.cut(g.delta,delta_edges,right=False),observed=True) if delta else g.groupby(axis)
            info=[(h.delta.mean() if delta else key,h[k].mean(),h[k].std(ddof=0),h['cv'+k].mean()) for key,h in groups if not delta or len(h)>=10]
            ax=pan[i][f,0];x,m,s,cv=np.array(info).T;ax.axhline(0,color='.55',lw=.7)
            ax.errorbar(x,m,yerr=s,fmt='o',color=colors[i],capsize=2.5,ms=3);ax.plot(x,cv,'--',lw=1,color='#48535c')
            ax.set(title=foilnames[f],ylabel='Residual (mm)',xlabel='Momentum deviation (%)' if delta else 'Sieve '+('column' if axis=='col' else 'row'))
            if delta:
                for v in [-15,22]:ax.axvline(v,color='.5',ls=':',lw=1)
            else:ax.set_xticks(range(-5,6))
        axs=[p[f,0] for p in pan];sync(axs,'y');sync(axs,'x')
    finish(fig,title)

def limits(vals,step):
    return np.arange(np.floor(min(v.min() for v in vals)/step)*step,(np.ceil(max(v.max() for v in vals)/step)+1)*step,step)
def heat(ax,v,xe,ye,norm):
    h=np.histogram2d(v[:,0],v[:,1],(xe,ye))[0]
    ax.pcolormesh(xe,ye,np.ma.masked_equal(h,0).T,norm=norm,cmap='GnBu',rasterized=True)
    return int(len(v)-h.sum())
def repro(mode):
    title,ks,xlabel,ylabel={'sieve':('Reprojection on the sieve',('sx','sy'),'Sieve X (cm)','Sieve Y (cm)'),
                            'foil':('Intersection with the foil plane',('fx','fy'),'Foil-plane X (cm)','Foil-plane Y (cm)'),
                            'beam':('Foil intersection relative to the event beam position',('dfx','dfy'),'X residual (mm)','Y residual (mm)')}[mode]
    caption={'sieve':'assigned nominal hole centers','foil':'the mean event beam position on each foil','beam':'zero displacement'}[mode]
    fig,pan=layout(title,f'Log color scale gives events per bin. Crosses mark {caption}.',height=11 if mode=='sieve' else 9)
    if mode=='sieve':xe,ye=np.linspace(-16,16,161),np.linspace(-11,11,111)
    else:
        xe,ye=[np.linspace(min(g[k].min() for g in usable)-1e-8,max(g[k].max() for g in usable)+1e-8,b) for k,b in zip(ks,[83,111])]
    vals={(i,f):byfoil(i,f)[list(ks)].to_numpy() for i in range(2) for f in range(3)}
    norm=LogNorm(1,max(np.histogram2d(*v.T,(xe,ye))[0].max() for v in vals.values()))
    for (i,f),v in vals.items():
        ax=pan[i][f,0];outside=heat(ax,v,xe,ye,norm);g=byfoil(i,f)
        if mode=='sieve':
            ax.set_aspect('equal')
            q=g[['hx','hy']].drop_duplicates();ax.scatter(q.hx,q.hy,marker='+',s=11,lw=.55,c='#ac511f')
        else:ax.plot(g.xv.mean() if mode=='foil' else 0,g.yv.mean() if mode=='foil' else 0,'+',color='#ac511f')
        ax.set(title=f'{foilnames[f]}  |  N = {len(v):,}  |  outside view = {outside}',xlabel=xlabel,ylabel=ylabel)
        fig.colorbar(ax.collections[0],ax=ax,pad=.015,label='Events / bin')
    finish(fig,title)

def individuals(f,axis):
    col=axis=='col';k='dx' if col else 'ds';name='column' if col else 'row'
    title=foilnames[f]+f'  |  Sieve residuals for each {name}'
    fig,pan=layout(title,'Horizontal axis: reconstructed minus nominal sieve coordinate (mm). Vertical axis: events per 0.5 mm.',3,4,
                   foot='Each panel reports the event count, mean and population SD in mm. Both sides use identical bins and axes for each corresponding index.')
    e=edges(k,.5)
    for n,idx in enumerate(range(-5,6)):
        for i in range(2):
            g=byfoil(i,f);v=g.loc[g[axis]==idx,k];ax=pan[i].flat[n]
            hist(ax,v,e,color=colors[i],lw=1.1);ax.axvline(0,color='.5',lw=.6);ax.tick_params(labelsize=7)
            ax.set_title(f'{name.title()} {idx:+d}  |  {len(v):,}',fontsize=8)
            if not len(v):
                ax.text(.5,.5,'No events',ha='center',transform=ax.transAxes,color='.45',fontsize=8)
                ax.set_yticks([])
            if len(v):ax.text(.03,.86,f'{v.mean():+.1f} / {v.std(ddof=0):.1f}',transform=ax.transAxes,fontsize=7)
            if n in [8,9,10]:ax.set_xlabel('Residual (mm)',fontsize=8)
        sync([p.flat[n] for p in pan],'y')
    for p in pan:p.flat[-1].axis('off')
    finish(fig,title)

def overlay():
    title='The three reconstructed target Y peaks'
    fig,pan=layout(title,'Each dashed line marks the mean geometric target Y for the corresponding foil.',1,1,height=6)
    e=edges('y',1)
    for i in range(2):
        ax=pan[i][0,0]
        for f,color in enumerate(['#1b7691','#b6632a','#496d45']):
            g=byfoil(i,f);hist(ax,g.y,e,color=color,label=foilnames[f]);ax.axvline(g.theory.mean(),color=color,ls='--',lw=1)
        ax.legend(frameon=False);ax.set(xlabel='Reconstructed Ytar (mm)',ylabel='Events / mm')
    sync([p[0,0] for p in pan],'y');finish(fig,title)

def edgecol():
    title='Foil 0, column -5  |  Momentum dependence of the sieve X residual'
    fig,pan=layout(title,'The fit window is -15% < delta < +22%. Histograms and trend points use all resolved events in this column.',2,1)
    data=[byfoil(i,0).query('col == -5') for i in range(2)]
    e=limits([g.dx for g in data],1)
    for i,g in enumerate(data):
        a,b=pan[i][:,0]
        for inside,color,label in [(True,colors[i],'In window'),(False,'#586574','Outside window')]:
            h=g[g.inside==inside];hist(a,h.dx,e,label=f'{label}  |  N = {len(h):,}',color=color)
            for _,s in h.groupby(pd.cut(h.delta,delta_edges,right=False),observed=True):
                b.errorbar(s.delta.mean(),s.dx.mean(),yerr=s.dx.std(ddof=0),color=color,fmt='o',ms=4,capsize=3)
        a.legend(frameon=False);a.set(xlabel='Sieve X residual (mm)',ylabel='Events / mm')
        b.axvline(-15,ls=':',color='.4');b.axhline(0,color='.5',lw=.8);b.set(xlabel='Momentum deviation (%)',ylabel='Mean +/- SD (mm)')
    for row in range(2):
        for axis in ['x','y']:sync([p[row,0] for p in pan],axis)
    finish(fig,title)

def clusters():
    results=[]
    for i,df in enumerate(data):
        names=index.set_index('cluster_id').cluster.to_dict() if i else dict(zip(index.previous_id,index.previous_cluster))
        items=[]
        for obj,g in df.groupby('obj'):
            pts=g[['sx','sy']].to_numpy()
            assert len(g[['foil','hx','hy']].drop_duplicates())==1
            loops,count=empirical_boundary(pts)
            items.append(dict(name=names[obj],parent=g.parent.iloc[0],foil=g.foil.iloc[0],med=np.median(pts,axis=0),
                              hole=g[['hx','hy']].iloc[0].to_numpy(),loops=loops,n=len(pts),inside=int(count),review=not g.ok.all()))
        assert sum(k['n'] for k in items)==106058 and len(items)==[219,222][i]
        results.append(items)
    return results

def marks(ax,items,labels=False):
    for k in items:
        x,y=k['med'];hx,hy=k['hole'];color='#ac3343' if k['review'] else plt.cm.tab10((int(k['name'][1:])*7)%10)
        for loop in k['loops']:ax.plot(*loop.T,color=color,lw=.85)
        ax.plot([x,hx],[y,hy],color='#6f463a',lw=.8)
        ax.scatter(hx,hy,marker='+',s=17,lw=.8,color='#28292c')
        ax.plot(x,y,'.',color=color,ms=3)
        if labels:
            ax.annotate(k['name']+('*' if k['review'] else ''),(x,y),xytext=(0,4),textcoords='offset points',
                        ha='center',fontsize=6.3,color='#182c3e',bbox=dict(facecolor='white',alpha=.63,edgecolor='none',pad=.25))
def clusterpage(mode,f=None):
    title={'density':'Sieve hit densities and cluster boundaries','vectors':'Cluster centers and their assigned holes',
           'labeled':f'Indexed cluster map  |  {foilnames[f] if f is not None else ""}'}[mode]
    fig,pan=layout(title,'Curves enclose approximately 90% of members. Dots are cluster medians; segments lead to assigned hole centers (+).',
                   1 if f is not None else 3,1,height=6.7 if f is not None else 11,
                   foot='All 106,058 events on both sides. Indexed-222 includes provisional C199 (147 events) here. Asterisk marks its label.')
    xe=np.linspace(-16,16,161);ye=np.linspace(-11,11,111)
    fs=range(3) if f is None else [f]
    norm=LogNorm(1,max(np.histogram2d(g.sx,g.sy,(xe,ye))[0].max() for d in data for _,g in d.groupby('foil')))
    for i in range(2):
        for j,foil in enumerate(fs):
            ax=pan[i][j,0];g=data[i].query('foil==@foil');items=[s for s in cluster_data[i] if s['foil']==foil]
            missing=int((~g.sx.between(-16,16)|~g.sy.between(-11,11)).sum())
            if mode=='density':heat(ax,g[['sx','sy']].to_numpy(),xe,ye,norm)
            marks(ax,items,mode=='labeled');ax.set(xlim=(-16,16),ylim=(-11,11),xlabel='Sieve X (cm)',ylabel='Sieve Y (cm)',
                title=f'{foilnames[foil]}  |  {len(g):,} events   {len(items)} clusters   {missing} outside view')
            ax.set_aspect('equal',adjustable='box')
            if mode=='density':fig.colorbar(ax.collections[0],ax=ax,pad=.02,label='Events / bin')
    finish(fig,title)

def splits():
    title='Three focal-plane cluster splits'
    fig,pan=layout(title,'The same member events are shown on each side, projected with their respective matrix and target assignment.',
                   foot='Each marker cloud uses the full parent membership. Lines connect cluster medians to assigned holes. Rows share the same X and Y limits across sides.')
    for j,parent in enumerate([500,800,2400]):
        clouds=[d[d.parent==parent] for d in data]
        assert np.array_equal(clouds[0].entry.values,clouds[1].entry.values)
        for i,g in enumerate(clouds):
            ax=pan[i][j,0];items=[s for s in cluster_data[i] if s['parent']==parent]
            for _,h in g.groupby('obj'):
                ax.scatter(h.sx,h.sy,s=3,alpha=.42,rasterized=True)
            marks(ax,items);ax.set(title='   +   '.join(k['name']+'  ('+str(k['n'])+')' for k in items),xlabel='Sieve X (cm)',ylabel='Sieve Y (cm)')
        for axis in ['x','y']:sync([p[j,0] for p in pan],axis)
    finish(fig,title)

if __name__=='__main__':
    left_npz,right_npz,left_matrix=map(Path,sys.argv[1:])
    right_matrix=repo/'run25521/unweighted/matrices/nps_shms_newfit_25521_unweighted.dat'
    with threadpool_limits(limits=2):
        data=[compute(left_npz,left_matrix,0),compute(right_npz,right_matrix,1)]
        assert np.array_equal(data[0].entry,data[1].entry) and np.array_equal(data[0].obj,data[1].parent)
        usable=[g[g.ok] for g in data]
        cluster_data=clusters()
    with PdfPages(out/'foil_diagnostics_comparison.pdf',metadata=dict(Title='Run 25521 | Foil and sieve diagnostics',Author='Run 25521 optics analysis')) as pdf:
        histpage(False);histpage(True)
        for mode in ['nominal','reconstructed','held']:maps(mode)
        for axis,k in [('row','ds'),('col','dx'),('row','dy'),('col','dy')]:trend(axis,k)
        for mode in ['sieve','foil','beam']:repro(mode)
        trend('delta','dy')
        for axis in ['row','col']:
            for f in range(3):individuals(f,axis)
        overlay();edgecol()
        clusterpage('density');splits();clusterpage('vectors')
        for f in range(3):clusterpage('labeled',f)
        assert pdf.get_pagecount()==len(audit['pages'])==27
    for i,g in enumerate(usable):
        m=g.groupby('foil').dy.agg(['size','mean',lambda x:x.std(ddof=0)])
        audit['matrix_checks'][i]['ytar_by_foil']=[dict(foil=int(f),events=int(n),bias_mm=b,sd_mm=s) for f,(n,b,s) in m.iterrows()]
    (out/'foil_diagnostics_index.json').write_bytes((json.dumps(audit,indent=2)+'\n').encode())
