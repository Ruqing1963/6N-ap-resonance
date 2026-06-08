import csv, numpy as np, matplotlib, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})


def rd(name):
    with open(os.path.join(_DATA,name)) as f: rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) if _isnum(r[k]) else r[k] for r in rows],dtype=object) for k in rows[0]}
def _isnum(s):
    try: float(s); return True
    except: return False

# spectrum
sp=rd("ap_spectrum_k3.csv"); D=sp["Delta"].astype(int); N3=sp["N3"].astype(float)
# resonance table
import collections
res=list(csv.DictReader(open(os.path.join(_DATA,"ap_resonance.csv"))))

# ============ FIGURE 1 : diffraction spectrum ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) stem spectrum 1..300
sel=D<=300
ax[0].vlines(D[sel],0,N3[sel]/1000,color="#1f77b4",lw=1.0)
for d,lab in [(6,"6"),(30,"30"),(210,"210")]:
    if d<=300: ax[0].annotate(lab,(d,N3[d-1]/1000),fontsize=7,xytext=(0,3),textcoords="offset points",ha="center",color="#d62728")
ax[0].set_xlabel(r"common difference $\Delta$"); ax[0].set_ylabel(r"$N_3(\Delta)$  (thousands)")
ax[0].set_title(r"(A) 3-AP diffraction spectrum: extinction off $6\mathbb{Z}$",fontsize=9)
ax[0].set_xlim(0,300)
# (B) admissible sub-spectrum vs HL singular series, multiples of 6
m6=(D%6==0)&(D<=1512); Dm=D[m6]; ratio=N3[m6]/N3[5]   # N3[5] is Delta=6 (index 5)
def relS3(d):
    r=1.0
    for q in [5,7,11,13,17,19,23,29,31,37,41,43,47]:
        if d%q==0: r*=(q-1)/(q-3)
    return r
pred=np.array([relS3(int(d)) for d in Dm])
ax[1].plot(Dm,ratio,"o",ms=2.5,color="#1f77b4",alpha=0.6,label="measured")
ax[1].plot(Dm,pred,".",ms=2.0,color="#d62728",label=r"HL $\mathfrak{S}_3(\Delta)/\mathfrak{S}_3(6)$")
ax[1].set_xlabel(r"common difference $\Delta$ (multiples of 6)")
ax[1].set_ylabel(r"$N_3(\Delta)/N_3(6)$")
ax[1].set_title(r"(B) resonance enhancement vs singular series",fontsize=9)
ax[1].legend(fontsize=7,loc="upper left")
fig.suptitle(r"Prime arithmetic progressions on the $6N$ skeleton ($X=10^8$)",fontsize=10)
fig.savefig(os.path.join(_FIG,"p30_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 : unified (q-1)/(q-k) law ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
col={3:"#1f77b4",4:"#2ca02c",5:"#d62728"}; mk={3:"o",4:"s",5:"^"}
# (A) measured enhancement vs predicted, all k
for k in (3,4,5):
    meas=[float(r["N_over_base_measured"]) for r in res if int(r["k"])==k]
    pr=[float(r["enhancement_predicted"]) for r in res if int(r["k"])==k]
    ax2[0].scatter(pr,meas,s=22,marker=mk[k],color=col[k],alpha=0.8,label="$k=%d$"%k,edgecolors="none")
lim=[0.8,12]; ax2[0].plot(lim,lim,"k-",lw=0.7,alpha=0.6)
ax2[0].set_xlim(lim); ax2[0].set_ylim(lim)
ax2[0].set_xlabel(r"predicted $\prod_{q\mid\Delta,\,q>k}(q-1)/(q-k)$")
ax2[0].set_ylabel(r"measured $N_k(\Delta)/N_k(\mathrm{primorial})$")
ax2[0].set_title("(A) unified resonance law, $k=3,4,5$",fontsize=9)
ax2[0].legend(fontsize=8,loc="upper left")
# (B) single-prime gain (q-1)/(q-k) vs q for each k, with measured single-prime points
qs=np.arange(5,48)
for k in (3,4,5):
    qg=qs[qs>k]
    ax2[1].plot(qg,(qg-1)/(qg-k),"-",color=col[k],lw=1.0,label="$(q-1)/(q-%d)$"%k)
    # measured single-prime points: Delta = primorial(k)*q
    prim={3:6,4:6,5:30}[k]
    for r in res:
        if int(r["k"])==k and r["prime_factors_gt_k"] not in ("-","") and "," not in r["prime_factors_gt_k"]:
            q=int(r["prime_factors_gt_k"]); 
            ax2[1].plot(q,float(r["N_over_base_measured"]),mk[k],color=col[k],ms=5,alpha=0.9)
ax2[1].set_xlabel(r"prime $q$ dividing $\Delta$ (single extra factor)")
ax2[1].set_ylabel("enhancement")
ax2[1].set_title(r"(B) per-prime gain $(q-1)/(q-k)$",fontsize=9)
ax2[1].legend(fontsize=8); ax2[1].set_xlim(4,48)
fig2.savefig(os.path.join(_FIG,"p30_fig2.pdf")); print("fig2 done")
