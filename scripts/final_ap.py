import numpy as np, csv, time, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
t0=time.time()

X=100_000_000
sieve=np.ones(X+1,dtype=bool); sieve[:2]=False
for p in range(2,int(X**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
pos=np.nonzero(sieve)[0]; pos=pos[pos>=5]
nprimes=int(sieve.sum())
print("X=%d primes=%d (%.1fs)"%(X,nprimes,time.time()-t0))

def Ncount(D,k):                     # # of k-term APs p,p+D,...,p+(k-1)D, p>=5, last<=X
    m=pos[pos+(k-1)*D<=X]
    acc=sieve[m+D]
    for j in range(2,k): acc=acc&sieve[m+j*D]
    return int(np.count_nonzero(acc))

PRIMES_GT3=[5,7,11,13,17,19,23,29,31,37,41,43,47]
def primorial(k): 
    r=1
    for p in [2,3,5,7,11,13]:
        if p<=k: r*=p
    return r
def enhance(D,k):                    # predicted enhancement vs primorial(k): prod_{q>k, q|D}(q-1)/(q-k)
    r=1.0
    for q in PRIMES_GT3:
        if q>k and D%q==0: r*=(q-1)/(q-k)
    return r

# ---------- (1) k=3 dense diffraction scan: ALL Delta 1..Dmax ----------
Dmax=1512
spec=np.array([Ncount(D,3) for D in range(1,Dmax+1)])
print("k=3 dense scan done %.1fs"%(time.time()-t0))
with open(os.path.join(_DATA,"ap_spectrum_k3.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["Delta","N3"])
    for D in range(1,Dmax+1): w.writerow([D,int(spec[D-1])])
nz=np.nonzero(spec)[0]+1
offlattice=[D for D in nz if D%6!=0]
print("  nonzero Delta that are NOT multiples of 6:",offlattice,
      "(counts:",[int(spec[D-1]) for D in offlattice],")")

# ---------- (2) cross-k resonance verification table ----------
rows=[]
sel={3:[6,12,18,24,30,36,42,60,66,78,90,210,330,390,420,462,2310],
     4:[6,12,30,42,210,330,420,2310],
     5:[30,60,90,150,210,330,420,2310]}
for k in (3,4,5):
    pr=primorial(k); base=Ncount(pr,k)
    for D in sel[k]:
        n=Ncount(D,k); fac=[q for q in PRIMES_GT3 if q>k and D%q==0]
        rows.append((k,pr,D,",".join(map(str,fac)) or "-",n,n/base,enhance(D,k)))
print("resonance table done %.1fs"%(time.time()-t0))
with open(os.path.join(_DATA,"ap_resonance.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["k","primorial_k","Delta","prime_factors_gt_k","N_k","N_over_base_measured","enhancement_predicted"])
    for r in rows: w.writerow([r[0],r[1],r[2],r[3],r[4],"%.5f"%r[5],"%.5f"%r[6]])

print("\n k   primorial  Delta  facts>k     N_k     meas    pred")
for r in rows:
    print(" %d   %6d   %6d  %-10s %7d  %.4f  %.4f"%(r[0],r[1],r[2],r[3],r[4],r[5],r[6]))

# ---------- (3) summary ----------
with open(os.path.join(_DATA,"ap_summary.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["parameter","value"])
    for kk,vv in [("X_max",X),("n_primes",nprimes),("primorial_3",primorial(3)),
        ("primorial_4",primorial(4)),("primorial_5",primorial(5)),
        ("N3_Delta6",Ncount(6,3)),("N4_Delta6",Ncount(6,4)),("N5_Delta30",Ncount(30,5)),
        ("Dmax_dense_scan",Dmax),
        ("offlattice_nonzero_count",len(offlattice))]:
        w.writerow([kk,vv])
print("\nTOTAL %.1fs"%(time.time()-t0))
