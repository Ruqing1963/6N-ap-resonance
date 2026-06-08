import numpy as np, time
t0=time.time()
X=100_000_000
sieve=np.ones(X+1,dtype=bool); sieve[:2]=False
for p in range(2,int(X**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
pos=np.nonzero(sieve)[0]; pos=pos[pos>=5]               # primes >=5
npr=int(sieve.sum())
print("X=%d  primes=%d  primes>=5=%d  (%.1fs)"%(X,npr,len(pos),time.time()-t0))

def N3(D):                                              # count 3-APs p,p+D,p+2D (p>=5) with p+2D<=X
    m=pos[pos+2*D<=X]
    return int(np.count_nonzero(sieve[m+D]&sieve[m+2*D]))

def N4(D):
    m=pos[pos+3*D<=X]
    return int(np.count_nonzero(sieve[m+D]&sieve[m+2*D]&sieve[m+3*D]))

# --- confirm the mod-2,3 obstruction: non-multiples of 6 should give ~0 ---
print("\n[k=3] non-admissible Delta (expect ~0, only sporadic small-prime APs):")
for D in [2,4,8,10,3,9,15,21]:
    print("  D=%-3d  N3=%d"%(D,N3(D)))

# --- the resonance ladder: multiples of 6 ---
print("\n[k=3] admissible Delta = multiples of 6:")
print("   D      factors>=5     N3(D)     N3/N3(6)   S3/S3(6)=prod (q-1)/(q-3)")
def relS3(D):                                           # singular-series enhancement vs D=6
    r=1.0; d=D
    for q in [5,7,11,13,17,19,23,29,31,37]:
        if d%q==0: r*=(q-1)/(q-3)
    return r
base=N3(6)
for D in [6,12,18,24,30,36,42,60,66,90,150,210,330,420,2310]:
    fac=[q for q in [5,7,11,13,17,19,23] if D%q==0]
    n=N3(D)
    print("  %5d  %-12s  %7d   %.4f     %.4f"%(D,",".join(map(str,fac)) or "-",n,n/base,relS3(D)))

print("\n[k=4] admissible Delta = multiples of 6:")
print("   D      factors>=5     N4(D)     N4/N4(6)")
b4=N4(6)
for D in [6,12,30,210,2310]:
    fac=[q for q in [5,7,11,13] if D%q==0]
    print("  %5d  %-10s  %7d   %.4f"%(D,",".join(map(str,fac)) or "-",N4(D),N4(D)/b4))

print("\nTOTAL %.1fs"%(time.time()-t0))
