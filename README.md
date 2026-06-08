# 6N-ap-resonance

**Part XXX — Prime Arithmetic Progressions as a Lattice Resonance: the Admissibility Selection Rule and the Singular-Series Enhancement Spectrum**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXX of *Arithmetic Geodynamics on the 6N Skeleton*.
We read the count of *k*-term arithmetic progressions of primes as a function of the common
difference Δ — and it behaves like a diffraction spectrum. **Everything here is a measured sieve
result — no fitted parameters, no fabricated numbers.**

## The question

A *k*-term AP of primes is `p, p+Δ, …, p+(k−1)Δ`, all prime. We sieve all primes to `X = 10⁸`
(5,761,455 primes) and count `N_k(Δ)` as a function of the common difference Δ.

## The result, in one paragraph

Two exact features of a lattice spectrum appear.

1. **Selection rule (extinction).** `N_k(Δ) = 0` unless Δ is divisible by the primorial
   `Π_{p≤k} p` (6 for k=3,4; 30 for k=5). This is admissibility: if a prime `p ≤ k` does not divide
   Δ, the *k* terms cover every residue mod *p*, so one is composite. We verify it as a **literal
   zero** — scanning *every* `Δ ≤ 1512` for k=3, the number of off-lattice (`6 ∤ Δ`) progressions
   found is exactly 0.

2. **Resonance enhancement.** On the admissible lattice, each extra prime `q > k` dividing Δ
   multiplies the count by the universal factor `(q−1)/(q−k)`:

   ```
   N_k(Δ·q) / N_k(Δ) = (q−1)/(q−k) ,    q > k, q ∤ Δ.
   ```

   This is the relative Hardy–Littlewood singular series
   `𝔖_k(Δ) = Π_p (1 − ν_p(Δ)/p)/(1 − 1/p)^k`, with `ν_p = 1` if `p|Δ` else `min(k,p)`. Verified to
   ≲1% across k = 3, 4, 5 (e.g. k=3 Δ=2310 → 3.768 vs 3.750; k=4 Δ=210 → 8.11 vs 8.00; k=5 Δ=2310
   → 4.955 vs 5.000). For k=3 the gain is `(q−1)/(q−3)`, for k=4 `(q−1)/(q−4)`, for k=5 `(q−1)/(q−5)`.

**Scope (honest).** Neither feature is a new theorem — both follow from the singular series. What this
adds is the spectral framing and the precision check at `X = 10⁸`, plus the unified per-prime
resonance law across *k*. The Green–Tao theorem (arbitrarily long progressions) is the untouchable
backdrop, not a target; no claim is made about the infinitude of any constellation.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_ap.py     # quick look: extinction off 6Z + the k=3 resonance ladder (console)
python3 final_ap.py       # dense k=3 spectrum (Δ=1..1512) + cross-k resonance table
                          #   -> data/ap_spectrum_k3.csv, ap_resonance.csv, ap_summary.csv
                          #   (~220 s, single core; the dense Δ-scan is the slow part)
python3 makefigs_ap.py    # reads the CSVs -> figures/p30_fig1.pdf, p30_fig2.pdf
```

Paths resolve relative to the script (outputs land in `../data` and `../figures`). The prime sieve to
`10⁸` needs ~0.1 GB RAM; single-threaded.

## Files

```
code/    explore_ap.py   final_ap.py   makefigs_ap.py
data/    ap_spectrum_k3.csv   Delta, N3            (Δ = 1..1512; the diffraction spectrum)
         ap_resonance.csv     k, primorial_k, Delta, prime_factors_gt_k, N_k,
                              N_over_base_measured, enhancement_predicted
         ap_summary.csv       parameter, value
figures/ p30_fig1.pdf  p30_fig2.pdf
paper/   paper30.tex   paper30.pdf
```

All data files are plain CSV — openable in any text editor or spreadsheet.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
