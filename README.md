# 🚀 Goldbach Checkpoint Verifier (4×10¹⁸ Scale)

High-performance **Goldbach verification engine** for extremely large ranges with:

- ✅ Checkpoint / Resume support  
- ⚡ QHot witness reuse (~99.9% hit rate)  
- 🔍 Anchor primes (bounded search)  
- 🧠 Segmented sieve fallback  
- 📦 Sieve-Q witness dump (MR-confirmed)  
- 🔁 Independent reconstruction checker  

---

## 📌 What this does

Verifies that every even number `N` in a given range satisfies:

```
N = p + q   (p and q are prime)
```

⚠️ This is **computational verification**, not a formal proof.

---

## 🧩 Core Ideas

**QHot Cache**
- Reuses recent `q` values  
- ~99.9% hit rate → most evens resolved instantly  

**Anchor Primes**
- Only small primes `p ≤ limit` are tested  
- Reduces search drastically  

**Segmented Sieve Fallback**
- Used when QHot fails (~0.08%)  
- Produces **Miller–Rabin verified q values**

---

## 📊 Example Run

```
Range: [4e18 .. 4e18 + 1e10]
Total evens: 5,000,000,001
QHot hits:   4,995,840,973
Sieve hits:  4,159,028
Missing:     0
Bad:         0

QHot ≈ 99.9168%
Sieve ≈ 0.0832%
```

---

## 📂 Repository Structure

```
src/
  goldbach_verifier_ckpt.cpp
  goldbach_sieve_q_dump.cpp

tools/
  reconstruct_sieve_q_checker.cpp

data/
  checkpoints/
  samples/

SHA256SUMS.txt
```

---

## ⚙️ Build

```bash
g++ -O3 -march=native -std=gnu++17 -fopenmp src/goldbach_verifier_ckpt.cpp -o goldbach

g++ -O3 -march=native -std=gnu++17 -fopenmp src/goldbach_sieve_q_dump.cpp -o goldbach_sieve_q_dump

g++ -O3 -march=native -std=gnu++17 tools/reconstruct_sieve_q_checker.cpp -o reconstruct_checker
```

---

## ▶️ Run Verifier

```bash
./goldbach <start> <end> <block_bits> <threads> [p_anchor_limit] [sample_limit]
```

Example:

```bash
./goldbach 4000000000000000000 4000010000000000000 24 12 100000 10000
```

---

## 🔄 Resume from Checkpoint

```bash
./goldbach  0 0 24 12 100000 10000 --resume
```

---

## 🧪 Sieve-Q Dump Mode

```bash
./goldbach_sieve_q_dump ... --dump-sieve-q sieve_q_witness.csv
```

Output:

```
# adjusted_start=...
# adjusted_end=...
# p_anchor_limit=...
# qhot_size=1024
N,p,q
```

- Only fallback (~0.08%) stored  
- Each `q` is **MR-confirmed**  
- QHot hits are not stored  

---

## 🔁 Reconstruction Checker

```bash
./reconstruct_checker sieve_q_witness.csv
```

- Reconstructs all pairs using QHot + sieve data  
- Validates:
  - `p + q = N`
  - `p` is anchor prime  
  - `q` is prime (MR)

---

## 🔐 Integrity

```bash
sha256sum -c SHA256SUMS.txt
```

---

## 📈 Current Status

- Verified: **up to 100 trillion evens (checkpointed)**  
- Misses: **0**  
- Reconstruction check: **PASS**

---

## 🧠 Key Insight

Instead of storing all pairs (petabytes):

```
Store ~0.08% sieve witnesses → reconstruct 100% pairs
```

---

## ⚠️ Notes

- Not a formal mathematical proof  
- Full output intentionally avoided (too large)  
- Designed for **scaling beyond 10^18**

---

## 📜 License

MIT
