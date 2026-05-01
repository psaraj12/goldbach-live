#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <signal.h>
#include <atomic>

using namespace std;

using u64 = uint64_t;
using u32 = uint32_t;
using u8  = uint8_t;

// ------------------ Checkpoint Config ------------------
static constexpr u64 CKPT_INTERVAL = 10000000000ULL;
static const char* CKPT_FILE = "goldbach_checkpoint.csv";

static atomic<int> g_shutdown{0};

static void signal_handler(int sig) {
    g_shutdown.store(1, memory_order_relaxed);
}

// ------------------ Deterministic 64-bit Miller-Rabin ------------------

static inline u64 mod_mul64(u64 a, u64 b, u64 mod) {
    return (u64)((__uint128_t)a * b % mod);
}

static inline u64 mod_pow64(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1ULL) r = mod_mul64(r, a, mod);
        a = mod_mul64(a, a, mod);
        d >>= 1ULL;
    }
    return r;
}

static bool is_prime64(u64 n) {
    if (n < 2) return false;
    static const u64 small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (u64 p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    static const u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a >= n) continue;
        u64 x = mod_pow64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul64(x, x, n);
            if (x == n - 1) {
                witness = false;
                break;
            }
        }
        if (witness) return false;
    }
    return true;
}

// ------------------ Prime Generation ------------------

static vector<u32> base_primes;

void build_base_primes(u32 limit) {
    vector<bool> is_prime((size_t)limit + 1, true);
    is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;

    for (u32 i = 2; (u64)i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (u64 j = (u64)i * i; j <= limit; j += i)
                is_prime[(size_t)j] = false;
        }
    }

    base_primes.clear();
    for (u32 i = 3; i <= limit; i += 2) {
        if (is_prime[i]) base_primes.push_back(i);
    }
}

// ------------------ Segmented Sieve ------------------

struct SegSieve {
    u64 lo = 0, hi = 0;
    vector<u8> bits;

    void build(u64 slo, u64 shi) {
        if (slo < 3) slo = 3;
        if ((slo & 1ULL) == 0) ++slo;

        lo = slo;
        hi = shi;

        u64 nbits  = (hi > lo) ? ((hi - lo) / 2 + 1) : 0;
        u64 nbytes = (nbits + 7) / 8;

        bits.assign((size_t)nbytes, 0xFF);

        if ((nbits & 7ULL) && nbytes > 0)
            bits[(size_t)nbytes - 1] &= (u8)((1u << (nbits & 7ULL)) - 1u);

        for (u32 p : base_primes) {
            u64 pp = (u64)p * p;
            if (pp >= hi) break;

            u64 start;
            if (pp >= lo) {
                start = pp;
            } else {
                u64 rem = lo % p;
                start = (rem == 0) ? lo : lo + (p - rem);
            }

            if ((start & 1ULL) == 0) start += p;
            if (start == p) start += 2ULL * p;

            u64 b = (start - lo) >> 1;
            for (; b < nbits; b += p)
                bits[(size_t)(b >> 3)] &= (u8)~(1u << (b & 7));
        }
    }

    inline bool is_prime(u64 q) const {
        if (q < 3 || (q & 1ULL) == 0) return q == 2;
        if (q < lo || q >= hi) return false;

        u64 b = (q - lo) >> 1;
        return (bits[(size_t)(b >> 3)] >> (b & 7)) & 1u;
    }
};

// ------------------ Simple QHot ------------------

static constexpr int QHOT_SIZE = 1024;

struct QHot {
    u64 buf[QHOT_SIZE];
    int pos = 0;
    int count = 0;

    inline void record(u64 q) {
        if (q < 3 || ((q & 1ULL) == 0)) return;

        if (count < QHOT_SIZE) {
            buf[pos] = q;
            ++pos;
            if (pos == QHOT_SIZE) pos = 0;
            ++count;
            return;
        }

        if (q > buf[pos]) {
            buf[pos] = q;
        }

        ++pos;
        if (pos == QHOT_SIZE) pos = 0;
    }

    inline bool hit(u64 N, u64 p_limit, const vector<u8>& anchor_lookup,
                    u64& p_out, u64& q_out) const {
        for (int i = 0; i < count; ++i) {
            int idx = pos - 1 - i;
            if (idx < 0) idx += QHOT_SIZE;

            u64 q = buf[idx];
            if (q >= N) continue;

            u64 p = N - q;
            if (p >= 3 && p <= p_limit && anchor_lookup[(size_t)p]) {
                p_out = p;
                q_out = q;
                return true;
            }
        }
        return false;
    }
};

// ------------------ Checkpoint (unified CSV) ------------------
// File format:
//   # Goldbach verification checkpoint + samples
//   # original_start=<val>
//   # original_end=<val>
//   # total_verified=<val>  ... etc
//   N,p,q,source
//   4000000000000000006,3,4000000000000000003,SIEVE
//   4000000000000000008,5,4000000000000000003,QHOT
//   ...
//
// Resume: last data row's N -> continue from N+2
// Stats: from # comment header
// Every row is a verifiable sample: p+q=N, both prime

struct CkptState {
    u64 original_start = 0;
    u64 original_end   = 0;
    u64 last_verified  = 0;    // actual last verified N (from header, not samples)
    u64 total_verified = 0;
    u64 total_total    = 0;
    u64 total_tries    = 0;
    u64 total_qhot     = 0;
    u64 total_misses   = 0;
    double elapsed_sec = 0.0;
    bool valid = false;
};

static CkptState load_checkpoint() {
    CkptState ck;
    ifstream f(CKPT_FILE);
    if (!f.is_open()) return ck;

    string line;

    while (getline(f, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            auto eq = line.find('=');
            if (eq == string::npos) continue;
            size_t ks = 1;
            while (ks < eq && line[ks] == ' ') ++ks;
            string key = line.substr(ks, eq - ks);
            string val = line.substr(eq + 1);

            try {
                if (key == "original_start")       ck.original_start  = stoull(val);
                else if (key == "original_end")    ck.original_end    = stoull(val);
                else if (key == "last_verified")   ck.last_verified   = stoull(val);
                else if (key == "total_verified")  ck.total_verified  = stoull(val);
                else if (key == "total_total")     ck.total_total     = stoull(val);
                else if (key == "total_tries")     ck.total_tries     = stoull(val);
                else if (key == "total_qhot")      ck.total_qhot      = stoull(val);
                else if (key == "total_misses")    ck.total_misses    = stoull(val);
                else if (key == "elapsed_sec")     ck.elapsed_sec     = stod(val);
            } catch (...) {}
            continue;
        }

        if (line.substr(0, 2) == "N,") continue;
        // data lines are samples, not needed for resume
    }

    if (ck.last_verified > 0) ck.valid = true;

    return ck;
}

static void save_checkpoint(u64 original_start, u64 original_end, u64 last_verified,
                            u64 total_verified, u64 total_total,
                            u64 total_tries, u64 total_qhot, u64 total_misses,
                            double elapsed_sec,
                            const vector<string>& all_samples) {
    string tmp = string(CKPT_FILE) + ".tmp";
    ofstream f(tmp, ios::out | ios::trunc);

    f << "# Goldbach verification checkpoint + samples\n";
    f << "# original_start=" << original_start << "\n";
    f << "# original_end=" << original_end << "\n";
    f << "# last_verified=" << last_verified << "\n";
    f << "# total_verified=" << total_verified << "\n";
    f << "# total_total=" << total_total << "\n";
    f << "# total_tries=" << total_tries << "\n";
    f << "# total_qhot=" << total_qhot << "\n";
    f << "# total_misses=" << total_misses << "\n";
    f << fixed << setprecision(6) << "# elapsed_sec=" << elapsed_sec << "\n";
    f << "N,p,q,source\n";

    for (auto& s : all_samples) f << s << "\n";

    f.close();
    rename(tmp.c_str(), CKPT_FILE);
}

// ------------------ Main ------------------

int main(int argc, char** argv) {
    if (argc < 5) {
        cout << "Usage: " << argv[0]
             << " <start> <end> <block_bits> <threads> [p_anchor_limit] [sample_limit]\n\n";
        cout << "Options:\n";
        cout << "  --resume    Resume from " << CKPT_FILE << "\n\n";
        cout << "Example:\n";
        cout << "  " << argv[0]
             << " 4000000000000000000 4001000000000000000 24 44 100000 10000\n";
        cout << "  " << argv[0]
             << " --resume 0 0 24 44 100000 10000\n";
        return 1;
    }

    signal(SIGINT,  signal_handler);
    signal(SIGTERM, signal_handler);

    bool resuming = false;
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "--resume") {
            resuming = true;
            break;
        }
    }

    u64 start = strtoull(argv[1], nullptr, 10);
    u64 end   = strtoull(argv[2], nullptr, 10);
    int block_bits = atoi(argv[3]);
    int threads    = atoi(argv[4]);
    u64 p_anchor_limit = (argc >= 6) ? strtoull(argv[5], nullptr, 10) : 1000000ULL;
    u64 sample_limit   = (argc >= 7) ? strtoull(argv[6], nullptr, 10) : 10000ULL;

    u64 cum_total = 0, cum_verified = 0, cum_tries = 0, cum_qhot = 0, cum_misses = 0;
    double cum_elapsed = 0.0;
    u64 original_start, original_end;

    // All sample lines (previous + new)
    vector<string> all_samples;

    if (resuming) {
        CkptState ck = load_checkpoint();
        if (!ck.valid) {
            cerr << "ERROR: --resume specified but no valid checkpoint found (" << CKPT_FILE << ")\n";
            return 1;
        }
        original_start = ck.original_start;
        original_end   = ck.original_end;
        start          = ck.last_verified + 2;
        end            = original_end;
        cum_total      = ck.total_total;
        cum_verified   = ck.total_verified;
        cum_tries      = ck.total_tries;
        cum_qhot       = ck.total_qhot;
        cum_misses     = ck.total_misses;
        cum_elapsed    = ck.elapsed_sec;

        // Read existing sample lines
        ifstream rf(CKPT_FILE);
        string line;
        while (getline(rf, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.substr(0, 2) == "N,") continue;
            all_samples.push_back(line);
        }

        cout << "RESUMING from checkpoint\n";
        cout << "  Original range: [" << original_start << " .. " << original_end << "]\n";
        cout << "  Last verified:  " << ck.last_verified << "\n";
        cout << "  Resume at:      " << start << "\n";
        cout << "  Samples kept:   " << all_samples.size() << "\n";
        cout << "  Cumulative:     " << cum_verified << " verified, "
             << fixed << setprecision(1) << cum_elapsed << " s elapsed\n\n";
    } else {
        original_start = start;
        original_end   = end;
    }

    if (start & 1ULL) ++start;
    if (end & 1ULL) --end;
    if (start < 6) start = 6;
    if (start > end) {
        if (resuming) {
            cout << "Checkpoint indicates run is already complete.\n";
            return 0;
        }
        cerr << "Invalid range after even adjustment.\n";
        return 1;
    }

    u64 BLOCK = 1ULL << block_bits;
    u64 SPAN  = 2ULL * BLOCK;

    u64 prime_limit64 = p_anchor_limit;
    if (prime_limit64 > numeric_limits<u32>::max()) {
        cerr << "p_anchor_limit too large for u32 base prime table: " << prime_limit64 << "\n";
        return 1;
    }

    cout << "goldbach_sieve_simple_qhot_verified\n";
    cout << "Range: [" << start << " .. " << end << "]\n";
    cout << "Block bits: " << block_bits << " | Threads: " << threads << "\n";
    cout << "P-anchor limit: " << p_anchor_limit << " | QHot size: " << QHOT_SIZE << "\n";
    cout << "Checkpoint interval: " << CKPT_INTERVAL << " even numbers\n";
    cout << "Sample limit: " << sample_limit << "\n";
    cout << "Building base primes up to " << prime_limit64 << "..." << flush;

    auto setup0 = chrono::steady_clock::now();
    build_base_primes((u32)prime_limit64);
    auto setup1 = chrono::steady_clock::now();
    double setup_sec = chrono::duration<double>(setup1 - setup0).count();

    cout << " done. Base primes: " << base_primes.size()
         << " | Setup: " << fixed << setprecision(3) << setup_sec << " s\n";

    if (p_anchor_limit > (u64)numeric_limits<size_t>::max() - 1) {
        cerr << "p_anchor_limit too large for lookup vector.\n";
        return 1;
    }

    vector<u8> anchor_lookup((size_t)p_anchor_limit + 1, 0);
    for (u32 p : base_primes) {
        if ((u64)p <= p_anchor_limit)
            anchor_lookup[(size_t)p] = 1;
        else
            break;
    }

    // ---- Process in mega-chunks for checkpointing ----
    u64 run_total = 0, run_verified = 0, run_tries = 0, run_qhot = 0, run_misses = 0;
    double run_elapsed = 0.0;
    u64 chunk_start = start;
    int ckpt_count = 0;

    while (chunk_start <= end) {
        u64 chunk_end = chunk_start + CKPT_INTERVAL * 2ULL - 2ULL;
        if (chunk_end > end) chunk_end = end;

        vector<pair<u64, u64>> blocks;
        for (u64 b = chunk_start; b <= chunk_end; ) {
            u64 blk_end = std::min<u64>(b + SPAN - 2ULL, chunk_end);
            blocks.push_back({b, blk_end});
            if (chunk_end - b < SPAN) break;
            b += SPAN;
        }

#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif

        u64 c_total = 0, c_verified = 0, c_tries = 0, c_qhot = 0, c_misses = 0;

        // Thread-0 sample buffer for this chunk
        vector<string> chunk_samples;

        auto t0 = chrono::steady_clock::now();

#pragma omp parallel
        {
            SegSieve seg;
            QHot qhot;

#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif

            vector<string> local_samples;
            u64 l_total = 0, l_verified = 0, l_tries = 0, l_qhot = 0, l_misses = 0;

#pragma omp for schedule(dynamic)
            for (size_t bi = 0; bi < blocks.size(); ++bi) {
                if (g_shutdown.load(memory_order_relaxed)) continue;

                auto [blk, blk_end] = blocks[bi];

                u64 lo = (blk > p_anchor_limit ? blk - p_anchor_limit : 3ULL);
                u64 hi = blk_end + 2ULL;

                seg.build(lo, hi);

                for (u64 N = blk; N <= blk_end; N += 2ULL) {
                    ++l_total;

                    u64 sp = 0, sq = 0;
                    if (qhot.hit(N, p_anchor_limit, anchor_lookup, sp, sq)) {
                        ++l_verified;
                        ++l_qhot;

                        if (tid == 0 && all_samples.size() + local_samples.size() < sample_limit) {
                            local_samples.push_back(
                                to_string(N) + "," + to_string(sp) + "," + to_string(sq) + ",QHOT");
                        }
                        continue;
                    }

                    bool ok = false;

                    for (u32 pr : base_primes) {
                        if ((u64)pr > p_anchor_limit) break;
                        if ((u64)pr >= N) break;

                        u64 q = N - pr;
                        ++l_tries;

                        if (seg.is_prime(q) && is_prime64(q)) {
                            ++l_verified;
                            qhot.record(q);

                            if (tid == 0 && all_samples.size() + local_samples.size() < sample_limit) {
                                local_samples.push_back(
                                    to_string(N) + "," + to_string(pr) + "," + to_string(q) + ",SIEVE");
                            }

                            ok = true;
                            break;
                        }
                    }

                    if (!ok) ++l_misses;
                }
            }

            // Merge thread-0 samples
            if (tid == 0) {
                chunk_samples = std::move(local_samples);
            }

#pragma omp atomic
            c_total += l_total;
#pragma omp atomic
            c_verified += l_verified;
#pragma omp atomic
            c_tries += l_tries;
#pragma omp atomic
            c_qhot += l_qhot;
#pragma omp atomic
            c_misses += l_misses;
        }

        auto t1 = chrono::steady_clock::now();
        double chunk_sec = chrono::duration<double>(t1 - t0).count();

        // If shutdown interrupted this chunk, don't save partial results
        // The previous checkpoint is still valid for resume
        if (g_shutdown.load(memory_order_relaxed)) {
            cout << "\n*** Shutdown signal received. ***\n";
            if (ckpt_count > 0) {
                cout << "Last valid checkpoint covers up to the previous chunk.\n";
            } else {
                cout << "No checkpoint was saved yet (interrupted during first chunk).\n";
            }
            cout << "Resume with: " << argv[0] << " --resume 0 0 "
                 << block_bits << " " << threads << " " << p_anchor_limit
                 << " " << sample_limit << "\n";
            break;
        }

        // Chunk completed fully — accumulate and save
        for (auto& s : chunk_samples) all_samples.push_back(s);

        run_total    += c_total;
        run_verified += c_verified;
        run_tries    += c_tries;
        run_qhot     += c_qhot;
        run_misses   += c_misses;
        run_elapsed  += chunk_sec;

        u64 gt  = cum_total + run_total;
        u64 gv  = cum_verified + run_verified;
        u64 gtr = cum_tries + run_tries;
        u64 gq  = cum_qhot + run_qhot;
        u64 gm  = cum_misses + run_misses;
        double ge = cum_elapsed + run_elapsed;

        save_checkpoint(original_start, original_end, chunk_end, gv, gt, gtr, gq, gm, ge, all_samples);
        ++ckpt_count;

        double overall_tp = ge > 0.0 ? (double)gt / ge / 1e6 : 0.0;
        double pct = 100.0 * (double)(chunk_end - original_start + 2) /
                     (double)(original_end - original_start + 2);

        cout << fixed << setprecision(2);
        cout << "[CKPT " << ckpt_count << "] "
             << "Up to " << chunk_end
             << " | " << pct << "% | "
             << gv << " verified | "
             << all_samples.size() << " samples | "
             << setprecision(1) << overall_tp << " M/s | "
             << setprecision(1) << ge << " s";
        if (c_misses > 0) cout << " | MISSES: " << c_misses;
        cout << "\n" << flush;

        if (chunk_end >= end) break;
        chunk_start = chunk_end + 2;
    }

    u64 gt  = cum_total + run_total;
    u64 gv  = cum_verified + run_verified;
    u64 gtr = cum_tries + run_tries;
    u64 gq  = cum_qhot + run_qhot;
    u64 gm  = cum_misses + run_misses;
    double ge = cum_elapsed + run_elapsed;

    cout << fixed << setprecision(6);
    cout << "\n===== FINAL SUMMARY =====\n";
    cout << "Full range:  [" << original_start << " .. " << original_end << "]\n";
    cout << "Coverage:    " << gv << " / " << gt << "\n";
    cout << "Misses:      " << gm << "\n";
    cout << "Tries/N:     " << (gt ? (double)gtr / (double)gt : 0.0) << "\n";
    cout << "QHot %:      " << (gt ? 100.0 * (double)gq / (double)gt : 0.0) << "%\n";
    cout << "Total time:  " << ge << " s\n";
    cout << "Throughput:  " << (ge > 0.0 ? (double)gt / ge / 1e6 : 0.0) << " M/s\n";
    cout << "Samples:     " << all_samples.size() << "\n";
    cout << "Checkpoints: " << ckpt_count << "\n";
    cout << "Output:      " << CKPT_FILE << "\n";

    if (gm != 0) return 2;
    if (g_shutdown.load(memory_order_relaxed)) return 3;
    return 0;
}
