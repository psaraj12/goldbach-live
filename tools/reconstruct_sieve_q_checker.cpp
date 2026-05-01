#include <bits/stdc++.h>
using namespace std;

using u64 = uint64_t;
using u32 = uint32_t;
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
static constexpr int QHOT_SIZE = 1024;

struct QHot {
    u64 buf[QHOT_SIZE]{};
    int pos = 0, count = 0;

    void record(u64 q) {
        if (q < 3 || (q & 1ULL) == 0) return;

        if (count < QHOT_SIZE) {
            buf[pos] = q;
            pos = (pos + 1) % QHOT_SIZE;
            ++count;
            return;
        }

        if (q > buf[pos]) buf[pos] = q;
        pos = (pos + 1) % QHOT_SIZE;
    }

    bool hit(u64 N, const vector<uint8_t>& anchor_lookup, u64& p_out, u64& q_out) const {
        for (int i = 0; i < count; ++i) {
            int idx = pos - 1 - i;
            if (idx < 0) idx += QHOT_SIZE;

            u64 q = buf[idx];
            if (q >= N) continue;

            u64 p = N - q;
            if (p < anchor_lookup.size() && anchor_lookup[(size_t)p]) {
                p_out = p;
                q_out = q;
                return true;
            }
        }
        return false;
    }
};

vector<u32> build_primes(u32 limit) {
    vector<uint8_t> is_prime(limit + 1, 1);
    if (limit >= 0) is_prime[0] = 0;
    if (limit >= 1) is_prime[1] = 0;

    for (u32 i = 2; (u64)i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (u64 j = (u64)i * i; j <= limit; j += i)
                is_prime[(size_t)j] = 0;
        }
    }

    vector<u32> primes;
    for (u32 i = 3; i <= limit; i += 2)
        if (is_prime[i]) primes.push_back(i);

    return primes;
}

static inline string trim(string s) {
    while (!s.empty() && isspace((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && isspace((unsigned char)s.back())) s.pop_back();
    return s;
}

struct Meta {
    u64 adjusted_start = 0;
    u64 adjusted_end = 0;
    u64 p_anchor_limit = 100000;
};

struct SieveRow {
    u64 N, p, q;
};

bool parse_csv_row(const string& line, SieveRow& r) {
    stringstream ss(line);
    string a, b, c;

    if (!getline(ss, a, ',')) return false;
    if (!getline(ss, b, ',')) return false;
    if (!getline(ss, c, ',')) return false;

    r.N = stoull(trim(a));
    r.p = stoull(trim(b));
    r.q = stoull(trim(c));
    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage:\n";
        cerr << "  " << argv[0] << " sieve_q_witness.csv [reconstructed_pairs.csv]\n\n";
        cerr << "Example:\n";
        cerr << "  " << argv[0] << " sieve_q_witness.csv reconstructed_pairs.csv\n";
        return 1;
    }

    string in_file = argv[1];
    string out_file = (argc >= 3) ? argv[2] : "reconstructed_pairs.csv";

    ifstream in(in_file);
    if (!in.is_open()) {
        cerr << "ERROR: cannot open input file: " << in_file << "\n";
        return 1;
    }

    Meta meta;
    vector<SieveRow> rows;
    string line;

    while (getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (line[0] == '#') {
            auto eq = line.find('=');
            if (eq == string::npos) continue;

            string key = trim(line.substr(1, eq - 1));
            string val = trim(line.substr(eq + 1));

            if (key == "adjusted_start") meta.adjusted_start = stoull(val);
            else if (key == "adjusted_end") meta.adjusted_end = stoull(val);
            else if (key == "p_anchor_limit") meta.p_anchor_limit = stoull(val);

            continue;
        }

        if (line == "N,p,q") continue;

        SieveRow r{};
        if (parse_csv_row(line, r))
            rows.push_back(r);
    }

    if (meta.adjusted_start == 0 || meta.adjusted_end == 0) {
        cerr << "ERROR: missing # adjusted_start / # adjusted_end in sieve_q file.\n";
        return 1;
    }

    sort(rows.begin(), rows.end(), [](const SieveRow& a, const SieveRow& b) {
        return a.N < b.N;
    });

    cout << "Range: [" << meta.adjusted_start << " .. " << meta.adjusted_end << "]\n";
    cout << "p_anchor_limit: " << meta.p_anchor_limit << "\n";
    cout << "sieve rows: " << rows.size() << "\n";

    if (meta.p_anchor_limit > numeric_limits<u32>::max()) {
        cerr << "ERROR: p_anchor_limit too large.\n";
        return 1;
    }

    auto primes = build_primes((u32)meta.p_anchor_limit);

    vector<uint8_t> anchor_lookup((size_t)meta.p_anchor_limit + 1, 0);
    for (u32 p : primes) anchor_lookup[p] = 1;

    ofstream out(out_file);
    if (!out.is_open()) {
        cerr << "ERROR: cannot open output file: " << out_file << "\n";
        return 1;
    }

    out << "N,p,q,source\n";

    QHot qhot;
    size_t row_idx = 0;

    u64 total = 0;
    u64 qhot_hits = 0;
    u64 sieve_hits = 0;
    u64 missing = 0;
    u64 bad = 0;

    u64 start = meta.adjusted_start;
    u64 end = meta.adjusted_end;

    if (start & 1ULL) ++start;
    if (end & 1ULL) --end;
    if (start < 6) start = 6;

    for (u64 N = start; N <= end; N += 2) {
        ++total;

        u64 p = 0, q = 0;

        if (qhot.hit(N, anchor_lookup, p, q)) {
            ++qhot_hits;
			if (qhot_hits <10000)
			{
            out << N << "," << p << "," << q << ",QHOT\n";
			}
            continue;
        }

        while (row_idx < rows.size() && rows[row_idx].N < N)
            ++row_idx;

        if (row_idx >= rows.size() || rows[row_idx].N != N) {
            ++missing;
            cerr << "MISSING sieve fallback for N=" << N << "\n";
            continue;
        }

        q = rows[row_idx].q;
        p = N - q;

        // sanity checks
if (p < 3 || (p & 1ULL) == 0) {
    cerr << "INVALID p (even or <3): N=" << N << " p=" << p << " q=" << q << "\n";
    ++bad;
    ++row_idx;
    continue;
}

// anchor check (fast)
if (p >= anchor_lookup.size() || !anchor_lookup[(size_t)p]) {
    cerr << "ANCHOR FAIL: N=" << N << " p=" << p << " q=" << q << "\n";
    ++bad;
    ++row_idx;
    continue;
}

// MR validation for credibility
if (!is_prime64(q) || !is_prime64(p)) {
    cerr << "MR FAIL: N=" << N << " p=" << p << " q=" << q << "\n";
    ++bad;
    ++row_idx;
    continue;
}

        ++sieve_hits;
        qhot.record(q);
		if (sieve_hits< 10000)
		{
        out << N << "," << p << "," << q << ",SIEVE\n";
		}
        ++row_idx;
    }

    cout << "\n===== SUMMARY =====\n";
    cout << "Total evens: " << total << "\n";
    cout << "QHot hits:   " << qhot_hits << "\n";
    cout << "Sieve hits:  " << sieve_hits << "\n";
    cout << "Missing:     " << missing << "\n";
    cout << "Bad:         " << bad << "\n";
    cout << "Output:      " << out_file << "\n";

    if (missing || bad) return 2;
    return 0;
}