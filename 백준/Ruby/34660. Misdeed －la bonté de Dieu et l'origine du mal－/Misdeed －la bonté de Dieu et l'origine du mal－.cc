#include <bits/stdc++.h>
using namespace std;

// GF(16) with primitive polynomial x^4 + x + 1 -> 0b1_0011 = 0x13
static inline uint8_t gf_add(uint8_t a, uint8_t b) { return a ^ b; }

static inline uint8_t gf_mul(uint8_t a, uint8_t b) {
    uint8_t res = 0;
    while (b) {
        if (b & 1) res ^= a;
        b >>= 1;
        uint8_t hi = a & 0x8; // bit for x^3 before shift
        a <<= 1;
        if (hi) a ^= 0x3;     // because (x^4) ≡ (x+1) -> after shift overflow bit reduces by 0b0011
        // Simpler/common form using carry at bit 4:
        // a <<= 1; if (a & 0x10) a ^= 0x13;
        // But since we masked only 4 bits above, the hi trick also works.
        a &= 0xF;
    }
    return res & 0xF;
}

static inline uint8_t gf_mul_poly(uint8_t a, uint8_t b) {
    // Safer standard version (avoid confusion), keep this and use it.
    uint8_t res = 0;
    while (b) {
        if (b & 1) res ^= a;
        b >>= 1;
        uint8_t carry = a & 0x8;
        a <<= 1;
        if (carry) a ^= 0x3; // reduction for x^4 -> x+1 (0x13 without the 1<<4 bit)
        a &= 0xF;
    }
    return res;
}

// Use the standard carry form to be crystal clear
static inline uint8_t gf_mul_std(uint8_t a, uint8_t b) {
    uint8_t res = 0;
    while (b) {
        if (b & 1) res ^= a;
        b >>= 1;
        uint8_t na = a << 1;
        if (a & 0x8) na ^= 0x13; // reduce when degree >= 4
        a = na & 0xF;
    }
    return res & 0xF;
}

static inline uint8_t gf_pow(uint8_t a, int e) {
    uint8_t res = 1;
    while (e > 0) {
        if (e & 1) res = gf_mul_std(res, a);
        a = gf_mul_std(a, a);
        e >>= 1;
    }
    return res;
}

static inline uint8_t gf_inv(uint8_t a) {
    // inverse in GF(16): a^(15-1) = a^14, for a != 0
    if (a == 0) return 0; // inverse unused for zero in valid paths
    return gf_pow(a, 14);
}

vector<vector<uint8_t>> mat_inv_gf16(vector<vector<uint8_t>> A) {
    int n = (int)A.size();
    vector<vector<uint8_t>> I(n, vector<uint8_t>(n, 0));
    for (int i = 0; i < n; ++i) I[i][i] = 1;

    for (int col = 0; col < n; ++col) {
        int pivot = -1;
        for (int r = col; r < n; ++r) {
            if (A[r][col] != 0) { pivot = r; break; }
        }
        if (pivot == -1) {
            // Singular (shouldn't happen with distinct nodes)
            // But to be safe, we still try to continue (undefined).
            // In contest, this shouldn't occur.
            continue;
        }
        if (pivot != col) {
            swap(A[pivot], A[col]);
            swap(I[pivot], I[col]);
        }
        uint8_t inv_p = gf_inv(A[col][col]);
        // scale pivot row
        for (int j = 0; j < n; ++j) {
            A[col][j] = gf_mul_std(A[col][j], inv_p);
            I[col][j] = gf_mul_std(I[col][j], inv_p);
        }
        // eliminate
        for (int r = 0; r < n; ++r) if (r != col) {
            uint8_t f = A[r][col];
            if (f == 0) continue;
            for (int j = 0; j < n; ++j) {
                A[r][j] = gf_add(A[r][j], gf_mul_std(f, A[col][j]));
                I[r][j] = gf_add(I[r][j], gf_mul_std(f, I[col][j]));
            }
        }
    }
    return I;
}

vector<uint8_t> mat_vec_mul_gf16(const vector<vector<uint8_t>>& M, const vector<uint8_t>& v) {
    int n = (int)M.size();
    vector<uint8_t> out(n, 0);
    for (int i = 0; i < n; ++i) {
        uint8_t acc = 0;
        for (int j = 0; j < n; ++j) {
            acc = gf_add(acc, gf_mul_std(M[i][j], v[j]));
        }
        out[i] = acc;
    }
    return out;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t;
    if (!(cin >> t)) return 0;

    const int D = 7;   // degree bound per variable: 0..6 (7 coefficients)
    const int N = 13;  // grid size per dimension

    auto toNibble = [](const string& s, int idx) -> uint8_t {
        // take bits s[4*idx .. 4*idx+3], MSB first
        uint8_t v = 0;
        for (int k = 0; k < 4; ++k) {
            v <<= 1;
            v |= (s[4*idx + k] == '1') ? 1 : 0;
        }
        return v & 0xF;
    };

    if (t == 0) {
        // Encode: S (196 bits) -> A (13x13)
        string S;
        cin >> S;

        // coefficients c[i][j] for X^i * Y^j, i,j in 0..6
        uint8_t coeff[D][D];
        int idx = 0;
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
                coeff[i][j] = toNibble(S, idx++);
            }
        }

        // precompute powers of x=1..13 and y=1..13 up to power 6
        uint8_t xpow[N+1][D], ypow[N+1][D];
        for (int x = 1; x <= N; ++x) {
            xpow[x][0] = 1;
            for (int p = 1; p < D; ++p) xpow[x][p] = gf_mul_std(xpow[x][p-1], (uint8_t)x);
        }
        for (int y = 1; y <= N; ++y) {
            ypow[y][0] = 1;
            for (int p = 1; p < D; ++p) ypow[y][p] = gf_mul_std(ypow[y][p-1], (uint8_t)y);
        }

        // evaluate P at (x=1..13, y=1..13)
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                uint8_t val = 0;
                for (int a = 0; a < D; ++a) {
                    for (int b = 0; b < D; ++b) {
                        uint8_t term = gf_mul_std(coeff[a][b], gf_mul_std(xpow[i][a], ypow[j][b]));
                        val ^= term;
                    }
                }
                cout << (int)val << (j == N ? '\n' : ' ');
            }
        }
    } else {
        // Decode: Given r[7], c[7], and B(7x7), recover S (196 bits)
        vector<int> r(D), c(D);
        for (int i = 0; i < D; ++i) cin >> r[i];
        for (int i = 0; i < D; ++i) cin >> c[i];

        uint8_t B[D][D];
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
                int x; cin >> x;
                B[i][j] = (uint8_t)(x & 0xF);
            }
        }

        // Build Vandermonde matrices for x and y
        vector<vector<uint8_t>> Mx(D, vector<uint8_t>(D, 0));
        for (int k = 0; k < D; ++k) {
            uint8_t x = (uint8_t)r[k];
            uint8_t p = 1;
            for (int d = 0; d < D; ++d) {
                Mx[k][d] = p;
                p = gf_mul_std(p, x);
            }
        }
        auto Mx_inv = mat_inv_gf16(Mx);

        vector<vector<uint8_t>> My(D, vector<uint8_t>(D, 0));
        for (int k = 0; k < D; ++k) {
            uint8_t y = (uint8_t)c[k];
            uint8_t p = 1;
            for (int d = 0; d < D; ++d) {
                My[k][d] = p;
                p = gf_mul_std(p, y);
            }
        }
        auto My_inv = mat_inv_gf16(My);

        // Step 1: For each fixed y = c[j], interpolate in x to get coeffs of X^i (i=0..6)
        uint8_t GX[D][D]; // GX[i][j] = coefficient of X^i evaluated at y = c[j]
        for (int j = 0; j < D; ++j) {
            vector<uint8_t> v(D);
            for (int k = 0; k < D; ++k) v[k] = B[k][j]; // f_j(x_k)
            auto coeff_x = mat_vec_mul_gf16(Mx_inv, v); // size D
            for (int i = 0; i < D; ++i) GX[i][j] = coeff_x[i];
        }

        // Step 2: For each i (coeff of X^i), interpolate in y to get full polynomial in Y
        uint8_t coeff[D][D]; // coeff[i][j] for X^i * Y^j
        for (int i = 0; i < D; ++i) {
            vector<uint8_t> v(D);
            for (int j = 0; j < D; ++j) v[j] = GX[i][j]; // g_i(y_j) at y=c[j]
            auto coeff_y = mat_vec_mul_gf16(My_inv, v); // size D
            for (int j = 0; j < D; ++j) coeff[i][j] = coeff_y[j];
        }

        // Output S: coefficients in order i=0..6, j=0..6, each as 4 bits (MSB first)
        string out;
        out.reserve(196);
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
                uint8_t v = coeff[i][j] & 0xF;
                for (int b = 3; b >= 0; --b) {
                    out.push_back(((v >> b) & 1) ? '1' : '0');
                }
            }
        }
        cout << out << '\n';
    }
    return 0;
}