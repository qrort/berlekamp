#include "Berlekamp.h"
#include "Polynomial.h"
#include "Matrix.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <cassert>
#include <unordered_map>

using namespace std;

std::vector<std::pair<Polynomial, int>> squarefree_decompose(const Polynomial& poly) {
    std::vector<std::pair<Polynomial, int>> result;

    Polynomial f = poly;
    Polynomial g;
    ll modp = poly.get_modp();
    ll m = 1;

    do {
        auto df = f.diff();
        g = Polynomial::gcd(f, df);
        auto t = Polynomial::div(f, g);
        int i = 1;
        while (!t.is_one()) {
            auto tt = Polynomial::gcd(t, g);
            auto qq = Polynomial::div(t, tt);
            if (!qq.is_one()) {
                result.emplace_back( qq, i * m );
            }
            t = tt;
            g = Polynomial::div(g, tt);
            i++;
        }
        if (!g.is_one()) {
            f = g.get_pth_root();
            m = m * modp;
        }
    } while (!g.is_one());

    return result;
}

Matrix calculate_Q(const Polynomial &poly, const ll &modp) {//O(qd^2)
    int sz = poly.get_degree();
    Matrix res(sz, modp);
    Polynomial p("1", modp);
    auto cf = p.get_coeffs(sz);
    for (int i = 0; i < sz; i++) {
        res.set(0, i, cf[i]);
    }
    vector <ll> tmp(modp + 1, 0);
    tmp[tmp.size() - 1] = 1;
    auto pn = Polynomial(tmp, modp);
    for (int i = 1; i < sz; i++) {
        p = p * pn;
        p = p % poly;
        cf = p.get_coeffs(sz);
        cout << "cf:\n";
        for (auto c : cf) {
            cout << c << " ";
        }
        cout << '\n';
        for (int j = 0; j < sz; j++) {
            res.set(i, j, cf[j]);
        }
    }
    return res;
}

Matrix rowEchelonForm(Matrix M) {
    int lead = 0;
    int n = M.get_size();
    for (int r = 0; r < n; r++) {
        if (lead >= n) {
            break;
        }
        int i = r;
        while (M.get(i, lead) == 0) {
            i++;
            if (i == n) {
                i = r;
                lead++;
                if (lead == n) {
                    break;
                }
            }
        }
        M.swap_rows(i, r);
        if (M.get(r, lead) != 0) {
            M.divide_row(r, M.get(r, lead));
        }
        for (i = 0; i < n; i++) {
            if (i != r) {
                M.sub_rows(i, r, M.get(i, lead));
            }
        }
        lead++;
    }
    return M;
}

// Returns a list of vectors in the null space of (Q-I)^T
std::vector<Polynomial> Q_eigenvectors(const Matrix &Q) {
    cout << "Q:\n" << Q << "\n";
    auto A = (Q - Matrix::identity(Q.get_size(), Q.get_modp())).get_transpose();
    // Gaussian elimination
    cout << "A:\n" << A << "\n";
    A = rowEchelonForm(A);
    cout << "AFTER GAUSS:\n" << A << "\n";
    // Now we'll solve Au = 0
    int n = A.get_size();
    vector <int> pivots(n, -1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A.get(i, j) != 0) {
                pivots[j] = i;
                break;
            }
        }
    }
    std::vector<Polynomial> basis;
    std::vector<vector<ll>> bas;
    auto modp = A.get_modp();
    for (int i = 0; i < n; i++) {
        if (pivots[i] == -1) {
            bas.emplace_back(n, -1);
            auto& vec = bas[bas.size() - 1];
            vec[i] = 1;
            for (int j = 0; j < n; j++) if (j != i && pivots[j] == -1) {
                vec[j] = 0;
            }
            int id = 0;
            for (int j = 0; j < n; j++) {
                if (vec[j] == -1) {
                    ll tmp = (modp - A.get(id, i)) % modp;
                    vec[j] = tmp;
                    assert(vec[id] >= 0);
                    id++;
                }
            }
            //reverse(vec.begin(), vec.end());
            basis.emplace_back(vec, modp);
        }
    }
    cout << "BASIS:\n";
    for (const auto& vec : bas) {
        for (auto el : vec) {
            cout << el << ' ';
        }
        cout << '\n';
    }
    cout << "as poly:\n";
    for (const auto& vec : basis) {
        cout << vec.to_string() << '\n';

    }
    return basis;
}

vector<Polynomial> factor(const Polynomial &poly, const ll &modp) {
    cout << "DECOMPOSING:\n";
    cout << poly.to_string() << '\n';
    if ( poly.get_degree() <= 1 ) {
        return std::vector<Polynomial>{poly};
    }
    auto Q = calculate_Q(poly, modp);
    auto basis = Q_eigenvectors(Q);
    std::vector<Polynomial> factors{poly};
    int k = 1;
    while (factors.size() < basis.size()) {
        std::vector<Polynomial> newfactors;
        for (int s = 0; s < modp; s++) {
            for (const auto& w : factors) {
                Polynomial ww = Polynomial::gcd(w, basis[k] - Polynomial(vector<ll>{s}, modp));
                if (!ww.is_one()) {
                    newfactors.push_back(ww);
                }
            }
        }
        swap(factors, newfactors);
        k += 1;
    }
    return factors;
}

vector<pair<Polynomial, int>> berlekamp_factor(const Polynomial& poly, const ll& modp) {
    vector<pair<Polynomial, int>> result;
    vector<pair<Polynomial, int>> sqrfree = squarefree_decompose(poly);
    for (auto const& value : sqrfree) {
        auto r1 = factor(value.first, modp);
        for (const auto& i : r1) {
            result.emplace_back(i, value.second);
        }
    }
    return result;

}
