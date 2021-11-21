#include "Polynomial.h"

#include <string>
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>
#include <cassert>

using namespace std;

namespace {
    std::vector<std::string> split_by(std::string s, std::string delimiter) {
        std::vector<std::string> result;

        auto start = 0U;
        auto end = s.find(delimiter);
        while (end != std::string::npos)
        {
            result.push_back(s.substr(start, end - start));
            start = end + delimiter.length();
            end = s.find(delimiter, start);
        }

        result.push_back(s.substr(start, end));

        return result;
    }
}

Polynomial::Polynomial(std::string str, ll modp) : modp(modp) {
    try {
        auto monoms_string = split_by(str, "+");
        vector<pair<ll, int>> monoms(monoms_string.size());
        for (size_t i = 0; i < monoms.size(); i++) {
            auto ss = split_by(monoms_string[i], "x");
            ll coeff;
            int pwr;

            if (ss.size() == 1) {
                pwr = 0;
            } else {
                if (ss[1] == "") {
                    pwr = 1;
                } else {
                    pwr = stoi(ss[1].substr(1, ss[1].size()));
                }
            }

            if (ss[0] == "") {
                coeff = 1;
            } else {
                coeff = stoull(ss[0]);
            }
            monoms[i] = { coeff, pwr };
        }


        int pwr = monoms[0].second;
        vector<ll> result(pwr + 1);

        for (auto pr : monoms)
        {
            result[pr.second] = pr.first;
        }

        coeff = result;
    } catch (...) {
        cout << "Error during polynomial parse" << endl;
    }
}


void Polynomial::prune()
{
    for (int i = coeff.size(); i > 0 && coeff.back() == 0; i--)
    {
        coeff.pop_back();
    }
}


int Polynomial::get_degree() const {
    if (coeff.empty()) return 0;
    return coeff.size() - 1;
}


std::string Polynomial::to_string(const std::string& default_variable_name) const
{
    if (coeff.empty()) {
        return "0";
    }

    int degree = get_degree();

    std::string result;

    for (int i = degree; i >= 0; i--)
    {
        if (coeff[i] == 0) continue;
        if (!result.empty()) {
            result += "+";
        }

        if (coeff[i] != 1 || i == 0) {
            result += std::to_string(coeff[i]);
        }
        if (i != 0) {
            result += default_variable_name;
            if (i != 1) {
                result += "^" + std::to_string(i);
            }
        }
    }

    return result;
}


Polynomial Polynomial::diff() const {
    vector<ll> v(get_degree());

    for (size_t i = 0; i < get_degree(); i++)
    {
        v[i] = (coeff[i + 1] * (i + 1)) % modp;
    }

    auto p = Polynomial(v, modp);
    p.prune();

    return p;
}

Polynomial Polynomial::get_pth_root() const {
    vector <ll> root_coeff(get_degree() / modp + 1);
    for (int i = 0; i <= get_degree(); i += modp) {
        root_coeff[i / modp] = coeff[i];
    }
    auto res = Polynomial(root_coeff, modp);
    res.prune();
    return res;
}


std::ostream& operator<<(std::ostream& strm, const Polynomial& poly) {
    return strm << poly.to_string();
}


Polynomial Polynomial::add(const Polynomial & a, const Polynomial & b) {
    assert(a.modp == b.modp);
    if (a.is_zero()) {
        return b;
    }
    if (b.is_zero()) {
        return a;
    }

    const Polynomial& a1 = a.get_degree() > b.get_degree() ? a : b;
    const Polynomial& b1 = a.get_degree() > b.get_degree() ? b : a;

    std::vector<ll> v(a1.get_degree() + 1);

    for (size_t i = 0; i <= b1.get_degree(); i++)
    {
        v[i] = (a1.coeff[i] + b1.coeff[i] + a.modp) % a.modp;
    }

    for (size_t i = b1.get_degree() + 1; i <= a1.get_degree(); i++)
    {
        v[i] = a1.coeff[i];
    }

    Polynomial p = Polynomial(v, a.modp);
    p.prune();
    return p;
}


Polynomial Polynomial::mul(const Polynomial & a, const Polynomial & b) {
    assert(a.modp == b.modp);
    if (a.is_zero() || b.is_zero()) {
        return Polynomial();
    }

    std::vector<ll> v(a.get_degree() + b.get_degree() + 1);

    for (int i = 0; i <= a.get_degree(); ++i) {
        for (int j = 0; j <= b.get_degree(); ++j) {
            v[i + j] = (v[i + j] + a.coeff[i] * b.coeff[j] % a.modp + a.modp) % a.modp;
        }
    }

    return Polynomial(v, a.modp);
}



Polynomial Polynomial::sub(const Polynomial & a, const Polynomial & b) {
    vector<ll> v = b.coeff;

    for (int i = 0; i < v.size(); i++) {
        v[i] = a.modp - v[i];
    }

    return add(a, Polynomial(v, a.modp));
}


std::pair< Polynomial, Polynomial> Polynomial::div_internal(const Polynomial & a, const Polynomial & b) {
    assert(a.modp == b.modp);
    ll bl = b.coeff.back();
    Polynomial b2 = b.normalize();

    int degree_of_result = a.get_degree() - b.get_degree() + 1;

    if (degree_of_result < 1) {
        return { Polynomial(), a };
    }

    std::vector<ll> coeff_result(degree_of_result);

    Polynomial at = a;

    for (int i = 0; i < degree_of_result; i++)
    {
        coeff_result[degree_of_result - 1 - i] = at.coeff[a.get_degree() - i];

        for (int j = 0; j <= b2.get_degree(); j++)
        {
            at.coeff[a.get_degree() - i - j] = (at.coeff[a.get_degree() - i - j] - coeff_result[degree_of_result - 1 - i] * b2.coeff[b2.get_degree() - j]) % a.modp;
            at.coeff[a.get_degree() - i - j] = (at.coeff[a.get_degree() - i - j] + a.modp) % a.modp;
        }
    }

    Polynomial quo = mul(Polynomial(coeff_result, a.modp), Polynomial(vector<ll>{ bl }, a.modp));
    
    at.prune();

    return { quo, at };
}


Polynomial Polynomial::div(const Polynomial & a, const Polynomial & b) {
    return Polynomial::div_internal(a, b).first;
}


Polynomial Polynomial::mod(const Polynomial & a, const Polynomial & b) {
    return Polynomial::div_internal(a, b).second;
}


Polynomial Polynomial::get_one(ll modp) {
    return Polynomial(vector<ll>{ 1 }, modp);
}


bool operator==(const Polynomial & poly1, const Polynomial & poly2) {
    return poly1.coeff == poly2.coeff;
}


bool operator!=(const Polynomial & poly1, const Polynomial & poly2) {
    return !(poly1 == poly2);
}


bool operator< (const Polynomial &poly1, const Polynomial &poly2) {
    if (poly1.coeff.size() == poly2.coeff.size()) {
        for (int i = poly1.coeff.size() - 1; i >= 0; i--) {
            if (poly1.coeff[i] != poly2.coeff[i]) {
                return poly1.coeff[i] < poly2.coeff[i];
            }
        }

        return false;
    }

    return poly1.coeff.size() < poly2.coeff.size();
}


Polynomial Polynomial::gcd(const Polynomial& a1, const Polynomial& b1) {
    assert(a1.modp == a1.modp);
    if (a1.is_zero()) {
        return b1;
    }
    if (b1.is_zero()) {
        return a1;
    }

    auto a = a1.normalize();
    auto b = b1.normalize();

    while (!b.is_zero()) {
        a = Polynomial::mod(a, b);
        auto z = a;
        a = b;
        b = z;
    }
    return a.normalize();
}

Polynomial Polynomial::powmod(const Polynomial &a, ll b, const Polynomial &mod){
    assert(a.modp == mod.modp);
    ll power = b;
    Polynomial rez = Polynomial::get_one(a.modp);
    Polynomial aa = a;
    while (power > 0) {
        if (power % 2 == 1) {
            rez = Polynomial::mul(rez, aa);
            rez = Polynomial::mod(rez, mod);
        }
        aa = Polynomial::mul(aa, aa);
        aa = Polynomial::mod(aa, mod);
        power /= 2;
    }
    return rez;
}

ll binpow(ll a, ll b, ll mod) {
    if (!b) return 1;
    if (b % 2) {
        return a * binpow(a, b - 1, mod) % mod;
    } else {
        ll res = binpow(a, b / 2, mod) % mod;
        return res * res;
    }
}


ll Polynomial::inverse(const ll& a, const ll& modp) {
    return binpow(a, modp - 2, modp);
}


Polynomial Polynomial::normalize() const
{
    ll bl = coeff.back();
    ll ib = inverse(bl, modp);

    vector<ll> v = coeff;

    for (int i = 0; i < v.size(); i++) {
        v[i] = (v[i] * ib) % modp;
    }

    return Polynomial(v, modp);
}


bool Polynomial::is_zero() const {
    for (auto i : coeff) {
        if (i != 0) {
            return false;
        }
    }
    return true;
}


bool Polynomial::is_one() const
{
    return coeff.size() == 1 && coeff[0] == 1;
}


Polynomial Polynomial::get_random_polynomial(int max_degree, const ll& modq)
{
    long long poly_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 rng(poly_seed);

    uniform_int_distribution<int> dist(0, modq - 1);
    uniform_int_distribution<int> dist_degree(0, max_degree - 1);

    auto degree = dist_degree(rng) + 1;

    std::vector<ll> vr(degree + 1);
    vr[degree] = 1;

    for (size_t i = 0; i < degree; i++)
    {
        vr[i] = dist(rng);
    }

    Polynomial result(vr, modq);

    result.prune();

    return result;
}
