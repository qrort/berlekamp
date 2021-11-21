#pragma once

#include <cstdlib>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <algorithm>

typedef long long ll;

class Polynomial
{
private:
    std::vector <ll> coeff;
    ll modp;
    void prune();

    static std::pair<Polynomial, Polynomial> div_internal(const Polynomial& a, const Polynomial& b);

public:
    Polynomial(std::vector<ll> coeff, ll modp) : coeff(std::move(coeff)), modp(modp) {
        prune();
    };

    explicit Polynomial() {}

    Polynomial(std::string s, ll modp);

    Polynomial(const Polynomial &polynomial) : coeff(polynomial.coeff), modp(polynomial.modp) {};

    Polynomial &operator= (const Polynomial &polynomial) = default;

    int get_degree() const;

    ll get_modp() const {
        return modp;
    }

    std::string to_string(const std::string& default_variable_name = "x") const;

    Polynomial diff() const;

    Polynomial get_pth_root() const;

    std::vector <ll> get_coeffs(int n) const {
        auto x = coeff;
        while (x.size() < n) {
            x.emplace_back(0);
        }
        return x;
    }

    static ll inverse(const ll& a, const ll& modp);

    static Polynomial add(const Polynomial& a, const Polynomial& b);

    static Polynomial sub(const Polynomial& a, const Polynomial& b);

    static Polynomial mul(const Polynomial& a, const Polynomial& b);

    static Polynomial div(const Polynomial& a, const Polynomial& b);

    static Polynomial mod(const Polynomial& a, const Polynomial& b);

    Polynomial operator+(const Polynomial &rhs) const {
        return add(*this, rhs);
    }

    Polynomial operator-(const Polynomial &rhs) const {
        return sub(*this, rhs);
    }
    Polynomial operator*(const Polynomial &rhs) const {
        return mul(*this, rhs);
    }

    Polynomial operator/(const Polynomial &rhs) const {
        return div(*this, rhs);
    }

    Polynomial operator%(const Polynomial &rhs) const {
        return mod(*this, rhs);
    }

    static Polynomial get_one(ll modp);

    static Polynomial get_random_polynomial(int max_degree, const ll& modq);

    Polynomial normalize() const;

    static Polynomial gcd(const Polynomial& a, const Polynomial& b);

    Polynomial powmod(const Polynomial& a, ll b, const Polynomial& mod);

    bool is_zero() const;

    bool is_one() const;

    friend bool operator== (const Polynomial &poly1, const Polynomial &poly2);
    friend bool operator< (const Polynomial &poly1, const Polynomial &poly2);
    friend bool operator!= (const Polynomial &poly1, const Polynomial &poly2);
};


std::ostream & operator<<(std::ostream & Str, Polynomial const & v);

