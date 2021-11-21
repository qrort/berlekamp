#include "Matrix.h"
#include <iostream>
#include <utility>

Matrix::Matrix(int si, ll modp) : size(si), entries(si*si), modp(std::move(modp)) {
}

void Matrix::set(int row, int column, const ll& value)
{
    entries[column + size * row] = value;
}

ll Matrix::get(int row, int column) const {
    return entries[column + size * row];
}

int Matrix::get_size() const {
    return size;
}

ll Matrix::get_modp() const {
    return modp;
}

Matrix Matrix::identity(int size, ll modp) {
    Matrix m{size, modp};
    for (int i=0; i < size; ++i) {
        m.set(i,i,1);
    }
    return m;
}

Matrix Matrix::operator+(const Matrix& rhs) const {
    int si = std::min(get_size(), rhs.get_size());
    Matrix m{si, modp};
    for (int row = 0; row < si; ++row) {
        for (int col = 0; col < si; ++col) {
            m.set(row, col, (get(row, col) + rhs.get(row, col) + modp) % modp);
        }
    }
    return m;
}

Matrix Matrix::operator-(const Matrix& rhs)const {
    int si = std::min(get_size(), rhs.get_size());
    Matrix m{si, modp};
    for (int row = 0; row < si; ++row) {
        for (int col = 0; col < si; ++col) {
            m.set(row, col, (modp + modp + get(row, col) - rhs.get(row, col)) % modp);
        }
    }
    return m;
}

Matrix Matrix::get_transpose()const
{
    Matrix m{get_size(), modp};
    for (int row = 0; row < get_size(); ++row) {
        for (int col = 0; col < get_size(); ++col) {
            m.set(col, row, get(row, col));
        }
    }
    return m;
}

void Matrix::swap_rows(int row1, int row2)
{
    for (int col = 0; col < get_size(); ++col) {
        ll a = get(row1, col);
        ll b = get(row2, col);
        set(row1, col, b);
        set(row2, col, a);
    }
}

void Matrix::sub_rows(int subfrom, int sub, ll multiplier) {
    for (int col = 0; col < get_size(); ++col) {
        set(subfrom, col, (get(subfrom, col) - (get(sub, col) * multiplier) % modp + modp + modp) % modp);
    }
}

std::ostream &operator<<(std::ostream &o, const Matrix &m) {
    for (int row = 0; row < m.get_size(); ++row) {
        for (int col = 0; col < m.get_size(); ++col) {
            o << m.get(row, col) << ' ';
        }
        o << '\n';
    }
    return o;
}

void Matrix::divide_row(int r, const ll& x) {
    ll inv = Polynomial::inverse(x, modp);
    for (int col = 0; col < get_size(); ++col) {
        auto value = get(r, col);
        value = (value * inv) % modp;
        set(r, col, value);
    }
}
