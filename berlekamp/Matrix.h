#pragma once
#include <vector>
#include "Polynomial.h"

class Matrix {
public:
    Matrix(int size, ll modp);

    static Matrix identity(int size, ll modp);

    void set(int row, int column, const ll& value);

    ll get(int row, int column) const;

    int get_size() const;

    ll get_modp() const;

    Matrix operator+(const Matrix &rhs) const;

    Matrix operator-(const Matrix &rhs) const;

    Matrix get_transpose() const;

    void swap_rows(int row1, int row2);

    void sub_rows(int subfrom, int sub, ll multiplier = 1);

    void divide_row(int r, const ll& x);

    friend std::ostream& operator<<(std::ostream& o, const Matrix& m);
private:
    int size;
    ll modp;
    std::vector<ll> entries;
};

