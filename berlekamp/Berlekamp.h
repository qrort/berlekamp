#pragma once

#include <vector>

#include "Polynomial.h"

std::vector<std::pair<Polynomial, int>> berlekamp_factor(const Polynomial& poly, const ll& modp);
