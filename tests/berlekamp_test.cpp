#include "gtest/gtest.h"

#include <set>
#include "Polynomial.h"
#include "Berlekamp.h"

using namespace std;

namespace {
    bool
    check_answer(std::vector<std::pair<Polynomial, int>> expected, std::vector<std::pair<Polynomial, int>> result) {
        set<std::pair<Polynomial, int>> s1(expected.begin(), expected.end());
        set<std::pair<Polynomial, int>> s2(result.begin(), result.end());

        return s1 == s2;
    };
}


TEST(Berlekamp, basic) {
    ll modp = 3;
    Polynomial poly = Polynomial("x^12+x^11+2x^9+2x^8+2x^6+x^5+2x^4+2x^3", modp);

    std::vector<std::pair<Polynomial, int>> expected = {
        {Polynomial("x^2+x+2", 3), 1},
        {Polynomial("x^3+2x+1", 3), 1},
        {Polynomial("x^2+1", 3), 2},
        {Polynomial("x", 3), 3},
    };

    auto result = berlekamp_factor(poly, 3);
    EXPECT_TRUE(check_answer(expected, result));
}


TEST(Berlekamp, basic_2) {
    Polynomial poly = Polynomial("x^11+x^10+x^9+x^8+x^6+x^5+1", 2);

    std::vector<std::pair<Polynomial, int>> expected = {
        {Polynomial("x^3+x^1+1", 2), 1},
        {Polynomial("x^2+x^1+1", 2), 1},
        {Polynomial("x^3+x^2+1", 2), 2},
    };

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(check_answer(expected, result));
}


TEST(Berlekamp, equal_degree_factors) {
    Polynomial poly = Polynomial("x^24+10x^23+32x^22+28x^21+13x^20+13x^19+17x^18+12x^17+29x^16+16x^15+13x^14+30x^13+31x^12+31x^11+2x^10+15x^9+5x^8+15x^7+3x^6+10x^5+18x^4+4x^3+6x^2+36x+19", 37);

    std::vector<std::pair<Polynomial, int>> expected = {
        {Polynomial("x^12+2x^11+30x^10+27x^9+23x^8+7x^7+4x^6+15x^5+8x^4+10x^3+16x^2+8x+15", 37), 1},
        {Polynomial("x^12+8x^11+23x^10+11x^9+24x^8+7x^7+25x^6+12x^5+30x^4+16x^3+34x^2+6x+21", 37), 1},
    };

    auto result = berlekamp_factor(poly, 37);

    EXPECT_TRUE(check_answer(expected, result));
}


TEST(Berlekamp, equal_degree_factors_2) {
    Polynomial poly = Polynomial("x^14+x^10+x^9+x^8+x^7+x^6+x^3+x^2+1", 2);

    std::vector<std::pair<Polynomial, int>> expected = {
        {Polynomial("x^7+x^4+1", 2), 1},
        {Polynomial("x^7+x^4+x^3+x^2+1", 2), 1},
    };

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(check_answer(expected, result));
}


TEST(Berlekamp, many_factors) {
    Polynomial poly = Polynomial("x^63+1", 2);

    std::vector<std::pair<Polynomial, int>> expected = {
        {Polynomial("x+1", 2), 1},
        {Polynomial("x^2+x+1", 2), 1},
        {Polynomial("x^3+x^2+1", 2), 1},
        {Polynomial("x^3+x+1", 2), 1},
        {Polynomial("x^6+x^4+x^3+x+1", 2), 1},
        {Polynomial("x^6+x^5+1", 2), 1},
        {Polynomial("x^6+x^5+x^3+x^2+1", 2), 1},
        {Polynomial("x^6+x^5+x^4+x+1", 2), 1},
        {Polynomial("x^6+x^4+x^2+x+1", 2), 1},
        {Polynomial("x^6+x^3+1", 2), 1},
        {Polynomial("x^6+x+1", 2), 1},
        {Polynomial("x^6+x^5+x^4+x^2+1", 2), 1},
        {Polynomial("x^6+x^5+x^2+x+1", 2), 1},
    };

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE(check_answer(expected, result));
}