#include <iostream>
#include <algorithm>
#include "BigFloat.h"
#include "catch2/catch_session.hpp"
#include "catch2/generators/catch_generators.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

TEST_CASE( "[BigFloat public methods]", "[All]" ) {
    BigFloat a = 1.96_bf;
    REQUIRE(a.toString(2) == "1.96");
    a.inverseSign();
    REQUIRE(a == -1.96_bf);
    REQUIRE(abs(a) == 1.96_bf);
    REQUIRE(a.toString(2) == "-1.96");
}

TEST_CASE("[BigFloat operators]", "[All]") {
    SECTION("+ whole numbers") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        BigFloat res = BigFloat(a) + BigFloat(b);
        REQUIRE(res == BigFloat(a + b));
    }
    SECTION("+ real numbers") {
        int a = GENERATE(take(5,random(-100, 100)));
        int b = GENERATE(take(5,random(-100, 100)));
        int delA = GENERATE(take(5,random(1, 100)));
        int delB = GENERATE(take(5,random(1, 100)));
        BigFloat x = BigFloat(a) / BigFloat(delA);
        BigFloat y = BigFloat(b) / BigFloat(delB);
        int lcm = delA * delB / std::gcd(delA, delB);
        BigFloat res = BigFloat(a * (lcm / delA) + b * (lcm / delB)) / BigFloat(lcm);
        REQUIRE(abs(x + y - res) < BigFloat("0.000000001"));
    }
    SECTION("*") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        REQUIRE(BigFloat(a) * BigFloat(b) == BigFloat(a * b));
    }
    SECTION("- whole numbers") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        REQUIRE(BigFloat(a) - BigFloat(b) == BigFloat(a - b));
    }
}

TEST_CASE("[BigFloat comparation operators]", "[All]") {
    SECTION("== and != operators") {
        int a = GENERATE(take(10,random(-100, 100)));
        REQUIRE(BigFloat(a) == BigFloat(a));
        REQUIRE(BigFloat(a) != BigFloat(a + 1));
    }
    SECTION("< and > operators") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        if (a < b) {
            INFO(a);
            INFO(b);
            REQUIRE(BigFloat(a) < BigFloat(b));
            REQUIRE(BigFloat(b) > BigFloat(a));
        } else if (a > b) {
            INFO(a);
            INFO(b);
            REQUIRE(BigFloat(a) > BigFloat(b));
            REQUIRE(BigFloat(b) < BigFloat(a));
        }
    }
}

