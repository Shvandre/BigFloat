#include <iostream>
#include "BigFloat.h"
#include "catch2/catch_session.hpp"
#include "catch2/generators/catch_generators.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

TEST_CASE( "All", "[BigFloat public methods]" ) {
    BigFloat a = 1.96_bf;
    REQUIRE(a.toString(2) == "1.96");
    a.inverseSign();
    REQUIRE(a == -1.96_bf);
    REQUIRE(abs(a) == 1.96_bf);
    REQUIRE(a.toString(2) == "-1.96");
}

TEST_CASE("All", "[BigFloat operators]") {
    SECTION("+ whole numbers") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        REQUIRE(BigFloat(a) + BigFloat(b) == BigFloat(a + b));
    }
    SECTION("+ real numbers") {
        int a = GENERATE(take(5,random(-100, 100)));
        int b = GENERATE(take(5,random(-100, 100)));
        int delA = GENERATE(10, 100, 1000, 1000);
        int delB = GENERATE(10, 100, 1000, 1000);
        BigFloat x = BigFloat(a) / BigFloat(delA);
        BigFloat y = BigFloat(b) / BigFloat(delB);
        REQUIRE(abs((x + y - BigFloat(a + b)) / BigFloat(delA * delB)) < BigFloat("0.000000001"));
    }
    SECTION("*") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        INFO(a);
        INFO(b);
        a = -60;
        b = 0;
        REQUIRE(BigFloat(a) * BigFloat(b) == BigFloat(a * b));
    }
    SECTION("- whole numbers") {
        int a = GENERATE(take(10,random(-100, 100)));
        int b = GENERATE(take(10,random(-100, 100)));
        REQUIRE(BigFloat(a) - BigFloat(b) == BigFloat(a - b));
    }
}

TEST_CASE("All", "[BigFloat comparation operators]") {
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


