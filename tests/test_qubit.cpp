#include "../Qubit.h"

// #define CATCH_CONFIG_MAIN
// #include "catch.hpp"

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace std;

TEST_CASE("Benchmark constructors", "[constructor]"){
    cout << "Constructor tests..." << endl;
    cout << "\tDefault constructor..." << endl;
    Qubit<double> default_constructor;
    cout << "\tCopy constructor..." << endl;
    Qubit<double> copy_constructor((const Qubit<double>&) default_constructor);
    cout << "\tMove constructor..." << endl;
    Qubit<double> move_constructor((Qubit<double>&) copy_constructor);
    cout << "Done!" << endl;
}

TEST_CASE("Benchmark assignment operators", "[assignment]"){
    cout << "Assignment tests..." << endl;
    Qubit<float> default_cons;
    cout << "\tCopy assignment..." << endl;
    Qubit<float> copy_assignment = (const Qubit<float>&) default_cons;
    cout << "\tMove assignment..." << endl;
    Qubit<float> move_assignment = (Qubit<float>&) copy_assignment;
    cout << "Done!" << endl;
}

TEST_CASE("Benchmark setting states", "[setstate]"){
    cout << "Setting states tests..." << endl;
    Qubit<double> test_qubit;
    cout << test_qubit << endl;
    cout << "\tTurn on..." << endl;
    test_qubit.turnOn();
    cout << test_qubit << endl;
    cout << "\tTurn off..." << endl;
    test_qubit.turnOff();
    cout << test_qubit << endl;
    cout << "\tSet state..." << endl;
    ComplexMatrix<double> test_value(2, 1);
    test_value.setElement(0, 0, sqrt(2) / 2);
    test_value.setElement(1, 0, sqrt(2) / 2);
    test_qubit.setState(test_value);
    cout << test_qubit << endl;
    cout << "Done!" << endl;
}

TEST_CASE("Benchmark normalization", "[normalization]"){
    cout << "Normalization tests..." << endl;
    ComplexMatrix<double> test_1(2, 1);
    test_1.setElement(0, 0, complex<double>(1, 2));
    test_1.setElement(1, 0, complex<double>(3, 4));
    Qubit<double> bit_1;
    bit_1.setState(test_1);
    REQUIRE(bit_1.magnitude() > 0.999999);
    REQUIRE(bit_1.magnitude() < 1.000001);
    ComplexMatrix<double> test_2(2, 1);
    test_2.setElement(0, 0, complex<double>(0, 0));
    test_2.setElement(1, 0, complex<double>(0, 2));
    Qubit<double> bit_2;
    bit_2.setState((const ComplexMatrix<double>) test_2);
    REQUIRE(bit_2.magnitude() > 0.999999);
    REQUIRE(bit_2.magnitude() < 1.000001);
    ComplexMatrix<double> test_3(2, 1);
    test_3.setElement(0, 0, complex<double>(3, 0));
    test_3.setElement(1, 0, complex<double>(4, 0));
    Qubit<double> bit_3;
    bit_3.setState((const ComplexMatrix<double>)test_3);
    REQUIRE(bit_3.magnitude() > 0.999999);
    REQUIRE(bit_3.magnitude() < 1.000001);
    cout << "Done!" << endl;
}

int main( int argc, char* argv[] ) {
    Catch::Session session;
    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 ) { return returnCode; }
    int numFailed = session.run();
    return numFailed;
}