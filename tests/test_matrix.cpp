#include "../Matrix.h"
#include <vector>
#include <cmath>
#include <math.h>

// #define CATCH_CONFIG_MAIN
// #include "catch.hpp"

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace std;

TEST_CASE("Benchmark constructors", "[constructor]"){
    std::cout << "Constructor tests..." << endl << "\tDefault constructor:" << endl;
    Matrix<double> default_val = Matrix<double>();
    ComplexMatrix<double> cdefault = ComplexMatrix<double>();
    std::cout << "\tCustom constructor..." << endl;
    Matrix<double> custom_constructor(5, 8);
    ComplexMatrix<double> ccustom(5, 8);
    std::cout << "\tCopy constructor..." << endl;
    Matrix<double> copy_constructor((const Matrix<double>&) custom_constructor);
    ComplexMatrix<double> ccopy((const ComplexMatrix<double>&) ccustom);
    std::cout << "\tMove constructor..." << endl;
    Matrix<double> move_constructor((Matrix<double>&) copy_constructor);
    ComplexMatrix<double> cmove((ComplexMatrix<double>&) move_constructor);
    std::cout << "Done!" << endl;
}

TEST_CASE("Benchmark assignment operators", "[assignment]"){
    std::cout << "Assignment tests..." << endl;
    Matrix<double> default_val = Matrix<double>();
    ComplexMatrix<double> cdefault = ComplexMatrix<double>();
    std::cout << "\tCustom constructor..." << endl;
    Matrix<double> custom_constructor(5, 8);
    ComplexMatrix<double> ccustom(5, 8);
    std::cout << "\tCopy assignments..." << endl;
    Matrix<double> copy_assignment = (const Matrix<double>&) custom_constructor;
    ComplexMatrix<double> ccopy = (const ComplexMatrix<double>&) ccustom;
    std::cout << "\tMove assignment..." << endl;
    Matrix<double> move_constructor = (Matrix<double>&) copy_assignment;
    ComplexMatrix<double> cmove = (ComplexMatrix<double>&) move_constructor;
    std::cout << "Done!" << endl;
}

TEST_CASE("Benchmark basic math tests...", "[math]"){
    std::cout << "Basic math tests..." << endl;
    Matrix<int> mat_1_int(3, 3);
    Matrix<int> mat_2_int(3, 3);
    ComplexMatrix<int> mat_1_c(3, 3);
    ComplexMatrix<int> mat_2_c(3, 3);
    for(unsigned long long int i = 0; i < 3; ++i){
        for(unsigned long long int j = 0; j < 3; ++j){
            mat_1_int.setElement(i, j, i + j);
            mat_2_int.setElement(i, j, i - j);
            mat_1_c.setElement(i, j, complex<int>(i, j));
            mat_2_c.setElement(i, j, complex<int>(j, i));
        }
    }
    std::cout << "\tAddition tests..." << endl;
    std::cout << mat_1_int << " + " << mat_2_int << " = " << (mat_1_int + mat_2_int) << endl;
    std::cout << mat_1_c << " + " << mat_2_c << " = " << (mat_1_c + mat_2_c) << endl;
    std::cout << "\tSubtraction tests..." << endl;
    std::cout << mat_1_int << " - " << mat_2_int << " = " << (mat_1_int - mat_2_int) << endl;
    std::cout << mat_1_c << " - " << mat_2_c << " = " << (mat_1_c - mat_2_c) << endl;
    std::cout << "\tMultiplication tests..." << endl;
    std::cout << mat_1_int << " * " << mat_2_int << " = " << (mat_1_int * mat_2_int) << endl;
    std::cout << mat_1_c << " * " << mat_2_c << " = " << (mat_1_c * mat_2_c) << endl;
    std::cout << "Done!" << endl;
}

TEST_CASE("Benchmark determinant tests...", "[determinant]"){
    std::cout << "Determinant tests..." << endl;
    std::cout << "\tInt tests..." << endl;
    Matrix<int> identity_2(2, 2);
    identity_2.setElement(0, 0, 1);
    identity_2.setElement(1, 1, 1);
    REQUIRE(identity_2.determinant() == 1);
    std::cout << "\tDouble tests..." << endl;
    Matrix<double> rand_double(4, 4);
    ComplexMatrix<double> rand_complex(4, 4);
    for(unsigned long long int i = 0; i < 4; ++i){
        for(unsigned long long int j = 0; j < 4; ++j){
            rand_double.setElement(i, j, i + j + sin(i) + cos(j));
            rand_complex.setElement(i, j, complex<double>(j - i + sin(i) - cos(j), i - j + sin(i) - cos(j)));
        }
    }
    std::cout << "\t\tRandom matrix: " << rand_double << endl;
    std::cout << "\t\tDeterminant: " << rand_double.determinant() << endl;
    std::cout << "\t\tRandom matrix 2: " << rand_complex << endl;
    std::cout << "\t\tDeterminant: " << rand_complex.determinant() << endl;
    std::cout << "Done!" << endl;
}

int main( int argc, char* argv[] ) {
    Catch::Session session;
    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 ) { return returnCode; }
    int numFailed = session.run();
    return numFailed;
}