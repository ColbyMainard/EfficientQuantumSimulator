#include "../Matrix.h"
#include <vector>
#include <cmath>
#include <math.h>

void constructor_tests(){
    cout << "Constructor tests..." << endl << "\tDefault constructor:" << endl;
    Matrix<double> default_val = Matrix<double>();
    ComplexMatrix<double> cdefault = ComplexMatrix<double>();
    cout << "\tCustom constructor..." << endl;
    Matrix<double> custom_constructor(5, 8);
    ComplexMatrix<double> ccustom(5, 8);
    cout << "\tCopy constructor..." << endl;
    Matrix<double> copy_constructor((const Matrix<double>&) custom_constructor);
    ComplexMatrix<double> ccopy((const ComplexMatrix<double>&) ccustom);
    cout << "\tMove constructor..." << endl;
    Matrix<double> move_constructor((Matrix<double>&) copy_constructor);
    ComplexMatrix<double> cmove((ComplexMatrix<double>&) move_constructor);
    cout << "Done!" << endl;
}

void assignment_tests(){
    cout << "Assignment tests..." << endl;
    Matrix<double> default_val = Matrix<double>();
    ComplexMatrix<double> cdefault = ComplexMatrix<double>();
    cout << "\tCustom constructor..." << endl;
    Matrix<double> custom_constructor(5, 8);
    ComplexMatrix<double> ccustom(5, 8);
    cout << "\tCopy assignments..." << endl;
    Matrix<double> copy_assignment = (const Matrix<double>&) custom_constructor;
    ComplexMatrix<double> ccopy = (const ComplexMatrix<double>&) ccustom;
    cout << "\tMove assignment..." << endl;
    Matrix<double> move_constructor = (Matrix<double>&) copy_assignment;
    ComplexMatrix<double> cmove = (ComplexMatrix<double>&) move_constructor;
    cout << "Done!" << endl;
}

void basic_math_tests(){
    cout << "Basic math tests..." << endl;
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
    cout << "\tAddition tests..." << endl;
    cout << mat_1_int << " + " << mat_2_int << " = " << (mat_1_int + mat_2_int) << endl;
    cout << mat_1_c << " + " << mat_2_c << " = " << (mat_1_c + mat_2_c) << endl;
    cout << "\tSubtraction tests..." << endl;
    cout << mat_1_int << " - " << mat_2_int << " = " << (mat_1_int - mat_2_int) << endl;
    cout << mat_1_c << " - " << mat_2_c << " = " << (mat_1_c - mat_2_c) << endl;
    cout << "\tMultiplication tests..." << endl;
    cout << mat_1_int << " * " << mat_2_int << " = " << (mat_1_int * mat_2_int) << endl;
    cout << mat_1_c << " * " << mat_2_c << " = " << (mat_1_c * mat_2_c) << endl;
    cout << "Done!" << endl;
}

void determinant_tests(){
    cout << "Determinant tests..." << endl;
    cout << "\tInt tests..." << endl;
    Matrix<int> identity_2(2, 2);
    identity_2.setElement(0, 0, 1);
    identity_2.setElement(1, 1, 1);
    if(identity_2.determinant() != 1){
        throw runtime_error("Answer should be 1.");
    }
    cout << "\tDouble tests..." << endl;
    Matrix<double> rand_double(4, 4);
    ComplexMatrix<double> rand_complex(4, 4);
    for(unsigned long long int i = 0; i < 4; ++i){
        for(unsigned long long int j = 0; j < 4; ++j){
            rand_double.setElement(i, j, i + j + sin(i) + cos(j));
            rand_complex.setElement(i, j, complex<double>(j - i + sin(i) - cos(j), i - j + sin(i) - cos(j)));
        }
    }
    cout << "\t\tRandom matrix: " << rand_double << endl;
    cout << "\t\tDeterminant: " << rand_double.determinant() << endl;
    cout << "\t\tRandom matrix 2: " << rand_complex << endl;
    cout << "\t\tDeterminant: " << rand_complex.determinant() << endl;
    cout << "Done!" << endl;
}


int main(){
    try {
        constructor_tests();
        assignment_tests();
        basic_math_tests();
        determinant_tests();
        return 0;
    } catch(exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Unknown exception." << endl;
        return 2;
    }
}