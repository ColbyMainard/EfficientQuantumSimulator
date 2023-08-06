#include "../Qubit.h"

void constructor_tests(){
    cout << "Constructor tests..." << endl;
    cout << "\tDefault constructor..." << endl;
    Qubit<double> default_constructor;
    cout << "\tCopy constructor..." << endl;
    Qubit<double> copy_constructor((const Qubit<double>&) default_constructor);
    cout << "\tMove constructor..." << endl;
    Qubit<double> move_constructor((Qubit<double>&) copy_constructor);
    cout << "Done!" << endl;
}

void assignment_tests(){
    cout << "Assignment tests..." << endl;
    Qubit<float> default_cons;
    cout << "\tCopy assignment..." << endl;
    Qubit<float> copy_assignment = (const Qubit<float>&) default_cons;
    cout << "\tMove assignment..." << endl;
    Qubit<float> move_assignment = (Qubit<float>&) copy_assignment;
    cout << "Done!" << endl;
}

void setting_states(){
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

void normalization_tests(){
    cout << "Normalization tests..." << endl;
    ComplexMatrix<double> test_1(2, 1);
    test_1.setElement(0, 0, complex<double>(1, 2));
    test_1.setElement(1, 0, complex<double>(3, 4));
    Qubit<double> bit_1;
    bit_1.setState(test_1);
    if(bit_1.magnitude() < 0.999999 || bit_1.magnitude() > 1.000001){
        cout << bit_1 << endl;
        throw runtime_error("Failed to normalize bit 1.");
    } else {
        cout << "\tBit 1 test succeeded..." << endl;
    }
    ComplexMatrix<double> test_2(2, 1);
    test_2.setElement(0, 0, complex<double>(0, 0));
    test_2.setElement(1, 0, complex<double>(0, 2));
    Qubit<double> bit_2;
    bit_2.setState((const ComplexMatrix<double>) test_2);
    if(bit_2.magnitude() < 0.999999 || bit_2.magnitude() > 1.000001){
        cout << bit_2 << endl;
        throw runtime_error("Failed to normalize bit 2.");
    } else {
        cout << "\tBit 2 test succeeded..." << endl;
    }
    ComplexMatrix<double> test_3(2, 1);
    test_3.setElement(0, 0, complex<double>(3, 0));
    test_3.setElement(1, 0, complex<double>(4, 0));
    Qubit<double> bit_3;
    bit_3.setState((const ComplexMatrix<double>)test_3);
    if(bit_3.magnitude() < 0.999999 || bit_3.magnitude() > 1.000001){
        cout << bit_3 << endl;
        throw runtime_error("Failed to normalize bit 3.");
    } else {
        cout << "\tBit 3 test succeeded..." << endl;
    }
    cout << "Done!" << endl;
}

int main(){
    try{
        constructor_tests();
        assignment_tests();
        setting_states();
        normalization_tests();
    } catch(exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Unknown exception." << endl;
        return 2;
    }
}