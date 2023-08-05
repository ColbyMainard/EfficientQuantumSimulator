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