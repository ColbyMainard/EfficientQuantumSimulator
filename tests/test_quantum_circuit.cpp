#include "../QuantumCircuit.h"

using namespace std;

void test_constructors(){
    cout << "Test constructors..." << endl;
    cout << "\tDefault constructor..." << endl;
    QuantumCircuit<float> default_con_float;
    QuantumCircuit<double> default_con_double;
    cout << "\tInt constructor..." << endl;
    QuantumCircuit<float> int_con_float(3);
    QuantumCircuit<double> int_con_double(3);
    cout << "\tAdding gates..." << endl;
    vector<vector<unsigned long long int>> indices = {{0}, {1}, {2}, {0}, {1}, {2}};
    int_con_float.addGate(X_gate<float>(), indices[0]);
    int_con_double.addGate(X_gate<double>(), indices[0]);
    int_con_float.addGate(Y_gate<float>(), indices[1]);
    int_con_double.addGate(Y_gate<double>(), indices[1]);
    int_con_float.addGate(Z_gate<float>(), indices[2]);
    int_con_double.addGate(Z_gate<double>(), indices[2]);
    int_con_float.addGate(Rx_gate<float>(.1), indices[3]);
    int_con_double.addGate(Rx_gate<double>(.1), indices[3]);
    int_con_float.addGate(Ry_gate<float>(.2), indices[4]);
    int_con_double.addGate(Ry_gate<double>(.2), indices[4]);
    int_con_float.addGate(S_gate<float>(), indices[5]);
    int_con_double.addGate(S_gate<double>(), indices[5]);
    cout << "\tCopy constructor..." << endl;
    QuantumCircuit<float> copy_float((const QuantumCircuit<float>&) int_con_float);
    QuantumCircuit<double> copy_double((const QuantumCircuit<double>&) int_con_double);
    cout << "\tMove constructor..." << endl;
    QuantumCircuit<float> move_float((QuantumCircuit<float>&) copy_float);
    QuantumCircuit<double> move_double((QuantumCircuit<double>&) copy_double);
    cout << "Done!" << endl;
}

void test_assignment_operators(){
    cout << "Test assignment operators..." << endl;
    cout << "\tDefault constructor..." << endl;
    QuantumCircuit<float> default_con_float;
    QuantumCircuit<double> default_con_double;
    cout << "\tInt constructor..." << endl;
    QuantumCircuit<float> int_con_float(3);
    QuantumCircuit<double> int_con_double(3);
    cout << "\tAdding gates..." << endl;
    vector<vector<unsigned long long int>> indices = {{0}, {1}, {2}, {0}, {1}, {2}};
    int_con_float.addGate(X_gate<float>(), indices[0]);
    int_con_double.addGate(X_gate<double>(), indices[0]);
    int_con_float.addGate(Y_gate<float>(), indices[1]);
    int_con_double.addGate(Y_gate<double>(), indices[1]);
    int_con_float.addGate(Z_gate<float>(), indices[2]);
    int_con_double.addGate(Z_gate<double>(), indices[2]);
    int_con_float.addGate(Rx_gate<float>(.1), indices[3]);
    int_con_double.addGate(Rx_gate<double>(.1), indices[3]);
    int_con_float.addGate(Ry_gate<float>(.2), indices[4]);
    int_con_double.addGate(Ry_gate<double>(.2), indices[4]);
    int_con_float.addGate(S_gate<float>(), indices[5]);
    int_con_double.addGate(S_gate<double>(), indices[5]);
    cout << "\tCopy constructor..." << endl;
    QuantumCircuit<float> copy_float = (const QuantumCircuit<float>&) int_con_float;
    QuantumCircuit<double> copy_double = (const QuantumCircuit<double>&) int_con_double;
    cout << "\tMove constructor..." << endl;
    QuantumCircuit<float> move_float = (QuantumCircuit<float>&) copy_float;
    QuantumCircuit<double> move_double = (QuantumCircuit<double>&) copy_double;
    cout << "Done!" << endl;
}

void x_test(){
    cout << "\tX Test..." << endl;
    auto gate = X_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void y_test(){
    cout << "\tY Test..." << endl;
    auto gate = Y_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void z_test(){
    cout << "\tZ Test..." << endl;
    auto gate = Z_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void h_test(){
    cout << "\tH Test..." << endl;
    auto gate = H_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void rx_test(){
    cout << "\tRx Test..." << endl;
    auto gate = Rx_gate<double>(.2);
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void ry_test(){
    cout << "\tRy Test..." << endl;
    auto gate = Ry_gate<double>(.2);
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void rz_test(){
    cout << "\tRz Test..." << endl;
    auto gate = Rz_gate<double>(.2);
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void s_test(){
    cout << "\tS Test..." << endl;
    auto gate = S_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void t_test(){
    cout << "\tT Test..." << endl;
    auto gate = T_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void cnot_test(){
    cout << "\tCNOT Test..." << endl;
    auto gate = CNOT_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void cz_test(){
    cout << "\tCZ Test..." << endl;
    auto gate = CZ_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void swqp_test(){
    cout << "\tSWAP Test..." << endl;
    auto gate = SWAP_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void toffoli_test(){
    cout << "\tToffoli Test..." << endl;
    auto gate = Toffoli_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void fredkin_test(){
    cout << "\tFredkin Test..." << endl;
    auto gate = Fredkin_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void bell_test(){
    cout << "\tBell Test..." << endl;
    auto gate = Bell_gate<double>();
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void test_gate_behaviors(){
    cout << "Gate tests..." << endl;
    x_test();
    y_test();
    z_test();
    h_test();
    rx_test();
    ry_test();
    rz_test();
    s_test();
    t_test();
    cnot_test();
    cz_test();
    swqp_test();
    toffoli_test();
    fredkin_test();
    bell_test();
    cout << "Done!" << endl;
}

void x_circuit_test(){
    cout << "\tX Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void y_circuit_test(){
    cout << "\tY Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void z_circuit_test(){
    cout << "\tZ Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void h_circuit_test(){
    cout << "\tH Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void rx_circuit_test(){
    cout << "\tRx Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void ry_circuit_test(){
    cout << "\tRy Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void rz_circuit_test(){
    cout << "\tRz Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void s_circuit_test(){
    cout << "\tS Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void t_circuit_test(){
    cout << "\tT Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void cnot_circuit_test(){
    cout << "\tCNOT Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void cz_circuit_test(){
    cout << "\tCZ Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void swqp_circuit_test(){
    cout << "\tSWAP Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void toffoli_circuit_test(){
    cout << "\tToffoli Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void fredkin_circuit_test(){
    cout << "\tFredkin Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void bell_circuit_test(){
    cout << "\tBell Test..." << endl;
    //TODO: Finish
    cout << "\tDone!" << endl;
}

void test_single_gate_circuits(){
    cout << "Single gate circuit tests..." << endl;
    x_circuit_test();
    y_circuit_test();
    z_circuit_test();
    h_circuit_test();
    rx_circuit_test();
    ry_circuit_test();
    rz_circuit_test();
    s_circuit_test();
    t_circuit_test();
    cnot_circuit_test();
    cz_circuit_test();
    swqp_circuit_test();
    toffoli_circuit_test();
    fredkin_circuit_test();
    bell_circuit_test();
    cout << "Done!" << endl;
}

void test_multi_gate_circuits(){
    cout << "Multi-gate circuit tests..." << endl;
    //TODO: Finish
    cout << "Done!" << endl;
}

int main(){
    try {
        test_constructors();
        test_assignment_operators();
        test_gate_behaviors();
        test_single_gate_circuits();
        test_multi_gate_circuits();
        return 0;
    } catch(exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Unknown exception." << endl;
        return 2;
    }
}