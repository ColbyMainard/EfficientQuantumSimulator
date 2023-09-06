#include "../QuantumCircuit.h"

// #define CATCH_CONFIG_MAIN
// #include "catch.hpp"

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <random>

TEST_CASE("Benchmark constructors", "[constructors]"){
    std::cout << "Test constructors..." << std::endl;
    std::cout << "\tDefault constructor..." << std::endl;
    QuantumCircuit<float> default_con_float;
    QuantumCircuit<double> default_con_double;
    std::cout << "\tInt constructor..." << std::endl;
    QuantumCircuit<float> int_con_float(3);
    QuantumCircuit<double> int_con_double(3);
    std::cout << "\tAdding gates..." << std::endl;
    std::vector<std::vector<unsigned long long int>> indices = {{0}, {1}, {2}, {0}, {1}, {2}};
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
    std::cout << "\tCopy constructor..." << std::endl;
    QuantumCircuit<float> copy_float((const QuantumCircuit<float>&) int_con_float);
    QuantumCircuit<double> copy_double((const QuantumCircuit<double>&) int_con_double);
    std::cout << "\tMove constructor..." << std::endl;
    QuantumCircuit<float> move_float((QuantumCircuit<float>&) copy_float);
    QuantumCircuit<double> move_double((QuantumCircuit<double>&) copy_double);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark assignment operators", "[assignment]"){
    std::cout << "Test assignment operators..." << std::endl;
    std::cout << "\tDefault constructor..." << std::endl;
    QuantumCircuit<float> default_con_float;
    QuantumCircuit<double> default_con_double;
    std::cout << "\tInt constructor..." << std::endl;
    QuantumCircuit<float> int_con_float(3);
    QuantumCircuit<double> int_con_double(3);
    std::cout << "\tAdding gates..." << std::endl;
    std::vector<std::vector<unsigned long long int>> indices = {{0}, {1}, {2}, {0}, {1}, {2}};
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
    std::cout << "\tCopy constructor..." << std::endl;
    QuantumCircuit<float> copy_float = (const QuantumCircuit<float>&) int_con_float;
    QuantumCircuit<double> copy_double = (const QuantumCircuit<double>&) int_con_double;
    std::cout << "\tMove constructor..." << std::endl;
    QuantumCircuit<float> move_float = (QuantumCircuit<float>&) copy_float;
    QuantumCircuit<double> move_double = (QuantumCircuit<double>&) copy_double;
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark X Matrix", "[x]"){
    std::cout << "X Test..." << std::endl;
    auto gate = X_gate<double>();
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(0, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Y Matrix", "[y]"){
    std::cout << "Y Test..." << std::endl;
    auto gate = Y_gate<double>();
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(1, 0) * std::complex<double>(0, 1)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(0, 0) * std::complex<double>(0, 1)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Z Matrix", "[z]"){
    std::cout << "Z Test..." << std::endl;
    auto gate = Z_gate<double>();
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) + random_val.getElement(1, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark H Matrix", "[h]"){
    std::cout << "H Test..." << std::endl;
    auto gate = H_gate<double>();
    std::cout << "H Gate:" << gate << std::endl;
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    auto val_1 = (random_val.getElement(0,0) + random_val.getElement(1, 0)) * std::complex<double>(1.0 / sqrt(2));
    REQUIRE(abs(ans.getElement(0, 0) - val_1) < 0.001);
    auto val_2 = (random_val.getElement(0,0) - random_val.getElement(1, 0)) * std::complex<double>(1.0 / sqrt(2));
    REQUIRE(abs(ans.getElement(1, 0) - val_2) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Rx Matrix", "[rx]"){
    std::cout << "Rx Test..." << std::endl;
    auto theta = 0.2;
    auto gate = Rx_gate<double>(0.2);
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    std::complex<double> val_1 = (random_val.getElement(0, 0) * std::complex<double>(cos(theta / 2))) - (random_val.getElement(1, 0) * std::complex<double>(sin(theta / 2)) * std::complex<double>(0, 1));
    if(abs(ans.getElement(0, 0) - val_1) > 0.001){
        throw std::runtime_error("Mismatch of expected value 0.");
    }
    std::complex<double> val_2 = (random_val.getElement(1, 0) * std::complex<double>(cos(theta / 2))) - (random_val.getElement(0, 0) * std::complex<double>(sin(theta / 2)) * std::complex<double>(0, 1));
    if(abs(ans.getElement(1, 0) - val_2) > 0.001){
        throw std::runtime_error("Mismatch of expected value 1.");
    }
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Ry Matrix", "[ry]"){
    std::cout << "Ry Test..." << std::endl;
    double theta = 0.2;
    auto gate = Ry_gate<double>(theta);
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    std::complex<double> val_1 = (random_val.getElement(0, 0) * std::complex<double>(cos(theta / 2))) - (random_val.getElement(1, 0) * std::complex<double>(sin(theta / 2)));
    REQUIRE(abs(ans.getElement(0, 0) - val_1) < 0.001);
    std::complex<double> val_2 = (random_val.getElement(0, 0) * std::complex<double>(sin(theta / 2))) + (random_val.getElement(1, 0) * std::complex<double>(cos(theta / 2)));
    REQUIRE(abs(ans.getElement(1, 0) - val_2) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Rz Matrix", "[rz]"){
    std::cout << "Rz Test..." << std::endl;
    double phi = 0.2;
    auto gate = Rz_gate<double>(phi);
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    std::complex<double> val = random_val.getElement(1, 0) * exp(std::complex<double>(0, phi));
    REQUIRE(abs(ans.getElement(1, 0) - val) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark S Matrix", "[s]"){
    std::cout << "S Test..." << std::endl;
    auto gate = S_gate<double>();
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - (random_val(1, 0) * std::complex<double>(0, 1))) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark T Matrix", "[t]"){
    std::cout << "T Test..." << std::endl;
    auto gate = T_gate<double>();
    ComplexMatrix<double> random_val(2, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    std::complex<double> val = random_val.getElement(1, 0) * exp(std::complex<double>(0, M_PI / 4));
    REQUIRE(abs(ans.getElement(1, 0) - val) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark CNOT Matrix", "[cnot]"){
    std::cout << "CNOT Test..." << std::endl;
    auto gate = CNOT_gate<double>();
    ComplexMatrix<double> random_val(4, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(2, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(3, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(2, 0) - random_val.getElement(3, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(3, 0) - random_val.getElement(2, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark CZ Matrix", "[cz]"){
    std::cout << "CZ Test..." << std::endl;
    auto gate = CZ_gate<double>();
    ComplexMatrix<double> random_val(4, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(2, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(3, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(2, 0) - random_val.getElement(2, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(3, 0) + random_val.getElement(3, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark SWAP Matrix", "[swap]"){
    std::cout << "SWAP Test..." << std::endl;
    auto gate = SWAP_gate<double>();
    ComplexMatrix<double> random_val(4, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(2, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(3, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(2, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(2, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(3, 0) - random_val.getElement(3, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Toffoli Matrix", "[toffoli]"){
    std::cout << "Toffoli Test..." << std::endl;
    auto gate = Toffoli_gate<double>();
    ComplexMatrix<double> random_val(8, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(2, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(3, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(4, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(5, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(6, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(7, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(2, 0) - random_val.getElement(2, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(3, 0) - random_val.getElement(3, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(4, 0) - random_val.getElement(4, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(5, 0) - random_val.getElement(5, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(6, 0) - random_val.getElement(7, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(7, 0) - random_val.getElement(6, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Fredkin Matrix", "[fredkin]"){
    std::cout << "Fredkin Test..." << std::endl;
    auto gate = Fredkin_gate<double>();
    ComplexMatrix<double> random_val(8, 1);
    random_val.setElement(0, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(1, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(2, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(3, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(4, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(5, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(6, 0, std::complex<double>(rand() % 5, rand() % 5));
    random_val.setElement(7, 0, std::complex<double>(rand() % 5, rand() % 5));
    auto ans = gate * random_val;
    std::cout << "\tBefore: " << random_val << std::endl << "\tAfter: " << ans << std::endl;
    REQUIRE(abs(ans.getElement(0, 0) - random_val.getElement(0, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(1, 0) - random_val.getElement(1, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(2, 0) - random_val.getElement(2, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(3, 0) - random_val.getElement(3, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(4, 0) - random_val.getElement(4, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(5, 0) - random_val.getElement(6, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(6, 0) - random_val.getElement(5, 0)) < 0.001);
    REQUIRE(abs(ans.getElement(7, 0) - random_val.getElement(7, 0)) < 0.001);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark X Circuit", "[x]"){
    std::cout << "X Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(X_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 0);
    REQUIRE(results["1"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Y Circuit", "[y]"){
    std::cout << "Y Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(Y_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 0);
    REQUIRE(results["1"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Z Circuit", "[z]"){
    std::cout << "\tZ Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(Z_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 100);
    REQUIRE(results["1"] == 0);
    std::cout << "\tDone!" << std::endl;
}

TEST_CASE("Benchmark H Circuit", "[h]"){
    std::cout << "H Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(H_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(48 < results["0"]);
    REQUIRE(52 > results["0"]);
    REQUIRE(48 < results["1"]);
    REQUIRE(52 > results["1"]);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Rx Circuit", "[rx]"){
    std::cout << "Rx Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(Rx_gate<double>(.7), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(85 < results["0"]);
    REQUIRE(90 > results["0"]);
    REQUIRE(10 < results["1"]);
    REQUIRE(15 > results["1"]);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Ry Circuit", "[ry]"){
    std::cout << "Ry Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(Ry_gate<double>(.7), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(85 < results["0"]);
    REQUIRE(90 > results["0"]);
    REQUIRE(10 < results["1"]);
    REQUIRE(15 > results["1"]);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Rz Circuit", "[rz]"){
    std::cout << "Rz Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(Rz_gate<double>(.7), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark S Circuit", "[s]"){
    std::cout << "S Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(S_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark T Circuit", "[t]"){
    std::cout << "T Test..." << std::endl;
    QuantumCircuit<double> test(1);
    test.addGate(T_gate<double>(), {0});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["0"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark CNOT Circuit", "[cnot]"){
    std::cout << "\tCNOT Test..." << std::endl;
    QuantumCircuit<double> test(2);
    test.addGate(CNOT_gate<double>(), {0, 1});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["00"] == 100);
    std::cout << "\tDone!" << std::endl;
}

TEST_CASE("Benchmark CZ Circuit", "[cz]"){
    std::cout << "CZ Test..." << std::endl;
    QuantumCircuit<double> test(2);
    test.addGate(CZ_gate<double>(), {0, 1});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["00"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark SWAP Circuit", "[swap]"){
    std::cout << "SWAP Test..." << std::endl;
    QuantumCircuit<double> test(2);
    test.addGate(SWAP_gate<double>(), {0, 1});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["00"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Toffoli Circuit", "[toffoli]"){
    std::cout << "Toffoli Test..." << std::endl;
    QuantumCircuit<double> test(3);
    test.addGate(Toffoli_gate<double>(), {0, 1, 2});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["000"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Fredkin Circuit", "[fredkin]"){
    std::cout << "Fredkin Test..." << std::endl;
    QuantumCircuit<double> test(3);
    test.addGate(Fredkin_gate<double>(), {0, 1, 2});
    auto results = test.sampleCircuit(100);
    REQUIRE(results["000"] == 100);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Bell Circuit", "[bell]"){
    std::cout << "Bell circuit test..." << std::endl;
    QuantumCircuit<double> test(2);
    test.addGate(H_gate<double>(), {0});
    test.addGate(CNOT_gate<double>(), {0, 1});
    auto results = test.sampleCircuit(100);
    REQUIRE(48 < results["00"]);
    REQUIRE(52 > results["00"]);
    REQUIRE(48 < results["11"]);
    REQUIRE(52 > results["11"]);
    std::cout << "Done!" << std::endl;
}

TEST_CASE("Benchmark Random Circuit", "[random]"){
    std::cout << "Nonsense circuit test..." << std::endl;
    QuantumCircuit<double> test(2);
    test.addGate(H_gate<double>(), {0});
    test.addGate(H_gate<double>(), {1});
    test.addGate(H_gate<double>(), {2});
    test.addGate(H_gate<double>(), {3});
    test.addGate(CNOT_gate<double>(), {0, 1});
    test.addGate(CNOT_gate<double>(), {2, 3});
    test.addGate(X_gate<double>(), {0});
    test.addGate(Y_gate<double>(), {1});
    test.addGate(Z_gate<double>(), {2});
    test.addGate(Rx_gate<double>(0.7), {3});
    test.addGate(Ry_gate<double>(0.7), {0});
    test.addGate(Rz_gate<double>(0.7), {1});
    test.addGate(S_gate<double>(), {2});
    test.addGate(T_gate<double>(), {3});
    test.addGate(CZ_gate<double>(), {0, 1});
    test.addGate(CZ_gate<double>(), {2, 3});
    test.addGate(Toffoli_gate<double>(), {0, 1, 2});
    test.addGate(Fredkin_gate<double>(), {1, 2, 3});
    auto results = test.sampleCircuit(100);
    REQUIRE(1 < results["0000"]);
    REQUIRE(4 > results["0000"]);
    REQUIRE(1 < results["0001"]);
    REQUIRE(4 > results["0001"]);
    REQUIRE(1 < results["0010"]);
    REQUIRE(4 > results["0010"]);
    REQUIRE(1 < results["0011"]);
    REQUIRE(4 > results["0011"]);
    REQUIRE(1 < results["0100"]);
    REQUIRE(4 > results["0100"]);
    REQUIRE(1 < results["0101"]);
    REQUIRE(4 > results["0101"]);
    REQUIRE(1 < results["0110"]);
    REQUIRE(4 > results["0110"]);
    REQUIRE(1 < results["0111"]);
    REQUIRE(4 > results["0111"]);
    REQUIRE(8 < results["1000"]);
    REQUIRE(12 > results["1000"]);
    REQUIRE(8 < results["1001"]);
    REQUIRE(12 > results["1001"]);
    REQUIRE(8 < results["1010"]);
    REQUIRE(12 > results["1010"]);
    REQUIRE(8 < results["1011"]);
    REQUIRE(12 > results["1011"]);
    REQUIRE(8 < results["1100"]);
    REQUIRE(12 > results["1100"]);
    REQUIRE(8 < results["1101"]);
    REQUIRE(12 > results["1101"]);
    REQUIRE(8 < results["1110"]);
    REQUIRE(12 > results["1110"]);
    REQUIRE(8 < results["1111"]);
    REQUIRE(12 > results["1111"]);
    std::cout << "Done!" << std::endl;
}

int main( int argc, char* argv[] ) {
    Catch::Session session;
    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 ) { return returnCode; }
    int numFailed = session.run();
    return numFailed;
}