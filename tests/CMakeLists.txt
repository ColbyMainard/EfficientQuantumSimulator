cmake_minimum_required(VERSION 3.16)
project(simple_example)
set(CMAKE_CXX_STANDARD 20)

SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fpermissive") # compiles with debug symbols
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fpermissive") # compiles without debug symbols

add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

enable_testing()

add_executable(matrix_tests.bin test_matrix.cpp)
add_test(NAME MatrixTests COMMAND matrix_tests.bin)
add_executable(qubit_tests.bin test_qubit.cpp)
add_test(NAME QubitTests COMMAND qubit_tests.bin)
add_executable(quantum_circuit_tests.bin test_quantum_circuit.cpp)
add_test(NAME QuantumCircuitTests COMMAND quantum_circuit_tests.bin)