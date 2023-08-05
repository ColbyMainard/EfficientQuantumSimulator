cmake_minimum_required(VERSION 3.16)
project(simple_example)
set(CMAKE_CXX_STANDARD 20)

SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fpermissive")

add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

add_executable(matrix_tests.bin test_matrix.cpp)
add_executable(qubit_tests.bin test_qubit.cpp)
add_executable(quantum_circuit_tests.bin test_quantum_circuit.cpp)