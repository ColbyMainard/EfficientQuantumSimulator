cmake_minimum_required(VERSION 3.16)
project(simple_example)
set(CMAKE_CXX_STANDARD 20)

find_package(Catch2 REQUIRED)

SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fpermissive") # compiles with debug symbols
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fpermissive") # compiles without debug symbols

add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

enable_testing()
include(CTest)
include(Catch)

add_executable(matrix_tests.bin test_matrix.cpp)
#target_link_libraries(matrix_tests.bin PRIVATE Catch2::Catch2WithMain)
target_link_libraries(matrix_tests.bin PRIVATE Catch2::Catch2)
catch_discover_tests(matrix_tests.bin)
add_executable(qubit_tests.bin test_qubit.cpp)
#target_link_libraries(qubit_tests.bin PRIVATE Catch2::Catch2WithMain)
target_link_libraries(qubit_tests.bin PRIVATE Catch2::Catch2)
catch_discover_tests(qubit_tests.bin)
add_executable(quantum_circuit_tests.bin test_quantum_circuit.cpp)
#target_link_libraries(quantum_circuit_tests.bin PRIVATE Catch2::Catch2WithMain)
target_link_libraries(quantum_circuit_tests.bin PRIVATE Catch2::Catch2)
catch_discover_tests(quantum_circuit_tests.bin)