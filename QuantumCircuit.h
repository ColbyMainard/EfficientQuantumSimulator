#pragma once

#include "Matrix.h"
#include "Qubit.h"
#include <exception>

template<typename T> ComplexMatrix<T> X_gate(){}

template<typename T> ComplexMatrix<T> Y_gate(){}

template<typename T> ComplexMatrix<T> Z_gate(){}

template<typename T> ComplexMatrix<T> Rx_gate(T theta){}

template<typename T> ComplexMatrix<T> Ry_gate(T theta){}

template<typename T> ComplexMatrix<T> Rz_gate(T lambda){}

template<typename T> ComplexMatrix<T> S_gate(){}

template<typename T> ComplexMatrix<T> T_gate(){}

template<typename T> ComplexMatrix<T> CNOT_gate(){}

template<typename T> ComplexMatrix<T> SWAP_gate(){}

template<typename T> ComplexMatrix<T> CZ_gate(){}

template<typename T> ComplexMatrix<T> Toffoli_gate(){}

template<typename T> ComplexMatrix<T> Fredkin_gate(){}

template<typename T> ComplexMatrix<T> Bell_gate(){}

template<typename T> struct QuantumCircuit{
    private:
    vector<UniqueComplexMatrixPtr<T>> qubits;
    vector<UniqueComplexMatrixPtr<T>> gate_matrices;
    vector<vector<unsigned long long int>> qubit_indices;
    protected:
    void runGate(SharedComplexMatrixPtr<T>& gate_matrix, vector<unsigned long long int> qubits_involved);
    public:
    QuantumCircuit(){}
    QuantumCircuit(const QuantumCircuit& other){}
    QuantumCircuit(QuantumCircuit& other){}
    QuantumCircuit operator=(const QuantumCircuit& other){}
    QuantumCircuit operator=(QuantumCircuit& other){}
    void addGate(){}
    void sampleCircuit(){}
};