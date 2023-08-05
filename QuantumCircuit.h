#pragma once

#include "Matrix.h"
#include "Qubit.h"
#include <exception>

template<typename T> struct QuantumCircuit{
    private:
    vector<SharedComplexMatrixPtr<T>> qubits;
    vector<SharedComplexMatrixPtr<T>> gate_matrices;
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