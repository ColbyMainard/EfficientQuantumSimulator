#pragma once

#include "Matrix.h"
#include "Qubit.h"
#include <exception>
#include <complex>
#include <unordered_map>
#include <string>

template<typename T> ComplexMatrix<T> X_gate(){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(0));
    ans.setElement(0, 1, std::complex<T>(1));
    ans.setElement(1, 0, std::complex<T>(1));
    ans.setElement(1, 1, std::complex<T>(0));
    return ans;
}

template<typename T> ComplexMatrix<T> Y_gate(){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(0));
    ans.setElement(0, 1, std::complex<T>(0, 1));
    ans.setElement(1, 0, std::complex<T>(0, 1));
    ans.setElement(1, 1, std::complex<T>(0));
    return ans;
}

template<typename T> ComplexMatrix<T> Z_gate(){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(0, 1, std::complex<T>(0));
    ans.setElement(1, 0, std::complex<T>(0));
    ans.setElement(1, 1, std::complex<T>(-1));
    return ans;
}

template<typename T> ComplexMatrix<T> H_gate(){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(1.0 / sqrt(2)));
    ans.setElement(0, 1, std::complex<T>(1.0 / sqrt(2)));
    ans.setElement(1, 0, std::complex<T>(1.0 / sqrt(2)));
    ans.setElement(1, 1, std::complex<T>(-1.0 / sqrt(2)));
    return ans;
}

template<typename T> ComplexMatrix<T> Rx_gate(T theta){
    ComplexMatrix<T> ans(2, 2);
    std::complex<T> cos_val(std::cos(theta / 2));
    std::complex<T> sin_val(0, -1 * std::sin(theta / 2));
    ans.setElement(0, 0, cos_val);
    ans.setElement(0, 1, sin_val);
    ans.setElement(1, 0, sin_val);
    ans.setElement(1, 1, cos_val);
    return ans;
}

template<typename T> ComplexMatrix<T> Ry_gate(T theta){
    ComplexMatrix<T> ans(2, 2);
    std::complex<T> cos_val(std::cos(theta / 2));
    std::complex<T> sin_val(std::sin(theta / 2));
    ans.setElement(0, 0, cos_val);
    ans.setElement(0, 1, std::complex<T>(-1) * sin_val);
    ans.setElement(1, 0, sin_val);
    ans.setElement(1, 1, cos_val);
    return ans;
}

template<typename T> ComplexMatrix<T> Rz_gate(T lambda){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(0, 1, std::complex<T>(0));
    ans.setElement(1, 0, std::complex<T>(0));
    auto power_val = std::complex<T>(0, lambda);
    ans.setElement(1, 1, std::exp(power_val));
    return ans;
}

template<typename T> ComplexMatrix<T> S_gate(){
    ComplexMatrix<T> ans(2, 2);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(0, 1, std::complex<T>(0));
    ans.setElement(1, 0, std::complex<T>(0));
    ans.setElement(1, 1, std::complex<T>(0, 1));
    return ans;
}

template<typename T> ComplexMatrix<T> T_gate(){
    return Rz_gate<T>(M_PI / 4);
}

template<typename T> ComplexMatrix<T> CNOT_gate(){
    ComplexMatrix<T> ans(4, 4);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(1, 1, std::complex<T>(1));
    ans.setElement(2, 3, std::complex<T>(1));
    ans.setElement(3, 2, std::complex<T>(1));
    return ans;
}

template<typename T> ComplexMatrix<T> CZ_gate(){
    ComplexMatrix<T> ans(4, 4);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(1, 2, std::complex<T>(1));
    ans.setElement(2, 1, std::complex<T>(1));
    ans.setElement(3, 3, std::complex<T>(-11));
    return ans;
}

template<typename T> ComplexMatrix<T> SWAP_gate(){
    ComplexMatrix<T> ans(4, 4);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(1, 2, std::complex<T>(1));
    ans.setElement(2, 1, std::complex<T>(1));
    ans.setElement(3, 3, std::complex<T>(1));
    return ans;
}

template<typename T> ComplexMatrix<T> Toffoli_gate(){
    ComplexMatrix<T> ans(8, 8);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(1, 1, std::complex<T>(1));
    ans.setElement(2, 2, std::complex<T>(1));
    ans.setElement(3, 3, std::complex<T>(1));
    ans.setElement(4, 4, std::complex<T>(1));
    ans.setElement(5, 5, std::complex<T>(1));
    ans.setElement(6, 7, std::complex<T>(1));
    ans.setElement(7, 6, std::complex<T>(1));
    return ans;
}

template<typename T> ComplexMatrix<T> Fredkin_gate(){
    ComplexMatrix<T> ans(8, 8);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(1, 1, std::complex<T>(1));
    ans.setElement(2, 2, std::complex<T>(1));
    ans.setElement(3, 3, std::complex<T>(1));
    ans.setElement(4, 4, std::complex<T>(1));
    ans.setElement(5, 6, std::complex<T>(1));
    ans.setElement(6, 5, std::complex<T>(1));
    ans.setElement(7, 7, std::complex<T>(1));
    return ans;
}

template<typename T> ComplexMatrix<T> Bell_gate(){
    ComplexMatrix<T> ans(4, 4);
    ans.setElement(0, 0, std::complex<T>(1));
    ans.setElement(0, 2, std::complex<T>(1));
    ans.setElement(1, 1, std::complex<T>(1));
    ans.setElement(1, 3, std::complex<T>(1));
    ans.setElement(2, 1, std::complex<T>(1));
    ans.setElement(2, 3, std::complex<T>(-1));
    ans.setElement(3, 0, std::complex<T>(1));
    ans.setElement(3, 2, std::complex<T>(-1));
    return ans;
}

template<typename T> struct QuantumCircuit{
    private:
    unsigned long long int num_qubits;
    std::vector<Qubit<T>> qubits;
    std::vector<UniqueComplexMatrixPtr<T>> gate_matrices;
    std::vector<std::vector<unsigned long long int>> qubit_indices;
    
    protected:
    void runGate(UniqueComplexMatrixPtr<T>& gate_matrix, std::vector<unsigned long long int> qubits_involved){}
    
    public:
    QuantumCircuit():
        num_qubits(0),
        qubits(),
        gate_matrices(),
        qubit_indices()
    {}
    
    QuantumCircuit(unsigned long long int num_qubits):
        num_qubits(num_qubits),
        qubits(),
        gate_matrices(),
        qubit_indices()
    {
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
    }
    
    QuantumCircuit(const QuantumCircuit& other){
        this->num_qubits = other.numQubits();
        this->qubits = std::vector<Qubit<T>>();
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
        gate_matrices = std::vector<UniqueComplexMatrixPtr<T>>();
        std::vector<ComplexMatrix<T>> other_gates;
        other.getGates(other_gates);
        for(unsigned long long int i = 0; i < other_gates.size(); ++i){
            gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(((const ComplexMatrix<T>&) other_gates[i]))));
        }
        auto val = getIndices();
        std::copy(val.begin(), val.end(), qubit_indices.begin());
    }
    
    QuantumCircuit(QuantumCircuit& other)
    {
        this->num_qubits = other.numQubits();
        this->qubits = std::vector<Qubit<T>>();
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
        gate_matrices = std::vector<UniqueComplexMatrixPtr<T>>();
        std::vector<ComplexMatrix<T>> other_gates;
        other.getGates(other_gates);
        for(unsigned long long int i = 0; i < other_gates.size(); ++i){
            gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(((const ComplexMatrix<T>&) other_gates[i]))));
        }
        auto val = getIndices();
        std::copy(val.begin(), val.end(), qubit_indices.begin());
    }

    ~QuantumCircuit(){
        while(!qubits.empty()){
            qubits.pop_back();
        }
        while(!gate_matrices.empty()){
            gate_matrices.pop_back();
        }
        while(!qubit_indices.empty()){
            qubit_indices.pop_back();
        }
        num_qubits = 0;
    }
    
    QuantumCircuit operator=(const QuantumCircuit& other){
        this->num_qubits = other.numQubits();
        this->qubits = std::vector<Qubit<T>>();
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
        gate_matrices = std::vector<UniqueComplexMatrixPtr<T>>();
        std::vector<ComplexMatrix<T>> other_gates;
        other.getGates(other_gates);
        for(unsigned long long int i = 0; i < other_gates.size(); ++i){
            gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(((const ComplexMatrix<T>&) other_gates[i]))));
        }
        auto val = getIndices();
        std::copy(val.begin(), val.end(), qubit_indices.begin());
    }
    
    QuantumCircuit operator=(QuantumCircuit& other){
        this->num_qubits = other.numQubits();
        this->qubits = std::vector<Qubit<T>>();
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
        gate_matrices = std::vector<UniqueComplexMatrixPtr<T>>();
        std::vector<ComplexMatrix<T>> other_gates;
        other.getGates(other_gates);
        for(unsigned long long int i = 0; i < other_gates.size(); ++i){
            gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(((const ComplexMatrix<T>&) other_gates[i]))));
        }
        auto val = getIndices();
        std::copy(val.begin(), val.end(), qubit_indices.begin());
        delete other;
    }
    
    void getGates(std::vector<ComplexMatrix<T>>& gates){
        for(unsigned long long int i = 0; i < gate_matrices.size(); ++i){
            auto gate_to_copy = gate_matrices[i].get();
            gates.push_back(*gate_to_copy);
        }
    }
    
    std::vector<std::vector<unsigned long long int>> getIndices(){
        return qubit_indices;
    }

    unsigned long long int numQubits(){
        return this->num_qubits;
    }
    
    void addGate(const ComplexMatrix<T>& mat, const std::vector<unsigned long long int>& indices){
        gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(mat)));
        qubit_indices.push_back(indices);
    }
    
    void sampleCircuit(){
        //TODO: Finish
    }
};

template<typename T> std::ostream& operator<<(std::ostream& os, const QuantumCircuit<T>& circ){
    std::vector<ComplexMatrix<T>> gates;
    circ.getGates(gates);
    auto indices = circ.getIndices();
    os << circ.numQubits() << std::endl;
    os << gates.size() << gates << std::endl;
    os << indices.size() << indices << std::endl;
    for(unsigned long long int i = 0; i < gates.size(); ++i){
        os << gates[i] << std::endl << indices[i] << std::endl;
    }
    return os;
}