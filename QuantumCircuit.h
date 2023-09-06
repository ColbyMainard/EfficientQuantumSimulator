#pragma once

#include "Matrix.h"
#include "Qubit.h"
#include <random>
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
    ans.setElement(1, 1, std::complex<T>(1));
    ans.setElement(2, 2, std::complex<T>(1));
    ans.setElement(3, 3, std::complex<T>(-1));
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

template<typename T> struct QuantumCircuit{
    private:
    unsigned long long int num_qubits; // the total number of qubits that should be in the circuit
    std::vector<Qubit<T>> qubits; // list qubits in martrix
    std::vector<UniqueComplexMatrixPtr<T>> gate_matrices; // matrices representing gates to be applied
    std::vector<std::vector<unsigned long long int>> qubit_indices; // indices associated with gates contained in gate_matrices
    
    protected:

    // Setting all qubits to an off state
    void turnAllQubitsOff(){
        for(unsigned long long int i = 0; i < num_qubits; ++i){
            qubits[i].turnOff();
        }
    }

    // Returns a matrix representing transformation, as well as updating qubit internal states for next gate
    ComplexMatrix<T> runGate(unsigned long long int gate_index){
        ComplexMatrix<T> gate_matrix((const ComplexMatrix<T>&) *(gate_matrices.at(gate_index)));
        std::vector<unsigned long long int> qubits_involved = qubit_indices[gate_index];
        // input validation
        if(qubits_involved.size() > 3){
            throw std::runtime_error("Gate is too large, and is not supported.");
        }
        // convert qubit state to larger Hilbert space (big-endian)
        ComplexMatrix<T> start(pow(2, qubits_involved.size()), 1);
        if(qubits_involved.size() == 1){
            start.setElement(0, 0, qubits[qubits_involved[0]].getState().getElement(0, 0));
            start.setElement(1, 0, qubits[qubits_involved[0]].getState().getElement(1, 0));
        }
        if(qubits_involved.size() == 2){
            start.setElement(0, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0));
            start.setElement(1, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0));
            start.setElement(2, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0));
            start.setElement(3, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0));
        }
        if(qubits_involved.size() == 3){
            start.setElement(0, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0) * qubits[qubits_involved[2]].getState().getElement(0, 0));
            start.setElement(1, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0) * qubits[qubits_involved[2]].getState().getElement(1, 0));
            start.setElement(2, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0) * qubits[qubits_involved[2]].getState().getElement(0, 0));
            start.setElement(3, 0, qubits[qubits_involved[0]].getState().getElement(0, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0) * qubits[qubits_involved[2]].getState().getElement(1, 0));
            start.setElement(4, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0) * qubits[qubits_involved[2]].getState().getElement(0, 0));
            start.setElement(5, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(0, 0) * qubits[qubits_involved[2]].getState().getElement(1, 0));
            start.setElement(6, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0) * qubits[qubits_involved[2]].getState().getElement(0, 0));
            start.setElement(7, 0, qubits[qubits_involved[0]].getState().getElement(1, 0) * qubits[qubits_involved[1]].getState().getElement(1, 0) * qubits[qubits_involved[2]].getState().getElement(1, 0));
        }
        //perform calculation of final Hilbert space
        ComplexMatrix<T> ans = gate_matrix * start;
        // convert final Hilbert space to qubit final states (big-endian)
        if(qubits_involved.size() == 1){
            ComplexMatrix<T> qubit_1_state(2, 1);
            qubit_1_state.setElement(0, 0, ans.getElement(0, 0));
            qubit_1_state.setElement(1, 0, ans.getElement(1, 0));
            this->qubits[qubits_involved[0]].setState(qubit_1_state);
        }
        if(qubits_involved.size() == 2){
            ComplexMatrix<T> qubit_1_state(2, 1);
            qubit_1_state.setElement(0, 0, ans.getElement(0, 0) + ans.getElement(1, 0));
            qubit_1_state.setElement(1, 0, ans.getElement(2, 0) + ans.getElement(3, 0));
            this->qubits[qubits_involved[0]].setState(qubit_1_state);
            ComplexMatrix<T> qubit_2_state(2, 1);
            qubit_2_state.setElement(0, 0, ans.getElement(0, 0) + ans.getElement(2, 0));
            qubit_2_state.setElement(1, 0, ans.getElement(1, 0) + ans.getElement(3, 0));
            this->qubits[qubits_involved[1]].setState(qubit_2_state);
        }
        if(qubits_involved.size() == 3){
            ComplexMatrix<T> qubit_1_state(2, 1);
            qubit_1_state.setElement(0, 0, ans.getElement(0, 0) + ans.getElement(1, 0) + ans.getElement(2, 0) + ans.getElement(3, 0));
            qubit_1_state.setElement(1, 0, ans.getElement(4, 0) + ans.getElement(5, 0) + ans.getElement(6, 0) + ans.getElement(7, 0));
            this->qubits[qubits_involved[0]].setState(qubit_1_state);
            ComplexMatrix<T> qubit_2_state(2, 1);
            qubit_2_state.setElement(0, 0, ans.getElement(0, 0) + ans.getElement(1, 0) + ans.getElement(4, 0) + ans.getElement(5, 0));
            qubit_2_state.setElement(1, 0, ans.getElement(2, 0) + ans.getElement(3, 0) + ans.getElement(6, 0) + ans.getElement(7, 0));
            this->qubits[qubits_involved[0]].setState(qubit_2_state);
            ComplexMatrix<T> qubit_3_state(2, 1);
            qubit_3_state.setElement(0, 0, ans.getElement(0, 0) + ans.getElement(2, 0) + ans.getElement(4, 0) + ans.getElement(6, 0));
            qubit_3_state.setElement(1, 0, ans.getElement(1, 0) + ans.getElement(3, 0) + ans.getElement(5, 0) + ans.getElement(7, 0));
            this->qubits[qubits_involved[0]].setState(qubit_3_state);
        }
        return ans;
    }

    //Helper method to int_to_string
    std::string unary_int_to_string(unsigned long long int val){
        switch(val){
            case 0:
                return "0";
            case 1:
                return "1";
            default:
                [[unlikely]]
                throw std::runtime_error("Value not recognized.  Please rewrite algorithm to use, unary, binary, and / or ternary gates.");
        }
    }

    //Helper method to int_to_string
    std::string binary_int_to_string(unsigned long long int val){
        switch(val){
            case 0:
                return "00";
            case 1:
                return "01";
            case 2:
                return "10";
            case 3:
                return "11";
            default:
                [[unlikely]]
                throw std::runtime_error("Value not recognized.  Please rewrite algorithm to use, unary, binary, and / or ternary gates.");
        }
    }

    //Helper method to int_to_string
    std::string trinary_int_to_string(unsigned long long int val){
        switch (val) {
            case 0:
                return "000";
            case 1:
                return "001";
            case 2:
                return "010";
            case 3:
                return "011";
            case 4:
                return "100";
            case 5:
                return "101";
            case 6:
                return "110";
            case 7:
                return "111";
            default:
                [[unlikely]]
                throw std::runtime_error("Value not recognized.  Please rewrite algorithm to use, unary, binary, and / or ternary gates.");
        }
    }

    //Convert integer to big-endian binary stored in string
    std::string int_to_string(unsigned long long int val, unsigned long long int gate_rows){
        switch(gate_rows){
            case 2:
                [[likely]]
                return unary_int_to_string(val);
            case 4:
                return binary_int_to_string(val);
            case 8:
                [[unlikely]]
                return trinary_int_to_string(val);
            default:
                [[unlikely]]
                throw std::runtime_error("Gate size not recognized.  Please rewrite algorithm to use, unary, binary, and / or ternary gates.");
        }
    }

    // Collect single sample of gate
    std::string sampleGate(const ComplexMatrix<T>& results){
        // set up boundaries between probabilities of outputs
        std::vector<double> boundaries = {0};
        double boundary = 0;
        for(unsigned long long int i = 0; i < results.getRows(); ++i){
            auto val = std::abs(results.getElement(0, i));
            boundary += val;
            boundaries.push_back(boundary);
        }
        // generate random value in range
        double lower_bound = 0;
        double upper_bound = boundaries[boundaries.size() - 1];
        std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
        std::default_random_engine re;
        double sample = unif(re);
        // compare random value in range to boundaries
        int val = 0;
        for(unsigned long long int i = 0; i < boundaries.size() - 1; ++i){
            if(sample > boundaries[i] && sample <= boundaries[i + 1]){
                return int_to_string(i, results.getRows());
            }
        }
        throw std::runtime_error("No value determined for gate sample.");
    }
    
    public:
    // Default constructor
    QuantumCircuit():
        num_qubits(0),
        qubits(),
        gate_matrices(),
        qubit_indices()
    {
        //all necessary operations handled by initializer lists
    }
    
    // Custom constructor
    QuantumCircuit(unsigned long long int num_qubits):
        num_qubits(num_qubits),
        qubits(),
        gate_matrices(),
        qubit_indices()
    {
        //initialize qubit list
        for(unsigned long long int i = 0; i < num_qubits; ++i) {
            qubits.push_back(Qubit<T>());
        }
    }
    
    // Copy constructor
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

    // Move constructor
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

    // Destructor
    ~QuantumCircuit(){
        // pop back elements from all vectors
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
    
    // Copy constructor
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
    
    // Move constructor
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
    
    // Pushes copies of gates onto gates
    void getGates(std::vector<ComplexMatrix<T>>& gates){
        for(unsigned long long int i = 0; i < gate_matrices.size(); ++i){
            const ComplexMatrix<T>& gate_to_copy = *gate_matrices[i].get();
            gates.push_back(gate_to_copy);
        }
    }
    
    // returns a copy of the indices corresponding to the gates in the circuit
    const std::vector<std::vector<unsigned long long int>>& getIndices(){
        return qubit_indices;
    }

    // for querying the size of the
    const unsigned long long int numQubits(){
        return this->num_qubits;
    }

    // allows users to add gates of preferred properties
    void addGate(const ComplexMatrix<T>& mat, const std::vector<unsigned long long int>& indices){
        gate_matrices.push_back(UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(mat)));
        qubit_indices.push_back(indices);
    }

    // collects samples and returns them in a map
    const std::unordered_map<std::string, unsigned long long int>& sampleCircuit(unsigned long long int samples){
        std::unordered_map<std::string, unsigned long long int> ans;
        std::vector<ComplexMatrix<T>> gate_results;
        //make sure to start in consistent state
        turnAllQubitsOff();
        // Calculate before and after of each gate
        for(unsigned long long int i = 0; i < qubit_indices.size(); ++i){
            gate_results.push_back(runGate(i));
        }
        // for each sample to be taken
        for(unsigned long long int i = 0; i < samples; ++i){
            //make sure to start in consistent state
            turnAllQubitsOff();
            // for each gate in the circuit
            for(unsigned long long int j = 0; j < gate_results.size(); ++j){
                //propagate results of gate through circuit, storing intermediate results in appropriate qubits
                std::string gate_sample = sampleGate(gate_results[j]);
                for(unsigned long long int bit = 0; bit < gate_sample.length(); ++bit){
                    if(gate_sample[bit] == '0'){
                        this->qubits[this->qubit_indices[j][bit]].turnOff();
                    } else if(gate_sample[bit] == '1'){
                        this->qubits[this->qubit_indices[j][bit]].turnOn();
                    } else{
                        [[unlikely]]
                        throw std::runtime_error("Could not identify value to set candidate qubit to.  Please check configuration of circuit.");
                    }
                }
            }
            // Turn end results into string and add it to map of results and count occurrences
            std::string sample = "";
            for(unsigned long long int index = 0; index < qubits.size(); ++index){
                ComplexMatrix<T> state = (const ComplexMatrix<T>&) qubits[i].getState();
                if(std::abs(state.getElement(0, 0)) >= 0.5){
                    sample += '0';
                } else if(std::abs(state.getElement(1, 0)) >= 0.5){
                    sample += '1';
                } else{
                    [[unlikely]]
                    throw std::runtime_error("Qubit sample failed.");
                }
            }
            ++ans[sample];
        }
        return ans;
    }
};