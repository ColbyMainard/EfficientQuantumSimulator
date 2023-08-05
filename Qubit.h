#pragma once

#include "Matrix.h"
#include <exception>

template<typename T> struct Qubit{
    private:
    UniqueComplexMatrixPtr<T> hilbert_space;
    protected:
    void normalize(){
        complex<T> total = hilbert_space.get()->getElement(0, 0) + hilbert_space.get()->getElement(1, 0);
        hilbert_space.get()->setElement(0, 0, hilbert_space.get()->getElement(0, 0) / total);
        hilbert_space.get()->setElement(1, 0, hilbert_space.get()->getElement(1, 0) / total);
    }
    public:
    Qubit(){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, complex<T>(1));
        hilbert_space.get()->setElement(1, 0, complex<T>(0));
    }
    
    Qubit(const Qubit& other){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, other.getState().getElement(0, 0));
        hilbert_space.get()->setElement(1, 0, other.getState().getElement(1, 0));
    }
    
    Qubit(Qubit& other){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, other.getState().getElement(0, 0));
        hilbert_space.get()->setElement(1, 0, other.getState().getElement(1, 0));
    }
    
    Qubit operator=(const Qubit& other){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, other.getState().getElement(0, 0));
        hilbert_space.get()->setElement(1, 0, other.getState().getElement(1, 0));
    }
    
    Qubit operator=(Qubit& other){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, other.getState().getElement(0, 0));
        hilbert_space.get()->setElement(1, 0, other.getState().getElement(0, 1));
    }

    ~Qubit(){
        hilbert_space.release();
    }
    
    void turnOn(){
        this->turnOff();
        complex<T> zero_val = this->hilbert_space.get()->getElement(0, 0);
        complex<T> one_val = this->hilbert_space.get()->getElement(1, 0);
        this->hilbert_space.get()->setElement(0, 0, complex<T>(0));
        this->hilbert_space.get()->setElement(1, 0, complex<T>(1));
    }
    
    void turnOff(){
        hilbert_space.get()->setElement(0, 0, complex<T>(1));
        hilbert_space.get()->setElement(1, 0, complex<T>(0));
    }
    
    const ComplexMatrix<T>& getState(){
        return *this->hilbert_space;
    }
    
    void setState(const ComplexMatrix<T>& state){
        if(state.getRows() != 2){
            throw runtime_error("Incorrect number of rows.");
        }
        if(state.getCols() != 1){
            throw runtime_error("Incorrect number of cols.");
        }
        this->hilbert_space.get()->setElement(0, 0, state.getElement(0, 0));
        this->hilbert_space.get()->setElement(1, 0, state.getElement(1, 0));
        normalize();
    }
};

template<typename T> ostream& operator<<(ostream& os, const Qubit<T>& bit){
    os << bit.getState();
    return os;
}

template<typename T> istream& operator<<(istream& is, Qubit<T>& bit){
    ComplexMatrix<T> matrix(2, 1);
    is >> matrix;
    bit.setState(matrix);
    return is;
}