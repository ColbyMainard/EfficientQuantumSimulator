#pragma once

#include <math.h>
#include <cmath>
#include "Matrix.h"
#include <exception>

template<typename T> struct Qubit{
    private:
    UniqueComplexMatrixPtr<T> hilbert_space;
    
    protected:
    void normalize(){
        std::complex<T> total = this->magnitude();
        hilbert_space.get()->setElement(0, 0, hilbert_space.get()->getElement(0, 0) / total);
        hilbert_space.get()->setElement(1, 0, hilbert_space.get()->getElement(1, 0) / total);
    }
    
    public:
    Qubit(){
        hilbert_space = UniqueComplexMatrixPtr<T>(new ComplexMatrix<T>(2, 1));
        hilbert_space.get()->setElement(0, 0, std::complex<T>(1));
        hilbert_space.get()->setElement(1, 0, std::complex<T>(0));
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
        std::complex<T> zero_val = this->hilbert_space.get()->getElement(0, 0);
        std::complex<T> one_val = this->hilbert_space.get()->getElement(1, 0);
        this->hilbert_space.get()->setElement(0, 0, std::complex<T>(0));
        this->hilbert_space.get()->setElement(1, 0, std::complex<T>(1));
    }
    
    void turnOff(){
        hilbert_space.get()->setElement(0, 0, std::complex<T>(1));
        hilbert_space.get()->setElement(1, 0, std::complex<T>(0));
    }
    
    const ComplexMatrix<T>& getState(){
        return *this->hilbert_space;
    }
    
    void setState(const ComplexMatrix<T>& state){
        if(state.getRows() != 2){
            throw std::runtime_error("Incorrect number of rows.");
        }
        if(state.getCols() != 1){
            throw std::runtime_error("Incorrect number of cols.");
        }
        this->hilbert_space.get()->setElement(0, 0, state.getElement(0, 0));
        this->hilbert_space.get()->setElement(1, 0, state.getElement(1, 0));
        normalize();
    }

    T magnitude(){
        std::complex<T> val_1 = hilbert_space.get()->getElement(0, 0);
        val_1 *= val_1;
        std::complex<T> val_2 = hilbert_space.get()->getElement(1, 0);
        val_2 *= val_2;
        return sqrt(abs(val_1 + val_2));
    }
};

template<typename T> std::ostream& operator<<(std::ostream& os, const Qubit<T>& bit){
    os << bit.getState();
    return os;
}

template<typename T> std::istream& operator<<(std::istream& is, Qubit<T>& bit){
    ComplexMatrix<T> matrix(2, 1);
    is >> matrix;
    bit.setState(matrix);
    return is;
}