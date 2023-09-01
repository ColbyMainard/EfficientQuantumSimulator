#pragma once

#include <vector>

#include <memory>

#include <algorithm>

#include <complex>

#include <ostream>
#include <istream>
#include <fstream>
#include <iostream>
#include <strstream>
#include <sstream>

#include <exception>

/*
MISCELLANEOUS COMPARISON OPERATORS
*/

template<typename T> bool operator>(const std::complex<T>& c1, const std::complex<T>& c2){
    return abs(c1) > abs(c2);
}

template<typename T> bool operator<(const std::complex<T>& c1, const std::complex<T>& c2){
    return abs(c1) < abs(c2);
}

template<typename T> bool operator<=(const std::complex<T>& c1, const std::complex<T>& c2){
    return abs(c1) <= abs(c2);
}

template<typename T> bool operator>=(const std::complex<T>& c1, const std::complex<T>& c2){
    return abs(c1) >= abs(c2);
}

template<typename T> bool operator>(const std::complex<T>& c1, const double& d){
    return abs(c1) > d;
}

template<typename T> bool operator<(const std::complex<T>& c1, const double& d){
    return abs(c1) < d;
}

template<typename T> bool operator<=(const std::complex<T>& c1, const double& d){
    return abs(c1) <= d;
}

template<typename T> bool operator>=(const std::complex<T>& c1, const double& d){
    return abs(c1) >= d;
}

template<typename T> bool operator==(const std::complex<T>& c1, const double& d){
    return abs(c1) == d;
}

template<typename T> bool operator>(const double& d, const std::complex<T>& c1){
    return abs(c1) > d;
}

template<typename T> bool operator<(const double& d, const std::complex<T>& c1){
    return abs(c1) < d;
}

template<typename T> bool operator<=(const double& d, const std::complex<T>& c1){
    return abs(c1) <= d;
}

template<typename T> bool operator>=(const double& d, const std::complex<T>& c1){
    return abs(c1) >= d;
}

template<typename T> bool operator==(const double& d, const std::complex<T>& c1){
    return abs(c1) == d;
}

template<typename T> struct Matrix{
    private:
        std::vector<std::vector<T>> elements;
        unsigned long long int rows;
        unsigned long long int cols;
    protected:
    public:
        //default constructor
        Matrix(){
            this->rows = 0;
            this->cols = 0;
            this->elements = std::vector<std::vector<T>>();
        }

        //custom constructor
        Matrix(unsigned long long int rownum, unsigned long long int colnum)
            :rows(rownum),
            cols(colnum),
            elements(std::vector<std::vector<T>>())
        {
            for(unsigned long long int i = 0; i < rownum; ++i){
                elements.push_back(std::vector<T>());
                for(unsigned long long int j = 0; j < colnum; ++j){
                    elements[i].push_back(T(0));
                }
            }
        }

        //copy constructor
        Matrix(const Matrix& other)
            :rows(other.rows),
            cols(other.cols),
            elements(std::vector<std::vector<T>>())
        {
            for(unsigned long long int i = 0; i < this->rows; ++i){
                this->elements.push_back(std::vector<T>());
                for(unsigned long long int j = 0; j < this->cols; ++j){
                    this->elements[i].push_back(T(0));
                }
                std::copy(other.elements[i].begin(), other.elements[i].end(), this->elements[i].begin());
            }
        }

        //move constructor
        Matrix(Matrix& other)
            :rows(other.rows),
            cols(other.cols),
            elements(std::vector<std::vector<T>>())
        {
            for(unsigned long long int i = 0; i < this->rows; ++i){
                this->elements.push_back(std::vector<T>());
                for(unsigned long long int j = 0; j < this->cols; ++j){
                    this->elements[i].push_back(T(0));
                }
                std::copy(other.elements[i].begin(), other.elements[i].end(), this->elements[i].begin());
            }
        }

        //copy assignment
        Matrix operator=(const Matrix& other){
            this->rows = other.getRows();
            this->cols = other.getCols();
            this->elements = std::vector<std::vector<T>>();
            
            for(unsigned long long int i = 0; i < this->rows; ++i){
                this->elements.push_back(std::vector<T>());
                for(unsigned long long int j = 0; j < this->cols; ++j){
                    this->elements[i].push_back(T(0));
                }
                std::copy(other.elements[i].begin(), other.elements[i].end(), this->elements[i].begin());
            }
        }

        //move assignment
        Matrix operator=(Matrix& other){
            this->rows = other.getRows();
            this->cols = other.getCols();
            this->elements = std::vector<std::vector<T>>();
            
            for(unsigned long long int i = 0; i < this->rows; ++i){
                this->elements.push_back(std::vector<T>());
                for(unsigned long long int j = 0; j < this->cols; ++j){
                    this->elements[i].push_back(T(0));
                }
                std::copy(other.elements[i].begin(), other.elements[i].end(), this->elements[i].begin());
            }

            delete other;
        }

        //destructor
        ~Matrix(){
            rows = 0;
            cols = 0;
        }

        //get row count
        const unsigned long long int getRows(){
            return this->rows;
        }

        //get column count
        const unsigned long long int getCols(){
            return this->cols;
        }

        //determinant
        T determinant(){
            if(this->rows != this->cols){
                throw std::runtime_error("Cannot get determinant.");
            }
            Matrix<T> mat_copy((const Matrix&) *this);
            // https://www.geeksforgeeks.org/determinant-of-a-matrix/
            T num1, num2, det(1), total(1);
            unsigned long long int index;
            std::vector<T> temp;
            for(unsigned int i = 0; i < rows; ++i){
                temp.push_back(0);
            }
            // loop for traversing the diagonal elements
            for (unsigned long long int i = 0; i < rows; i++){
                index = i; // initialize the index
                // finding the index which has non zero value
                while (index < rows && mat_copy(index, i) == T(0)){
                    index++;
                }
                if (index == rows) {
                    // if there is non zero element
                    // the determinant of matrix as zero
                    continue;
                }
                if (index != i) {
                    // loop for swapping the diagonal element row and
                    // index row
                    for (unsigned long long int j = 0; j < rows; j++) {
                        T tmp = mat_copy(index, j);
                        mat_copy.setElement(index, j, mat_copy.getElement(i, j));
                        mat_copy.setElement(i, j, tmp);
                    }
                    // determinant sign changes when we shift rows
                    // go through determinant properties
                    det = det * pow(-1, index - i);
                }
                // storing the values of diagonal row elements
                for (unsigned long long int j = 0; j < rows; j++) {
                    temp[j] = mat_copy(i, j);
                }
                // traversing every row below the diagonal element
                for (unsigned long long int j = i + 1; j < rows; j++) {
                    num1 = temp[i]; // value of diagonal element
                    num2 = mat_copy(j, i); // value of next row element

                    // traversing every column of row
                    // and multiplying to every row
                    for (unsigned long long int k = 0; k < rows; k++){
                        // multiplying to make the diagonal
                        // element and next row element equal
                        mat_copy.setElement(j, k, (num1 * mat_copy.getElement(j, k)) - (num2 * temp[k]));
                    }
                    total = total * num1; // Det(kA)=kDet(A);
                }
            }
 
            // multiplying the diagonal elements to get determinant
            for (unsigned long long int i = 0; i < rows; i++) {
                det = det * mat_copy(i, i);
            }
            return (det / total); // Det(kA)/k=Det(A);
        }

        Matrix transpose(){
            Matrix<T> ans(this->cols, this->rows);
            for(unsigned long long int i = 0; i < this->rows; ++i){
                for(unsigned long long int j = 0; j < this->cols; ++j){
                    ans.setElement(j, i, this->elements[i][j]);
                }
            }
        }

        //swap rows r1 and r2
        void swapRows(unsigned long long int r1, unsigned long long int r2){
            for(unsigned long long int i = 0; i < this->cols; ++i){
                T temp = elements[r1][i];
                elements[r1][i] = elements[r2][i];
                elements[r2][i] = temp;
            }
        }

        //equivalent to r[elements] * T = r[elements]
        void scaleRow(unsigned long long int r, T scalar){
            for(unsigned long long int i = 0; i < this->cols; ++i){
                this->elements[r][i] *= scalar;
            }
        }

        //add multiple of row 1 to row 2
        void scaleAndAddRows(unsigned long long int r1, unsigned long long int r2, T scalar){
            for(unsigned long long int i = 0; i < cols; ++i){
                this->elements[r1][i] += this->elements[r2][i];
            }
        }

        //get an element
        T& operator()(unsigned long long int i, unsigned long long int j){
            return this->elements[i][j];
        }

        void setElement(unsigned long long int r, unsigned long long int c, T val){
            this->elements[r][c] = val;
        }

        T getElement(unsigned long long int r, unsigned long long int c){
            return this->elements[r][c];
        }
};

/*
Shorthand for matrix-based subclasses and pointers
*/

template<typename T> using UniqueMatrixPtr = std::unique_ptr<Matrix<T>>;

template<typename T> using SharedMatrixPtr = std::shared_ptr<Matrix<T>>;

template<typename T> using ComplexMatrix = Matrix<std::complex<T>>;

template<typename T> using UniqueComplexMatrixPtr = std::unique_ptr<ComplexMatrix<T>>;

template<typename T> using SharedComplexMatrixPtr = std::shared_ptr<ComplexMatrix<T>>;

/*
Math operators
*/

template<typename T> Matrix<T> operator+(const Matrix<T>& mat1, const Matrix<T>& mat2){
    auto mat_1_rows = mat1.getRows();
    auto mat_2_rows = mat2.getRows();
    if(mat_1_rows != mat_2_rows){
        throw std::runtime_error("Rows not compatible.");
    }
    if(mat1.getCols() != mat2.getCols()){
        throw std::runtime_error("Cols not compatible.");
    }
    Matrix<T> ans = Matrix<T>(mat1.getRows(), mat1.getCols());
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            ans.setElement(i, j, mat1.getElement(i, j) + mat2.getElement(i, j));
        }
    }
    return ans;
}

template<typename T> Matrix<T> operator-(const Matrix<T>& mat1, const Matrix<T>& mat2){
    if(mat1.getRows() != mat2.getRows()){
        throw std::runtime_error("Rows not compatible.");
    }
    if(mat1.getCols() != mat2.getCols()){
        throw std::runtime_error("Cols not compatible.");
    }
    Matrix<T> ans = Matrix<T>(mat1.getRows(), mat1.getCols());
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            ans.setElement(i, j, mat1.getElement(i, j) - mat2.getElement(i, j));
        }
    }
    return ans;
}

template<typename T> Matrix<T> operator*(const Matrix<T>& mat1, const Matrix<T>& mat2){
    if (mat1.getCols() != mat2.getRows()) {
		throw std::runtime_error("Matrices to multiply incompatible");
	}
	Matrix<T> answer = Matrix<T>(mat1.getRows(), mat2.getCols());

	for (int i = 0; i < mat1.getRows(); ++i){
		for (int j = 0; j < mat2.getCols(); ++j){
			for (int k = 0; k < mat1.getCols(); ++k){
				answer.setElement(i, j, answer.getElement(i, j) + mat1.getElement(i, k) * mat2.getElement(k,j));
			}
		}
	}
	return answer;
}

template<typename T> Matrix<T> operator*(const Matrix<T>& mat1, T scalar){
    Matrix<T> answer(mat1);
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            answer(i, j) *= scalar;
        }
    }
    return answer;
}

template<typename T> Matrix<T> operator*(T scalar, const Matrix<T>& mat1){
    Matrix<T> answer(mat1);
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            answer(i, j) *= scalar;
        }
    }
    return answer;
}

template<typename T> Matrix<T> operator*(Matrix<T>& mat1, T scalar){
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            mat1(i, j) *= scalar;
        }
    }
}

template<typename T> Matrix<T> operator*(T scalar, Matrix<T>& mat1){
    for(unsigned long long int i = 0; i < mat1.getRows(); ++i){
        for(unsigned long long int j = 0; j < mat1.getCols(); ++j){
            mat1(i, j) *= scalar;
        }
    }
}

/*
Input and output operators
*/

template<typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat){
    unsigned long long int rows = mat.getRows();
    unsigned long long int cols = mat.getCols();
    os << rows << " " << cols << std::endl;
    for(unsigned long long int i = 0; i < rows; ++i){
        for(unsigned long long int j = 0; j < cols; ++j){
            os << mat.getElement(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

template<typename T> std::istream& operator<<(std::istream& is, Matrix<T>& mat){
    unsigned long long int rows, cols;
    is >> rows;
    is >> cols;
    for(unsigned long long int i; i < rows; ++i){
        for(unsigned long long int j = 0; j < cols; ++j){
            T tmp;
            is >> tmp;
        }
    }
    return is;
}