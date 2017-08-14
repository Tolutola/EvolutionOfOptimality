#include<vector>
#include<iostream>
#include<tuple>
#include<cmath>

using namespace std;
#ifndef MATRIX_H_
#define MATRIX_H_



template<class T>
class Matrix {
public:
	vector<vector<T> > data;  // Store data in 2D vector
	//	int numrows; // number of rows
	//	int numcols; // number of columns
	Matrix(int, int);     // Constructor for matrix class
	Matrix();
	int numCols() {
		return data.size();
	} // Returns number of columns
	int numRows() {
		return data[0].size();
	} // Returns number of rows
	void setMatrix(T); // sets all elements of the matrix to the same number
	void printMatrix(); // prints out all elements of the matrix

	// add row or column to the matrix
	void addRow(T myvalue); // adds a new row of 'myvalue'
	void addCol(T myvalue); // adds a new col of 'myvalue'

	// scalar operators for my matrix class
	void operator +=(T);
	void operator -=(T);
	void operator *=(T);
	void operator /=(T);

	T& operator ()(int, int); // access individual elements

	// matrix operators (addition and subtraction)
	Matrix<T> operator +(const Matrix&);
	Matrix<T> operator -(const Matrix&);

	// overloading the << operator
	friend ostream& operator <<(ostream &os, Matrix<T> &M) {
		//		os << "Dimensions m x n = ";
		//		os << M.numRows() << " x " << M.numCols() << ".\n";
		for (int i=0;i<M.numRows();i++){
			for (int j=0;j<M.numCols();j++){
				os<<M(i,j)<<" \t ";
			}
			os<<endl;
		}
		return os;
	}

	// transpose function
	void transpose(); // this tranposes the matrix
	//	Matrix<T> transpose(Matrix<T>);// this outputs the result of the transpose
	// into a new matrix

};

template<class T>
Matrix<T> transpose(Matrix<T>);

// this function multiplies two matrices
template<class T>
Matrix<T> multiply(Matrix<T>, Matrix<T>);

// this function  finds the dot product of two matrices
// the result is a vector of dot products between corresponding
// columns of the two matrices
// (inspired by a similar functionality in MATLAB)

template<class T>
Matrix<T> dotProduct(Matrix<T>, Matrix<T>);

// Additional functionality
// Reshaping a matrix

// reshape a matrix; convert from vector to matrix and vice versa
// B=reshape(A,m,n) reshapes matrix A to m rows and n cols. The number of elements of
// A must be equal to m*n otherwise the matrix is returned unchanged
// reshape(A) converts a matrix  to a regular vector (i.e converts from
// matrix object to vector object)
// B=reshape(A,m,n), where A is a regular vector object, produces B which is
// matrix object

template<class T>
Matrix<T> reshape(Matrix<T>, int, int);

template<class T>
vector<T> reshape(Matrix<T>);

template<class T>
Matrix<T> reshape(vector<T>, int, int);

template<class T>
Matrix<T> IdentityMatrix(T );


template<class T>
Matrix<T> backSubstitution(Matrix<T>);


template<class T>
tuple <Matrix<T>,Matrix<T>,T> pivot(Matrix<T>,Matrix<T>, int RowIndex);

template<class T>
tuple <Matrix<T>,Matrix<T>,Matrix<T>,T > forwardElimination(Matrix<T>);


template<class T>
tuple <Matrix<T>,Matrix<T>,Matrix<T>, T > LUPdecomposition(Matrix<T>);

template<class T>
T Determinant(Matrix <T>);

#include "BasicMatrix.tpp"
#include "AdvancedMatrix.tpp"


#endif
