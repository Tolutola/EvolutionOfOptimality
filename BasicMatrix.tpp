#include "Matrix.h"

template<class T>
Matrix<T>::Matrix(int rows, int cols) {
	// Allocate memory according for the required number of rows and columns
	data.resize(cols);
	// the matrix will be initialized to all zeros
	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {
			data[i].push_back(0);
		}
	}

}

template<class T>
Matrix<T>::Matrix() {
	// Allocate memory according for the required number of rows and columns
	int cols=1;
	int rows=1;
	data.resize(cols);
	// the matrix will be initialized to all zeros
	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {
			data[i].push_back(0);
		}
	}

}


template<class T>
void Matrix<T>::setMatrix(T myvalue) {
	// sets all elements of the matrix to the same number
	for (int i = 0; i < data.size(); i++) {
		for (auto &j : data[i]) {
			j = myvalue;
		}
	}

}
template<class T>
void Matrix<T>::printMatrix() {
//	cout<<data[0].size()<<endl;
//	cout<<data.size()<<endl;
	for (int i = 0; i < data[0].size(); i++) {
		for (int j = 0; j < data.size(); j++) {
			cout << data[j][i] << "\t";

		}
		cout << endl;
	}
	cout << endl;

}

template<class T>
void Matrix<T>::addCol(T myvalue) {

	int oldRowSize = data[0].size();
	data.resize(data.size() + 1); // add column
	// set right number of rows in new column
	data[data.size() - 1].resize(oldRowSize);

	// set all elements in the new column to myvalue
	for (auto &q : data[data.size() - 1]) {
		q = myvalue;
	}
}

template<class T>
void Matrix<T>::addRow(T myvalue) {
	// set all elements in the new row to myvalue
	for (auto &q : data) {
		q.push_back(myvalue);
	}
}

template<class T>
void Matrix<T>::operator +=(T addValue) {
	// Add addValue to each element in the matrix
	for (int i = 0; i < data.size(); i++) {

		for (auto &j : data[i])
			j += addValue;
	}
}

template<class T>
void Matrix<T>::operator -=(T subtractValue) {
	// Add addValue to each element in the matrix
	Matrix<T>::operator +=(-subtractValue);

}

template<class T>
void Matrix<T>::operator *=(T multiplyValue) {
	// multiply each element by multiplyValue
	for (int i = 0; i < data.size(); i++) {

		for (auto &j : data[i])
			j *= multiplyValue;
	}
}

template<class T>
void Matrix<T>::operator /=(T divideValue) {
	// divide each element by divideValue
	Matrix<T>::operator*=(1 / divideValue);
}

template<class T>
T& Matrix<T>::operator ()(int rowIndex, int colIndex) {
	return data[colIndex][rowIndex];
}

template<class T>
Matrix<T> Matrix<T>::operator +(const Matrix& myMatrix) {
	Matrix<T> C(data[0].size(), data.size());
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[0].size(); j++) {
			C.data[i][j] = data[i][j] + myMatrix.data[i][j];

		}
	}
	return C;

}

template<class T>
Matrix<T> Matrix<T>::operator -(const Matrix& myMatrix) {
	Matrix<T> C(data[0].size(), data.size());
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[0].size(); j++) {
			C.data[i][j] = data[i][j] - myMatrix.data[i][j];

		}
	}

	return C;

}

template<class T>
void Matrix<T>::transpose() {
	Matrix<T> temp(this->numRows(), this->numCols());
	temp.data = data; // temporarily store the data

	int newRowSize = data.size();
	int newColSize = data[0].size();

	data.resize(newColSize); // make the new matrix the right size
	// set right number of rows in new column
	for (auto &q : data) {
		q.resize(newRowSize);
	}

	// fill up newly sized matrix with old data
	for (int i = 0; i < newColSize; i++) {
		for (int j = 0; j < newRowSize; j++) {
			data[i][j] = temp.data[j][i];

		}
	}
}

template<class T>
Matrix<T> transpose(Matrix<T> A) {
	Matrix<T> temp(A.numRows(), A.numCols());
	temp.data = A.data;
	temp.transpose();
	return temp;
}

// this function multiplies two matrices
template<class T>
Matrix<T> multiply(Matrix<T> A, Matrix<T> B) {
	// first check if A and B are compatible
	if (A.numCols() != B.numRows()) {
		cout << "Matrices are incompatible" << endl;
		//return a matrix of zeros of size A
		Matrix<T> errorMatrix(A.numRows(), A.numCols());
		return errorMatrix;
	} else { // compute matrix product
		int colSize = B.numCols(); // number of columns of new matrix
		int rowSize = A.numRows(); // number of rows of new matrix
		int commonDim = A.numCols(); // common dimension

		Matrix<T> matrixProduct(A.numRows(), B.numCols());

		for (int i = 0; i < rowSize; i++) {
			for (int j = 0; j < colSize; j++) {
				for (int q = 0; q < commonDim; q++) {
					matrixProduct(i, j) += A(i, q) * B(q, j);
				}
			}
		}

		return matrixProduct;
	}

}

template<class T>
Matrix<T> dotProduct(Matrix<T> A, Matrix<T> B) {
	// first check if A and B are compatible
	if ((A.numCols() == B.numCols()) & (A.numRows() == B.numRows())) {

		// compute matrix product
		int commonDim = A.numCols(); // common dimension
		int rowSize = 1; // a one dimensional matrix

		Matrix<T> myVector(1, commonDim);

		for (int i = 0; i < A.numCols(); i++) {

			for (int j = 0; j < A.numRows(); j++) {
				myVector.data[i][0] += A.data[i][j] * B.data[i][j];
			}
		}

		return myVector;
	} else {
		cout << "Matrices are incompatible" << endl;
		//return a matrix of zeros of size A
		Matrix<T> errorMatrix(A.numRows(), A.numCols());
		return errorMatrix;
	}
}

template<class T>
Matrix<T> reshape(Matrix<T> A, int m, int n) {
	// first check if the matrix of required dimensions can be formed
	if ((A.numCols() * A.numRows()) != (m * n)) {
		cout << "The required Matrix cannot be formed" << endl;
		return A;
	} else {
		// empty contents of matrix into one vector
		vector<T> mydata;
		for (auto &q : A.data) {
			for (auto &r : q) {
				mydata.push_back(r);
			}
		}

		// create empty matrix of required dimensions
		Matrix<T> resultMatrix(m, n);
		// fill up th matrix with data, column by column
		int counter = 0;
		for (auto &q : resultMatrix.data) {
			for (auto &r : q) {
				r = mydata[counter];
				counter++;

			}
		}
		return resultMatrix;
	}
}

template<class T>
vector<T> reshape(Matrix<T> A) {

	// convert matrix to vector

	// empty contents of matrix into one vector
	vector<T> mydata;
	for (auto &q : A.data) {
		for (auto &r : q) {
			mydata.push_back(r);
		}
	}

	return mydata;

}

template<class T>
Matrix<T> reshape(vector<T> A, int m, int n) {
	// create empty matrix of required dimensions
	Matrix<T> resultMatrix(m, n);
	// first check if the matrix of required dimensions can be formed
	if ((A.size()) != (m * n)) {
		cout << "The required Matrix cannot be formed" << endl;
		// return a matrix of zeros of size m
		return resultMatrix;
	} else {
		// fill up th matrix with data, column by column
		int counter = 0;
		for (auto &q : resultMatrix.data) {
			for (auto &r : q) {
				r = A[counter];
				counter++;
			}
		}
	}
	return resultMatrix;
}


template<class T>
Matrix<T> IdentityMatrix(T n) {
	Matrix<T> A(n,n);

	for (int i=0;i<n;i++){
		A(i,i)=1;
	}

	return A;
}
