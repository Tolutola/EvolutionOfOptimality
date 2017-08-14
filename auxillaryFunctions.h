#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<tuple>
#include<cmath>
#include "Matrix.h"
#include "mosek.h"

using namespace std;
#ifndef AUX_H_
#define AUX_H_

template<class T>
Matrix<T> readfile (string filename, T numOfCols);

template<class T>
vector<T> readfile (string filename);

#include "readfile.tpp"


void writeFluxDistributions(Matrix<double> F1,Matrix<double> F2,Matrix<double>F3,Matrix<double>F4) {
	// saving the data in a text file
	// delete any exiting file with the names i'm going to use
	ofstream myfileD1, myfileD2, myfileD3,myfileD4;
	myfileD1.open("FluxMatrix1.txt", ofstream::out | ofstream::trunc);
	myfileD1.close();

	myfileD2.open("FluxMatrix2.txt", ofstream::out | ofstream::trunc);
	myfileD2.close();

	myfileD3.open("FluxMatrix3.txt", ofstream::out | ofstream::trunc);
	myfileD3.close();

	myfileD4.open("FluxMatrix4.txt", ofstream::out | ofstream::trunc);
	myfileD4.close();

	ofstream myfile1("FluxMatrix1.txt",ios::app);
	if (myfile1.is_open()) {
		myfile1<<F1<<endl;
		myfile1.close();
	}

	ofstream myfile2("FluxMatrix2.txt",ios::app);
	if (myfile2.is_open()) {
		myfile2<<F2<<endl;
		myfile2.close();
	}

	ofstream myfile3("FluxMatrix3.txt",ios::app);
	if (myfile3.is_open()) {
		myfile3<<F3<<endl;
		myfile3.close();
	}

	ofstream myfile4("FluxMatrix4.txt",ios::app);
	if (myfile4.is_open()) {
		myfile4<<F4<<endl;
		myfile4.close();
	}
}

#endif
