#include <sstream>
#include <fstream>

#include "auxillaryFunctions.h"

template<class T>

Matrix<T>  readfile (string filename, T numOfCols){
	// create an input filestream
	ifstream input (filename.c_str());
	Matrix<T> data(1,numOfCols);

	string line; // string to read in a line from the file
	getline(input,line); // get first line from file
	istringstream inputss(line); //input string stream

	// get each column from file
	T num;
	for(int k=1;k<numOfCols+1;++k){
		T rownum=0;
		while(getline(input,line)) {
				inputss.clear(); //need to clear the string stream before reuse
				inputss.str(line); // add the new line into the string stream
				//skip to the required column
				for (int i=0;i<k;i++) {
					inputss>>num;
				}
				data(rownum,k-1)=num;
				++rownum;
				if (k==1) data.addRow(0);
				
		}
		
	}
	
	return data;
}


template<class T>

vector<T>  readfile (string filename,T numOfCols,int dummy){
	// create an input filestream
	ifstream input (filename.c_str());
	vector<T> data;

	string line; // string to read in a line from the file
	getline(input,line); // get first line from file
	istringstream inputss(line); //input string stream

	// get specified column from file
	float num;
	while(getline(input,line)) {
			inputss.clear(); //need to clear the string stream before reuse
			inputss.str(line); // add the new line into the string stream
			//skip to the required column
			for (int i=0;i<1;++i) {
				inputss>>num;
				//std::cout<<num<<std::endl;
			}
			data.push_back(num);
	}
	return(data);
}



