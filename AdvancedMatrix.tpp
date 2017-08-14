#include "Matrix.h"

template<class T>
tuple <Matrix<T> , Matrix<T>,T> pivot(Matrix<T> A, Matrix<T> P,int RowIndex ){
	T rowSwap=0; // assume no row swapping to start
	T MaxNumber=A(RowIndex,RowIndex);
	T MaxIndex=RowIndex;

	int n=A.numRows();
	for (int i=RowIndex+1;i<n;i++){
		if (A(i,RowIndex)>MaxNumber){
			MaxNumber=A(i,RowIndex);
			MaxIndex=i;
		}

	}

	if (MaxIndex!=RowIndex){
		rowSwap=MaxIndex;
		T tempVarA; // to help with the swapping
		T tempVarP;

		for (int j=0;j<n+1;j++){
			tempVarA=A(RowIndex,j);
			A(RowIndex,j)=A(MaxIndex,j);
			A(MaxIndex,j)=tempVarA;

			if (j<n){
				tempVarP=P(RowIndex,j);
				P(RowIndex,j)=P(MaxIndex,j);
			    P(MaxIndex,j)=tempVarP;

			}

		}

	}



	return make_tuple(A,P, rowSwap);
}

template<class T>
tuple <Matrix<T>,Matrix<T>,Matrix<T>,T > forwardElimination(Matrix<T> A){
	int n=A.numRows();
	// first create identity matrices for P and L
	Matrix<T> P =IdentityMatrix(T(n));
	Matrix<T> L =IdentityMatrix(T(n));
	T rowSwapped=0;
	T factor;
	T tempVar;

	// counter to track number of swaps
	T counter =0;

	for (int k=0;k<n-1;k++){
		// partial pivoting
		tie(A,P,rowSwapped)=pivot(A,P,k);
		if (rowSwapped!=0) {
			counter+=1;
			for (int r=0;r<k;r++){
				tempVar=L(k,r);
				L(k,r)=L(rowSwapped,r);
				L(rowSwapped,r)=tempVar;
			}
		}


		for (int i=k+1;i<n;i++){
			factor=A(i,k)/A(k,k);
			L(i,k)=factor;
			for(int j=k;j<n+1;j++){
				A(i,j)=A(i,j)-factor*A(k,j);
			}


		}


	}



	return make_tuple(A,P, L,counter);



}



template<class T>
Matrix<T> backSubstitution(Matrix<T> A) {
	int n=A.numRows();
	int m=A.numCols();

	Matrix <T> X(n,1);

	X(n-1,0)=A(n-1,m-1)/A(n-1,n-1);


	T theSum=0;

	for(int i=n-2;i>-1;i--) {
		theSum=0;
		for (int j=i+1;j<n;j++){
			theSum=theSum+ A(i,j)*X(j,0);
		}

		X(i,0)=(A(i,m-1)-theSum)/A(i,i);
	}

	return X;
}


template<class T>
tuple <Matrix<T>,Matrix<T>,Matrix<T>,T > LUPdecomposition(Matrix<T> A){
	int n=A.numRows();
	Matrix<T> L(n,n);
	Matrix<T> U(n,n);
	Matrix<T> P(n,n);
	T counter;

	// forward elimination expects an extra column so we'll create a dummy one
	A.addCol(T(0));
	tie(A,P,L,counter)=forwardElimination(A);
	// we now remove the extra column added
	for(int i=0;i<n;i++){
		for (int j=0;j<n;j++) {
			U(i,j)=A(i,j);
		}
	}

	return make_tuple (L,U,P,counter);


}

template<class T>
T Determinant(Matrix<T> A){
	int n=A.numRows();
	Matrix<T> L(n,n);
	Matrix<T> U(n,n);
	Matrix<T> P(n,n);
	T counter;
	T detA=1;

	tie(L,U,P,counter)=LUPdecomposition(A);

	for (int i=0;i<n;i++){
		detA=detA*L(i,i)*U(i,i);
	}
	detA=detA*pow(-1,counter);

	return detA;


}