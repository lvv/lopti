
#ifndef  CONVERT_GSl_H
#define  CONVERT_GSl_H

	#include <gsl/gsl_vector.h>
	#include <lvv/array.h>
		using lvv::array;
	
			template <typename T, int N, int B> array<T,N,B>&
	 operator<<=  (array<T,N,B>& A, const gsl_vector* gV)  {	// operator= should be member, so we are using operator<<
		assert(A.size()==gV->size);  assert(A.ibegin()==0);  
		for (int i=A.ibegin(); i<A.iend(); i++)  A[i] = gsl_vector_get(gV, i);
		return A;
	 };

			template <typename T, int N, int B> array<T,N,B>&
	 operator<<  (array<T,N,B>& A, const gsl_vector* gV)  {	// operator= should be member, so we are using operator<<
		assert(A.size()==gV->size); 
		for (int i=0; i<N; i++)  A[i+B] = gsl_vector_get(gV, i);
		return A;
	 };

			template <typename T, int N, int B> gsl_vector*
	 operator<<=  (gsl_vector* gV, array<T,N,B>& A)  {
		assert(A.size()==gV->size);  assert(A.ibegin()==0);  
		for (int i=A.ibegin(); i<A.iend(); i++)  gsl_vector_set(gV, i, A[i]);
		return gV;
	 };

			template <typename T, int N, int B> gsl_vector*
	 operator<<   (gsl_vector* gV, array<T,N,B>& A)  {
		assert(A.size()==gV->size);
		for (int i=0; i<N; i++)    gsl_vector_set(gV, i, A[i+B]);
		return gV;
	 };

 
  // double[]  to  gsl_vector  wrapper
  struct mk_gsl_vector: public gsl_vector {
  	 mk_gsl_vector(double *a, int n) {
	 	size = n;
		stride = 1;
		data = a;
		block = NULL;
		owner = 0;
	 };
  };
#endif
