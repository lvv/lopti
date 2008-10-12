#include <newuoa-wrap.h>
			template<typename V>  typename V::value_type
of_chebyquad		(V& X) 		{			// The Chebyquad test problem (Fletcher, 1965) 
//of_chebyquad		(V& X, void* var)   {			// The Chebyquad test problem (Fletcher, 1965) 
	
	const int N = V::size();
	typedef  typename V::value_type fp_t;

	array<array<fp_t,V::sz,1>, V::sz+1,1> Y;

     	for (int J=1; J<=N; J++)  {
		Y[1][J] = 1.0;
		Y[2][J] = 2.0*X[J]-1.0;
	}

     	for (int I=2; I<=N; I++) 
		for (int J=1; J<=N; J++)
			Y[I+1][J]=2.0*Y[2][J]*Y[I][J]-Y[I-1][J];

     	fp_t 	F  = 0.0;
     	int	NP = N+1;
     	int	IW = 1;

     	for (int I=1; I<=NP; I++)  {
		fp_t  SUM=0.0;
		for (int J=1; J<=N; J++) 	SUM += Y[I][J];
		SUM = SUM/N;
		if (IW > 0)  SUM += 1.0/(I*I-2*I);
		IW =-IW;
	   	F += SUM*SUM;
	}

	return F;
}

			template<typename V>  typename V::value_type
//of_rosenberg		(V& X, void* var=NULL)   {  return 100*pow2(X[2]-pow2(X[1]))+pow2(1-X[1]);  };
of_rosenberg		(V& X)   {  return 100*pow2(X[2]-pow2(X[1]))+pow2(1-X[1]);  };



int main() {

	/*
	int    IPRINT = 2;
	int    MAXFUN = 5000;
	double RHOEND = 1.0e-6;
	{
		const int N = 2; const int NPT = 2*N+1;
		typedef array<double,N, 1> vector;
		vector X;
		for (int I=1; I<=N; I++) X[I] = I/double(N+1);
		double RHOBEG = 0.2 * X[1];
		newuoa<vector NPT>(&of_chebyquad, X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	{
		const int N = 4; const int NPT = 2*N+1;
		typedef array<double,N, 1> vector;
		vector X;
		for (int I=1; I<=N; I++) X[I] = I/double(N+1);
		double RHOBEG = 0.2 * X[1];
		newuoa<vector, NPT>(&of_chebyquad, X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	{
		const int N = 6; const int NPT = 2*N+1;
		typedef array<double,N, 1> vector;
		vector X;
		for (int I=1; I<=N; I++) X[I] = I/double(N+1);
		double RHOBEG = 0.2 * X[1];
		newuoa<vector, NPT>(&of_chebyquad, X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}
	{
		const int N = 8; const int NPT = 2*N+1;
		typedef array<double,N, 1> vector;
		vector X;
		for (int I=1; I<=N; I++) X[I] = I/double(N+1);
		double RHOBEG = 0.2 * X[1];
		newuoa<vector, NPT>(&of_chebyquad, X, RHOBEG, RHOEND, IPRINT, MAXFUN);
	}

	*/
	{
		const int N=2;
		//typedef lvv::array<double,N>		array_t;	

		typedef array<double,N, 1> vector;
		//array_t		X ; // *= {{ -1.2, 1 }};
		vector		X;
		for (int I=1; I<=N; I++)	X[I] = I/double(N+1);

		newuoa_wrap<vector, 2*N+1>	mzr(of_chebyquad<vector>, X);
		mzr.rho_begin		(0.2*X[1]);
		mzr.rho_end		(1e-4);
		mzr.verbose		(true);
		vector	Xmin = mzr.argmin();
		
		MSG("# Result: Xmin%.10g   y=%.10g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());
	}

	{
		const int N=2;
		typedef array<double,N,1>		array_t;	

		array_t			X= {{ -1.2, 1 }};
												// good too
												//function<double(array_t&)>	fct;
												//fct = of_rosenberg<array_t>;
												//newuoa_wrap<array_t>	mzr(fct,  X);

		newuoa_wrap<array_t>	mzr(of_rosenberg<array_t>,  X);
		mzr.rho_begin		(0.5);
		mzr.rho_end		(1e-4);
		mzr.verbose		(true);
		array_t			Xmin = mzr.argmin();
		
		MSG("# Result: Xmin%.10g   y=%.10g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());
	}
return 0 ;}
