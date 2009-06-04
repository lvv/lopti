	
	#ifndef LOPTI_OF_H
	#define LOPTI_OF_H

	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr; // used as ofstream ptr

	//#include <gzstream.h>
	#include <iostream>                                                                                                                                         
	#include <fstream>
		using std::ios;
		using std::ofstream;

	#include <lvv/timer.h>
		using lvv::Timer;

	#include <lvv/array.h>
		using lvv::array;
		using lvv::matrix;

	#include <boost/function.hpp>
		using boost::function;

 namespace lopti {

						#define CLONER(CLASS)		\
							virtual CLASS& 	clone()  const		{  return  *new CLASS(*this); }


						 #define         OBJECTIVE_MEMBERS    \
							using		objective0<V>::iter_; \
							using		objective0<V>::wrapped_objective_v; \
							using		objective0<V>::name_; \
							using		objective0<V>::name; \
							using		objective0<V>::X_opt_; \
							static const int B = V::ibg; \
							static const int N = V::sz;

						
						#define 	OBJECTIVE_TYPES	\
							typedef		shared_ptr<objective0<V> >	objective_p_t; \
							typedef 	const objective0<V>&		objective_cref_t;	\
							typedef		typename V::value_type		T;

						
						#define		NaN	numeric_limits<T>::quiet_NaN()
			template<typename V>
struct	objective0		{									// Lopti Object FuncTor

				OBJECTIVE_TYPES;

			string			name_;
			V			X_opt_;
			int			iter_;
			objective_p_t		wrapped_objective_v;

	// CTOR
	explicit		objective0	()			:  	iter_(0)				{ X_opt_ = -1; };
	explicit		objective0	(const string& s)	:  	iter_(0),	name_(s)		{ X_opt_ = -1; };
	virtual			~objective0	()	= 0;
	virtual objective0<V>&	clone		() const = 0; 

	// set-ters
	void 			known_optimum	(const V& X_answ)	{ X_opt_ = X_answ; };		// known optimum, used for testing optimizers
	virtual void		objective	(objective_cref_t ref)	{ wrapped_objective_v = objective_p_t(&ref.clone());	name_ = ref.name(); };

	// get-ers
	virtual const string	name		()	const		{ return  name_; };
	virtual int 		iter		()	const		{ return  ! wrapped_objective_v  ?  iter_  :  wrapped_objective_v->iter(); };
	virtual	T 		opt_distance	(V& X)	const		{ return  distance_norm2(X_opt_, X);  };
	virtual	bool 		empty		()	const		{ return  ! wrapped_objective_v;  };

	// do-ers
	virtual T		operator()	(V&  X)			{
					assert( wrapped_objective_v != 0 );
		T   y = (*wrapped_objective_v)(X);
		return y;
	}

	//virtual void		reset		()		{ iter_ = 0; };
 };

	template<typename V>	objective0<V>::~objective0 () {}; 				// we need this (see http://www.devx.com/tips/Tip/12729)

			template<typename V>
struct	objective1: objective0<V>	{									// Lopti Object FuncTor
	// do-ers
	virtual T		eval0	(V&  X)			{
					assert( wrapped_objective_v != 0 );
		T   y = (*wrapped_objective_v).eval0(X);
		return y;
	}
	virtual V		eval1	(V&  X)			{
					assert( wrapped_objective_v != 0 );
		V   G = (*wrapped_objective_v).eval1(X);
		return G;
	}
}

 /////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBROCK
template<typename V>	struct	rosenbrock  : objective0<V> { 				OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(rosenbrock)
	rosenbrock	()	: objective0<V>("rosenbrock") 	{ V const  X_answ = {{ 1.0, 1.0 }};   known_optimum(X_answ); };
	T	operator() 		(V& X)  {  iter_++;      return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
	T	eval0	 		(V& X)  {  operator() (X) };
	V&&	eval1	 		(V& X)  {
				V G; 
				G[0+B] = -400 * X[0+B] * (X[1+B]-pow2(X[0+B])) ;
				G[1+B] =  200 *          (X[1+B]-pow2(X[0+B])) ;
		};
									// (%o3) rb(x0,x1):=100*(x1-x0^2)^2+(1-x0)^2
									// (%i5) diff(rb(x0,x1),x0);
									// (%o5) 	-400*x0*(x1-x0^2)-2*(1-x0)
									// (%i6) diff(rb(x0,x1),x1);
									// (%o6) 	200*(x1-x0^2)
 };

 template<typename V>   typename V::value_type    plain_fn_rosenbrock  (V& X) { const int B = V::ibg;   return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: BAD SCALE ROSENBROCK
template<typename V, int FACTOR>	struct	bad_scale_rosenbrock	 : objective0<V> { 	OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(bad_scale_rosenbrock);
	bad_scale_rosenbrock() : objective0<V>("bad_scale_rosenbrock") { V const  X_answ = {{ 1.0, 1.0*FACTOR }};   known_optimum(X_answ); };
	T	operator() 		(V& X)   {  iter_++; return  100 * pow2(X[1+B]/FACTOR-pow2(X[0+B])) + pow2(1-X[0+B]); };
 };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: CHEBYQUAD

			template<typename V> 
struct  chebyquad: objective0<V>		{	 OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(chebyquad)		// The Chebyquad test problem (Fletcher, 1965) 

			chebyquad		()	: objective0<V>("chebyquad") 	{};

	T		operator() 		(V& X)  {

			iter_++;

			matrix <T, V::sz, V::sz+1, 0+B, 0+B>		Y;

			for (int J=0+B; J<N+B; J++)  {
				Y[0+B][J] = 1.0;
				Y[1+B][J] = 2.0*X[J]-1.0;
			}

			for (int I=1+B; I<N+B; I++) 
				for (int J=0+B; J<N+B; J++)
					Y[I+1][J]=2.0*Y[1+B][J]*Y[I][J]-Y[I-1][J];

			T 	F  = 0.0;
			int	IW = 1;

			for (int I=0+B; I<N+1+B; I++)  {
				T  SUM = 0.0;
				for (int J=0+B; J<N+B; J++) 	SUM += Y[I][J];
				SUM = SUM/N;
				if (IW > 0)  SUM += 1.0/((I+1-B)*(I+1-B)-2*(I+1-B));
				IW =-IW;
				F += SUM*SUM;
			}

			return F;
	}
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER: RESCALE
template<typename V>	struct	rescale :  objective0<V> 	{				OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(rescale)
				V R;
			rescale		(objective_cref_t objective_v, V& _R)  : 	R(_R)  {  X_opt_ = objective_v.X_opt_; X_opt_ /= _R;   objective(objective_v); /*cout << "rescaled opt: " << X_opt_ << endl <<  "opt: " << objective_v.X_opt_ << endl;*/  };
	const string	name		()	const	{ return  name_ + " rescaled"; };
	T		operator()	(V&  X)		{ iter_++;  	V XR = X;  XR *= R;    return  (*wrapped_objective_v)(XR); };
	T 		opt_distance	(V& X)	const	{ V XR = X;   XR*=R;    return  distance_norm2(wrapped_objective_v->X_opt_, XR); }; // distance in normalized coord
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER:  XG_LOG (xgraphic)
 template<typename V> class minimizer;

template<typename V>	struct	xg_log : objective0<V> 		{				OBJECTIVE_TYPES;   OBJECTIVE_MEMBERS;  CLONER(xg_log)
				shared_ptr<ofstream>	log_file;  // need smart ptr becase xg_log dtor-ed on coping
	xg_log	 (objective_cref_t _objective_v, minimizer<V>& mzr)	{
		objective(_objective_v);											assert(log_file == 0 );
		// FIXME: test if there is a "log" dir
		log_file = shared_ptr<ofstream>(new ofstream(("log/" +  (&mzr)->name() + "(" + name() + ")" ).c_str()));	assert(log_file->good());
	};

	T		operator()	(V&  X)			{
		T   y = (*wrapped_objective_v)(X); 			assert(log_file->good());
		// 1 - iter no(gp: splot label);  2 - hight (y);  3 - |X-opt| ignored;  4,5 - X coord;  
		//*log_file << format("%d \t %22.18g \t %22.18g \t %22.18g \n")   % (wrapped_objective_v->iter())   %y   % wrapped_objective_v->opt_distance(X)  %X  << flush;  
		*log_file << setprecision(18) << wrapped_objective_v->iter() << " 	" <<  y  << " " <<  wrapped_objective_v->opt_distance(X)  << " " <<  X  << endl;  
		return  y;
	};
 };

template<typename V>	struct	trace : objective0<V> 		{				  CLONER(trace); OBJECTIVE_TYPES;   OBJECTIVE_MEMBERS;
	Timer	timer;
	trace	 (objective_cref_t _objective_v)	{ objective(_objective_v);	};

	T		operator()	(V&  X)			{

		if (wrapped_objective_v->iter()==0)   cout << "# (iter)              X[*]               ==F-value" << endl;
		//cout << format("(%d) \t % 11.8g")   % (wrapped_objective_v->iter())  %X;
		printf("(%d) 	 ",  wrapped_objective_v->iter());  cout <<  X;

		//MSG("(%d/%.1fs") % cnt++  %timer();                                                                                                                
		//for(int i=0; i<param_size; ++i)      MSG("%=10.6g")   %param[i];

		T   y = (*wrapped_objective_v)(X); 			
		//cout << format("\t(%.1fs)   ==%18.13g\n")  %timer()    %y<< flush;  
		printf("	(%.1fs)   ==%18.13g\n", timer(), y);	cout << flush;
		return  y;
	};
 };

 /////////////////////////////////////////////////////////////////////////////////////////  ADAPTER: PLAIN_FN
template<typename V>	struct	make_objective : objective0<V>		{  				OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(make_objective)
						function<T(V&)>	of;
	explicit	make_objective	(function<T(V&)> _of, const string& _name)	 : objective0<V>(_name) , of(_of)   {};
	T	operator()	(V&  X)	  { assert(!of.empty() && ">> NOT DEFINED OBJ FUNC <<");  iter_++;   T y = (of)(X); return y; }
 };
 }; // namespace lopti
 #endif // LOPTI_OF_H
