	
	#ifndef LOPTI_OF_H
	#define LOPTI_OF_H

		using std::shared_ptr; // used as ofstream ptr

	//#include <gzstream.h>
	#include <iostream>                                                                                                                                         
	#include <fstream>
		using std::ios;
		using std::ofstream;

	#include <lvv/timer.h>
		using lvv::timer_t;

	#include <lvv/array.h>
		using lvv::array;
		using lvv::matrix;

	#include	<lvv/math.h>
		using lvv::pow2;


 namespace lopti {

/* TODO: replace CLONER macro with CRTP: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template<typename Derived> struct cloneable { Derived* clone() const {  return new Derived(static_cast<Derived const&>(*this)); } };
used as:
	obj: cloneable<obj> { ...};
	obj* o = other_o.clone();
*/

						#define CLONER(CLASS)		\
							CLASS& 	clone()  const		{  return  *new CLASS(*this); }

						
						 #define         OBJECTIVE_MEMBERS    \
							using		objective0<V>::iter_; \
							using		objective0<V>::X_opt_; \
							static const int B = V::ibg; \
							static const int N = V::sz;


						
						#define 	OBJECTIVE_TYPES	\
							typedef		std::shared_ptr<objective0<V> >	objective_p_t; \
							typedef 	const objective0<V>&		objective_cref_t;	\
							typedef		typename V::value_type		T;

						
						#define		NaN	numeric_limits<T>::quiet_NaN()

			template<typename V>
struct	objective0		{									

				OBJECTIVE_TYPES;

			V			X_opt_;
			int			iter_;

	// CTOR
				objective0	()			:  	iter_(0)				{ X_opt_ = -1; };
	virtual			~objective0	()	 = 0;
	virtual objective0<V>&	clone		() const = 0; 

	// set-ters
	void 			known_optimum	(const V& X_answ)	{ X_opt_ = X_answ; };		// known optimum, used for testing optimizers

	// get-ers
	virtual const string	name		()	const		{ return  "not-defined"; };
	virtual int 		iter		()	const		{ return  iter_;  };
	virtual	T 		opt_distance	(const V& X) const	{ return  distance_norm2(X_opt_, X);  };

	// do-ers
	virtual T		operator()	(const V&  X) 		{ return T(); }
	virtual T		eval0		(const V&  X) 		{ return  operator()(X); };
	virtual V&&		eval1		(const V&  X) 		{ assert(false);      return std::move(V()); };
 };

	template<typename V>	objective0<V>::~objective0 () {}; 				// we need this (see http://www.devx.com/tips/Tip/12729)

			template<typename V>
struct	wrapper: objective0<V>		{							
				OBJECTIVE_TYPES;
			objective_p_t		objective_v;

	virtual void		objective	(objective_cref_t ref)	{ objective_v = objective_p_t(&ref.clone()); };
	virtual int 		iter		()	const		{ return   objective_v->iter(); };
	        const string	name		()	const		{ return   objective_v->name(); };
	virtual T		operator()	(const V&  X)		{
										assert( objective_v != 0 );
		T   y = (*objective_v)(X);
		return y;
	}
	virtual T		eval0		(const V&  X) 		{ return  operator()(X); };
	virtual V&&		eval1		(const V&  X) 		{
						assert( objective_v != 0 );
		return objective_v->eval1(X);
	};
 };

//
 /////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBROCK
template<typename V>	struct	rosenbrock  : objective0<V> { 
			OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(rosenbrock)
				rosenbrock		()		{ V const  X_answ = {{ 1.0, 1.0 }};   this->known_optimum(X_answ); };
									// unset view; set surface;  set isosamples 150,150;  set contour base; set cntrparam levels 20; splot  [-3:4] [-2:8]  log10 (100*(y-x**2)**2 + (1-x)**2)
									// set view map ; unset surface;  set grid ; set samples 500,500;  set contour base; set cntrparam levels 20; splot  [-3:4] [-2:8]  log10 (100*(y-x**2)**2 + (1-x)**2)
	virtual const string	name() const  		{  return  "rosenbrock";  };
	T		operator() 		(const V& X)	{  iter_++;      return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
	T		eval0	 		(const V& X)	{  return operator() (X); };
	V&&		eval1	 		(const V& X)	{
				V G; 
				G[0+B] = -400 * X[0+B] * (X[1+B] - pow2(X[0+B]))  -  2*(1-X[0+B]) ;
				G[1+B] =  200 *          (X[1+B] - pow2(X[0+B])) ;
									// (%o3) rb(x0,x1):=100*(x1-x0^2)^2+(1-x0)^2
									// (%i5) diff(rb(x0,x1),x0);
									// (%o5) 	-400*x0*(x1-x0^2)-2*(1-x0)
									// (%i6) diff(rb(x0,x1),x1);
									// (%o6) 	200*(x1-x0^2)
				return std::move(G);
		};

	/*			typedef   matrix<T,V::sz, V::sz>   M; 

	M&&	eval2	 		(V& X)  {
				M H; 
				H[1,1] = 1200*pow2(X[0+B]) − 400*X[1+B] + 2; 		H[1,2] = −400 * X[0+B];
				H[2,1] = -400*X[0+B]; 					H[2,2] = 200;


				return H;
		}; */
 };

 template<typename V>   typename V::value_type    plain_fn_rosenbrock  (V& X) { const int B = V::ibg;   return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: BAD SCALE ROSENBROCK
template<typename V, int FACTOR>	struct	bad_scale_rosenbrock	 : objective0<V> { 	OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(bad_scale_rosenbrock);
	bad_scale_rosenbrock() { V const  X_answ = {{ 1.0, 1.0*FACTOR }};   known_optimum(X_answ); };
	T	operator() 		(V& X)   {  iter_++; return  100 * pow2(X[1+B]/FACTOR-pow2(X[0+B])) + pow2(1-X[0+B]); };
 };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: CHEBYQUAD

template<typename V>	struct	chebyquad: objective0<V>		{	 OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(chebyquad)		// The Chebyquad test problem (Fletcher, 1965) 


	T		operator() 		(const V& X)  {

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

	virtual const string 	name() const		{ return "chebyquad"; };
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER: RESCALE
template<typename V>	struct	rescale :  wrapper<V> 	{		OBJECTIVE_TYPES;   OBJECTIVE_MEMBERS;  CLONER(rescale)
				V R;
			rescale		(objective_cref_t ref, V& _R)  : 	R(_R)  {  X_opt_ = ref.X_opt_;  X_opt_ /= _R;   objective(ref); /*cout << "rescaled opt: " << X_opt_ << endl <<  "opt: " << objective_v.X_opt_ << endl;*/  };
	virtual const string	name		()	const	{ return  wrapper<V>::objective_v->name() + " rescaled"; };
	T		operator()	(const V&  X)		{ iter_++;  	V XR = X;  XR *= R;    return  (*wrapper<V>::objective_v)(XR); };
	T 		opt_distance	(const V& X)	const	{ V XR = X;   XR*=R;    return  distance_norm2(wrapper<V>::objective_v->X_opt_, XR); }; // distance in normalized coord
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER:  XG_LOG (xgraphic)
 template<typename V> class minimizer;

template<typename V>	struct	xg_log : wrapper<V> 		{		OBJECTIVE_TYPES;   OBJECTIVE_MEMBERS;  CLONER(xg_log)
				std::shared_ptr<ofstream>	log_file;  // need smart ptr becase xg_log dtor-ed on coping
	xg_log	 (objective_cref_t _objective_v, minimizer<V>& mzr)	{
		this->objective(_objective_v);								assert(log_file == 0 );
		log_file = shared_ptr<ofstream>(new ofstream(("log/" +  (&mzr)->name() + "(" + wrapper<V>::name() + ")" ).c_str()));	assert(log_file->good());
	};

	T		operator()	(const V&  X)			{
		T   y = (*wrapper<V>::objective_v)(X); 			assert(log_file->good());
		*log_file << setprecision(18) << wrapper<V>::objective_v->iter() << " 	" <<  y  << " " <<  wrapper<V>::objective_v->opt_distance(X)  << " " <<  X  << endl;  
		return  y;
	};
 };


template<typename V>	struct	trace :  wrapper<V> 		{		  CLONER(trace); OBJECTIVE_TYPES;   OBJECTIVE_MEMBERS;
						lvv::timer_t	 timer;
	explicit 	trace	 	(objective_cref_t ref)	{ this->objective(ref); };

	T		operator()	(const V&  X)			{

		if (wrapper<V>::objective_v->iter()==0)   cout << "# (iter)              X[*]               ==F-value" << endl;
		printf("(%d) 	 ",  wrapper<V>::objective_v->iter());  cout <<  X;

		T   y = (*wrapper<V>::objective_v)(X); 			
		cout << "	(" << timer.interval_ticks() << " ticks)  " << y << endl;
		return  y;
	};
 };


 /////////////////////////////////////////////////////////////////////////////////////////  ADAPTER: PLAIN_FN
template<typename V>	struct	make_objective : objective0<V>	{OBJECTIVE_TYPES;  OBJECTIVE_MEMBERS;  CLONER(make_objective)
						std::function<T(V&)>	of;
						const string 		name_;
	explicit	make_objective	(std::function<T(V&)> _of, const string& _name)	 : name_(_name), objective0<V>() , of(_of)   {};
	T		operator()	(V&  X)	  	{ assert(!of.empty() && ">> NOT DEFINED OBJ FUNC <<");  iter_++;   T y = (of)(X); return y; }
	const string 	name		() const	{ return std::move(name_); };
 };
 }; // namespace lopti
 #endif // LOPTI_OF_H
