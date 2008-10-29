
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr; // used as ofstream ptr

	//#include <gzstream.h>
	#include <iostream>                                                                                                                                         
	#include <fstream>
		using std::ios;
		using std::ofstream;


 namespace lopti {

	#define CLONER(T)	 virtual T&  clone()  const   {  return  *new T(*this); }

 /////////////////////////////////////////////////////////////////////////////////////////  OF: CHEBYQUAD

			template<typename V> 
struct  chebyquad: loft_base<V>		{	 LOFT_TYPES;  LOFT_MEMBERS;  CLONER(chebyquad)		// The Chebyquad test problem (Fletcher, 1965) 

			chebyquad		()	: loft_base<V>("chebyquad") 	{};

	fp_t		operator() 		(V& X)  {

			iter_++;

			matrix <fp_t, V::sz, V::sz+1, 0+B, 0+B>		Y;

			for (int J=0+B; J<N+B; J++)  {
				Y[0+B][J] = 1.0;
				Y[1+B][J] = 2.0*X[J]-1.0;
			}

			for (int I=0+1+B; I<N+B; I++) 
				for (int J=0+B; J<N+B; J++)
					Y[I+1][J]=2.0*Y[1+B][J]*Y[I][J]-Y[I-1][J];

			fp_t 	F  = 0.0;
			int	IW = 1;

			for (int I=0+B; I<N+1+B; I++)  {
				fp_t  SUM = 0.0;
				for (int J=0+B; J<N+B; J++) 	SUM += Y[I][J];
				SUM = SUM/N;
				if (IW > 0)  SUM += 1.0/((I+1-B)*(I+1-B)-2*(I+1-B));
				IW =-IW;
				F += SUM*SUM;
			}

			return F;
	}
 };
 /////////////////////////////////////////////////////////////////////////////////////////  ADAPTER: PLAIN_FN
template<typename V>	struct	plain_fn : loft_base<V>		{  				LOFT_TYPES;  LOFT_MEMBERS;  CLONER(plain_fn)
						function<fp_t(V&)>	of;
	explicit	plain_fn	(function<fp_t(V&)> _of, const string& _name)	 : loft_base<V>(_name) , of(_of)   {};
	fp_t		operator()	(V&  X)	  { assert(!of.empty() && ">> NOT DEFINED OBJ FUNC <<");  iter_++;   fp_t y = (of)(X); return y; }
 };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBERG
template<typename V>	struct	rosenberg  : loft_base<V> { 				LOFT_TYPES;  LOFT_MEMBERS;  CLONER(rosenberg)
			rosenberg		()	: loft_base<V>("rosenberg") 	{ V const  X_answ = {{ 1.0, 1.0 }};   opt(X_answ); };
	fp_t		operator() 		(V& X)  {  iter_++;      return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
 };

 template<typename V>   typename V::value_type    plain_fn_rosenberg  (V& X) { const int B = V::ibg;   return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
 /////////////////////////////////////////////////////////////////////////////////////////  OF: BAD SCALE ROSENBERG
template<typename V, int FACTOR>	struct	bad_scale_rosenberg	 : loft_base<V> { 	LOFT_TYPES;  LOFT_MEMBERS;  CLONER(bad_scale_rosenberg);
			bad_scale_rosenberg	() : loft_base<V>("bad_scale_rosenberg") { V const  X_answ = {{ 1.0, 1.0*FACTOR }};   opt(X_answ); };
	fp_t		operator() 		(V& X)   {  iter_++; return  100 * pow2(X[1+B]/FACTOR-pow2(X[0+B])) + pow2(1-X[0+B]); };
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER: RESCALE
template<typename V>	struct	rescale :  loft_base<V> 	{				LOFT_TYPES;  LOFT_MEMBERS;  CLONER(rescale)
				V R;
			rescale		(loft_cref_t loft_v, V& _R)  : 	R(_R)  {  X_opt_ = loft_v.X_opt_; X_opt_ /= _R;   loft(loft_v); cout << "rescaled opt: " << X_opt_ << endl <<  "opt: " << loft_v.X_opt_ << endl;  };
	const string	name		()	const	{ return  name_ + " rescaled"; };
	fp_t		operator()	(V&  X)		{ iter_++;  	V XR = X;  XR *= R;    return  (*wrapped_loft_v)(XR); };
	fp_t 		opt_distance	(V& X)	const	{ V XR = X;   XR*=R;    return  distance_norm2(wrapped_loft_v->X_opt_, XR); }; // distance in normalized coord
 };
 /////////////////////////////////////////////////////////////////////////////////////////  WRAPPER:  XG_LOG (xgraphic)
template<typename V>	struct	xg_log : loft_base<V> 		{				LOFT_TYPES;   LOFT_MEMBERS;  CLONER(xg_log)
				shared_ptr<ofstream>	log_file;  // need smart ptr becase xg_log dtor-ed on coping
	xg_log	 (loft_cref_t _loft_v, minimizer<V>& mzr)	{
		loft(_loft_v);													assert(log_file == 0 );
		log_file = shared_ptr<ofstream>(new ofstream(("log/" +  (&mzr)->name() + "(" + name() + ")" ).c_str()));	assert(log_file->good());
	};

	fp_t		operator()	(V&  X)			{
		fp_t   y = (*wrapped_loft_v)(X); 			assert(log_file->good());
		// 1 - iter no(gp: splot label);  2 - hight (y);  3 - |X-opt| ignored;  4,5 - X coord;  
		*log_file << format("%d \t %g \t %g \t %22.15g \n")   % (wrapped_loft_v->iter())   %y   % wrapped_loft_v->opt_distance(X)  %X  << flush;  
		return  y;
	};
 };

 }; // namespace lopti
