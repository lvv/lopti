	#ifndef LVV_LOPTI_H
	#define LVV_LOPTI_H

	#include <lvv/lvv.h>
	#include <lvv/array.h>
		using lvv::array;
	#include <limits>
		using 	std::numeric_limits;

	// functional 
	#include <functional>
		using std::binder1st;
		using std::unary_function;
		using std::mem_fun;
	#include <boost/function.hpp>
		using boost::function;
	#include <boost/bind.hpp>
		//using boost::bind;
		

			template<typename V>
class	loft		{  public:				// Lopti Object FuncTor

				typedef		loft<V>*			loft_v_t;	// virt_ptr to  loft
				typedef		typename V::value_type		fp_t;
			string			name_;
			V			X_opt_;
			int			iter_;
			loft_v_t		wrapped_loft_v;

	// CTOR
	loft		(loft_v_t loft_v)	:  wrapped_loft_v(loft_v),	name_(loft_v->name()), iter_(0)	{ X_opt_ = 0; };
	loft		(const string& s)	:  wrapped_loft_v(0),		name_(s),              iter_(0)	{ X_opt_ = 0; };

	// set-ters
	void 			opt		(const V& X_answ)	{ X_opt_ = X_answ; };

	// get-ers
	virtual const string&	name		()	const	{ return  name_; };
	virtual int 		size		()	const	{ return  V::size(); };
	virtual int 		iter		()	const	{ return  wrapped_loft_v ? wrapped_loft_v->iter(): iter_; };
	virtual	fp_t 		opt_distance	(V& X)	const	{ return  distance_norm2(X_opt_, X); };

	// do-ers
	virtual fp_t		operator()	(V&  X)			{
					assert(this->wrapped_loft_v);
		fp_t   y = (*this->wrapped_loft_v)(X); 
		return y;
	}
	virtual void		reset		()		{ iter_ = 0;};
 };

			template<typename V>
class	plain_fn : public loft<V>		{  public:				// Lopti Object FuncTor
				typedef		typename V::value_type		fp_t;
			function<fp_t(V&)>	of;
	explicit		plain_fn	(function<fp_t(V&)> _of, const string& _name="unknown")	 :  of(_of), loft<V>(_name)  {};
	virtual fp_t		operator()	(V&  X)			{ assert(this->of); fp_t   y = (this->of)(X); return y; }
 };


                 template<typename V>
class	minimizer { public:
				typedef 	loft<V>				loft_t;
				typedef 	loft<V>*			loft_v_t;
				typedef		typename V::value_type		fp_t;

			loft_v_t			loft_v;	
			int				max_iter_;
			bool				verbose_;
			V				X;
			V				Xmin_;
			fp_t				ymin_;
			int				iter_;
			bool				found_;
			string				name_;

				explicit 		
	minimizer		(V& _X, const string& _name = "unknown")       
	:
		max_iter_	(500),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(0),
		X        	(_X),
		verbose_ 	(false),
		found_ 		(false),
		name_		(_name),
		loft_v		(0)		// Lopti Object FuncTor  ptr
	{};

	virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters
	virtual minimizer<V>&		object_functOR		(loft_v_t  p)	{ loft_v = p; return *this; };
	virtual const string		name			() 	const	{  return (format("%s-%d")  %name_  %(V::size()) ).str(); };

	virtual minimizer<V>&		max_iter		(int mx)	{  max_iter_   = mx;	return *this;  };
	virtual minimizer<V>&		verbose			(bool flag)	{  verbose_ = flag;	return *this;  };
	virtual minimizer<V>&		rho_begin		(fp_t rho)	{  return *this;  };
	virtual minimizer<V>&		rho_end			(fp_t rho)	{  return *this;  };
	virtual minimizer<V>&		step			(V&)		{  return *this;  };
	virtual minimizer<V>&		characteristic_size	(fp_t)		{  return *this;  };

	// get-ters
	virtual fp_t 	 		ymin			()	const	{  return ymin_; };
	virtual V 	 		Xmin			()	const	{  return Xmin_; };
	virtual fp_t 	 		iter			()	const	{  return iter_; };
	virtual bool			found			() 	const	{  return found_; };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_; };
	//virtual void			print			()		{ MSG("%s(%s)  %25t  iter=%d  \t ymin=%g \t Xmin=%20.12g \n") %name() %loft_v->name() %loft_v->iter()  %ymin()  %Xmin();};
	virtual void			print			()		{ MSG("%s(%s)  %25t  iter=%d  \t ymin=%g \t Xmin=%20.12g \n") %name() %loft_v->name() %loft_v->iter()  %ymin()  %Xmin();};
 };

				 template<typename V>
class	trust_region_minimizer : public minimizer<V>    { public:

				typedef  typename minimizer<V>::fp_t			fp_t;
				using minimizer<V>::name;

			fp_t 				rho_begin_;	// r(rho) start
			fp_t 				rho_end_;	// r end

					explicit
	trust_region_minimizer		(V& _X, const char* _name= "unknown (trust region type)"):  
		minimizer<V>	(_X, _name),
		rho_begin_ 	(numeric_limits<fp_t>::quiet_NaN ()),
		rho_end_   	(numeric_limits<fp_t>::quiet_NaN ())
	{};

	virtual minimizer<V>&	rho_begin		(fp_t rho)	{ rho_begin_ = rho;  return *this; };
	virtual minimizer<V>&	rho_end			(fp_t rho)	{ rho_end_   = rho;  return *this; };
 };


	// of prog inlude <lopti.h> then include all
	#ifndef MINIMIZER
		#include	<lopti/condor-wrap.h>
		#include   	<lopti/newuoa-wrap.h>
		#include	<lopti/gsl-nelder-mead-wrap.h>
	#endif 

	#endif 
