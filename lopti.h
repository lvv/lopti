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
		
		//using namespace boost;
		//using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////   OF_BASE

			template<typename V>
class	loft		{  public:

			typedef		loft<V>*	loft_v_t;	// virt_ptr to  loft
			typedef		typename V::value_type		fp_t;
		string			name_;
		V			X_opt_;
		int			iter_;
		//typename loft<typename V>*	wrapped_loft_v;
		loft_v_t		wrapped_loft_v;

	// CTOR
	loft		() 			: name_("unknown"), wrapped_loft_v(0)		{ X_opt_ = 0; };
	loft		(loft_v_t loft_v)	: name_(loft_v->name()), wrapped_loft_v(loft_v)		{ X_opt_ = 0; };
	loft		(const string& s)	: name_(s), wrapped_loft_v(0)		{ X_opt_ = 0; };
	//loft		(loft<V&>* loft_v, const string& s) : name_(s), wrapped_loft_v(loft_v)		{ X_opt_ = 0; };

	// set-ters
	void 			opt		(const V& X_answ)	{ X_opt_ = X_answ; };

	// get-ers
	virtual const string&	name		()	const	{ return  name_; };
	virtual int 		size		()	const	{ return  V::size(); };
	virtual int 		iter		()	const	{ return  iter_; };
	virtual	fp_t 		opt_distance	(V& X)	const	{ return  distance_norm2(X_opt_, X); };

	// do-ers
	virtual fp_t		operator()	(V&  X)			{
		assert(this->wrapped_loft_v);
		fp_t   y = (*this->wrapped_loft_v)(X); 
		return y;
	}
	//virtual fp_t		operator()	(V& X) 		{ assert(false); return X[0]; };
	virtual void		reset		()		{ iter_ = 0;};
 };


                 template<typename V>
class	minimizer { public:

				typedef 	loft<V>				loft_t;
				typedef 	loft<V>*			loft_v_t;
				typedef		typename V::value_type		fp_t;

			loft_v_t			loft_v;	
			function<fp_t(V&)>		oco;
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
		oco		(0),		// Object Callable Object (boost::functor)
		loft_v		(0)		// Lopti Object FuncTor  ptr
	{};

	virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters
	virtual minimizer<V>&		object_functOR		(loft_t* oft_p)	{
		loft_v = oft_p;
			// oco = *loft_v;    - compiles but operator() of base class called
			// oco = boost::bind(&loft<V&>::operator(),_1)(*loft_v);   //  error: ‘operator()’ is not a member of ‘loft
			//oco = boost::bind(mem_fun(&loft_t::operator()),_1)(loft_v);   //   ‘* f’ cannot be used as a function
		name_ += "-" + loft_v->name();
		return *this; 
	};

	virtual minimizer<V>&		object_functION		(function<fp_t(V&)> of)	{  oco = of;		return *this;  };

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
	virtual const string		name			() 	const	{  return (format("%s-%d") %name_ %(V::size())).str(); };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_;  };
	virtual void			print			()		{ MSG("%s  %25t iter=%d  \t ymin=%g \t Xmin=%20.12g") %name() %iter()  %ymin()  %Xmin();};
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
