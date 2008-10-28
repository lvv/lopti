// TODO 
//	rename loft_base to loft
//	macro for loft member
//
	
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
struct	loft_base		{									// Lopti Object FuncTor

				typedef		loft_base<V>*			loft_p_t;	// virt_ptr to  loft_base
				typedef 	const loft_base<V>&		loft_cref_t;	// virt ref
				typedef		typename V::value_type		fp_t;
			string			name_;
			V			X_opt_;
			int			iter_;
			loft_p_t		wrapped_loft_v;

	// CTOR
	explicit		loft_base	()			:  wrapped_loft_v(0),	iter_(0)				{ X_opt_ = 0; };
	// like copy-CTOR  
	//explicit		loft_base	(loft_cref_t loft_cref)	:  wrapped_loft_v(&loft_cref.clone()),	iter_(0), name_(loft_cref.name())	{cout << "loft_base::loft_base(cref) \n"; X_opt_ = 0; };
	explicit		loft_base	(const string& s)	:  wrapped_loft_v(0),	iter_(0), name_(s)			{ X_opt_ = 0; };
	virtual loft_base<V>&	clone		()	const		{ cout << "base.clone\n"; assert(false); return  *new loft_base<V>(*this); }


	// set-ters
	void 			opt		(const V& X_answ)	{ X_opt_ = X_answ; };
	virtual void		loft		(loft_cref_t loft_cref)	{ cout << "loft_base::loft clone: \n"; wrapped_loft_v = &loft_cref.clone();	name_ = loft_cref.name(); };

	// get-ers
	virtual const string	name		()	const		{ return  name_; };
	virtual int 		size		()	const		{ return  V::size(); }; // TODO to delete
	virtual int 		iter		()	const		{ return  wrapped_loft_v == 0  ?  iter_  :  wrapped_loft_v->iter(); };
	virtual	fp_t 		opt_distance	(V& X)	const		{ return  distance_norm2(X_opt_, X); };
	virtual	bool 		empty		()	const		{ return  wrapped_loft_v == 0; };

	// do-ers
	virtual fp_t		operator()	(V&  X)			{
		cout << "base.name(): " << name() << "\n";
					assert( this->wrapped_loft_v != 0 );
		fp_t   y = (const_cast<loft_base<V>&> (*this->wrapped_loft_v))(X); 
		return y;
	}

	virtual void		reset		()		{ iter_ = 0; };
 };


                 template<typename V>
class	minimizer { public:
				typedef		typename V::value_type		fp_t;
				typedef		loft_base<V>*			loft_p_t;	// virt_ptr to  loft_base
				typedef 	const loft_base<V>&		loft_cref_t;	// virt ref

			loft_p_t			loft_v;	
			int				max_iter_;
			bool				verbose_;
			V				X;
			V				Xmin_;
			fp_t				ymin_;
			int				iter_;
			bool				found_;
			string				name_;

				/*explicit 		
	minimizer		(loft_v_t  ref, const string& _name = "unknown")  
	:
		loft_v		(ref.clone()),
		max_iter_	(500),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(0),
		verbose_ 	(false),
		found_ 		(false),
		name_		(_name)
	{};*/

				explicit 		
	minimizer		(const string& _name = "unknown")  
	:
		name_		(_name),
		max_iter_	(500),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(0),
		verbose_ 	(false),
		found_ 		(false)
	{};

	virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters
	virtual minimizer<V>&		loft			(loft_cref_t  ref){ cout << "minimizer::loft()\n"; loft_v = &ref.clone(); return *this; };

	virtual minimizer<V>&		X0			(V& _X) 	{  X  = _X;	return *this;  };
	virtual minimizer<V>&		max_iter		(int mx)	{  max_iter_   = mx;	return *this;  };
	virtual minimizer<V>&		verbose			(bool flag)	{  verbose_ = flag;	return *this;  };
	virtual minimizer<V>&		rho_begin		(fp_t rho)	{  assert(false); return *this;  }; // NOP, will be defined in derived
	virtual minimizer<V>&		rho_end			(fp_t rho)	{  assert(false); return *this;  }; // (can not use pure virtual)
	virtual minimizer<V>&		step			(V&)		{  assert(false); return *this;  };  // TODO remove this from base
	virtual minimizer<V>&		characteristic_size	(fp_t)		{  assert(false); return *this;  };  // TODO remove this from base

	// get-ters
	virtual const string		name			() 	const	{  return (format("%s-%d")  %name_  %(V::size()) ).str(); };
	virtual fp_t 	 		ymin			()	const	{  return ymin_; };
	virtual V 	 		Xmin			()	const	{  return Xmin_; };
	virtual fp_t 	 		iter			()	const	{  return iter_; };
	virtual bool			found			() 	const	{  return found_; };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_; };
	virtual void			print			()		{ MSG("%s(%s)  %35t  iter=%d  \t ymin=%g \t Xmin=%22.15g \n") %name() %loft_v->name() %loft_v->iter()  %ymin()  %Xmin();};
 };

				 template<typename V>
class	trust_region_minimizer : public minimizer<V>    { public:

				typedef  typename minimizer<V>::fp_t			fp_t;
				using minimizer<V>::name;
				typedef 	const loft_base<V>&		loft_v_t;

			fp_t 				rho_begin_;	// r(rho) start
			fp_t 				rho_end_;	// r end

					explicit
	trust_region_minimizer		(const char* _name= "unknown (trust region type)")
	: 
		minimizer<V>(_name),
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
