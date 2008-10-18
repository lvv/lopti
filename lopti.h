	#ifndef LVV_LOPTI_H
	#define LVV_LOPTI_H


	// this should be included 1st so that array.h will know how to convert
	#ifndef MINIMIZER
		#include	<gsl/gsl_vector.h>
		#include	<condor/Vector.h>
	#endif

	#include <lvv/lvv.h>
	#include <limits>
		using 	std::numeric_limits;

	// functional 
	#include <functional>
		using std::binder1st;
		using std::unary_function;
	#include <boost/function.hpp>
		using boost::function;
	#include <boost/bind.hpp>
		//using boost::bind;
		
		//using namespace boost;
		//using namespace std;

                 template<typename V>
class	minimizer { 

	public:
		typedef  typename V::value_type fp_t;
		typedef  function<fp_t(V&)>	of_ptr_t;
		of_ptr_t			of_;
		int				max_iter_;
		bool				verbose_;
		V				X;
		V				Xmin_;
		fp_t				ymin_;
		int				iter_;
		bool				found_;
		const char*			name_;

	explicit 		minimizer		(V& _X, const char* _name = "unknown")       
	:
		max_iter_	(500),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(-1),
		X        	(_X),
		verbose_ 	(false),
		found_ 		(false),
		name_		(_name)
	{};

	virtual 		~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 

	// set-ters
	virtual minimizer<V>&		object_function		(of_ptr_t of)	{  of_ = of;		return *this;  };
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
	virtual V&			argmin			() 		{  return Xmin_;  };
	virtual bool			found			() 	const	{  return found_; };
	virtual const char*		name			() 	const	{  return name_;  };
	virtual void			print			()		{ MSG("%s %25t iter=%d  \t ymin=%g \t Xmin=%g \n") %name()  %iter()  %ymin()  %Xmin();};
};

                 template<typename V>
class	trust_region_minimizer : public minimizer<V>    { public:

		typedef  typename minimizer<V>::fp_t		fp_t;
		typedef  typename minimizer<V>::of_ptr_t	of_ptr_t;
		using minimizer<V>::name;

		fp_t 				rho_begin_;	// r(rho) start
		fp_t 				rho_end_;	// r end


	explicit trust_region_minimizer		(V& _X, const char* _name= "unknown (trust region type)"):  
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
