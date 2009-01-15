	
	#ifndef LOPTI_H
	#define LOPTI_H

	#include <lvv/lvv.h>
	#include <lvv/array.h>
		using lvv::array;

	#include <cstdio>
		// sprintf

	#include <limits>
		using 	std::numeric_limits;

	// functional  		// TODO: try to remove 
	//#include <functional>
		//using std::binder1st;
		//using std::unary_function;
		//using std::mem_fun;
	#include <boost/function.hpp>
		using boost::function;
	#include <boost/bind.hpp>
		//using boost::bind;
	
	#include <lopti/object_function.h>

namespace lopti {
				
						 #define         MINIMIZER_MEMBERS    \
							using minimizer<V>::X;  		\
							using minimizer<V>::x0;			\
							using minimizer<V>::max_iter_;		\
							using minimizer<V>::iter_;		\
							using minimizer<V>::verbose_;		\
							using minimizer<V>::verbose_;		\
							using minimizer<V>::ymin_;		\
							using minimizer<V>::Xmin_;		\
							using minimizer<V>::name;		\
							using minimizer<V>::name_;		\
							using minimizer<V>::found_;		\
							using minimizer<V>::loft_v;

						 #define         TR_MINIMIZER_MEMBERS    \
							using trust_region_minimizer<V>::rho_begin_;	\
							using trust_region_minimizer<V>::rho_end_;


					 template<typename V>
struct	minimizer {
					LOFT_TYPES;

			loft_p_t			loft_v;	
			int				max_iter_;
			bool				verbose_;
			V				X;
			V				Xmin_;
			fp_t				ymin_;
			int				iter_;
			bool				found_;
			string				name_;

				explicit 		
	minimizer		(const string& _name = "unknown")  
	:
		name_		(_name),
		max_iter_	(10000),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(0),
		verbose_ 	(false),
		found_ 		(false)
	{};

	//virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters
	virtual minimizer<V>&		loft			(loft_cref_t  ref){ loft_v = &ref.clone(); return *this; };

	virtual minimizer<V>&		x0			(V& _X) 	{  X  = _X;	return *this;  };
	virtual minimizer<V>&		max_iter		(int mx)	{  max_iter_   = mx;	return *this;  };
	virtual minimizer<V>&		verbose			(bool flag)	{  verbose_ = flag;	return *this;  };


	// get-ters
	//virtual const string		name			() 	const	{  return (format("%s-%d")  %name_  %(V::size()) ).str(); };
	virtual const string		name			() 	const	{ 
			char	buf[100];
			sprintf(buf, "%s-%d",  name_.c_str(), V::size());
			return string(buf);
	};
	virtual fp_t 	 		ymin			()	const	{  return ymin_; };
	virtual V 	 		Xmin			()	const	{  return Xmin_; };
	virtual fp_t 	 		iter			()	const	{  return iter_; };
	virtual bool			found			() 	const	{  return found_; };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_; };
	//virtual void			print			()		{ MSG("%s(%s)  %35t  iter=%d  \t ymin=%g \t Xmin=%22.15g \n") %name() %loft_v->name() %loft_v->iter()  %ymin()  %Xmin();};
	virtual void			print			()		{ printf("%s(%s)  	  iter=%d  	 ymin=%g 	 Xmin: ", name().c_str(),  loft_v->name().c_str(), loft_v->iter(),  ymin()); cout <<  Xmin() << endl; };
 };

				 template<typename V>
struct	trust_region_minimizer : minimizer<V>    { 
					LOFT_TYPES;  MINIMIZER_MEMBERS;
			fp_t 				rho_begin_, rho_end_;

	explicit		trust_region_minimizer		(const char* _name= "unknown (trust region type)")
	:	minimizer<V>(_name),
		rho_begin_ 	(1.0),
		//rho_end_   	(0.1)
		rho_end_   	(numeric_limits<fp_t>::min()*1000)
									//rho_begin_ 	(numeric_limits<fp_t>::quiet_NaN ()),
									//rho_end_   	(numeric_limits<fp_t>::quiet_NaN ())
	{};

	virtual minimizer<V>&		rho_begin		(fp_t rho)	{ rho_begin_ = rho;  return *this; };
	virtual minimizer<V>&		rho_end			(fp_t rho)	{ rho_end_   = rho;  return *this; };
 };


	// if app inludes <lopti.h>,  then include all
//	#ifndef MINIMIZER
//		#include	<lopti/condor-wrap.h>
//		#include   	<lopti/newuoa-wrap.h>
//		#include	<lopti/gsl-nelder-mead-wrap.h>
//		#include	<lopti/hook-jeevs.h>
//	#endif 

	} // namespace lopti
	#endif  // LOPTI_H
