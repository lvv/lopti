	
	#ifndef LOPTI_H
	#define LOPTI_H

	#include <lvv/lvv.h>
	#include <lvv/array.h>
		using lvv::array;

	#include <limits>
		using 	std::numeric_limits;

	#include <lopti/object_function.h>

namespace lopti {
				
						 #define         MINIMIZER_MEMBERS    \
							using minimizer<V>::X;  		\
							using minimizer<V>::x0;			\
							using minimizer<V>::max_iter_;		\
							using minimizer<V>::iter_;		\
							using minimizer<V>::verbose_;		\
							using minimizer<V>::ymin_;		\
							using minimizer<V>::Xmin_;		\
							using minimizer<V>::name;		\
							using minimizer<V>::name_base;		\
							using minimizer<V>::found_;		\
							using minimizer<V>::objective_v;

						 #define         TR_MINIMIZER_MEMBERS    \
							using trust_region_minimizer<V>::rho_begin_;	\
							using trust_region_minimizer<V>::rho_end_;


					 template<typename V>
struct	minimizer {
					OBJECTIVE_TYPES;

			objective_p_t			objective_v;	
			int				max_iter_;
			bool				verbose_;
			V				X;
			V				Xmin_;
			T				ymin_;
			int				iter_;
			bool				found_;
			string				name_base;

				explicit 		
	minimizer		()
	:
		name_base	(""),	
		max_iter_	(10000),
		ymin_    	(numeric_limits<T>::quiet_NaN ()),
		iter_    	(0),
		verbose_ 	(false),
		found_ 		(false)
	{};

	//virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters

        virtual minimizer<V>&		objective		(objective_cref_t ref)	{  objective_v = objective_p_t(&ref.clone());   return *this; };

	virtual minimizer<V>&		x0			(V& _X) 		{  X  = _X;		return *this;  };
	virtual minimizer<V>&		max_iter		(int mx)		{  max_iter_   = mx;	return *this;  };
	virtual minimizer<V>&		verbose			(bool flag)		{  verbose_ = flag;	return *this;  };


	// get-ters
	virtual const string		name			() 	const	{  return mk_name(""); };
		const string		mk_name			(const string&& name_base) 	const	{  return ({ std::ostringstream oss; oss << name_base << V::size(); oss.str();}); };
	virtual T 	 		ymin			()	const	{  return ymin_; };
	virtual V 	 		Xmin			()	const	{  return Xmin_; };
	virtual T 	 		iter			()	const	{  return iter_; };
	virtual bool			found			() 	const	{  return found_; };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_; };
	virtual void			print			()		{  cout << name() << "(" << objective_v->name() << ") \t iter/max=" << objective_v->iter()  << "/" << max_iter_  << "\t ymin=" << ymin()  << "\t Xmin=" << Xmin(); };
 };

				 template<typename V>
struct	trust_region_minimizer : minimizer<V>    { 
					OBJECTIVE_TYPES;  MINIMIZER_MEMBERS;
			T 				rho_begin_, rho_end_;

	explicit		trust_region_minimizer		(const char* _name= "unknown (trust region type)")
	:	//minimizer<V>(_name),
		rho_begin_ 	(1.0),
		//rho_end_   	(0.1)
		rho_end_   	(numeric_limits<T>::min()*1000)
									//rho_begin_ 	(numeric_limits<T>::quiet_NaN ()),
									//rho_end_   	(numeric_limits<T>::quiet_NaN ())
	{};

	virtual minimizer<V>&		rho_begin		(T rho)	{ rho_begin_ = rho;  return *this; };
	virtual minimizer<V>&		rho_end			(T rho)	{ rho_end_   = rho;  return *this; };
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
