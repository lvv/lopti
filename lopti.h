	
	#ifndef LOPTI_H
	#define LOPTI_H

	#include <lvv/lvv.h>
	#include <lvv/array.h>
		using lvv::array;

	#include <boost/format.hpp>
		using boost::format;
	//#include <cstdio>
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
							using minimizer<V>::ymin_;		\
							using minimizer<V>::Xmin_;		\
							using minimizer<V>::name;		\
							using minimizer<V>::name_;		\
							using minimizer<V>::found_;		\
							using minimizer<V>::objective_v;

						 #define         TR_MINIMIZER_MEMBERS    \
							using trust_region_minimizer<V>::rho_begin_;	\
							using trust_region_minimizer<V>::rho_end_;


					 template<typename V>
struct	minimizer {
					OBJECTIVE_TYPES;

			//objective_p_t			objective_v;	
			shared_ptr<objective_base<V>>	objective_v;	
			int				max_iter_;
			bool				verbose_;
			V				X;
			V				Xmin_;
			T				ymin_;
			int				iter_;
			bool				found_;
			string				name_;

				explicit 		
	minimizer		(const string& _name = "unknown")  
	:
		name_		(_name),
		max_iter_	(10000),
		ymin_    	(numeric_limits<T>::quiet_NaN ()),
		iter_    	(0),
		verbose_ 	(false),
		found_ 		(false)
	{};

	//virtual 			~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 


	// set-ters
	//virtual minimizer<V>&		objective			(objective_cref_t  ref)	{  objective_v = &ref.clone();   return *this;  };
	virtual minimizer<V>&		objective			(objective_cref_t  ref)	{  objective_v = shared_ptr<objective_base<V>>(&ref.clone());   return *this;  };
	//virtual void			objective			(objective_cref_t objective_cref)	{ wrapped_objective_v = shared_ptr<objective_base<V>>(&objective_cref.clone());	name_ = objective_cref.name(); };
	virtual minimizer<V>&		x0			(V& _X) 		{  X  = _X;	return *this;  };
	virtual minimizer<V>&		max_iter		(int mx)		{  max_iter_   = mx;	return *this;  };
	virtual minimizer<V>&		verbose			(bool flag)		{  verbose_ = flag;	return *this;  };


	// get-ters
	virtual const string		name			() 	const	{  return (format("%s-%d")  %name_  %(V::size()) ).str(); };
	/*virtual const string		name			() 	const	{ 
			char	buf[100];
			sprintf(buf, "%s-%d",  name_.c_str(), V::size());
			return string(buf);
	};*/
	virtual T 	 		ymin			()	const	{  return ymin_; };
	virtual V 	 		Xmin			()	const	{  return Xmin_; };
	virtual T 	 		iter			()	const	{  return iter_; };
	virtual bool			found			() 	const	{  return found_; };

	// do-ers
	virtual V&			argmin			() 		{  return Xmin_; };
	virtual void			print			()		{  cout << format("%s(%s)  %35t  iter=%d  \t ymin=%g \t Xmin=%22.15g \n") %name() %objective_v->name() %objective_v->iter()  %ymin()  %Xmin();};
	/*virtual void			print			()		{
		printf("%s(%s)  	  iter=%d  	 ymin=%g 	 Xmin: ",
			name().c_str(), 
			objective_v->name().c_str(),
			iter_, 
			ymin()
		);
		cout <<  Xmin() << endl; }; */
 };

				 template<typename V>
struct	trust_region_minimizer : minimizer<V>    { 
					OBJECTIVE_TYPES;  MINIMIZER_MEMBERS;
			T 				rho_begin_, rho_end_;

	explicit		trust_region_minimizer		(const char* _name= "unknown (trust region type)")
	:	minimizer<V>(_name),
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
