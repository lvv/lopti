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

	explicit 		minimizer		(of_ptr_t of,  V& _X)       
	:
		max_iter_	(500),
		ymin_    	(numeric_limits<fp_t>::quiet_NaN ()),
		iter_    	(-1),
		X        	(_X),
		verbose_ 	(false),
		of_      	(of),
		found_ 		(false)
	{};

	virtual 		~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 

	virtual void		max_iter		(int mx)	{  max_iter_   = mx;  };
	virtual fp_t 	 	ymin			()	const	{  return ymin_; };
	virtual V 	 	Xmin			()	const	{  return Xmin_; };
	virtual fp_t 	 	iter			()	const	{  return iter_; };

	virtual void 	 verbose		(bool flag)	{
		#define  GP_F "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";

		verbose_ = flag;
	};

	virtual V&		argmin			() 		= 0;
	virtual bool		found			() 	const	{ return found_;};
	virtual const char*	name			() 	const	{ return "n/a";};
	virtual void		print			()		{ MSG("%s %25t iter=%d  \t ymin=%g \t Xmin%g \n") %name()  %iter()  %ymin()  %Xmin();};
};

                 template<typename V>
class	trust_region_minimizer : public minimizer<V>    { public:

		typedef  typename minimizer<V>::fp_t		fp_t;
		typedef  typename minimizer<V>::of_ptr_t	of_ptr_t;

		fp_t 				rho_begin_;	// r(rho) start
		fp_t 				rho_end_;	// r end

	trust_region_minimizer		(of_ptr_t of, V& _X):  
		minimizer<V>	(of, _X),
		rho_begin_ 	(numeric_limits<fp_t>::quiet_NaN ()),
		rho_end_   	(numeric_limits<fp_t>::quiet_NaN ())
	{};

	virtual void		rho_begin		(fp_t rho)	{ rho_begin_ = rho; };
	virtual void		rho_end			(fp_t rho)	{ rho_end_   = rho; };
	virtual const char*	name			() 	const	{ return "trust region type"; };
};

	// of prog inlude <lopti.h> then include all
	#ifndef MINIMIZER
		#include	<lopti/condor-wrap.h>
		#include   	<lopti/newuoa-wrap.h>
		#include	<lopti/gsl-nelder-mead-wrap.h>
	#endif 

	#endif 
