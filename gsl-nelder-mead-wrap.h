
	#undef  NM
	#define	LOPTI_NM
	#define	MINIMIZER	nelder_mead_minimizer

	#include <lopti/lopti.h>
	#include <stdlib.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_multimin.h>

	#include <lvv/array.h>
		using lvv::array;
	#include <lvv/lvv.h>
		using std::cerr;

 
                 template<typename V>
class gsl_of_wrap {  public:
		typedef  typename V::value_type fp_t;
		typedef  function<fp_t(V&)>	of_ptr_t;

		static		function<typename V::value_type(V&)>  of_ptr;
		static void	init	(of_ptr_t of) 				{ of_ptr = of; }	
		static double	eval	(const gsl_vector* gX, void * var)	{ V X;    X << gX;    return  of_ptr(X); }	
 };

template<typename V>  typename boost::function<typename V::value_type(V&)>   gsl_of_wrap<V>::of_ptr; // this is in gsl_ow_wrap class, but we need to decl it 1 more time for compiler

//template<typename V>  typename V::value_type  (*gsl_of_wrap<V>::of_ptr)(V&); // this is in gsl_ow_wrap class, but we need to decl it 1 more time for compiler
//template<typename V>  typename V::value_type  (*gsl_of_wrap<V>::of_ptr)(V&, void*); // this is in gsl_ow_wrap class, but we need to decl it 1 more time for compiler

                 template<typename V>
class	nelder_mead_minimizer  :  public minimizer<V> { public:
		typedef  typename minimizer<V>::fp_t		fp_t;
		typedef  typename minimizer<V>::of_ptr_t	of_ptr_t;

	using minimizer<V>::X;  		// without this we woun't see minimizer members
	using minimizer<V>::max_iter_;
	using minimizer<V>::verbose_;
	using minimizer<V>::of_;
	using minimizer<V>::verbose_;
	using minimizer<V>::ymin_;
	using minimizer<V>::Xmin;
	using minimizer<V>::iter_;
	using minimizer<V>::found_;

	virtual const char*	name			()	const 		{ return "nelder-mead (aka simplex)"; };


		gsl_vector*			gX;	
		gsl_vector* 			gS;	
		const gsl_multimin_fminimizer_type *gsl_minimizer_type_ptr;
		gsl_multimin_fminimizer*	gsl_minimizer;
		gsl_multimin_function		minex_func;
		double 				characteristic_size;



	explicit 		nelder_mead_minimizer	(of_ptr_t of, V& _X)
	:		minimizer<V>(of, _X),
			gsl_minimizer_type_ptr (gsl_multimin_fminimizer_nmsimplex)
	{
			gsl_of_wrap<V>::init(of);
			//minex_func.f = &gsl_of_wrap<V>::eval2;
			minex_func.f = &gsl_of_wrap<V>::eval;
			minex_func.n = X.size();
			//minex_func.params = var;
			minex_func.params = NULL;
			gX  = gsl_vector_alloc(X.size());
			gX << X;
	};

	~nelder_mead_minimizer () { gsl_multimin_fminimizer_free(gsl_minimizer);  gsl_vector_free(gX);   gsl_vector_free(gS);  };

	void		step			(V S)		{
		gS = gsl_vector_alloc(V::size()); 
		gS << S; 
		gsl_minimizer = gsl_multimin_fminimizer_alloc(gsl_minimizer_type_ptr, V::size());
		//xmin = minimizer->x->data;
		if (verbose_)   cerr << "#  minimizer: "  <<  gsl_multimin_fminimizer_name (gsl_minimizer)  <<  endl;
		gsl_multimin_fminimizer_set(gsl_minimizer, &minex_func, gX  , gS);
	};

	void 		gsl_var			(void* var )	{ minex_func.params = var; };
	void   		gsl_characteristic_size	(double cs)	{ characteristic_size = cs; };

	fp_t 	 	ymin			()		{  return gsl_minimizer->fval; };
	fp_t 	 	iter			()		{  return iter_; };

	V&		argmin () {
		int	test_status=GSL_CONTINUE;			// test_status:  GSL_SUCCESS==0; GSL_CONTINUE==-2; 

		if (verbose_)  FMT("# Itr  %10t Y   %20t  X[0..]   Step\n");
	
		for  ( iter_=0;  iter_ < max_iter_  &&  (test_status==GSL_CONTINUE);   ++iter_ )   {

			int  iterate_status = gsl_multimin_fminimizer_iterate(gsl_minimizer);

			if (iterate_status) {
				cerr << "error: FMinimizer: in gsl_multimin_fminimizer_iterate()\n";
				break;
			}

			double size = gsl_multimin_fminimizer_size(gsl_minimizer);
			test_status = gsl_multimin_test_size(size, characteristic_size);

			if (test_status == GSL_SUCCESS)  {
				found_=true;
				//fmin=gsl_minimizer->fval;
				if (verbose_ )    cerr  << "# converged to minimum at\n";
			}
			
			if (verbose_) { 
				FMT("%5d \t %18.10g \t ") %iter_   %(gsl_minimizer->fval);
				for (int i=0; i < gsl_minimizer->x->size; ++i) FMT("%18.10d") %gsl_vector_get(gsl_minimizer->x, i);
				FMT("%18.10g")  %size;
				cout << endl;
			}
		}

		Xmin << gsl_minimizer->x;
		return Xmin;
	}


	void 	 verbose		(bool flag)	{
		verbose_ =flag;
		#define  GP_F  "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";
	};

 };
