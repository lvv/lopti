
#include <stdlib.h>
//#include <iostream>
using std::cerr;
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <lvv/lvv.h>

 
                 template<typename V>
class gsl_of_wrap {  public:
	typedef  typename V::value_type v;
	typedef  v (*of_ptr_t)(V ) ;	
	static of_ptr_t  of_ptr;

	gsl_of_wrap(){};


	static void  init(v (*F)(V)) { of_ptr = F; }	

	static double  eval(const gsl_vector* gX, void * var) {
		V	X;
		for (int i = X.ibegin();  i<X.iend();  i++)     X[i] = gsl_vector_get(gX, i);
		return  (*of_ptr)(X);
	}	
};

template<typename V>  typename V::value_type  (*gsl_of_wrap<V>::of_ptr)(V);
//template<typename V>  typename gsl_of_wrap<V>::of_ptr_t  gsl_of_wrap<V>::of_ptr;

                 template<typename V>
class	minimizer { public:
		typedef  typename V::value_type v;
		int				max_iter_;
		bool				gnuplot_print_;
		V				Xmin;
		int				iter_;

		gsl_vector*			gX;	
		gsl_vector* 			gS;	
		const gsl_multimin_fminimizer_type *T;
		gsl_multimin_fminimizer*	gsl_minimizer;
		gsl_multimin_function		minex_func;
		double 				characteristic_size;
		bool				found_;
		//gsl_of_wrap<V>::of_ptr;
		//double*			xmin;
		//double			fmin;

	minimizer		(v (*F)(V), V X)     
	:	
		max_iter_(300),
		T (gsl_multimin_fminimizer_nmsimplex)
	{
		gsl_of_wrap<V>::init(F);
		minex_func.f = &gsl_of_wrap<V>::eval;
		minex_func.n = X.size();
		gX  = gsl_vector_alloc(X.size());	for (int i=0; i<X.size(); i++)  gsl_vector_set(gX, i, X[i]);
	

	};

	~minimizer () { gsl_multimin_fminimizer_free(gsl_minimizer);  gsl_vector_free(gX);   gsl_vector_free(gS);  };

	void		step			(V S)		{
		gS = gsl_vector_alloc(V::size());  for (int i=0; i<V::size(); i++)  gsl_vector_set(gS, i, S[i]); 
		gsl_minimizer = gsl_multimin_fminimizer_alloc(T, V::size());
		//xmin = minimizer->x->data;
		if (gnuplot_print_)   cerr << "#  minimizer: "  <<  gsl_multimin_fminimizer_name (gsl_minimizer)  <<  endl;
		gsl_multimin_fminimizer_set(gsl_minimizer, &minex_func, gX  , gS);
	};

	void 		max_iter		(int mx)	{ max_iter_ 	= mx;	 };
	void 		gsl_var			(void* var )	{ minex_func.params = var; };
	void   		gsl_characteristic_size	(double cs)	{ characteristic_size = cs; };


	v 	 	ymin			()		{  return gsl_minimizer->fval; };
	v 	 	iter			()		{  return iter_; };


	V&		argmin () {

		int	test_status=GSL_CONTINUE;			// test_status:  GSL_SUCCESS==0; GSL_CONTINUE==-2; 

		if (gnuplot_print_)  FMT("# Itr  %10t Y   %20t  X[0..]   Step\n");
	
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
				if (gnuplot_print_ )    cerr  << "# converged to minimum at\n";
			}
			
			if (gnuplot_print_) { 
				FMT("%5d \t %18.10g \t ") %iter_   %(gsl_minimizer->fval);
				for (int i=0; i < gsl_minimizer->x->size; ++i) FMT("%18.10d") %gsl_vector_get(gsl_minimizer->x, i);
				FMT("%18.10g")  %size;
				cout << endl;
			}
		}

		for (int i = Xmin.ibegin();  i<Xmin.iend();  i++)     Xmin[i] = gsl_vector_get(gsl_minimizer->x, i);
		return Xmin;
	}


	void 	 gnuplot_print		(bool flag)	{
		gnuplot_print_ =flag;
		#define  GP_F  "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";
	};

 };
