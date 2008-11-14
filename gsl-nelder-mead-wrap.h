
	#ifndef LOPTI_NELDER_MEAD_H
	#define LOPTI_NELDER_MEAD_H

	#ifndef		MINIMIZER
		#define	MINIMIZER	gsl_nelder_mead_minimizer
	#endif

	#ifndef		V_IBEGIN
		#define	V_IBEGIN 0
	#endif
	

	#include <lvv/convert-gsl.h>
		using lvv::array;
	#include <lopti/lopti.h>
	#include <stdlib.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_multimin.h>

	#include <lvv/lvv.h>
		using std::cerr;

	
	namespace lopti {
 
 template<typename V>		struct gsl_of_wrap {
				LOFT_TYPES;
				static		loft_p_t	 loft_v;
		static void	init	(loft_p_t  _loft_v)			{ loft_v = _loft_v; }	
		static double	eval	(const gsl_vector* gX, void * var)	{ V X;    X << gX;    return  (*loft_v)(X); }	
 };

template<typename V>  loft_base<V>* gsl_of_wrap<V>::loft_v; // this is in gsl_of_wrap class, but we need to decl it 1 more time for compiler

                 template<typename V>
struct	gsl_nelder_mead_minimizer  :  minimizer<V> {
						MINIMIZER_MEMBERS;   LOFT_TYPES;
		gsl_vector*			gX;	
		gsl_vector* 			gS;	
		const gsl_multimin_fminimizer_type *gsl_minimizer_type_ptr;
		gsl_multimin_fminimizer*	gsl_minimizer;
		gsl_multimin_function		minex_func;
		double 				characteristic_size_;

	explicit 		gsl_nelder_mead_minimizer	() : minimizer<V>("nelder-mead"), gsl_minimizer_type_ptr (gsl_multimin_fminimizer_nmsimplex) {};

	~gsl_nelder_mead_minimizer () { gsl_multimin_fminimizer_free(gsl_minimizer);  gsl_vector_free(gX);   gsl_vector_free(gS);  };

	virtual minimizer<V>&		step0			(V& S)		{
		gS = gsl_vector_alloc(V::size()); 
		gS << S; 
		gsl_minimizer = gsl_multimin_fminimizer_alloc(gsl_minimizer_type_ptr, V::size());
		if (verbose_)   cerr << "#  minimizer: "  <<  gsl_multimin_fminimizer_name (gsl_minimizer)  <<  endl;
		return *this;
	};

	//void 		gsl_var			(void* var )	{ minex_func.params = var; };
	virtual  minimizer<V>&		characteristic_size	(fp_t cs)	{ characteristic_size_ = cs;  return *this; };

	virtual V&		argmin () {
		
		////  gsl init
		gsl_of_wrap<V>::init(loft_v);
		minex_func.f = &gsl_of_wrap<V>::eval;
		minex_func.n = X.size();
		//minex_func.params = var;
		minex_func.params = NULL;
		gX  = gsl_vector_alloc(X.size());
		gX << X;
		gsl_multimin_fminimizer_set(gsl_minimizer, &minex_func, gX  , gS);

		///// main cicle
		int	test_status=GSL_CONTINUE;			// test_status:  GSL_SUCCESS==0; GSL_CONTINUE==-2; 

		if (verbose_)  FMT("# Itr  %10t Y   %20t  X[0..]   Step\n");
	
		for  ( iter_=0;  iter_ < max_iter_  &&  (test_status==GSL_CONTINUE);   ++iter_ )   {

			int  iterate_status = gsl_multimin_fminimizer_iterate(gsl_minimizer);

			if (iterate_status) {
				cerr << "error: FMinimizer: in gsl_multimin_fminimizer_iterate()\n";
				break;
			}

			double size = gsl_multimin_fminimizer_size(gsl_minimizer);
			test_status = gsl_multimin_test_size(size, characteristic_size_);

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

		ymin_ = gsl_minimizer->fval;
		Xmin_  <<  gsl_minimizer->x;
		return Xmin_;
	}

 };
	}	// namespace  lopti
	#endif // LOPTI_NELDER_MEAD_H
