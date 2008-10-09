
#include <lvv/lvv.h>
#include <limits>
	using 	std::numeric_limits;

                 template<typename V>
class	minimizer { 

	public:
		typedef  typename V::value_type fp_t;
		typedef  fp_t (*of_ptr_t)(V&, void*);
		int				max_iter_;
		bool				verbose_;
		V				X;
		V				Xmin;
		fp_t				ymin_;
		int				iter_;
		of_ptr_t			of_;

	minimizer		(of_ptr_t of,  V& _X)       
	:
		max_iter_  (500),
		ymin_      (numeric_limits<fp_t>::quiet_NaN ()),
		iter_      (-1),
		X          (_X),
		verbose_   (false),
		of_        (of)
	{};

	virtual 		~minimizer		()		{};  // it it here so that approprite polimorfic DTOR called 

	virtual void		max_iter		(int mx)	{ max_iter_   = mx;  };
	virtual fp_t 	 	ymin			()		{  return ymin_; };
	virtual fp_t 	 	iter			()		{  return iter_; };

	virtual void 	 verbose		(bool flag)	{
		#define  GP_F "splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		cout << "# :gnuplot: set view 0,0,1.7;   set font \"arial,6\"; set dgrid3d;  set key off;"
			"  set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 3:4:2:1 with labels; \n";

		verbose_ = flag;
	};

	virtual V&		 argmin() = 0;
};

                 template<typename V>
class	trust_region_minimizer: public minimizer<V> { 
		typedef  typename V::value_type fp_t;
		typedef  fp_t (*of_ptr_t)(V&, void*);

	public:
		fp_t 				rho_begin_;	// r(rho) start
		fp_t 				rho_end_;	// r end

	trust_region_minimizer		(of_ptr_t of, V& _X):  
		minimizer<V>	(of, _X),
		rho_begin_ 	(numeric_limits<fp_t>::quiet_NaN ()),
		rho_end_   	(numeric_limits<fp_t>::quiet_NaN ())
	{};

	virtual void		rho_begin		(fp_t rho)		{ rho_begin_ = rho; };
	virtual void		rho_end			(fp_t rho)		{ rho_end_   = rho; };
};

