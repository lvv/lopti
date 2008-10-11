#include <lvv/lvv.h>
#include <string>

#include <functional>
using namespace std;
using namespace boost;
#include <boost/function.hpp>
using boost::function;

//////////////////////////////////////////////////////////////////////////////////////////////////// TEST CALLABLE OBJECT
int				plain_f		(string s)	{ cout  << "function: " << s << " =  "; return 0;}

struct	functor_t {
				functor_t	()		: value(0) {};
	void 			init		(int  i)	{ value = i; };
	int			operator()	(string s)	{ cout << s;  return value; };
	int value;
};

struct empty {}; 
template< typename RET=empty, typename ARG=empty> 
struct eval_count;

template<typename RET, typename ARG>        struct	eval_count;

	   template<typename RET, typename ARG>                                                                                                          
struct	eval_count<RET(ARG)> : public unary_function<RET,ARG>		{ 
	eval_count		(function<RET (ARG a)>    _of)	: of(_of), eval_cnt(0)	{ cout << "\neval_cnt_wrap CTor\n"; };
	//~eval_count		()		{ cout << "\neval_cnt_wrap dtOR:  object function called  total times: " << eval_cnt << endl;};
	RET	operator()	(ARG    arg)	{ FMT("eval_cnt_wrap: cnt=%d   arg=%s.   Executing:\n ") %eval_cnt %arg;  eval_cnt++; return of(arg);};
	int eval_cnt; 
	function<RET (ARG)> of;
}; 

class	obj_t 	{ public:
	int	static		static_mem_f	(string s)	{ cout << s; return 0; };
	int			mem_f		(string s)	{ cout << s << " ObjF::mem_f() = "; return 0; };
};



int main() {
	
	// FUNCTOR 
	functor_t	fct1;							cout << fct1("stand alone  fct1() ") << endl;

	// FUNCTOR  with init fun
	functor_t	fct2;	 						cout << fct2("stand alone  fct2 ") << endl;
			fct2.init(33);						cout << fct2("stand alone  fct2.init(33) ") << endl;
	
	// FUNCTOR  WRAP
	eval_count<int(string)>		ecw(fct2);				cout << ecw("boost: ecw(fct2)-1 ") << endl;
										cout << ecw("boost: ecw(fct2)-2 ") << endl;

	function<int(string s)> 	ecw_p = eval_count<int(string)>(fct2);	cout << ecw_p("boost: ecw ") << endl;

	// F POINTERS
	function<int(string s)>   		fct1_p = functor_t();		cout << fct1_p("boost: functor_t ") << endl;
	function<int(string s)>   		fct2_p = fct2;			cout << fct2_p("boost: fct2 ") << endl;

										cout << ecw_p("boost: ecw_p ") << endl;
										cout << ecw_p("boost: ecw_p ") << endl;

	// PLAIN F()
	function<int(string s)>       		bf_f = &plain_f;		cout << bf_f("boost: plain f()") << endl;
	
	// PLAIN static mem_f()
	function<int(string s)>			bf_smf = &obj_t::static_mem_f;	cout << bf_smf("boost: ObjF::static_mem_f = ") << endl;

	// MEMEBER-F
		//boost::function<int(string s)>       	bf_mf = &obj_t::mem_f;
	obj_t o;
	boost::function<int(string s)>       		bf_mfa_mf = bind1st( mem_fun(&obj_t::mem_f), &o);
	boost::function<int(obj_t* o, string s)>  	bf_mfa_mf2 = &obj_t::mem_f;
	cout << bf_mfa_mf  (    "boost:  bind1st( mem_fun(&obj_t::mem_f), &o)  ") << endl;
	cout << bf_mfa_mf2 (&o, "boost:  &obj_t::mem_f                         ") << endl;


	return 0;

}
