#include <lvv/lvv.h>
#include <lvv/array.h>
	using lvv::array;
#include <string>

#include <functional>
	using std::binder1st;
	using std::unary_function;
#include <boost/function.hpp>
	using boost::function;
#include <boost/bind.hpp>
	using boost::bind;

	//using namespace std;
	//using namespace boost;
	//
//////////////////////////////////////////////////////////////////////////////////////////////////// TEST CALLABLE OBJECT
typedef array<double,2,1> array_t;

int				plain_f		(string s)		{ cout  << "plain_f: " << s << " =  "; return 0;}
int				plain_f2	(string s, int i)	{ cout  << "plain_f2: " << s << " =  "; return i;}
int				plain_fa	(array_t A)		{ cout  << "plain_fA: " << A << " =  "; return 0;}
template<typename T> T		plain_fa	(array_t A)		{ cout  << "plain_fA: " << A << " =  "; return 0;}

struct	functor_t {
				functor_t	()		: value(0) {};
	void 			init		(int  i)	{ value = i; };
	int			operator()	(string s)	{ cout << s;  return value; };
	int value;
};


class	obj_t 	{ public:
	int	static		static_mem_f	(string s)	{ cout << s; return 0; };
	int			mem_f		(string s)	{ cout << s << " ObjF::mem_f() = "; return 0; };
};


////////////////////////////////////////////////////////////////////////////////////////////////  TEST WRAPPER
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

int main() {
	
	// FUNCTOR 
	functor_t	fct1;							cout << fct1("stand alone  fct1() ") << endl;

	// FUNCTOR  with init fun
	functor_t	fct2;	 						cout << fct2("stand alone  fct2 ") << endl;
			fct2.init(33);						cout << fct2("stand alone  fct2.init(33) ") << endl;
	function<int(string s)> 	fct2p = fct2;				cout << fct2p("fct2p ") << endl;
					fct2p->init(111);			cout << fct2p("fct2p 2") << endl;
	
	// FUNCTOR  WRAP
	eval_count<int(string)>		ecw(fct2);				cout << ecw("boost: ecw(fct2)-1 ") << endl;
										cout << ecw("boost: ecw(fct2)-2 ") << endl;

	function<int(string s)> 	ecw_p = eval_count<int(string)>(fct2);	cout << ecw_p("ecw_p ") << endl;
										cout << ecw_p("ecw_p ") << endl;
										cout << ecw_p("ecw_p ") << endl;


	// F POINTERS
	function<int(string s)>   		fct1_p = functor_t();		cout << fct1_p("fct1_p = functor_t() ") << endl;
	function<int(string s)>   		fct2_p = fct2;			cout << fct2_p("fct2_p = fct2 ") << endl;

	// PLAIN F() with arrays
	array_t  X = {{11, 22}};
	function<int(array_t&)>       		bf_fa;
	bf_fa = plain_fa<int>;						cout << "boost: plain fa()" << bf_fa(X) << endl;
	//bf_f = bind2nd(function<int(string,int)>(&plain_f2),100);	cout << bf_f("sdt::bind2nd plain-f2()") << endl;
	//bf_f = bind (plain_f2, _1, 100);				cout << bf_f("boost::bind  plain-f2()") << endl;

	// PLAIN F()
	function<int(string s)>       		bf_f;
	bf_f = plain_f;						cout << bf_f("boost: plain f()") << endl;
	bf_f = bind2nd(function<int(string,int)>(&plain_f2),100);	cout << bf_f("sdt::bind2nd plain-f2()") << endl;
	bf_f = bind (plain_f2, _1, 100);				cout << bf_f("boost::bind  plain-f2()") << endl;
	
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
