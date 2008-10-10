#include <lvv/lvv.h>
#include <string>

//////////////////////////////////////////////////////////////////////////////////////////////////// TEST FUNCTION
int				plain_f		(string s)	{ cout  << "function: " << s << " =  "; return 0;}
struct	functor_t { int		operator()	(string s)	{ cout << s;  return 77; }; }; // can not be inside main()

class	obj_t 	{ public:
	int	static		static_mem_f	(string s)	{ cout << s; return 0; };
	int			mem_f		(string s)	{ cout << s << " ObjF::mem_f()"; return 0; };
};



////////////////////////////////////////////////////////////////////////////////////////////////////  MY WRAP
// my wrap for condor
				template <typename FuncT> // <- What would go here
struct		fwrap {
		fwrap		(FuncT func) : func_(func) 	{}
	int 	eval		() 				{ return func_(" test "); }
		FuncT func_;
};


#include <boost/function.hpp>
#include <functional>
using namespace std;
using namespace boost;

int main() {
	
	//////////////////////////////////////////////////////////////////////////////////////////// MY WRAP TEST
	//fwrap<int (*)(string s)>        	t1	(&plain_f);		t1.eval("lvv fwrap: ");
					//fwrap<int (*)(string s)> 		t2	(&obj_t::static_mem_f); t2.eval("lvv fwrap: ");
					//fwrap<typeof(&obj_t::static_mem_f)>	t3	(&obj_t::static_mem_f);	t3.eval("lvv fwrap: ");
					
					//fwrap<int (ObjF::operator())(double, double)>         t4(Func2);      t4.eval();
					
					//fwrap<int(ObjF::*)(double, double)>           t4(objf.*operator());
					// no err:  fwrap<int(ObjF::*)(double, double)>         t4(&ObjF::operator());
	//fwrap<int(*)(string s)>          	t4(&obj_t::static_mem_f);		t4.eval("lvv fwrap: ");
	//fwrap<typeof(&obj_t::eval)>		t5(&obj_t::static_mem_f);		t5.eval("lvv fwrap: " );

	//////////////////////////////////////////////////////////////////////////////////////////   BOOST  TEST
	
	// FUNCTOR
	boost::function<int(string s)>   	bf_fct = functor_t();		cout << bf_fct("boost: functor ") << endl;


	// PLAIN F()
	function<int(string s)>       		bf_f = &plain_f;		cout << bf_f("boost: plain f()") << endl;
	
	// PLAIN static mem_f()
	function<int(string s)>			bf_smf = &obj_t::static_mem_f;	cout << bf_smf("boost: ObjF::static_mem_f = ") << endl;

	// MEMEBER-F
	//function<int(string s)>       	bf_mf = &ObjF::mem_f;
	function<int(string s)>       	bf_mfa_mf = mem_fun_ref(&obj_t::mem_f);
	obj_t o;
	//cout << bf_mfa_mf(o, "boost: o/mem_fun_ref ") << endl;

	return 0;

}
