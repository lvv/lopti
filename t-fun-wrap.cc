#include <lvv/lvv.h>
#include <string>

//////////////////////////////////////////////////////////////////////////////////////////////////// TEST FUNCTION
int				plain_f		(string s)	{ cout  << "function: " << s << " =  "; return 0;}
struct	functor_t { int		operator()	(string s)	{ cout << s;  return 77; }; }; // can not be inside main()

class	obj_t 	{ public:
	int	static		static_mem_f	(string s)	{ cout << s; return 0; };
	int			mem_f		(string s)	{ cout << s << " ObjF::mem_f()"; return 0; };
};



#include <boost/function.hpp>
#include <functional>
//using namespace std;
//using namespace boost;

int main() {
	

	//////////////////////////////////////////////////////////////////////////////////////////   BOOST  TEST
	
	// FUNCTOR
	boost::function<int(string s)>   	bf_fct = functor_t();		cout << bf_fct("boost: functor ") << endl;


	// PLAIN F()
	boost::function<int(string s)>       		bf_f = &plain_f;		cout << bf_f("boost: plain f()") << endl;
	
	// PLAIN static mem_f()
	boost::function<int(string s)>			bf_smf = &obj_t::static_mem_f;	cout << bf_smf("boost: ObjF::static_mem_f = ") << endl;

	// MEMEBER-F
		//boost::function<int(string s)>       	bf_mf = &obj_t::mem_f;
	obj_t o;
	boost::function<int(string s)>       		bf_mfa_mf = bind1st( mem_fun(&obj_t::mem_f), &o);
	cout << bf_mfa_mf("boost: o/mem_fun_ref ") << endl;
		//cout << bf_mfa_mf(o, "boost: o/mem_fun_ref ") << endl;

	//////////////////////////////////////////////////////////////////////////////////////////   BOOST  TEST
	
	/*
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
	*/

	return 0;

}
