#include <lvv/lvv.h>
#include <string>

int		f(string s)	{ cout  << "function: " << s << " =  "; return 0;}

class ObjF { public:
	//int		static operator()(double x, double y)	{ cout << x+y << " objf.operator() \n"; return 0; };	// error: must be non static
	int		static static_mem_f(string s)	{ cout << s; return 0; };
	int		mem_f(string s)	{ cout << s << " ObjF::mem_f()"; return 0; };
};


				template <typename FuncT> // <- What would go here
class		fwrap { public:
		fwrap		(FuncT func) : func_(func) 	{}
	int 	eval		() 				{ return func_(1, 2); }
		FuncT func_;
};

struct fobj_t { int operator()(int i) { cout << "fobj_t: ";  return i+100; }; }; // can not be inside main()

#include <boost/function.hpp>
#include <functional>
using namespace std;
using namespace boost;

int main() {
	
	//////////////////////////////////////////////////////////////////////////////////////////   BOOST
	
	// FUNCTOR
	boost::function<int(int i)>   	bf_fct = fobj_t();
	cout << bf_fct(11) << endl;


	// PLAIN F()
	function<int(string s)>       	bf_f = &f;
	bf_f = &f;
	cout << bf_f("plain f()") << endl;
	
	// PLAIN static mem_f()
	function<int(string s)>       	bf_smf = &ObjF::static_mem_f;
	cout << bf_smf("ObjF::static_mem_f = ") << endl;

	// MEMEBER-F
	//function<int(string s)>       	bf_mf = &ObjF::mem_f;
	//ObjF objf;
	//cout << bf_mf(&objf, "objf") << endl;

	return 0;

}
