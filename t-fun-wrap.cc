#include <lvv/lvv.h>
#include <functional>
using namespace std;

int		f1(int x, int y)	{ cout << x+y << " f1 \n"; return 0; }
int		f2(double x, double y)	{ cout << x+y << " f2 \n"; return 0;}

class ObjF { public:
	//int		static operator()(double x, double y)	{ cout << x+y << " objf.operator() \n"; return 0; };
	int		static eval(double x, double y)	{ cout << x+y << " ObjF::eval() \n"; return 0; };
};

ObjF objf;

				template <typename FuncT> // <- What would go here
class		fwrap { public:
		fwrap		(FuncT func) : func_(func) 	{}
	int 	eval		() 				{ return func_(1, 2); }
		FuncT func_;
};

int main() {
	fwrap<int (*)(int, int)>	t1(f1);	t1.eval();
	fwrap<int (*)(double, double)>	t2(f2);	t2.eval();
	fwrap<typeof(&f2)>		t3(f2);	t3.eval();

	//fwrap<int (ObjF::operator())(double, double)>		t4(Func2);	t4.eval();

	//fwrap<int(ObjF::*)(double, double)>		t4(objf.*operator());
	// no err:  fwrap<int(ObjF::*)(double, double)>		t4(&ObjF::operator());
	fwrap<int(*)(double, double)>		t4(ObjF::eval); t4.eval();
	fwrap<typeof(&ObjF::eval)>		t5(ObjF::eval); t5.eval();

	return 0;

}
