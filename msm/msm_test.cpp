#include"msm_pip.hpp"
//#include"curve.hpp"
using namespace std;
int main()
{
	epoint a;
	epoint b;
	epoint res=UniPadd(a,b,bls377_k,bls377_q,bls377_i);
    return 0;
}
