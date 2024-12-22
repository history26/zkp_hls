#include"arithmetic.hpp"
using namespace std;
/* 预计算,模数为381位
uint768 barrett_precompute(uint384 q) {
    // μ = ⌊2^762 / q⌋
    uint768 mu = (uint768(1) << 762) / q;
    return mu;
}*/

uint384 div(uint384 a, uint384 b)
{
	uint384 temp=a/b;
	return temp;
}
// 巴雷特约简
uint384 barrett_reduce(uint768 a, uint384 q, uint768 mu)
{
    uint768 t = karatsuba_768(a, mu) >> 762;//big_mul a*mu
	//uint768 t = (a*mu) >> 762;
    uint384 r = uint384(a - t * q);
    if (r >= q) {
        r = r - q;
    }
    return r;
}
ap_uint<1536> karatsuba_768(ap_uint<768> x, ap_uint<768> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<384> x_low = x(383, 0);
    ap_uint<384> x_high = x(767, 384);
    ap_uint<384> y_low = y(383, 0);
    ap_uint<384> y_high = y(767, 384);

    ap_uint<385> z1_left = x_low + x_high;
    ap_uint<385> z1_right = y_low + y_high;
    ap_uint<384> low1 = z1_left(383,0);
    ap_uint<384> low2 = z1_right(383,0);
    ap_uint<1536> z0 = karatsuba_384(x_low, y_low);
    ap_uint<1536> z2 = karatsuba_384(x_high, y_high);
    ap_uint<768> z1t = karatsuba_384(low1, low2);
    ap_uint<1536> temp=z1_left[384]*z1_right[384];
    ap_uint<1536> z1 = z1t+(ap_uint<768>(z1_left[384]*low2)<<384)+(ap_uint<768>(z1_right[384]*low1)<<384)+(temp<<768);
    ap_uint<1536> res = (z2 << 768) + ((z1-z0-z2) << 384) + z0;
    return res;
}
uint384 ADD(uint384 a, uint384 b, uint384 q)
{
	uint384 temp = a + b;//a<2^381 & b<2^381
	if(temp>=q)
		temp=temp-q;
	return temp;
}
uint384 SUB(uint384 a, uint384 b, uint384 q)
{
	uint384 temp;
	if(a>=b){
		temp = a - b;
	} else {
		temp = q - ( b - a );
	}
	return temp;
}
//uint384 MUL(uint384 a, uint384 b, uint384 q, uint768 u)
//{
//	uint768 temp=karatsuba_384(a,b);
// 	uint768 temp=a*b;
//	uint384 ans=barrett_reduce(temp, q, u);
//	return ans;
//}
uint384 montgomery_reduce(uint768 t, const uint384& m, const uint384& inv)
{
	uint384 t_low = t(383, 0);
	uint768 temp=karatsuba_384(t_low, inv);
//	uint768 temp=t_low*inv;
	uint384 temp_low=temp(383,0);
	uint768 x;
//	x=t+karatsuba_384(temp_low,m); //overflow warning!
	x=t+temp_low*m;
	uint384 res=x(767,384);
	if(res>=m)
	{
		res=res-m;
	}
	return res;
}
uint384 MUL(uint384 a, uint384 b, uint384 q, uint384 INV)
{
  uint768 temp=karatsuba_384(a,b);
//  uint768 temp=a*b;
  uint384 res=montgomery_reduce(temp,q,INV);
  return res;
}
//uint384 MUL(uint384 a, uint384 b, uint384 q, uint384 INV)
//{
//	uint768 temp=karatsuba_384(a,b);
////  uint768 temp=a*b;
//  uint384 res=temp % q;
//  return res;
//}

ap_uint<768> karatsuba_384(ap_uint<384> x, ap_uint<384> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<192> x_low = x(191, 0);
    ap_uint<192> x_high = x(383, 192);
    ap_uint<192> y_low = y(191, 0);
    ap_uint<192> y_high = y(383, 192);

    ap_uint<193> z1_left = x_low + x_high;
    ap_uint<193> z1_right = y_low + y_high;
    ap_uint<192> low1 = z1_left(191,0);
    ap_uint<192> low2 = z1_right(191,0);
    ap_uint<768> z0 = karatsuba_192(x_low, y_low);
    ap_uint<768> z2 = karatsuba_192(x_high, y_high);
    ap_uint<384> z1t = karatsuba_192(low1, low2);
//    ap_uint<768> z0 = x_low*y_low;
//    ap_uint<768> z2 = x_high*y_high;
//    ap_uint<384> z1t = low1*low2;
    ap_uint<768> temp=z1_left[192]*z1_right[192];
    ap_uint<768> z1 = z1t+(ap_uint<384>(z1_left[192]*low2)<<192)+(ap_uint<384>(z1_right[192]*low1)<<192)+(temp<<384);
    ap_uint<768> res = (z2 << 384) + ((z1-z0-z2) << 192) + z0;
    return res;
}

ap_uint<384> karatsuba_192(ap_uint<192> x, ap_uint<192> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<96> x_low = x(95, 0);
    ap_uint<96> x_high = x(191, 96);
    ap_uint<96> y_low = y(95, 0);
    ap_uint<96> y_high = y(191, 96);

    ap_uint<97> z1_left = x_low + x_high;
    ap_uint<97> z1_right = y_low + y_high;
    ap_uint<96> low1 = z1_left(95,0);
    ap_uint<96> low2 = z1_right(95,0);
    ap_uint<384> z0 = karatsuba_96(x_low, y_low);
    ap_uint<384> z2 = karatsuba_96(x_high, y_high);
    ap_uint<192> z1t = karatsuba_96(low1, low2);
//    ap_uint<384> z0 = x_low*y_low;
//    ap_uint<384> z2 = x_high*y_high;
//    ap_uint<192> z1t = low1*low2;
    ap_uint<384> temp=z1_left[96]*z1_right[96];
    ap_uint<384> z1 = z1t+(ap_uint<192>(z1_left[96]*low2)<<96)+(ap_uint<192>(z1_right[96]*low1)<<96)+(temp<<192);
    ap_uint<384> res = (z2 << 192) + ((z1-z0-z2) << 96) + z0;
    return res;
}

ap_uint<192> karatsuba_96(ap_uint<96> x, ap_uint<96> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<48> x_low = x(47, 0);
    ap_uint<48> x_high = x(95, 48);
    ap_uint<48> y_low = y(47, 0);
    ap_uint<48> y_high = y(95, 48);

    ap_uint<49> z1_left = x_low + x_high;
    ap_uint<49> z1_right = y_low + y_high;
    ap_uint<48> low1 = z1_left(47,0);
    ap_uint<48> low2 = z1_right(47,0);
    ap_uint<192> z0 = karatsuba_48(x_low, y_low);
    ap_uint<192> z2 = karatsuba_48(x_high, y_high);
    ap_uint<96> z1t = karatsuba_48(low1, low2);
//    ap_uint<192> z0 = x_low*y_low;
//    ap_uint<192> z2 = x_high*y_high;
//    ap_uint<96> z1t = low1*low2;
    ap_uint<192> temp=z1_left[48]*z1_right[48];
    ap_uint<192> z1 = z1t+(ap_uint<96>(z1_left[48]*low2)<<48)+(ap_uint<96>(z1_right[48]*low1)<<48)+(temp<<96);
    ap_uint<192> res = (z2 << 96) + ((z1-z0-z2) << 48) + z0;
    return res;
}
ap_uint<96> karatsuba_48(ap_uint<48> x, ap_uint<48> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<24> x_low = x(23, 0);
    ap_uint<24> x_high = x(47, 24);
    ap_uint<24> y_low = y(23, 0);
    ap_uint<24> y_high = y(47, 24);

    ap_uint<25> z1_left = x_low + x_high;
    ap_uint<25> z1_right = y_low + y_high;
    ap_uint<24> low1 = z1_left(23,0);
    ap_uint<24> low2 = z1_right(23,0);
//    ap_uint<96> z0 = karatsuba_24(x_low, y_low);
//    ap_uint<96> z2 = karatsuba_24(x_high, y_high);
//    ap_uint<48> z1t = karatsuba_24(low1, low2);
    ap_uint<96> z0 = x_low*y_low;
    ap_uint<96> z2 = x_high*y_high;
    ap_uint<48> z1t = low1*low2;
    ap_uint<96> temp=z1_left[24]*z1_right[24];
    ap_uint<96> z1 = z1t+(ap_uint<48>(z1_left[24]*low2)<<24)+(ap_uint<48>(z1_right[24]*low1)<<24)+(temp<<48);
    ap_uint<96> res = (z2 << 48) + ((z1-z0-z2) << 24) + z0;
    return res;
}
ap_uint<48> karatsuba_24(ap_uint<24> x, ap_uint<24> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<12> x_low = x(11, 0);
    ap_uint<12> x_high = x(23, 12);
    ap_uint<12> y_low = y(11, 0);
    ap_uint<12> y_high = y(23, 12);

    ap_uint<48> z0 = x_low*y_low;
    ap_uint<48> z2 = x_high*y_high;
    ap_uint<48> z1 = (x_low + x_high)*(y_low + y_high);
    ap_uint<48> res = (z2 << 24) + ((z1-z0-z2) << 12) + z0;
    return res;
}
