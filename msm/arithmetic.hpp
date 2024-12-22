#include<ap_int.h>
#include<iostream>
typedef ap_uint<256> uint256;
typedef ap_uint<384> uint384;
typedef ap_uint<768> uint768;
const uint64_t INV_r= 0xc2e1f593efffffff;
// Bls381 G1 curve equation: y^2 = x^3 + 4
// Bls377 G1 curve equation: y^2 = x^3 + 1


const uint384 bls377_r("8444461749428370424248824938781546531375899335154063827935233455917409239041");
const uint768 bls377_u("93782108216853554240664012119956165350543297590574812882380072409752409455923467528765225298139683367851683600137419");//for barrett_precompute

uint384 div(uint384 a, uint384 b);
uint384 barrett_reduce(uint768 a, uint384 q, uint768 mu);
uint384 montgomery_reduce(uint768 t, const uint384& m, const uint384& inv);
uint384 ADD(uint384 a, uint384 b, uint384 q);
uint384 SUB(uint384 a, uint384 b, uint384 q);
//uint384 MUL(uint384 a, uint384 b, uint384 q, uint768 u);
uint384 MUL(uint384 a, uint384 b, uint384 q, uint384 INV);
ap_uint<1536> karatsuba_768(ap_uint<768> x, ap_uint<768> y);
ap_uint<768> karatsuba_384(ap_uint<384> x, ap_uint<384> y);
ap_uint<384> karatsuba_192(ap_uint<192> x, ap_uint<192> y);
ap_uint<192> karatsuba_96(ap_uint<96> x, ap_uint<96> y);
ap_uint<96> karatsuba_48(ap_uint<48> x, ap_uint<48> y);
ap_uint<48> karatsuba_24(ap_uint<24> x, ap_uint<24> y);
