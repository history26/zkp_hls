#include"arithmetic.hpp"
struct epoint {
    uint384 x=0;
    uint384 y=1;
    uint384 z=0;
    uint384 t=0;
};


epoint UniPadd(const epoint P, const epoint Q,uint384 k, uint384 q, uint384 INV);
