#include"curve.hpp"

epoint UniPadd(const epoint P, const epoint Q, uint384 k, uint384 q, uint384 INV) {
#pragma HLS PIPELINE II=3
//k=d'/2
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	// Step 1: Calculate (Y1 - X1) and (Y2 - X2)
    uint384 Y1_minus_X1 = SUB(P.y, P.x, q);
    uint384 Y2_minus_X2 = SUB(Q.y, Q.x, q);
    // Step 2: Calculate A = (Y1 - X1) * (Y2 - X2)
    uint384 A = MUL(Y1_minus_X1, Y2_minus_X2, q, INV);
    // Step 3: Calculate (Y1 + X1) and (Y2 + X2)
    uint384 Y1_plus_X1 = ADD(P.y, P.x, q);
    uint384 Y2_plus_X2 = ADD(Q.y, Q.x, q);
    // Step 4: Calculate B = (Y1 + X1) * (Y2 + X2)
    uint384 B = MUL(Y1_plus_X1, Y2_plus_X2, q, INV);
    // Step 5: Calculate C = k * T1 * T2
    uint384 T1_times_T2 = MUL(P.t, Q.t, q, INV);
    uint384 C = MUL(k, T1_times_T2, q, INV);
    // Step 6: Calculate D = 2 * Z1 * Z2
    uint384 Z1_times_Z2 = MUL(P.z, Q.z, q, INV);
    uint384 D = ADD(Z1_times_Z2, Z1_times_Z2, q);  // D = 2 * Z1 * Z2
    // Step 7: Calculate E = B - A
    uint384 E = SUB(B, A, q);
    // Step 8: Calculate F = D - C
    uint384 F = SUB(D, C, q);
    // Step 9: Calculate G = D + C
    uint384 G = ADD(D, C, q);
    // Step 10: Calculate H = B + A
    uint384 H = ADD(B, A, q);
    // Step 11: Calculate the result coordinates
    epoint R;
    R.x = MUL(E, F, q, INV);            // X3 = E * F
    R.y = MUL(G, H, q, INV);            // Y3 = G * H
    R.t = MUL(E, H, q, INV);            // T3 = E * H
    R.z = MUL(F, G, q, INV);            // Z3 = F * G
    return R;
}





