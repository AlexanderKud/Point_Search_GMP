#include <iostream>
#include <gmpxx.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "secp256k1.h"
#include "../util/util.h"

Point::Point() {}

Point::Point(mpz_class cx, mpz_class cy) {
  x = cx;
  y = cy;
}

Point::Point(const Point &p) {
  x = p.x;
  y = p.y;
}

Elliptic_Curve::Elliptic_Curve() {}

Secp256k1::Secp256k1() {}

void Secp256k1::Init() {
    EC.a = 0; EC.b = 7;
    EC.p = "115792089237316195423570985008687907853269984665640564039457584007908834671663";
    EC.n = "115792089237316195423570985008687907852837564279074904382605163141518161494337";
    G.x = "55066263022277343669578718895168534326250603453777594175500187360389116729240";
    G.y = "32670510020758816978083085130507043184471273380659243275938904335757337482424";
    return;
}

mpz_class Secp256k1::ModNeg(const mpz_class &y) {
    mpz_class ret;
    mpz_sub(ret.get_mpz_t(), EC.p.get_mpz_t(), y.get_mpz_t());
    return ret;
}

std::string Secp256k1::GetPublicKeyHex(Point &a) {
    char cpub[68];
    if(mpz_tstbit(a.y.get_mpz_t(), 0) == 0) {
        gmp_snprintf(cpub, 67, "02%0.64Zx", a.x.get_mpz_t());
    } else {
        gmp_snprintf(cpub, 67, "03%0.64Zx", a.x.get_mpz_t());
    }
    return std::string(cpub);
}

std::string Secp256k1::GetXHex(const mpz_class &x) {
    return x.get_str(16);
}

Point Secp256k1::ParsePublicKeyHex(std::string pubkey) {
    Point A;
    char X_coord[65];
    if(startsWith("02", pubkey.c_str())) {
        mpz_class ysquared, mpz_aux1, mpz_aux2, Y1, Y2;
        substr(X_coord, const_cast<char*>(pubkey.c_str()), 3, 64);
        A.x = mpz_class(X_coord, 16);
        mpz_pow_ui(mpz_aux1.get_mpz_t(), A.x.get_mpz_t(), 3);
        mpz_add(mpz_aux1.get_mpz_t(), mpz_aux1.get_mpz_t(), EC.b.get_mpz_t());
        mpz_mod(ysquared.get_mpz_t(), mpz_aux1.get_mpz_t(), EC.p.get_mpz_t());
        mpz_add_ui(mpz_aux2.get_mpz_t(), EC.p.get_mpz_t(), 1);
        mpz_fdiv_q_ui(mpz_aux2.get_mpz_t(), mpz_aux2.get_mpz_t(), 4);
        mpz_powm(Y1.get_mpz_t(), ysquared.get_mpz_t(), mpz_aux2.get_mpz_t(), EC.p.get_mpz_t());
        mpz_sub(Y2.get_mpz_t(), EC.p.get_mpz_t(), Y1.get_mpz_t());
        if(mpz_tstbit(Y1.get_mpz_t(), 0) == 0) {
            A.y = Y1;
        } else {
            A.y = Y2;
        }
    }
    else {
        mpz_class ysquared, mpz_aux1, mpz_aux2, Y1, Y2;
        substr(X_coord, const_cast<char*>(pubkey.c_str()), 3, 64);
        A.x = mpz_class(X_coord, 16);
        mpz_pow_ui(mpz_aux1.get_mpz_t(), A.x.get_mpz_t(), 3);
        mpz_add(mpz_aux1.get_mpz_t(), mpz_aux1.get_mpz_t(), EC.b.get_mpz_t());
        mpz_mod(ysquared.get_mpz_t(), mpz_aux1.get_mpz_t(), EC.p.get_mpz_t());
        mpz_add_ui(mpz_aux2.get_mpz_t(), EC.p.get_mpz_t(), 1);
        mpz_fdiv_q_ui(mpz_aux2.get_mpz_t(), mpz_aux2.get_mpz_t(), 4);
        mpz_powm(Y1.get_mpz_t(), ysquared.get_mpz_t(), mpz_aux2.get_mpz_t(), EC.p.get_mpz_t());
        mpz_sub(Y2.get_mpz_t(), EC.p.get_mpz_t(), Y1.get_mpz_t());
        if(mpz_tstbit(Y1.get_mpz_t(), 0) == 0) {
            A.y = Y2;
        } else {
            A.y = Y1;
        }
    }
    return A;
}

Point Secp256k1::DoublePoint(Point &p) {
    Point A;
    mpz_class temp, slope;

    mpz_mul_ui(temp.get_mpz_t(), p.y.get_mpz_t(), 2);
    mpz_invert(temp.get_mpz_t(), temp.get_mpz_t(), EC.p.get_mpz_t());
    mpz_mul(slope.get_mpz_t(), p.x.get_mpz_t(), p.x.get_mpz_t());
    mpz_mul_ui(slope.get_mpz_t(), slope.get_mpz_t(), 3);
    mpz_add(slope.get_mpz_t(), slope.get_mpz_t(), EC.a.get_mpz_t());
    mpz_mul(slope.get_mpz_t(), slope.get_mpz_t(), temp.get_mpz_t());
    mpz_mod(slope.get_mpz_t(), slope.get_mpz_t(), EC.p.get_mpz_t());

    mpz_mul(A.x.get_mpz_t(), slope.get_mpz_t(), slope.get_mpz_t());
    mpz_sub(A.x.get_mpz_t(), A.x.get_mpz_t(), p.x.get_mpz_t());
    mpz_sub(A.x.get_mpz_t(), A.x.get_mpz_t(), p.x.get_mpz_t());
    mpz_mod(A.x.get_mpz_t(), A.x.get_mpz_t(), EC.p.get_mpz_t());

    mpz_sub(temp.get_mpz_t(), p.x.get_mpz_t(), A.x.get_mpz_t());
    mpz_mul(A.y.get_mpz_t(), slope.get_mpz_t(), temp.get_mpz_t());
    mpz_sub(A.y.get_mpz_t(), A.y.get_mpz_t(), p.y.get_mpz_t());
    mpz_mod(A.y.get_mpz_t(), A.y.get_mpz_t(), EC.p.get_mpz_t());
    
    return A;
}

Point Secp256k1::AddPoints(Point &a, Point &b) {
    Point A;
    mpz_class u, v, slope;

    mpz_sub(u.get_mpz_t(), a.y.get_mpz_t(), b.y.get_mpz_t());
    mpz_sub(v.get_mpz_t(), a.x.get_mpz_t(), b.x.get_mpz_t());
    mpz_invert(v.get_mpz_t(), v.get_mpz_t(), EC.p.get_mpz_t());
    mpz_mul(slope.get_mpz_t(), u.get_mpz_t(), v.get_mpz_t());
    mpz_mod(slope.get_mpz_t(), slope.get_mpz_t(), EC.p.get_mpz_t());

    mpz_mul(A.x.get_mpz_t(), slope.get_mpz_t(), slope.get_mpz_t());
    mpz_sub(A.x.get_mpz_t(), A.x.get_mpz_t(), a.x.get_mpz_t());
    mpz_sub(A.x.get_mpz_t(), A.x.get_mpz_t(), b.x.get_mpz_t());
    mpz_mod(A.x.get_mpz_t(), A.x.get_mpz_t(), EC.p.get_mpz_t());

    mpz_sub(u.get_mpz_t(), a.x.get_mpz_t(), A.x.get_mpz_t());
    mpz_mul(A.y.get_mpz_t(), slope.get_mpz_t(), u.get_mpz_t());
    mpz_sub(A.y.get_mpz_t(), A.y.get_mpz_t(), a.y.get_mpz_t());
    mpz_mod(A.y.get_mpz_t(), A.y.get_mpz_t(), EC.p.get_mpz_t());
    
    return A;
}


Point Secp256k1::SubtractPoints(Point &a, Point &b) {
    Point A, ATemp;
    mpz_set(ATemp.x.get_mpz_t(), b.x.get_mpz_t());
    mpz_set(ATemp.y.get_mpz_t(), b.y.get_mpz_t());
    mpz_sub(ATemp.y.get_mpz_t(), EC.p.get_mpz_t(), ATemp.y.get_mpz_t());
    A = AddPoints(a, ATemp);
    return A;
}

Point Secp256k1::ScalarMultiplication(mpz_class pk) {
    Point A;
    int no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(pk.get_mpz_t(), 2);
    mpz_set(A.x.get_mpz_t(), G.x.get_mpz_t());
    mpz_set(A.y.get_mpz_t(), G.y.get_mpz_t());
    for(loop = no_of_bits - 2; loop >= 0 ; loop--) {
        A = DoublePoint(A);
        if(mpz_tstbit(pk.get_mpz_t(), loop)) {
            A = AddPoints(A, G);
        }
    }
    return A;
}

Point Secp256k1::PointMultiplication(Point &a, mpz_class pk) {
    Point A, ATemp;
    mpz_set(ATemp.x.get_mpz_t(), a.x.get_mpz_t());
    mpz_set(ATemp.y.get_mpz_t(), a.y.get_mpz_t());
    int no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(pk.get_mpz_t(), 2);
    mpz_set(A.x.get_mpz_t(), a.x.get_mpz_t());
    mpz_set(A.y.get_mpz_t(), a.y.get_mpz_t());
    for(loop = no_of_bits - 2; loop >= 0 ; loop--) {
        A = DoublePoint(A);
        if(mpz_tstbit(pk.get_mpz_t(), loop)) {
            A = AddPoints(A, ATemp);
        }
    }
    return A;
}

Point Secp256k1::PointDivision(Point &a, mpz_class sc) {
    Point A;
    mpz_class div_s;
    mpz_invert(div_s.get_mpz_t(), sc.get_mpz_t(), EC.n.get_mpz_t());
    A = PointMultiplication(a, div_s);
    return A;
}

IntGroup::IntGroup(int size) {
    this->size = size;
    subp = (mpz_class*)malloc(size * sizeof(mpz_class));
}

IntGroup::~IntGroup() {
    free(subp);
}

void IntGroup::ModInv(mpz_class *ints) { // GMP library takes mpz_class array as input directly
    mpz_class newValue, inverse;         // and modifies it directly because no const specifier in front of it.
    subp[0] = ints[0];
    for (int i = 1; i < size; i++){
        mpz_mul(subp[i].get_mpz_t(), subp[i - 1].get_mpz_t(), ints[i].get_mpz_t());
        mpz_mod(subp[i].get_mpz_t(), subp[i].get_mpz_t(), P.get_mpz_t());
    }
    inverse = subp[size - 1];
    mpz_invert(inverse.get_mpz_t(), inverse.get_mpz_t(), P.get_mpz_t());
    for(int i = size - 1; i > 0; i--) {
        mpz_mul(newValue.get_mpz_t(), subp[i - 1].get_mpz_t(), inverse.get_mpz_t());
        mpz_mod(newValue.get_mpz_t(), newValue.get_mpz_t(), P.get_mpz_t());
        mpz_mul(inverse.get_mpz_t(), inverse.get_mpz_t(), ints[i].get_mpz_t());
        mpz_mod(inverse.get_mpz_t(), inverse.get_mpz_t(), P.get_mpz_t());
        ints[i] = newValue;
    }
    ints[0] = inverse;
}
