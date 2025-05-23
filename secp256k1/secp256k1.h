#ifndef SECP256K1_H
#define SECP256K1_H

class Point {
public:
    mpz_class x;
    mpz_class y;
    Point();
    Point(mpz_class x, mpz_class y);
    Point(const Point &p);
};

class Elliptic_Curve {
public:
    mpz_class a; 
    mpz_class b; 
    mpz_class p;
    mpz_class n;
    Elliptic_Curve(); 
};

class Secp256k1 {
public:
    Point G;
    Elliptic_Curve EC;
    Secp256k1();
    mpz_class ModNeg(const mpz_class &y);
    std::string GetPublicKeyHex(Point &p);
    Point ParsePublicKeyHex(std::string pubkey);
    Point DoublePoint(Point &p);
    Point AddPoints(Point &a, Point &b);
    Point SubtractPoints(Point &a, Point &b);
    Point ScalarMultiplication(mpz_class pk);
    Point PointMultiplication(Point &a, mpz_class pk);
    Point PointDivision(Point &a, mpz_class sc);
    void Init();
};

class IntGroup { // Batch Inversion Class and logic taken from JLP IntGroup.h/.cpp 
private:         // and with minor modifications adapted for GMP library.
    int size;
    mpz_class *subp;
    mpz_class P = mpz_class("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10);
public:
    IntGroup(int size);
    ~IntGroup();
    void ModInv(mpz_class *ints);
};

#endif
