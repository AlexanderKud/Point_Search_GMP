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
    //Point Infinity;
    Elliptic_Curve EC;
    Secp256k1();
    Point ParsePublicKeyHex(std::string pubkey);
    Point DoublePoint(Point &p);
    Point AddPoints(Point &a, Point &b);
    Point ScalarMultiplication(mpz_class pk);
    Point PointMultiplication(Point &a, mpz_class pk);
    Point NegatePoint(Point &p);
    Point SubtractPoints(Point &a, Point &b);
    Point PointDivision(Point &a, mpz_class sc);
    std::string GetPublicKeyHex(Point &p);
    void Context_Init();
};

#endif
