#include <cmath>

struct Complex {
	double r, i;
	
#define PI 3.141592653589793

	Complex(double x = 0, double y = 0) : r(x), i(y) {}
    Complex(int n) : r(cos(2 * PI / n)), i(sin(2 * PI / n)) {}

	Complex operator+ (const Complex& b)const {
		return Complex(r + b.r, i + b.i);
	}
	Complex operator- (const Complex& b)const {
		return Complex(r - b.r, i - b.i);
	}
	Complex operator* (const Complex& b)const {
		return Complex(r * b.r - i * b.i, r * b.i + i * b.r);
	}
	Complex operator/ (const Complex& b)const {
		double x = (r * b.r + i * b.i) / (b.r * b.r + b.i * b.i);
		double y = (i * b.r - r * b.i) / (b.r * b.r + b.i * b.i);
		return Complex(x, y);
	}
    friend Complex& operator+= (Complex& a, const Complex& b) {
        return (a = a + b);
    }
    friend Complex& operator-= (Complex& a, const Complex& b) {
        return (a = a - b);
    }
    friend Complex& operator*= (Complex& a, const Complex& b) {
        return (a = a * b);
    }
    friend Complex& operator/= (Complex& a, const Complex& b) {
        return (a = a / b);
    }
};