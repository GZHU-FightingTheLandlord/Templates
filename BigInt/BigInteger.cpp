#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <vector>
using namespace std;

#define pb      push_back
#define mp      make_pair
#define endl    '\n'
#define szz(a)  (int)a.size()
#define all(a)  a.begin(),a.end()

typedef long long ll;

struct BigInteger {
	// --------------------------- Constructor ----------------------------------
	BigInteger() {}
	BigInteger(char *s, int len) : n(len), Vu(len + 5) {
		for (int i = 0; i < len; i++) Vu[i] = s[len - i - 1] - '0';
	}
	// --------------------------------------------------------------------------
	int n;          // BitNum
	vector<int> Vu; // Val
	// ---------------------------- Add, Sub ------------------------------------
	friend BigInteger operator + (const BigInteger& a, const BigInteger& b) {
		// 
	}
	friend BigInteger& operator += (BigInteger& a, const BigInteger& b) {
		return a = a + b;
	}

	// friend BigInterger operator - (const BigInterger& a, const BigInterger& b);
	// friend BigInterger& operator -= (BigInterger& a, const BigInterger& b);
	// --------------------------------------------------------------------------
	// ------------------------------ FFT ---------------------------------------
#define PI 3.141592653589
	// Complex
	struct Complex {
		double r, i;
		Complex(double x = 0, double y = 0) : r(x), i(y) {}
		Complex(int n) : r(cos(2 * PI / n)), i(sin(2 * PI / n)) {}
		Complex operator + (const Complex& b)const {
			return Complex(r + b.r, i + b.i);
		}
		Complex operator - (const Complex& b)const {
			return Complex(r - b.r, i - b.i);
		}
		Complex operator * (const Complex& b)const {
			return Complex(r * b.r - i * b.i, r * b.i + i * b.r);
		}
		friend Complex& operator *= (Complex& a, const Complex& b) {
			return a = a * b;
		}
	};
	// bit reverse
	static void rev(vector<Complex>& a, int n)
	{
		for (int i = 1, j = n >> 1, k; i < n - 1; i++) {
			if (i < j) swap(a[i], a[j]);
			for (k = n >> 1; j >= k; j -= k, k >>= 1);
			j += k;
		}
	}
	// Discrete Fourier transform
	// t ->  1, DFT
	// t -> -1, IDFT
	static void DFT(vector<Complex>& a, int n, int t)
	{
		rev(a, n);
		for (int i = 2; i <= n; i <<= 1) {
			Complex wi(i * t);
			for (int j = 0; j < n; j += i) {
				Complex w(1, 0);
				for (int k = j, h = i >> 1; k < j + h; k++) {
					Complex t = w * a[k + h], u = a[k];
					a[k] = u + t;
					a[k + h] = u - t;
					w *= wi;
				}
			}
		}
		if (t == -1) for (int i = 0; i < n; i++) a[i].r /= n;
	}
	// Get FFT_Len
	// min(2 ^ p) which (2 ^ p) > x
	static int trans(int x)
	{
		int i = 0;
		for (; x >(1 << i); i++);
		return 1 << i;
	}
	// --------------------------------------------------------------------------
	// ---------------------------- Multiply ------------------------------------
	friend BigInteger operator * (const BigInteger& a, const BigInteger& b) {
		int len = trans(a.n + b.n - 1);

		vector<Complex> c1(len + 3), c2(len + 3);
		for (int i = 0; i < a.n; i++) c1[i] = Complex(a.Vu[i], 0);
		for (int i = a.n; i < len; i++) c1[i] = Complex();
		for (int i = 0; i < b.n; i++) c2[i] = Complex(b.Vu[i], 0);
		for (int i = b.n; i < len; i++) c2[i] = Complex();

		DFT(c1, len, 1); DFT(c2, len, 1);
		for (int i = 0; i < len; i++) c1[i] *= c2[i];
		DFT(c1, len, -1);

		BigInteger ret; ret.n = len; ret.Vu = vector<int>(len + 2, 0);
		for (int i = 0; i < len; i++) ret.Vu[i] = (int)(c1[i].r + 0.5); ret.Vu[len] = 0;
		for (int i = 0; i < len; i++) {
			ret.Vu[i + 1] += ret.Vu[i] / 10;
			ret.Vu[i] %= 10;
		}
		while (ret.Vu[ret.n] == 0 && ret.n > 0) {
			ret.n--;
		}
		ret.n++;
		return ret;
	}
	friend BigInteger& operator *= (BigInteger& a, const BigInteger& b) {
		return a = a * b;
	}
	// --------------------------------------------------------------------------
};