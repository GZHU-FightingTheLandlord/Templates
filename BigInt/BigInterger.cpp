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

struct BigInterger {

// --------------------------- Constructor ----------------------------------
    BigInterger(ll num = 0);
    BigInterger(char *s, int len);
// --------------------------------------------------------------------------

    int n;          // BitNum
    vector<int> Vu; // Val

// ---------------------------- Add, Sub ------------------------------------
    friend BigInterger operator + (const BigInterger& a, const BigInterger& b);
    friend BigInterger& operator += (BigInterger& a, const BigInterger& b);
    friend BigInterger operator - (const BigInterger& a, const BigInterger& b);
    friend BigInterger& operator -= (BigInterger& a, const BigInterger& b);
// --------------------------------------------------------------------------

// ------------------------------ FFT ---------------------------------------
    const double PI = acos(-1.0);
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
    void rev(vector<Complex>& a, int n)
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
    void DFT(vector<Complex>& a, int n, int t)
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
    int trans(int x)
    {
        int i = 0;
        for (; x > (1 << i); i++);
        return 1 << i;
    }
// --------------------------------------------------------------------------

// ---------------------------- Multiply ------------------------------------
    friend BigInterger operator * (const BigInterger& a, const BigInterger& b);
    friend BigInterger& operator *= (BigInterger& a, const BigInterger& b);
// --------------------------------------------------------------------------
};