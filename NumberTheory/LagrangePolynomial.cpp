// Author: ConanYu
// type: double
class poly {
  private:
    vector<double>v;
  public:
    template <typename ...T>
    poly(const T &... tail) {
        v.clear();
        init(tail...);
    }

    void print() const {
        const int N = v.size();
        for(int i = N - 1; i >= 0; i--) {
            cout << v[i] << "x^" << i << (i == 0 ? "" : " + ");
        }
    }

    friend poly operator + (const poly & x, const poly & y) {
        const bool f = x.v.size() > y.v.size();
        const poly *a = (f ? &x : &y), *b = (f ? &y : &x);
        poly z(0);
        z.v.clear();
        const int n1 = b->v.size(), n2 = a->v.size();
        for(int i = 0; i < n1; i++) {
            z.v.push_back(a->v[i] + b->v[i]);
        }
        for(int i = n1; i < n2; i++) {
            z.v.push_back(a->v[i]);
        }
        return z;
    }

    friend poly operator * (const poly & x, const poly & y) {
        poly z(0);
        z.v = vector<double>(x.v.size() + y.v.size() - 1, 0);
        for(int i = 0; i < int(x.v.size()); i++) {
            for(int j = 0; j < int(y.v.size()); j++) {
                z.v[i + j] += x.v[i] * y.v[j];
            }
        }
        return z;
    }

    friend ostream & operator << (ostream & out, const poly & x) {
        x.print();
        return out;
    }
  private:
    void init() {}
    template <typename ...T>
    void init(const double & head, const T &... tail) {
        init(tail...);
        v.push_back(head);
    }

};

typedef pair<double, double> pdd;

poly lagrange(const vector<pdd> & x) {
    poly ans(0);
    for(int i = 0; i < int(x.size()); i++) {
        poly cur(x[i].second);
        for(int j = 0; j < int(x.size()); j++) {
            if(i == j)continue;
            cur = cur * poly(1.0 / (x[i].first - x[j].first), x[j].first / (x[j].first - x[i].first));
        }
        ans = ans + cur;
    }
    return ans;
}

void solve() {
    double x, y;
    vector<pdd>v;
    int n;
    cin >> n;
    for(int i = 0; i < n; i++) {
        cin >> x >> y;
        v.push_back(make_pair(x, y));
    }
    cout << lagrange(v) << endl;
}
