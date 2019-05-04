## Snippet

### C++

```cpp
#define pb push_back
#define sz(s) ((int)s.size())
#define all(vec) vec.begin(), vec.end()

typedef long long ll;
typedef vector<ll> VL;
typedef vector<int> VI;
typedef pair<int, int> pii;

#ifdef local
#define debug(x...) do { cout << "[ "#x" ] -> "; err(x); } while (0)
template <class T>
inline void _E(T x) { cout << x; }
template <class L, class R>
inline void _E(pair<L, R> arg) {
  cout << "("; _E(arg.first), _E(','), _E(' '), _E(arg.second); cout << ")";
}
template <template <class...> class T, class t>
inline void _E(T<t> arr) {
  cout << "[ ";
  for (auto it = begin(arr), en = end(arr); it != en; it++) {
    if (it != begin(arr)) cout << ", "; _E(*it);
  }
  cout << " ]";
}
inline void _E(string s) { cout << "\"" + s + "\""; }

inline void err() { cout << std::endl; }
template <class T, class... U>
inline void err(T arg, U... args) {
  _E(arg); if (sizeof...(args)) cout << ", "; err(args...);
}
#else
#define debug(...) do {} while (0)
#endif
```

### Java

```java
import java.io.*;
import java.util.*;
import java.math.BigInteger;

public class Main {
  public static void main(String[] args) {
    InputReader in = new InputReader(System.in);
    PrintWriter out = new PrintWriter(System.out);
    Task solver = new Task();
    int taskNum = 1;
    // int taskNum = in.nextInt();
    solver.solve(taskNum, in, out);
    out.close();
  }

  public static class Task {
    void solve(int t, InputReader in, PrintWriter out) {

    }
  }

  static class InputReader {
    public BufferedReader reader;
    public StringTokenizer tokenizer;

    public InputReader(InputStream stream) {
      reader = new BufferedReader(new InputStreamReader(stream), 32768);
      tokenizer = null;
    }

    public String next() {
      while (tokenizer == null || !tokenizer.hasMoreTokens()) {
        try {
          tokenizer = new StringTokenizer(reader.readLine());
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
      return tokenizer.nextToken();
    }
    public int nextInt() {
      return Integer.parseInt(next());
    }
    public BigInteger nextBigInteger() {
      return new BigInteger(next());
    }
  }
}
```

### unordered_map

```cpp
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().\
                                          time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};

unordered_map<long long, int, custom_hash> safe_map;
```