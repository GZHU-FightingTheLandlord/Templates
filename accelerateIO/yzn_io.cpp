namespace io {
namespace detail {
  const size_t buflen = (1 << 21) + 1;
  char buf[buflen], *st = nullptr, *ed = nullptr;
  inline char gc() { return ((st == ed) ? (st = buf, ed = st + fread(buf, 1, buflen, stdin), (st == ed) ? EOF : *st++) : *st++); }
  inline bool blank(char c) { return c == ' ' || c == '-' || c == '\n' || c == '\t' || c == '\r'; }
  template <typename T> inline bool Re(T &x) {
    char c; int f = 1;
    while (blank(c = gc())) if (c == '-') f = -1;
    if (c == EOF) return false;
    for (x = c - '0'; (c = gc()) >= '0' && c <= '9'; x = x * 10 + c - '0');
    x *= f; return true;
  }
  char obuf[buflen], *_ptr = obuf, *_oed = obuf + buflen - 1, _stk[55];
  int _top;
  void flush() { fwrite(obuf, 1, _ptr - obuf, stdout); _ptr = obuf;}
  void pc(char c) { *(_ptr++) = c; if (_ptr == _oed) flush(); }
  template <typename T> inline void Pr(T x) {
    if (x == 0) { pc('0'); return; }
    if (x < 0) { pc('-'); x = -x; }
    for (_top = 0; x > 0; _top++, x /= 10) _stk[_top] = (x % 10) + '0';
    while (_top > 0) pc(_stk[--_top]);
  }
  template <typename T> inline void Prln(T x) { Pr(x); pc('\n'); }
  struct _IOflusher_ { ~_IOflusher_() { flush(); } } __flusher__;
}  // namespace detail
struct instream {
  instream &operator>>(char &__n) { __n = detail::gc(); return *this; }
  instream &operator>>(unsigned char &__n) { __n = detail::gc(); return *this; }
  instream &operator>>(short &__n) { detail::Re(__n); return *this; }
  instream &operator>>(unsigned short &__n) { detail::Re(__n); return *this; }
  instream &operator>>(int &__n) { detail::Re(__n); return *this; }
  instream &operator>>(unsigned int &__n) { detail::Re(__n); return *this; }
  instream &operator>>(long &__n) { detail::Re(__n); return *this; }
  instream &operator>>(unsigned long &__n) { detail::Re(__n); return *this; }
  instream &operator>>(long long &__n) { detail::Re(__n); return *this; }
  instream &operator>>(unsigned long long &__n) { detail::Re(__n); return *this; }
#ifdef __SIZEOF_INT128__
  instream &operator>>(__int128 &__n) { detail::Re(__n); return *this;}
  instream &operator>>(__uint128_t &__n) { detail::Re(__n); return *this; }
#endif
  instream &operator>>(std::string &__n) {
    __n.clear();
    char c = detail::gc();
    if(c == EOF) return *this;
    while(c != EOF && detail::blank(c)) c = detail::gc();
    if(c == EOF) return *this;
    while(c != EOF && !detail::blank(c)) {
      __n.push_back(c);
      c = detail::gc();
    }
    if(c != EOF) detail::st--;
    return *this;
  }
} in;
struct outstream {
  outstream &operator<<(const char &__n) { detail::pc(__n); return *this; }
  outstream &operator<<(const unsigned char &__n) { detail::pc(__n); return *this; }
  outstream &operator<<(const short &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const unsigned short &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const int &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const unsigned int &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const long &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const unsigned long &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const long long &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const unsigned long long &__n) { detail::Pr(__n); return *this; }
#ifdef __SIZEOF_INT128__
  outstream &operator<<(const __int128 &__n) { detail::Pr(__n); return *this; }
  outstream &operator<<(const __uint128_t &__n) { detail::Pr(__n); return *this; }
#endif
  outstream &operator<<(const std::string &__n) {
    for(const char &__c: __n) detail::pc(__c);
    return *this;
  }
} out;
}  // namespace io
using io::in;
using io::out;
