#include <stdio.h>
#include <string.h>
namespace io {
    const size_t buflen = (1 << 21) + 1;
    char buf[buflen], *st, *ed;
    inline char gc() { return ((st == ed) ? (st = buf, ed = st + fread(buf, 1, buflen, stdin), (st == ed) ? EOF : *st++) : *st++); }
    inline bool blank(char c) { return c == ' ' || c == '-' || c == '\n' || c == '\t' || c == '\r'; }
    template <typename T> inline bool Re(T& x) {
        register char c; int f = 1;
        while (blank(c = gc())) if (c == '-') f = -1;
        if (c == EOF) return false;
        for (x = c - '0'; (c = gc()) >= '0' && c <= '9'; x = x * 10 + c - '0'); x *= f;
        return true;
    }
    char obuf[buflen], *_ptr = obuf, *_oed = obuf + buflen - 1, _stk[55], _top;
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
}
using io::Re;
using io::Prln;