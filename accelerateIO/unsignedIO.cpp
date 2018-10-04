// Need C++11 support(maybe)
#include <stdio.h>
using namespace std;

namespace io {
    const int BUFLEN = (1 << 21) + 1;
    char buf[BUFLEN], *st, *ed;
#define gc() ((st == ed) ? ((ed = (st = buf) + fread(buf, 1, BUFLEN, stdin)), ((st == ed) ? -1 : *st++)) : *st++)
    inline bool blank(char x) { return x == '\n' || x == ' ' || x == '\r' || x == '\t'; }
    template <typename I> inline bool read(I& x) {
        register char c;
        while (blank(c = gc()));
        if (c == -1)  return false;
        for (x = c - '0'; (c = gc()) >= '0' && c <= '9'; x = x * 10 + (c & 15));
        return true;
    }
    inline bool read() { return true; }
    template <typename T, typename... U> inline bool read(T& head, U&... tail) {
        return read(head) && read(tail...);
    }
    char obuf[BUFLEN], *ost = obuf, *oed = obuf + BUFLEN - 1, stk[55], top;
    inline void flush() { fwrite(obuf, 1, ost - obuf, stdout); ost = obuf; }
    inline void pc(char x) { *ost++ = x; if (ost == oed) flush(); }
    template <typename I> inline void print(I x) {
        if (!x) pc('0');
        while (x) stk[++top] = x % 10 + '0', x /= 10;
        while (top) pc(stk[top--]);
    }
    template <typename I> inline void print(I x, int b) { print(x), pc(" \n"[b]); }
    template <typename I> inline void println(I x) { print(x), pc('\n'); }
    struct IOFLUSHER { ~IOFLUSHER() { flush(); } } _ioflusher_;
}
using io::read;
using io::flush;
using io::print;
using io::println;
