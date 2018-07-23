#include <stdio.h>
#include <string.h>
namespace io {
#define BUFLEN (1 << 21)
    bool EOFError;
    inline char gc() {
        static char buf[BUFLEN], *st = nullptr, *ed = nullptr;
        return (st == ed) ? ((ed = (st = buf) + fread(buf, 1, BUFLEN, stdin)), ((st == ed) ? -1 : *st++)) : *st++; 
    }
    inline bool check(char& x) { return x == '-' || x == '\n' || x == ' ' || x == '\r' || x == '\t'; }
    template <class I> inline void read(I& x) {
        char c; int f = 1;
        while (check(c = gc())) if (c == '-') f = -1;
        if (c == -1)  { EOFError = 1; return; }
        for (x = c - '0'; (c = gc()) >= '0' && c <= '9'; x = x * 10 + (c & 15)); x *= f;
    }
#define Isalpha(a) ((a >= 'a' && a <= 'z') || (a >= 'A' && a <= 'Z'))
    inline void gstr(char *s, int& len) {
        char c;
        for (c = gc(); Isalpha(c); c = gc());
        if (c == -1) return;
        for (len = 0; Isalpha(c); c = gc()) s[len++] = c;
        s[len] = 0;
    }
    char obuf[BUFLEN], *ost = obuf, *oed = obuf + BUFLEN - 1, Stack[55], Top;
    inline void flush() { fwrite(obuf, 1, ost - obuf, stdout); ost = obuf; }
    inline void pc(char& x) { *ost++ = x; if (ost == oed) flush(); }
    template <class I> inline void print(I x) {
        if (!x) pc('0');
        if (x < 0) pc('-'), x = -x;
        while (x) Stack[++Top] = x % 10 + '0', x /= 10;
        while (Top) pc(Stack[Top--]);
    }
    template <class I> inline void println(I x) { print(x), pc('\n'); }
    inline void pstr(char *s) { for (int i = 0, n = strlen(s); i < n; i++) pc(s[i]); }
    struct IOFLUSHER { ~IOFLUSHER() { flush(); } } _ioflusher_;
#undef BUFLEN
}
using io::EOFError;
using io::flush;
using io::gc;
using io::gstr;
using io::read;
using io::pc;
using io::print;
using io::println;
using io::pstr;