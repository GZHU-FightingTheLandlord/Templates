public static BigInteger rootM(BigInteger n, final int m) {
    final String tmp = n.toString();
    BigDecimal x = new BigDecimal(tmp.substring(0, tmp.length() / m + (m == 1 ? 0 : 1)));
    BigDecimal l = BigDecimal.ZERO;
    final BigDecimal M = BigDecimal.valueOf(m);
    final BigDecimal N = new BigDecimal(n), eps = BigDecimal.valueOf(1e-6);
    while(x.subtract(l).abs().compareTo(eps) > 0) {
        l = x;
        x = x.subtract(x.pow(m).subtract(N).divide(M.multiply(x.pow(m - 1)), 50, BigDecimal.ROUND_HALF_EVEN));
    }
    return x.toBigInteger();
}
// https://www.lydsy.com/JudgeOnline/problem.php?id=1213
