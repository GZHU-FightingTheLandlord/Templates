import java.io.BufferedInputStream;
import java.io.PrintWriter;
import java.util.NoSuchElementException;
import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        Scanner in = new Scanner(new BufferedInputStream(System.in));
        PrintWriter out = new PrintWriter(System.out);
        Task solver = new Task();
        try {
            solver.solve(in, out);
        } catch (NoSuchElementException e) {
            out.close();
        }
    }
    private static BigInteger[] a = new BigInteger[1005];
    public static class Task {
        public static void solve(Scanner in, PrintWriter out) {
            // Your code start here.
            
        }
    }
}