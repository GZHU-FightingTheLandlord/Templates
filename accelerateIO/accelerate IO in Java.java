import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.math.*;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

public class Main {
	public static void main(String[] args) {
		Scanner in = new Scanner(new BufferedInputStream(System.in));
		PrintWriter out = new PrintWriter(Sysytem.out);
		Task solver = new Task();
		try {
			solver.solve(in, out);
		} catch (NoSuchElementException e) {
			out.close();
		}
	}
	public static class Task {
		public static void solve(Scanner in, PrintWriter out) {
			// Your code start here.
			
		}
	}
}
