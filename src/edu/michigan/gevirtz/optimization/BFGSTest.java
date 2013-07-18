package edu.michigan.gevirtz.optimization;

import org.jblas.DoubleMatrix;

public class BFGSTest {

	class QuadraticProblem implements OptimizationProblem {

		public QuadraticProblem() {
		}

		@Override
		public double Function(DoubleMatrix x) {
			// TODO Auto-generated method stub
			double total = 0;
			for (int i = 0; i < x.rows; i++)
				total += x.get(i) * x.get(i);

			return total;
		}


		@Override
		public DoubleMatrix Gradient(DoubleMatrix x) {

			DoubleMatrix gradient = new DoubleMatrix(x.rows);
			for (int i = 0; i < x.rows; i++)
				gradient.put(i, 2 * x.get(i));

			return gradient;
		}

	}

	public BFGSTest() {
		QuadraticProblem quadratic = new QuadraticProblem();

		BFGSOptimizer optimizer = new BFGSOptimizer();
		
		DoubleMatrix x0 = new DoubleMatrix(10);
		x0.put(0, 1.0);
		x0.put(1, 1.0);
		x0.put(2, 1.0);
		x0.put(3, 1.0);
		x0.put(4, 1.0);
		x0.put(5, 1.0);
		x0.put(6, 1.0);
		x0.put(7, 1.0);
		x0.put(8, 1.0);
		x0.put(9, 1.0);
		
		
		optimizer.Minimize(quadratic, x0);
	}

	public static void main(String[] argv) {
		new BFGSTest();
	}
}
