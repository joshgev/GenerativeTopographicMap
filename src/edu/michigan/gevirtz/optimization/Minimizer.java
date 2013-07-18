package edu.michigan.gevirtz.optimization;

import org.jblas.DoubleMatrix;

public interface Minimizer {
	
	public DoubleMatrix Minimize(OptimizationProblem problem, DoubleMatrix x0);

}
