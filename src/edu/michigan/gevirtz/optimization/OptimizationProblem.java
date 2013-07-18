package edu.michigan.gevirtz.optimization;

import org.jblas.DoubleMatrix;

public interface OptimizationProblem {
	
	//Returns parameters to be optimized
	public double Function(DoubleMatrix x);
	public DoubleMatrix Gradient(DoubleMatrix x);

}
