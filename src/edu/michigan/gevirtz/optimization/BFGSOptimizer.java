package edu.michigan.gevirtz.optimization;

import org.jblas.DoubleMatrix;
import org.jblas.Solve;

public class BFGSOptimizer implements Minimizer{
	
	private OptimizationProblem problem = null;
	
	private double alpha = .1	;
	

	int k = 100; 	//number of steps in to search while doing line search for new parameter vector
	
	public BFGSOptimizer(){
		
	}
	
	//Setters
	/**
	 * Sets the largest (in magnitude) alpha to be used in updating the current parameter vector 
	 * @param alpha
	 */
	public void SetMaximumAlpha(double alpha){
		this.alpha = alpha;
	}
	
	public void SetK(int k){
		this.k = k;
	}
	
	private DoubleMatrix GetDirection(DoubleMatrix B, DoubleMatrix x){
		DoubleMatrix grad = problem.Gradient(x).mul(-1.0);
		return Solve.solve(B, grad);
	}
	

	private double PerformLineSearch(DoubleMatrix x, DoubleMatrix direction){
		double currentAlpha = 0;
		double currentMinimum = Double.MAX_VALUE;
		double step = alpha / k;
		for(int i = 1; i <= k; i++ ){
			DoubleMatrix new_x = x.add(direction.mul(step * i));
			double val = problem.Function(new_x);
			if(val < currentMinimum){
				currentAlpha = val;
				currentAlpha = step * i;
			}
		}

		return currentAlpha;

	}
	@Override
	public DoubleMatrix Minimize(OptimizationProblem problem, DoubleMatrix x0) {
		
		this.problem = problem;
		
		int N = x0.rows;
		
		DoubleMatrix B = DoubleMatrix.eye(N);
		DoubleMatrix x = x0.dup();
		
//		while(problem.Gradient(x).norm2() >= .05){
		while(problem.Function(x) >= .05){
			DoubleMatrix direction = GetDirection(B, x);
			double norm = direction.norm2();
			direction.muli(1/norm);
			
			double alpha = PerformLineSearch(x, direction);
			
			DoubleMatrix s = direction.mul(alpha);
			DoubleMatrix x_new = x.add( s );
			
			
			DoubleMatrix y = problem.Gradient(x_new).sub(  problem.Gradient(x) );
			
			double denominator1 = y.transpose().mmul(s).get(0);

			double denomintaor2 = s.transpose().mmul(B.mmul(s)).get(0);
			
			DoubleMatrix B_new = B.add(y.mmul(y.transpose()).mul( 1 / denominator1 ));
			B_new.subi(B.mmul(s.mmul(s.transpose().mmul(B))).mul( 1 / denomintaor2));
			 
			B = B_new;
			x = x_new;
			
			System.out.println("Grad: "+problem.Gradient(x).norm2());
			System.out.println("Func: "+problem.Function(x));
		}
		
		return x;
	}

}
