package edu.michigan.gevirtz.latent;

import java.awt.image.SampleModel;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.xml.crypto.Data;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jblas.Singular;
import org.jblas.Solve;

import edu.michigan.gevirtz.data.DataVector;
import edu.michigan.gevirtz.optimization.BFGSOptimizer;
import edu.michigan.gevirtz.optimization.OptimizationProblem;

public class GTM {

	private List<DataVector> data = null;
	private List<DoubleMatrix> dataVectors = null;		//holds the training data in jBLAS vectors
	private DoubleMatrix T = null;			//Holds the data in a matrix;

	private int L = 0; // Dimensionality of latent space
	private int M = 0; // Number of basis functions
	private int K = 0; // Number of points in latent spaces
	private int D = 0; // Dimensionality of observable space
	private int N = 0;

	private int K_per_latent_dimension = 0;
	private int M_per_latent_dimension = 0;

	private double sigma = 0; // Width of gaussian basis functions

	private DoubleMatrix phiMatrix = null; // Matrix containing evaluation of
											// all M
	private DoubleMatrix W = null; // Matrix defining y(x->t)
	private double beta = 0;		// variance of distribution of t
	// basis functions on all K latent space
	// sample points

	private DoubleMatrix[] basisFunctionCenters = null; // Centers for the basis
														// functions

	private DoubleMatrix[] latentSpaceSamplePoints = null; // Center for
																	// each x_i
																	// in the
																	// latent
																	// space
	


	/**
	 * 
	 * @param data
	 *            Data to train the latent model on
	 * @param L
	 *            Dimensionality of the latent space
	 * 
	 */
	
	private void SubtractMean(List<DataVector> data){
		double [] totals = new double[data.get(0).getValues().length];
		for(int i = 0; i < data.get(0).getValues().length; i++)
			totals[i] = 0;
		for(DataVector v : data)
			for(int i = 0; i < v.getValues().length;i++)
				totals[i] += v.getValues()[i];
		for(int i = 0; i < data.get(0).getValues().length; i++)
			totals[i] /= data.size();
		for(int i = 0; i < data.size(); i++){
			double [] values = data.get(i).getValues();
			for(int j = 0; j < values.length; j++)
				values[j] -= totals[j];
			
//			data.remove(i);
//			data.add(i, new DataVector(values));
			data.get(i).SetValues(values);
		}
		
	}
	public GTM(List<DataVector> data, int L) {

		SubtractMean(data);
		this.data = data;
		
		
		this.L = L;
		
		this.N = data.size();
		
		this.dataVectors = new ArrayList<DoubleMatrix>();
		for(DataVector dv : data)
			this.dataVectors.add( new DoubleMatrix(dv.getValues()) );
		
		this.T = dataVectors.get(0).transpose();
		for(int i = 1; i < dataVectors.size(); i++)
			T = DoubleMatrix.concatVertically(T, dataVectors.get(i).transpose());
		
		System.out.println("dataVectors.size(): "+dataVectors.size());
		System.out.println("T: "+T.rows+"x"+T.columns);
		this.D = data.get(0).getValues().length;
	}

	// Setters
	public void setSigma(double sigma) {
		this.sigma = sigma;
	}

	public void setM(int M_per_latent_dimension) {

		this.M_per_latent_dimension = M_per_latent_dimension;
		this.M = (int) Math.pow(M_per_latent_dimension, L);

	}

	public void setK(int K_per_latent_dimension) {

		this.K_per_latent_dimension = K_per_latent_dimension;
		this.K = (int) Math.pow(K_per_latent_dimension, L);

	}

	private void CalculatePhiMatrix() {
		
		phiMatrix = new DoubleMatrix(K,M);//CalculatePhiVector(latentSpaceSamplePoints[0]).transpose();
		for(int i = 0; i < K; i++)
			phiMatrix.putRow(i,CalculatePhiVector(latentSpaceSamplePoints[i]).transpose());
		
		
	}
	
	public List<DataVector> ProjectDataToLatentSpace(){
		
		
	
		return ProjectDataVectorToLatentSpace(this.data);
		
	}
	public List<DataVector> ProjectDataVectorToLatentSpace(List<DataVector> vectors){
		System.out.println("ProjectDataVectorToLatentSpace");
		List<DataVector> result = new ArrayList<DataVector>();
		System.out.print("\tCalculating R...");
		DoubleMatrix R = CalculateR(W, beta, vectors);
		System.out.println("done");
		int N = vectors.size();
		System.out.print("\tAdding contributions...");
		for(int i = 0; i < N; i++){
			DoubleMatrix vector = new DoubleMatrix(L);
			for(int j = 0 ; j < K; j++)
				vector.addi(latentSpaceSamplePoints[j].mul(R.get(j,i)) );
			DataVector dv = new DataVector(vector.toArray());
			dv.SetLabel(vectors.get(i).GetLabel());
			result.add(dv);
		}
		System.out.println("done");
		return result;
	}
	public DoubleMatrix ProjectToDataSpace(DoubleMatrix t){
		
		return W.mmul(CalculatePhiVector(t));
	}

	private DoubleMatrix CalculatePhiVector(DoubleMatrix x) {

		DoubleMatrix phi = new DoubleMatrix(M);
		for (int i = 0; i < M; i++)
			phi.put(i, calculateGaussian(x, basisFunctionCenters[i], sigma, L));
		return phi;
	}

	// Getters

	// Internal

	/**
	 * Calls GenerateBasisFunctionCenters, setups up W matrix, calls
	 * CalculatePhi
	 */

	public void InitializeModel() {
		
		System.out.println("Initializing model...");

		if (sigma == 0)
			throw new RuntimeException(
					"Attempted to initialize model withouth setting basis function width, sigma");

		System.out.print("\tGenerating basis function centers...");
		GenerateBasisFunctionCenters();
		System.out.println("done");
		System.out.print("\tGenerating latent space smaple points...");
		GenerateLatentSpaceSamplePointCenters();
		System.out.println("done");
		System.out.print("\tCalculating phi matrix...");
		CalculatePhiMatrix();
		System.out.println("done");
		
		//Find the initial W parameters
		W = DoubleMatrix.rand(D, M);
		
		ParameterInitializer paramInitializer = new ParameterInitializer(this, 2000);
		BFGSOptimizer optimizer = new BFGSOptimizer();
		DoubleMatrix W_new = optimizer.Minimize(paramInitializer, new DoubleMatrix(W.toArray()));
		//Put the column-vector representation of W_new into W 
		for(int i = 0; i < D; i++)
			for(int j = 0; j < M ; j ++)
				W.put(i ,j , W_new.get(i * M + j));

		
		//find initial beta parameter
		PCA pca = new PCA(data.subList(0, 1000));
		
		double beta_inv =pca.GetLargestEigenvalue();
		
		this.beta = 1.0 / beta_inv;
		System.out.println("done");

	}

	public double CalculateDataSpaceDistribution(DoubleMatrix t, DoubleMatrix W, double beta){
		
		double result = 0;
		for(int x_i = 0 ; x_i < K; x_i++){
		
			result += CalculateDataSpaceDistribution(t, x_i, W, beta);
					
		}
		
		return result / K;
		
	}
	public double CalculateDataSpaceDistribution(DoubleMatrix t, DoubleMatrix x, DoubleMatrix W, double beta){
		

		return calculateGaussian(t, W.mmul(CalculatePhiVector(x)), beta, this.D);
	}
	public double CalculateDataSpaceDistribution(DoubleMatrix t, int x_i, DoubleMatrix W, double beta){
		

		return calculateGaussian(t, W.mmul(phiMatrix.getRow(x_i).transpose()), beta, this.D); //Changed from calculating phi(x_i) every thime
	}
	public double CalculateDataSpaceDistribution(DoubleMatrix t){
		
		return CalculateDataSpaceDistribution(t, this.W, this.beta);
		
	}
	
	private DoubleMatrix CalculateR(DoubleMatrix W, double beta){
		
		return CalculateR(W, beta, this.data);
	}
	private DoubleMatrix CalculateR(DoubleMatrix W, double beta, List<DataVector> data){

		int N = data.size();
		DoubleMatrix R = new DoubleMatrix(K, N);
		
		
		double [] denominators = new double[N];
		
		
		for(int i = 0; i < N; i++){
			DoubleMatrix vector = dataVectors.get(i);	//Changed from creating a new thing based on the data
			denominators[i] = CalculateDataSpaceDistribution(vector, W, beta) * K;
		}
		
		for(int i = 0; i < K; i++){
			for(int j = 0; j < N; j++){
//				DoubleMatrix vector = new DoubleMatrix(data.get(j).getValues());
				DoubleMatrix vector = dataVectors.get(j);	//Changed from creating a new thing based on the data
				R.put(i,j, CalculateDataSpaceDistribution(vector, latentSpaceSamplePoints[i], W, beta) / denominators[j]);
			}
		}
		
		return R;
	}
	private DoubleMatrix CalculateG(DoubleMatrix W, double beta){
		
		DoubleMatrix R = CalculateR(W, beta);
		DoubleMatrix G = new DoubleMatrix(K,K);

		for(int i = 0; i < K; i++){
			double total = 0.0;
			for(int j = 0; j < N; j++)
				total += R.get(i, j);
			G.put(i,i,total);
		}
		
		return G;
	}
	/**
	 * Sets the lattice of basis functions
	 */
	private void GenerateBasisFunctionCenters() {
		if (M == 0)
			throw new RuntimeException(
					"Attempted to generate basis function centers without first setting M");
		if (L == 0)
			throw new RuntimeException(
					"Attempted to generate basis function centers without first setting L");

		int[] currentIndices = new int[L];
		for (int i = 0; i < L; i++)
			currentIndices[i] = 0;
		double spacing = 2.0 / (M_per_latent_dimension - 1);

		basisFunctionCenters = new DoubleMatrix[M];

		for (int i = 0; i < M; i++) {

			// update the indices
			for (int j = 0; j < L; j++) {
				int num = (int) Math.pow(M_per_latent_dimension, j);
				currentIndices[j] = (i / num) % M_per_latent_dimension;
			}
			DoubleMatrix center = new DoubleMatrix(L);
			// loop through each dimension, set the position of this center for
			// each dimension based on the current values of "currentIndex"
			for (int j = 0; j < L; j++)

				center.put(j, spacing * currentIndices[j] - 1);

//			System.out.println("Basis function center: "+center);
			basisFunctionCenters[i] = center;
			
		}
	}

	private void GenerateLatentSpaceSamplePointCenters() {
		if (K == 0)
			throw new RuntimeException(
					"Attempted to generate latent space sample point centers without first setting K");
		if (L == 0)
			throw new RuntimeException(
					"Attempted to generate latent space sample point centers without first setting L");

		int[] currentIndices = new int[L];
		for (int i = 0; i < L; i++)
			currentIndices[i] = 0;
		double spacing = 2.0 / (K_per_latent_dimension - 1);

		latentSpaceSamplePoints = new DoubleMatrix[K];

		for (int i = 0; i < K; i++) {
			// update the indices
			for (int j = 0; j < L; j++) {
				int num = (int) Math.pow(K_per_latent_dimension, j);
				currentIndices[j] = (i / num) % K_per_latent_dimension;
			}
			DoubleMatrix center = new DoubleMatrix(L);
			// loop through each dimension, set the position of this center for
			// each dimension based on the current values of "currentIndex"
			for (int j = 0; j < L; j++)

				center.put(j, spacing * currentIndices[j] - 1);
			
//			System.out.println("latent space sample point: "+center);

			latentSpaceSamplePoints[i] = center;
	
		}
	}

	private double calculateGaussian(DoubleMatrix point, DoubleMatrix center,
			double width, int D) {

		if (point.columns != 1)
			throw new RuntimeException("Point not a column vector");
		if (center.columns != 1)
			throw new RuntimeException("Center not a column vector!");
		if (point.columns != center.columns)
			throw new RuntimeException("Dimensions are not equal");

		return Math.pow(width / (2 * 3.14159), D / 2.0)
				* Math.exp(-1*center.squaredDistance(point));

	}

	private DoubleMatrix SVDInvert(DoubleMatrix matrix){
		
		DoubleMatrix [] svd = Singular.fullSVD(matrix);
		
		DoubleMatrix U = svd[0];
		DoubleMatrix W = svd[1];
		DoubleMatrix V = svd[2];
		
		
		DoubleMatrix W_inverse = new DoubleMatrix(W.rows, W.rows);
		for(int i = 0; i < W.rows; i++)
			W_inverse.put(i,i,1.0 / W.get(i,0));
		return V.mmul(W_inverse.mmul(U.transpose()));
	}
	
	private void PerformEMStep(){
		System.out.println("Performing EM Step");
		
		System.out.print("\tCalculating R using old parameters...");
		DoubleMatrix R_old = CalculateR(W, beta);
		System.out.println("done");
		System.out.print("\tCalculating A using old parameters...");
		DoubleMatrix A = phiMatrix.transpose().mmul(CalculateG(W, beta).mmul(phiMatrix));
		System.out.println("done");
		System.out.print("\tCalculating B using old parameters...");
		DoubleMatrix B = phiMatrix.transpose().mmul(R_old.mmul(T));
		System.out.println("done");

		System.out.print("\tInverting to get new W...");
		W =  SVDInvert(A).mmul(B).transpose();
		System.out.println("done");

		System.out.print("\tUpdating beta...");
		double beta_inv = 0;
		for(int i = 0; i < N; i++)
			for(int j = 0; j < K; j++){
				DoubleMatrix phi = phiMatrix.getRow(j).transpose();
				beta_inv += R_old.get(j,i) * W.mmul(phi).squaredDistance(dataVectors.get(i));
			}
		beta_inv /= (D * N);
		
		beta = 1.0 / beta_inv;
		System.out.println("done");
	}
	private void WriteDataVectors(List<DataVector> dataVectors, String filename){
		try{
			BufferedWriter writer = new BufferedWriter( new FileWriter(filename) );
			for(DataVector vector : dataVectors){
				writer.write(vector.GetLabel() + " ");
				double [] values = vector.getValues();
				for(int i = 0; i < values.length; i++)
					writer.write(values[i]+" ");
				writer.write("\n");
					
			}
			writer.close();
			
		}catch( IOException e){
			e.printStackTrace();
		}
	}
	
	
	private void WritePosteriorDistribution(List<DataVector> data, String filename){
		
		DoubleMatrix R = CalculateR(W, beta, data);
		Map<Integer, Double> values = new HashMap<Integer, Double>();
		BufferedWriter writer = null;
		try{
			writer = new BufferedWriter(new FileWriter(filename));
		}catch(IOException e){
			e.printStackTrace();
			return;
		}
		for(int i = 0; i < K; i++){
			values.clear();
			for(int n = 0; n < data.size(); n++){
				
				int label = data.get(n).GetLabel();
				if(!values.containsKey(label))
					values.put(label, 0.0);
				
				double value = values.get(label);
				value += R.get(i,n);
				
				values.put(label, value );
				
				
			}
			for(int label : values.keySet()){
				DoubleMatrix center = latentSpaceSamplePoints[i];
				try{
					
					writer.write(label+" ");
					writer.write(values.get(label) / data.size()+" ");
					for(int l = 0; l < L; l++)
						writer.write(center.get(l)+" ");
					writer.write("\n");
					
				}catch(IOException e){
					e.printStackTrace();
					return;
				}
			}
		}
		try{
			writer.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	private class ParameterInitializer implements OptimizationProblem {

		private GTM gtm = null;
		private int n = 0;
		private DoubleMatrix inversePCA = null;
		private List<DoubleMatrix> transformedLatentSamplePoints = null;

		/**
		 * 
		 * @param gtm
		 * @param n
		 *            number of the subset of training vectors to use for
		 *            parameter initialization. Not all training vectors are
		 *            needed; using a subsample increases speed
		 */
		public ParameterInitializer(GTM gtm, int n) {

			this.gtm = gtm;
			this.n = n;
			
			PCA pca = new PCA(gtm.data.subList(0, n));

			inversePCA = pca.getTransform(L).transpose();
			transformedLatentSamplePoints = new ArrayList<DoubleMatrix>();
			for(DoubleMatrix point : gtm.latentSpaceSamplePoints)
				transformedLatentSamplePoints.add(inversePCA.mmul(point));
			
//			NormalizeTransformedLatentSamplePoints();
			
		}
		
		private void NormalizeTransformedLatentSamplePoints(){
			double [] max = new double[gtm.D];
			double [] min = new double[gtm.D];
			double [] mean = new double[gtm.D];
			for(int i = 0; i < max.length; i++){
				max[i] = Double.MIN_VALUE;
				min[i] = Double.MAX_VALUE;
			}
			for(DoubleMatrix point : transformedLatentSamplePoints)
				for(int i = 0; i < gtm.D; i++){
					mean[i] += point.get(i);
					if(Math.abs(point.get(i)) > max[i])
						max[i] = Math.abs(point.get(i));
					if(Math.abs(point.get(i)) < min[i])
						min[i] = Math.abs(point.get(i));
				}
			double [] norms = new double[gtm.D];
			for(int i = 0; i < gtm.D; i++){
				norms[i] = max[i] - min [i];
				mean[i] /= transformedLatentSamplePoints.size();
			}
			
			
			for(DoubleMatrix point : transformedLatentSamplePoints)
				for(int i = 0; i < gtm.D; i++){
					double value = point.get(i);
					point.put(i, (value - mean[i])/ (norms[i] / 2) + mean[i]) ;
				}
					
				
		}

		/**
		 * x contains the elements of W
		 */
		@Override
		public double Function(DoubleMatrix x) {

			DoubleMatrix W = new DoubleMatrix(D,M);
			for(int i = 0; i < D; i++)
				for(int j = 0; j < M; j++)
					W.put(i, j, x.get(i * M + j));
			
			double total = 0.0;
			//Loop over all latent space sample points
			for(int i = 0; i < K; i++){
				DoubleMatrix transformedLatentSpaceSamplePoint = transformedLatentSamplePoints.get(i);
				total += W.mmul(gtm.phiMatrix.getRow(i).transpose()).distance2(transformedLatentSpaceSamplePoint );
			}
			return .5 * total / gtm.K; // / gtm.M;
		}

		@Override
		public DoubleMatrix Gradient(DoubleMatrix x) {
			DoubleMatrix W = new DoubleMatrix(D,M);
			for(int i = 0; i < D; i++)
				for(int j = 0; j < M; j++)
					W.put(i, j, x.get(i * M + j));
			DoubleMatrix gradient = new DoubleMatrix(D*M);
			DoubleMatrix vector = new DoubleMatrix(D*M);
			for(int i = 0; i < K; i++){
				DoubleMatrix phi = gtm.phiMatrix.getRow(i).transpose();
				DoubleMatrix Ux = transformedLatentSamplePoints.get(i);
				DoubleMatrix Wphi = W.mmul(phi);
				double denominator = Wphi.distance2(Ux);
				//D x M
				for(int j = 0 ; j < D * M; j++){
					int a = j / M;
					int b = j % M;
					
					vector.put(j, .5*phi.get(b)*( Wphi.get(a) - Ux.get(a) ) / denominator);
					
				}
				gradient.addi(vector);
			}
			
			return gradient;
		}

	}
	
	
	private static void test1(){
		List<DataVector> data = new ArrayList<DataVector>();
		Random rand = new Random();
		int N = 10000;
		for (int i = 0; i < N; i++){
			double x = rand.nextDouble();
			double y = x + ( rand.nextDouble() - .5 ) / .5 * .01;
			double [] values = {x,y};
			data.add(new DataVector(values));
			
		}
		//Initialize the GTM
		GTM gtm = new GTM(data, 1);
		System.out.println("df");
		gtm.setSigma(.01);
		gtm.setK(200);
		gtm.setM(2);
		gtm.InitializeModel();

		List<DoubleMatrix> testingData = new ArrayList<DoubleMatrix>();
		for (int i = 0; i < 1000; i++){
			double x = rand.nextDouble();
			double y = x + ( rand.nextDouble() - .5 ) / .5 * .01;
			double [] values = {x};
			DoubleMatrix latentPoint = new DoubleMatrix(values);
			
			testingData.add(latentPoint);
		}
		
		
		//Create a PCA transform from latent to data space:
		PCA pca = new PCA(data.subList(0, 1000));
		DoubleMatrix inversePCA = pca.getTransform(2).transpose();

		//First test PCA.  Write stuff to "pca.txt"
		try{
			BufferedWriter pcaWriter = new BufferedWriter(new FileWriter("pca.txt"));

			for(int i = 0; i < 1000; i++){
				DoubleMatrix dataSpaceProjectionPCA = inversePCA.mmul(testingData.get(i));
				for(int j = 0; j < dataSpaceProjectionPCA.rows; j++)
					pcaWriter.write(dataSpaceProjectionPCA.get(j)+" ");
				pcaWriter.write('\n');
				
		}
			pcaWriter.close();
			BufferedWriter gtmWriter = new BufferedWriter(new FileWriter("gtm.txt"));
			for(int i = 0; i < 1000; i++){
				DoubleMatrix dataSpaceProjectionGTM = gtm.ProjectToDataSpace(testingData.get(i));
				for(int j = 0; j < dataSpaceProjectionGTM.rows; j++)
					gtmWriter.write(dataSpaceProjectionGTM.get(j)+" ");
				gtmWriter.write('\n');
			}
			gtmWriter.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		//Now perform a single step of the EM aglorithm
		
		gtm.PerformEMStep();
		
	}
	private static void test2(String [] argv){
		
		
		String dataFile = argv[0];
		List<DataVector> data = new ArrayList<DataVector>();
		//Read data in
		try{
			BufferedReader reader = new BufferedReader(new FileReader(dataFile));
			String line = "";
			while((line = reader.readLine()) != null){
				String [] elements = line.split("\\s+");
				double [] values = new double[elements.length - 1];
				for(int i = 1; i < elements.length; i++)
					values[i-1] = Double.valueOf(elements[i]);
				DataVector vector = new DataVector(values);
				vector.SetLabel(Integer.valueOf(elements[0]));
				data.add(vector);
			}
				
			reader.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		GTM gtm = new GTM(data,2);
		gtm.setSigma(Double.valueOf(argv[1])); //.5
		gtm.setK(Integer.valueOf(argv[2])); //50
		gtm.setM(Integer.valueOf(argv[3]));	//4
		gtm.InitializeModel();
		System.out.print("Writing pca transformed latent space sample points...");
		PCA pca = new PCA(data);
		DoubleMatrix inversePCA = pca.getTransform(2).transpose();
		List<DataVector> transformedSamplePoints = new ArrayList<DataVector>();
		for(DoubleMatrix vector : gtm.latentSpaceSamplePoints){
			DataVector dv = new DataVector( inversePCA.mmul(vector).toArray() );
			
			transformedSamplePoints.add(dv);
		}
		String baseName = "gtm_output/";
		
		gtm.WriteDataVectors(transformedSamplePoints, baseName+"transformed_sample_points.txt");
		System.out.println("done");
		int N = 1000;
		
		for(int i = 0; i < N; i++){
//			System.out.println(gtm.W);
			List<DataVector> projection = gtm.ProjectDataToLatentSpace();
			
			
			gtm.WriteDataVectors(projection, baseName+String.format("%03d", i)+".txt");
			gtm.WritePosteriorDistribution(projection, baseName+String.format("prior/%03d", i)+".txt");
			
			//Look at how latent sample points transform
			transformedSamplePoints.clear();
			System.out.print("Projecting latent space sample points to data space...");
			for(DoubleMatrix vector : gtm.latentSpaceSamplePoints){
				DoubleMatrix phi = gtm.CalculatePhiVector( vector );
				DataVector dv = new DataVector( gtm.W.mmul(phi).toArray() );
				
				transformedSamplePoints.add(dv);
			}
			System.out.println("done.");
			
			gtm.WriteDataVectors(transformedSamplePoints, baseName+"latent/"+String.format("%03d", i)+".txt");
			System.out.println("Performing EM step!");
			System.out.println(gtm.beta);
			gtm.PerformEMStep();
		}
		
		
	}
	public static void main(String [] argv){
		
		test2(argv);
		
		
	}

}
