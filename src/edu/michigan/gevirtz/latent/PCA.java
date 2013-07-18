package edu.michigan.gevirtz.latent;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import edu.michigan.gevirtz.data.DataVector;

public class PCA {


	private DoubleMatrix dataMatrix = null;

	private DoubleMatrix covarianceMatrix = null;
	
	private DoubleMatrix eigenVectors = null;
	
	private DoubleMatrix eigenValues = null;
	
	private void BuildCovarianceMatrix(){
		
		//dataMatrix.columns is the number of variables there are per data vector
		double [][] result = new double[dataMatrix.columns][dataMatrix.columns];
		
		int N = dataMatrix.rows;
		
		for(int i = 0; i < dataMatrix.columns; i++)
			for(int j = 0; j < dataMatrix.columns; j++){
				
				double total = 0;
				for(int k = 0; k < dataMatrix.rows; k++)
					total += dataMatrix.get(k,i) * dataMatrix.get(k,j);  

				result[i][j] = total / ( N - 1);
				result[j][i] = total / ( N - 1);
				
			}
		
		covarianceMatrix = new DoubleMatrix(result);
		
		
	}
	
	
	public void CalculatePCA(){
		
		DoubleMatrix values = Eigen.eigenvectors(covarianceMatrix)[1].real();
		DoubleMatrix  vectors = Eigen.eigenvectors(covarianceMatrix)[0].real();
		

//		System.out.println(covarianceMatrix);
		List<EigenVectorEigenValuePair> pairs = new ArrayList<PCA.EigenVectorEigenValuePair>();
		for(int i = 0; i < vectors.columns; i++){
			pairs.add( new EigenVectorEigenValuePair(vectors.getColumn(i), values.get(i,i)));
//			System.out.println("Eigenvector: "+vectors.getColumn(i));
		}
		
		Collections.sort(pairs, pairs.get(0).comparator);
		
		//Create the column vector that will contain the eigenvalues
		eigenValues = new DoubleMatrix(pairs.size());
		//Fill the structures holding eigenvectors and eigenvalues
		int index = 0;
		for(EigenVectorEigenValuePair pair : pairs){
			
			//if eigenVectors is null, we need to add the first of the eigenvectors.
			if(eigenVectors == null){
				eigenVectors = pair.eigenvector.dup();
				continue;
			}
			//Otherwise, stick the next eigenvector at the end of the matrix
			eigenVectors = DoubleMatrix.concatHorizontally(eigenVectors, pair.eigenvector);
			eigenValues.put(index,0,pair.eigenvalue);
			index++;
			
			
			
		}
//		System.out.println("eigenValues: "+eigenValues);
		
	}
	
	public double GetLargestEigenvalue(){
		return eigenValues.get(0,0);
	}
		
	private void BuildDataMatrix(List<DataVector> data){

		int N = data.size();
		
		double [][] dataArray = new double[N][];
		
		for(int i = 0; i < data.size(); i++)
			dataArray[i] = data.get(i).getValues();
			
		
		dataMatrix = new DoubleMatrix(dataArray);
		
		//Subtract means
		DoubleMatrix columnMeans = dataMatrix.columnMeans();
		dataMatrix.subiRowVector(columnMeans);
	}
	
	public List<DataVector> transform(){
		
		List<DataVector> result = new ArrayList<DataVector>();
		DoubleMatrix product = eigenVectors.transpose().mmul( dataMatrix.transpose() );

		for(int i = 0; i < product.columns; i++){
			result.add( new DataVector(product.getColumn(i).toArray()) );
		}
		
		return result;
		
	}
	public List<DataVector> transform(List<DataVector> data, int L){
		
		List<DataVector> result = new ArrayList<DataVector>();
		
		DoubleMatrix transformation = getTransform(L);
		
		for(DataVector dataVector : data){
			DoubleMatrix vector = new DoubleMatrix(dataVector.getValues());
			DataVector newDataVector = new DataVector(transformation.mmul(vector).toArray());
			newDataVector.SetLabel( dataVector.GetLabel() );
			
			result.add(newDataVector);
			
		}
		
		return result;
		
	}
	
	public PCA(List<DataVector> data){
		
		BuildDataMatrix(data);
		
		BuildCovarianceMatrix();
		
		CalculatePCA();
		
	}
	
	public DoubleMatrix getTransform(int L){
		
		int rows = eigenVectors.rows;
		
		return eigenVectors.getRange(0,rows,0,L).transpose();
		
	}
	
	private class EigenVectorEigenValuePair{
		private Comparator<EigenVectorEigenValuePair> comparator = new Comparator<PCA.EigenVectorEigenValuePair>() {
			
			@Override
			public int compare(EigenVectorEigenValuePair o1,
					EigenVectorEigenValuePair o2) {
				// TODO Auto-generated method stub
				double val = (o2.eigenvalue - o1.eigenvalue); 
				return (int)(val / Math.abs(val));
			}
		};
		
		public double eigenvalue = 0;
		public DoubleMatrix eigenvector = null;
		
		public EigenVectorEigenValuePair(DoubleMatrix eigenvector, double eigenvalue){
			this.eigenvalue = eigenvalue;
			this.eigenvector = eigenvector;
		}
		
	}
	
	
	//Simple test to be run from commandline.  Loads file supplied as arg[0], reads data vectors, spits out PCA-transformed vectors
	public static void main(String [] args){
		
		List<DataVector> data = new ArrayList<DataVector>();
		
		try{
			BufferedReader reader = new BufferedReader( new FileReader(args[0]) );
			String line;
			while((line = reader.readLine()) != null){
				String [] elements  = line.split("\\s+");
				double [] values = new double[elements.length];;
				for(int i = 1; i < elements.length; i++)
					values[i] = Double.valueOf(elements[i]);
				DataVector dataVector = new DataVector( values ) ;
				dataVector.SetLabel(Integer.valueOf(elements[0]));
				data.add( dataVector );
			}
			
			reader.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(1);
		}
		
		PCA pca = new PCA(data);
		for(DataVector dataVector : pca.transform(data, 2)){
			double [] values = dataVector.getValues();
			System.out.print(dataVector.GetLabel()+" ");
			for(double value : values)
				System.out.print( value +" ");
			System.out.println();
		}
			
		
	}
}
