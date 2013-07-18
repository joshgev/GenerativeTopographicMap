package edu.michigan.gevirtz.data;

import java.util.Arrays;

public class DataVector {

	private double [] values;
	private int label = -1;
	
	public DataVector(double [] values){
		
		this.values = values;
		
	}
	
	public void SetValues(double [] values){
		this.values = values;
	}
	
	public double [] getValues(){
		return values;
	}
	public void SetLabel(int label){
		this.label = label;
	}

	public int GetLabel(){
		return this.label;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(values);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		DataVector other = (DataVector) obj;
		if (!Arrays.equals(values, other.values))
			return false;
		return true;
	}
	
	
}
