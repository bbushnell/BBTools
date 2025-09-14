package dna;
import java.util.HashMap;
import java.util.Set;

import fileIO.MatrixFile;


public class Matrix {
	
	
	public Matrix(float[][] g, int pre, String nm){
		grid=g;
		prefix=pre;
		name=nm;
	}
	
	public float[][] subGrid(int prefixLength, int length){
		float[][] r=new float[length][];
		int start=prefix-prefixLength;
		for(int i=0; i<length; i++){
			r[i]=grid[i+start].clone();
		}
		return r;
	}
	
	public float[][] grid;
	public int prefix;
	public String name;
	
	
	
	
	private static HashMap<String, Matrix> table=null;
	
	public static Set<?> keys(){return table.keySet();}
	
	public static Matrix get(String s){
		if(table==null){
			table=new HashMap<String, Matrix>(64);
//			fillTable("matrices.txt");
//			fillTable("matrices2.txt");

//			fillTable("matrixN1.txt");
//			fillTable("matrixN2.txt");
//			fillTable("matrixN3.txt");
//			fillTable("matrixN4.txt");
			
			fillTable("matrix_build37_N1.txt");
			fillTable("matrix_build37_N2.txt");
			fillTable("matrix_build37_N3.txt");
//			fillTable("matrix_build37_N4.txt");
			
			

//			fillTable("asmGstart_sept9.txt");
//			fillTable("asmEstart_sept9.txt");
//			fillTable("asmTRstart_sept9.txt");
//			fillTable("asmGstop_sept9.txt");
//			fillTable("asmEstop_sept9.txt");
//			fillTable("asmTRstop_sept9.txt");
//			fillTable("asmEstop_sept16.txt");

//			fillTable("SplicePercentiles_b37_Sept16.txt");
			fillTable("SplicePercentiles_b37_Nov24.txt");
			
		}
		Matrix m=table.get(s);
		
//		assert(table.containsKey(s)) : "\nCan't find "+s+" in\n\n"+table.keySet()+"\n";
//		assert(m!=null) : "\nValue for "+s+" is null\n";
		
		if(!table.containsKey(s) || m==null){
			if(!table.containsKey(s)){throw new RuntimeException("Can't find "+s+" in\n\n"+table.keySet()+"\n");}
			if(m==null){throw new RuntimeException("\nValue for "+s+" is null");}
		}
		
		
		return m;
	}
	
	private static void fillTable(String fname){
		MatrixFile mf=new MatrixFile(fname);
		Matrix mat=mf.nextMatrix();
		while(mat!=null){
//			System.out.println("Adding "+mat.name);
			table.put(mat.name, mat);
			table.put(mat.name.toLowerCase(), mat);
			mat=mf.nextMatrix();
		}
		mf.close();
	}
	
}
