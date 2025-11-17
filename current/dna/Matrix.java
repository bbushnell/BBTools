package dna;
import java.util.HashMap;
import java.util.Set;

import fileIO.MatrixFile;


/**
 * Genomic matrix management class for storing and retrieving float-based scoring matrices.
 * Provides static access to matrices loaded from text files with efficient HashMap-based lookup.
 * Supports extraction of matrix sub-grids for genomic analysis operations.
 * @author Brian Bushnell
 */
public class Matrix {
	
	
	/**
	 * Constructs a Matrix with the specified grid data, prefix length, and name.
	 * @param g The 2D float array containing matrix values
	 * @param pre The prefix length for matrix operations
	 * @param nm The name identifier for this matrix
	 */
	public Matrix(float[][] g, int pre, String nm){
		grid=g;
		prefix=pre;
		name=nm;
	}
	
	/**
	 * Extracts a sub-grid from the matrix with specified prefix length and dimensions.
	 * Creates a new 2D array containing cloned rows from the original grid.
	 *
	 * @param prefixLength The prefix length to account for when determining start position
	 * @param length The number of rows to extract from the matrix
	 * @return A new 2D float array containing the extracted sub-grid
	 */
	public float[][] subGrid(int prefixLength, int length){
		float[][] r=new float[length][];
		int start=prefix-prefixLength;
		for(int i=0; i<length; i++){
			r[i]=grid[i+start].clone();
		}
		return r;
	}
	
	/** The 2D float array containing the matrix values */
	public float[][] grid;
	/** The prefix length used for matrix operations and sub-grid extraction */
	public int prefix;
	/** The name identifier for this matrix */
	public String name;
	
	
	
	
	/** Static HashMap storing all loaded matrices with name-based lookup */
	private static HashMap<String, Matrix> table=null;
	
	/** Returns the set of matrix names currently loaded in the static table.
	 * @return Set of matrix names as keys from the internal HashMap */
	public static Set<?> keys(){return table.keySet();}
	
	/**
	 * Retrieves a matrix by name from the static table, loading default matrices if needed.
	 * Lazy-loads standard genomic matrices from build37 files and splice percentile data.
	 * Supports case-insensitive lookup by storing both original and lowercase keys.
	 *
	 * @param s The name of the matrix to retrieve
	 * @return The Matrix object associated with the given name
	 * @throws RuntimeException If the matrix name is not found or the value is null
	 */
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
	
	/**
	 * Loads matrices from a file and adds them to the static table.
	 * Reads matrices using MatrixFile and stores both original and lowercase name keys.
	 * @param fname The filename containing matrix data to load
	 */
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
