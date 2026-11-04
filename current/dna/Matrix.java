package dna;

import java.util.HashMap;
import java.util.Set;

import fileIO.MatrixFile;

/**
 * Mutable scoring matrix and legacy registry of named genomic matrices.
 * Instances retain caller-owned rows; subGrid makes independent row copies.
 * Registry initialization and later mutation require external synchronization.
 * The registry is retained only after every default file loads and closes successfully;
 * after a catchable load failure, the next lookup retries from the first file.
 * @author Brian Bushnell
 */
public class Matrix{

	/**
	 * Retains the supplied data, prefix offset and name without copying or validation.
	 * @param g Borrowed two-dimensional array; rows may have different lengths
	 * @param pre Row offset corresponding to the reference position
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
	 * @param prefixLength Places the reference position this many rows into the result;
	 * the first copied row is prefix-prefixLength
	 * @param length Number of rows to copy; the selected range must exist and contain
	 * nonnull rows. A zero length copies nothing.
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

	/** Shared mutable rows supplied to the constructor. */
	public float[][] grid;
	/** The prefix length used for matrix operations and sub-grid extraction */
	public int prefix;
	/** The name identifier for this matrix */
	public String name;

	/** Lazily published after complete loading; original and lowercase aliases share values. */
	private static HashMap<String, Matrix> table=null;

	/** Returns a live key-set view; removing keys or clearing it modifies the registry.
	 * Does not initialize the registry; calling before get throws NullPointerException.
	 * @return Original names and lowercase aliases currently present */
	public static Set<?> keys(){return table.keySet();}

	/**
	 * Retrieves a matrix by name from the static table, loading default matrices if needed.
	 * Lazy-loads standard genomic matrices from build37 files and splice percentile data.
	 * Lookup is exact: each matrix is registered under its original spelling and the
	 * lowercase spelling produced with the default locale. Arbitrary mixed case is
	 * not normalized. The returned mutable object is shared by both aliases.
	 * A catchable initialization failure leaves the registry unset for a later retry.
	 *
	 * @param s The name of the matrix to retrieve
	 * @return The Matrix object associated with the given name
	 * @throws RuntimeException If loading/closure fails, or the requested matrix is absent
	 */
	public static Matrix get(String s){
		if(table==null){
			final HashMap<String, Matrix> loaded=new HashMap<String, Matrix>(64);
//			fillTable("matrices.txt");
//			fillTable("matrices2.txt");

//			fillTable("matrixN1.txt");
//			fillTable("matrixN2.txt");
//			fillTable("matrixN3.txt");
//			fillTable("matrixN4.txt");

			fillTable("matrix_build37_N1.txt", loaded);
			fillTable("matrix_build37_N2.txt", loaded);
			fillTable("matrix_build37_N3.txt", loaded);
//			fillTable("matrix_build37_N4.txt");

//			fillTable("asmGstart_sept9.txt");
//			fillTable("asmEstart_sept9.txt");
//			fillTable("asmTRstart_sept9.txt");
//			fillTable("asmGstop_sept9.txt");
//			fillTable("asmEstop_sept9.txt");
//			fillTable("asmTRstop_sept9.txt");
//			fillTable("asmEstop_sept16.txt");

//			fillTable("SplicePercentiles_b37_Sept16.txt");
			fillTable("SplicePercentiles_b37_Nov24.txt", loaded);
			table=loaded;

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
	 * Loads matrices from a file into a private, not-yet-published registry.
	 * Reads matrices using MatrixFile and stores original and default-locale lowercase
	 * name keys. Later duplicate names overwrite earlier values. Always attempts reader
	 * closure after successful opening; cleanup failures are suppressed onto a primary
	 * parsing failure, otherwise propagated to prevent publishing incomplete input.
	 * @param fname The filename containing matrix data to load
	 * @param loaded Mutable destination shared by this initialization's file loads
	 */
	private static void fillTable(String fname, HashMap<String, Matrix> loaded){
		final MatrixFile mf=new MatrixFile(fname);
		Throwable failure=null;
		try{
			Matrix mat=mf.nextMatrix();
			while(mat!=null){
//			System.out.println("Adding "+mat.name);
				loaded.put(mat.name, mat);
				loaded.put(mat.name.toLowerCase(), mat);
				mat=mf.nextMatrix();
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			try{
				if(mf.close()){throw new RuntimeException("Error closing matrix input "+fname);}
			}catch(RuntimeException | Error e){
				if(failure==null){throw e;}
				if(failure!=e){failure.addSuppressed(e);}
			}
		}
	}

}
