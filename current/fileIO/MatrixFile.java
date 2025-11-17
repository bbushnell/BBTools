package fileIO;
import dna.Matrix;



/**
 * File reader for matrix data files in BBTools format.
 * Parses text files containing formatted matrix data with metadata headers.
 * Extends TextFile to provide specialized matrix parsing capabilities.
 * @author Brian Bushnell
 */
public class MatrixFile extends TextFile{
	
	/**
	 * Program entry point for testing matrix file reading.
	 * Reads and prints all lines from the specified matrix file.
	 * @param args Command-line arguments; first argument should be matrix file name
	 * @throws RuntimeException If file cannot be read or other error occurs
	 */
	public static void main(String[] args){
		
		try {
			//Name of mat file
			String name=args[0];
			
			MatrixFile mat=new MatrixFile(name);
			
			String s=null;
			
			for(s=mat.readLine(); s!=null; s=mat.readLine()){
				System.out.println(s);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}
	
	
	/** Constructs a MatrixFile reader for the specified file.
	 * @param name Path to the matrix file to read */
	public MatrixFile(String name){super(name, false);}
	
	@Override
	public String nextLine(){
		String line=readLine();
		
		while(line!=null && line.charAt(0)!='{' && line.charAt(0)!='/'){
			line=readLine();
		}
		return line;
	}
	
	/**
	 * Parses the next complete matrix from the file.
	 * Reads metadata headers (name, size, prefix, count) followed by matrix data rows.
	 * Matrix data is expected in bracket-enclosed, comma-separated format.
	 *
	 * @return Parsed Matrix object or null if end of file or "//end" marker reached
	 * @throws AssertionError If matrix format is invalid or expected headers missing
	 * @throws NumberFormatException If numeric values cannot be parsed
	 */
	public Matrix nextMatrix(){
		String line;
		String[] split;
		
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ","").trim();
		
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ","");
		split=line.split("x");
		int length=Integer.parseInt(split[0]);
		int width=Integer.parseInt(split[1]);
		
		line=nextLine();
		assert(line.startsWith("//prefix: ")) : line;
		line=line.replace("//prefix: ","");
		int prefix=Integer.parseInt(line);
		
		line=nextLine();
		assert(line.startsWith("//count: ")) : line;
		line=line.replace("//count: ","");
		int count=Integer.parseInt(line);
		
		
		float[][] grid=new float[length][width];
		for(int i=0; i<length; i++){
			line=nextLine();
			
			while(line.startsWith("//")){line=nextLine();}
			
			assert(line.startsWith("{"));
			if(line.endsWith(",")){line=line.substring(0, line.length()-1);}
			assert(line.endsWith("}"));
			line=line.replace("{", "").replace("}", "").replace(" ", "");
			split=line.split(",");
			assert(split.length==width);
			for(int j=0; j<split.length; j++){
				grid[i][j]=Float.parseFloat(split[j]);
			}
		}
		
		return new Matrix(grid, prefix, name);
	}
	
}
