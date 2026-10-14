package fileIO;

import dna.Matrix;

/**
 * Reads float matrices preceded by name, dimensions, prefix and count headers.
 * Returns the grid, prefix and name in a Matrix; the count is parsed but discarded.
 * Inherits TextFile's mutable reader state and blank-line filtering. Callers own
 * closure: neither EOF nor a //end marker closes the reader.
 * @author Brian Bushnell
 */
public class MatrixFile extends TextFile{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints the input's retained TextFile lines without parsing matrices.
	 * Uses readLine(), bypassing MatrixFile's slash/brace filter. This diagnostic
	 * entry point wraps caught exceptions and does not explicitly close its reader.
	 * @param args Input name as the first argument
	 * @throws RuntimeException if an exception is caught while opening or reading input
	 */
	public static void main(String[] args){
		try{
			String name=args[0];
			MatrixFile mat=new MatrixFile(name);
			String s=null;
			for(s=mat.readLine(); s!=null; s=mat.readLine()){
				System.out.println(s);
			}
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens the input immediately with subprocess permission disabled.
	 * ReadWrite still selects the actual backend, as for TextFile.
	 * @param name Input name interpreted by ReadWrite */
	public MatrixFile(String name){super(name, false);}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Skips logical lines until the first character is slash or opening brace.
	 * TextFile.readLine() removes its defined blank lines first. Retained lines are
	 * not trimmed, so leading whitespace prevents selection. A slash selects any
	 * such line; this method does not validate metadata or recognize end markers.
	 * @return Next selected line, or null at EOF or when the inherited reader is closed
	 */
	@Override
	public String nextLine(){
		String line=readLine();
		while(line!=null && line.charAt(0)!='{' && line.charAt(0)!='/'){
			line=readLine();
		}
		return line;
	}

	/**
	 * Parses four selected header lines followed by the declared number of data rows:
	 * <pre>
	 * //name: example
	 * //size: 2x2
	 * //prefix: 1
	 * //count: 25
	 * {1.0,2.0}
	 * {3.0,4.0}
	 * </pre>
	 * Trims the name; dimensions, prefix and count are parsed without trimming.
	 * Dimensions use a lowercase x separator. The count is discarded; grid, prefix
	 * and name are passed unchanged to the Matrix constructor.
	 * Each data row occupies one selected line. Lines starting // are skipped while
	 * seeking rows, but not between headers. ASCII spaces and braces are removed from
	 * each row before splitting; a comma after the closing brace is also tolerated.
	 * Header prefixes, boundary braces and row widths are checked with assertions;
	 * malformed or truncated input can also produce other unchecked exceptions.
	 * A first selected line beginning //end returns null without latching an end state;
	 * later calls can continue after that marker. Does not close the reader.
	 * @return Parsed matrix, or null when no first line is available or it begins //end
	 * @throws AssertionError if an enabled header, brace or row-width assertion fails
	 * @throws NumberFormatException if an integer header or float value cannot be parsed
	 */
	public Matrix nextMatrix(){
		String line;
		String[] split;
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ", "").trim();
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ", "");
		split=line.split("x");
		int length=Integer.parseInt(split[0]);
		int width=Integer.parseInt(split[1]);
		line=nextLine();
		assert(line.startsWith("//prefix: ")) : line;
		line=line.replace("//prefix: ", "");
		int prefix=Integer.parseInt(line);
		line=nextLine();
		assert(line.startsWith("//count: ")) : line;
		line=line.replace("//count: ", "");
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
