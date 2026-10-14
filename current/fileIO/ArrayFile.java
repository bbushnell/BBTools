package fileIO;

/**
 * Reads single-line float arrays preceded by name and size headers.
 * Array names are parsed but discarded; only the float values are returned.
 * Inherits TextFile's mutable reader state and blank-line filtering. Callers own
 * closure: neither EOF nor a //end marker closes the reader.
 * @author Brian Bushnell
 */
public class ArrayFile extends TextFile{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints the input's retained TextFile lines without parsing arrays.
	 * Uses readLine(), so ArrayFile's slash/brace filter is bypassed. This diagnostic
	 * entry point wraps caught exceptions and does not explicitly close its reader.
	 * @param args Input name as the first argument
	 */
	public static void main(String[] args){
		try{
			String name=args[0];
			ArrayFile mat=new ArrayFile(name);
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
	public ArrayFile(String name){super(name, false);}

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
	 * Parses one array from three successive selected lines, for example:
	 * <pre>
	 * //name: example
	 * //size: 2
	 * {1.0,2.0}
	 * </pre>
	 * Discards the name and allocates the declared number of floats. Header prefixes,
	 * boundary braces and value count are checked with assertions. The size is not
	 * trimmed; ASCII spaces and braces are removed from the value line before splitting.
	 * A comma after the closing brace is also tolerated. Values use Float.parseFloat.
	 * A line beginning //end returns null for this call without latching an end state;
	 * a later call can continue after that marker. Does not close the reader.
	 * @return Parsed float values, or null when the first selected line is EOF or //end
	 */
	public float[] nextArray(){
		String line;
		String[] split;
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ", "").trim();
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ", "");
		int length=Integer.parseInt(line);
		float[] grid=new float[length];
		line=nextLine();
		assert(line.startsWith("{"));
		if(line.endsWith(",")){line=line.substring(0, line.length()-1);}
		assert(line.endsWith("}"));
		line=line.replace("{", "").replace("}", "").replace(" ", "");
		split=line.split(",");
		assert(split.length==length);
		for(int i=0; i<split.length; i++){
			grid[i]=Float.parseFloat(split[i]);
		}
		return grid;
	}
}
