package fileIO;


/**
 * Specialized file reader for parsing array data files with metadata extraction.
 * Extends TextFile to handle files containing named float arrays.
 * Expects array files with metadata headers followed by brace-enclosed float values.
 * @author Brian Bushnell
 */
public class ArrayFile extends TextFile{
	
	/**
	 * Program entry point that demonstrates basic ArrayFile usage.
	 * Reads an array file specified as the first command-line argument and prints all lines to stdout.
	 * @param args Command-line arguments; expects array file path as first argument
	 */
	public static void main(String[] args){
		
		try {
			//Name of mat file
			String name=args[0];
			
			ArrayFile mat=new ArrayFile(name);
			
			String s=null;
			
			for(s=mat.readLine(); s!=null; s=mat.readLine()){
				System.out.println(s);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}
	
	
	/** Constructs an ArrayFile reader for the specified file.
	 * @param name Path to the array file to read */
	public ArrayFile(String name){super(name, false);}
	
	@Override
	public String nextLine(){
		String line=readLine();
		char c=line.charAt(0);
		
		while(line!=null && c!='{' && c!='/'){
			line=readLine();
			c=line.charAt(0);
		}
		return line;
	}
	
	/**
	 * Parses and returns the next float array from the file with its metadata.
	 * Expects a specific format with metadata headers followed by comma-separated values in braces.
	 * Format expected:
	 * //name: [array_name]
	 * //size: [array_length]
	 * {value1,value2,value3,...}
	 *
	 * @return Float array containing parsed values, null if end marker or EOF reached
	 */
	public float[] nextArray(){
		String line;
		String[] split;
		
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ","").trim();
		
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ","");
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
