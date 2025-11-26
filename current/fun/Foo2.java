package fun;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.regex.Pattern;

import shared.Parse;

/**
 * Efficient file processing utility for parsing delimited text files with
 * size-based filtering. Processes pipe-delimited text files extracting and
 * aggregating file size information using two parsing strategies for
 * performance comparison.
 *
 * @author Brian Bushnell
 */
public class Foo2 {
	
	/** Compiled regex pattern for splitting on pipe characters */
	static final Pattern pipePattern=Pattern.compile("\\|");

	/**
	 * Main entry point for file processing. Reads a pipe-delimited file and
	 * processes each line to extract file size information for files marked
	 * with 'F' in the seventh field.
	 *
	 * @param args Command-line arguments where args[0] is the input file path
	 * and optional args[1] specifies processing mode ("slow" for
	 * regex-based parsing, anything else for fast manual parsing)
	 * @throws Exception if file I/O operations fail
	 */
	public static void main(String[] args) throws Exception{
		final boolean slow=args.length<2 ? true : "slow".equalsIgnoreCase(args[1]);
		final BufferedReader reader=new BufferedReader(new FileReader(args[0]));
		long sum=0, lines=0, chars=0;

//		ArrayList<Pair> list=new ArrayList<Pair>(10000000);
		for(String line=reader.readLine(); line!=null; line=reader.readLine()) {
			final long size=(slow ? processSlow(line) : processFast(line));
			if(size>=0) {
				sum+=size;
				lines++;
				chars+=line.length();
			}
		}

		reader.close();
		System.out.println("sum="+sum);
		System.out.println("lines="+lines);
		System.out.println("chars="+chars);
	}
	
	/**
	 * Processes a pipe-delimited line using regex pattern splitting.
	 * Extracts file size from the fourth field if the seventh field starts
	 * with 'F'. This method demonstrates regex-based parsing but is slower
	 * than manual parsing for large datasets.
	 *
	 * @param line The pipe-delimited input line to process
	 * @return The file size from field 4 if field 7 starts with 'F',
	 * otherwise -1 to indicate filtering
	 */
	static long processSlow(String line) {
		String[] split=pipePattern.split(line);
		if(split[6].charAt(0)!='F') {return -1;}
		long size=Long.parseLong(split[3]);
		assert(size>=0) : line;
		return size;
	}

	/**
	 * Processes a pipe-delimited line using manual character-by-character
	 * parsing for optimal performance. Manually advances through pipe
	 * delimiters to locate the fourth field (file size) and seventh field
	 * (file type flag), filtering for entries starting with 'F'.
	 *
	 * @param line The pipe-delimited input line to process
	 * @return The file size from field 4 if field 7 starts with 'F',
	 * otherwise -1 to indicate filtering
	 */
	static long processFast(String line) {
		final int delimiter='|';
		final int len=line.length();
		int a=0, b=0;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		long size=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		if(line.charAt(a)!='F') {return -1;}
		b++;
		a=b;

		return size;
	}

}
