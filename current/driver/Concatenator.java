package driver;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * Utility class for concatenating text files.
 * Provides methods to merge multiple text files into a single output stream
 * or combine their contents into a StringBuilder.
 * @author Brian Bushnell
 */
public class Concatenator {
	
	
	/**
	 * Program entry point for concatenating comma-separated list of files.
	 * Takes two arguments: comma-separated input file names and output file path.
	 * Writes all input files sequentially to the specified output file.
	 * @param args Command-line arguments: [input_files_comma_separated, output_file]
	 */
	public static void main(String args[]){
		
		assert(args.length==2 && !args[1].contains(","));
		TextStreamWriter tsw=new TextStreamWriter(args[1], false, false, true);
		tsw.start();
		for(String s : args[0].split(",")){
			writeFile(s, tsw);
		}
		tsw.poison();
	}
	
	/**
	 * Writes the contents of a text file to output stream or standard output.
	 * Reads the specified file line by line and outputs each line.
	 * If writer is null, prints to standard output; otherwise writes to the stream.
	 *
	 * @param fname Path to the input file to read
	 * @param tsw Output stream writer, or null to print to standard output
	 */
	public static void writeFile(String fname, TextStreamWriter tsw){
		TextFile tf=new TextFile(fname, false);
		if(tsw==null){
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				System.out.println(s);
			}
		}else{
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				tsw.println(s);
			}
		}
		tf.close();
	}
	
	
	/**
	 * Merges multiple text files into a single StringBuilder.
	 * Reads each file completely into memory and appends all lines
	 * to the StringBuilder with newline separators.
	 *
	 * @param fnames Array of file paths to merge
	 * @return StringBuilder containing all file contents concatenated
	 */
	public static StringBuilder merge(String[] fnames){
		StringBuilder sb=new StringBuilder();
		
		for(int i=0; i<fnames.length; i++){
			String fname=fnames[i];
			if(fname!=null){
				TextFile tf=new TextFile(fname, false);
				String[] lines=tf.toStringLines();
				tf.close();
				for(int j=0; j<lines.length; j++){
					String s=lines[j];
					lines[j]=null;
//					if(i<2 || !s.startsWith("#")){
//						sb.append(s);
//						sb.append('\n');
//					}
					sb.append(s);
					sb.append('\n');
				}
			}
		}
		return sb;
	}
	
	
}
