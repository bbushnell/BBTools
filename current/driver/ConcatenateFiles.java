package driver;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Timer;

/**
 * Utility for concatenating multiple files into a single output file.
 * Supports both directory-based concatenation and pattern-based sequence file merging.
 * Provides binary file copying for general use and specialized FASTA processing.
 * @author Brian Bushnell
 */
public class ConcatenateFiles {
	
	/**
	 * Program entry point for file concatenation operations.
	 * Routes to directory concatenation if input path is a directory,
	 * otherwise performs pattern-based concatenation for sequence files.
	 *
	 * @param args Command-line arguments where args[0] is input path/pattern
	 * and args[1] is optional output path
	 */
	public static void main(String[] args){
		
		Timer t=new Timer();
		final String in=args[0];
		final String out=(args.length>1 ? args[1] : null);
		if(new File(in).isDirectory()){
			try {
				concatenateDirectory(in, out);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else{
			concatenatePattern(in, out);
		}
		t.stop();
		System.err.println(t);
		
	}
	
	/**
	 * Concatenates chromosome-specific sequence files using a numbered pattern.
	 * Replaces '#' placeholder in basename with chromosome numbers 1-25.
	 * Outputs FASTA format with chromosome headers, filtering out existing headers.
	 * Note: Contains human-specific logic as indicated by assertion.
	 *
	 * @param basename File pattern with '#' placeholder for chromosome number
	 * @param out Output filename, or null to auto-generate by replacing '#' with 'ALL'
	 */
	public static void concatenatePattern(final String basename, final String out){
		assert(false) : "This is human-specific.";
		String outname=(out==null ? basename.replace("#", "ALL") : out);
		
		TextStreamWriter tsw=new TextStreamWriter(outname, true, false, true);
		tsw.start();
		
		for(int chrom=1; chrom<26; chrom++){
			String fname=basename.replace("#", ""+chrom);
			TextFile tf=new TextFile(fname, false);
			
			tsw.print(">chr"+chrom+"\n");
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				char c=s.charAt(0);
				if(c!='>' && c!=';'){
					tsw.println(s);
				}
			}
			System.err.print(".");
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Concatenates all files in a directory using binary copying.
	 * Reads directory contents, sorts files alphabetically, and copies each file's
	 * content to the output stream using a 32KB buffer for efficiency.
	 * Skips the output file itself if it exists in the same directory.
	 *
	 * @param in Input directory path containing files to concatenate
	 * @param out Output file path, or null/"stdout" for standard output
	 * @throws IOException If directory cannot be read or file I/O operations fail
	 */
	public static void concatenateDirectory(final String in, String out) throws IOException{
		if(out==null){out="stdout";}
		
		final byte[] buf=new byte[32768];
		
		final File dir=new File(in);
		final File[] files=dir.listFiles();
		Arrays.sort(files);
		
		final File outfile=new File(out);
		final OutputStream os=ReadWrite.getOutputStream(out, false, true, true);
		
		for(File f : files){
			if(f!=null && f.isFile() && !f.equals(outfile)){
				String fname=f.getAbsolutePath();
				System.err.println("Processing "+fname);
				
				InputStream is=ReadWrite.getInputStream(fname, false, false, true);
				
				for(int lim=is.read(buf); lim>0; lim=is.read(buf)){
					os.write(buf, 0, lim);
				}
				
				is.close();
				System.err.print(".");
			}
		}
		ReadWrite.close(os);
	}
	
	
}
