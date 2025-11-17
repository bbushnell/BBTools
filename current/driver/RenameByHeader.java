package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;

/**
 * Renames files based on their headers
 * @author Brian Bushnell
 * @date May 19, 2016
 *
 */
public class RenameByHeader {
	
	/** Program entry point for file renaming utility.
	 * @param args Command-line arguments specifying input files or directories */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameByHeader x=new RenameByHeader(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs RenameByHeader instance and parses command-line arguments.
	 * Processes input arguments to identify files and directories containing
	 * FASTA/FASTQ files. Supports verbose mode and directory traversal.
	 * @param args Command-line arguments including file paths and options
	 */
	public RenameByHeader(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_UNBGZIP=false;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			File f=(b==null ? new File(arg) : null);

			if(f!=null && f.exists()){
				if(f.isDirectory()){
					for(File f2 : f.listFiles()){
						String name=f2.getAbsolutePath();
						if(f2.isFile() && FileFormat.hasFastqOrFastqExtension(name)){
							list.add(name);
						}
					}
				}else{
					list.add(f.getAbsolutePath());
				}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
	}
	
	/**
	 * Processes all files in the input list for renaming.
	 * Iterates through collected file paths and applies header-based renaming.
	 * @param t Timer for tracking execution duration
	 */
	void process(Timer t){
		for(String s : list){
			processFile(s);
		}
	}
	
	/**
	 * Processes a single file for header-based renaming.
	 * Reads the first header line, extracts taxonomic information, and renames
	 * the file using genus and species identifiers. Handles species abbreviations
	 * and constructs new filename with original extension preserved.
	 *
	 * @param path Absolute path to the file to be renamed
	 */
	void processFile(String path){
		TextFile tf=new TextFile(path);
		String line=tf.nextLine();
		tf.close();
		if(line==null){return;}
		
		StringBuilder sb=new StringBuilder();
		File f=new File(path);
		String dir=f.getParent();
		if(dir!=null){sb.append(dir).append('/');}
		try {
			String[] split=line.substring(1).replace(",", "").split(" ");
			sb.append(split[1]);
			sb.append('_');
			sb.append(split[2]);
			sb.append('_');
			if(split[2].equals("sp.")){
				sb.append(split[3]);
				sb.append('_');
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.err.println(path);
			e.printStackTrace();
			return;
		}
		if(sb.length()>0){
			String name=f.getName();
			sb.append(name);
			f.renameTo(new File(sb.toString()));
		}
	}
	
	/*--------------------------------------------------------------*/

	/** List of file paths to process for renaming */
	private ArrayList<String> list=new ArrayList<String>();
	/** Output stream for messages and error reporting */
	private PrintStream outstream=System.err;
	/** Controls verbose output during file processing */
	private static boolean verbose=false;
	
}
