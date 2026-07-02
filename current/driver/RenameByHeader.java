package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;

/**
 * Renames files based on their headers by parsing biological sequence files.
 * Extracts organism names from FASTA/FASTQ headers and renames files with
 * taxonomic identifiers. Processes individual files or entire directories
 * containing sequence files.
 *
 * @author Brian Bushnell
 */
public class RenameByHeader {
	
	/** Entry point: parses args, processes files, and closes redirected streams.
	 * @param args Input files or directories to rename (supports fastq/fasta) */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameByHeader x=new RenameByHeader(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses command-line arguments, collects FASTA/FASTQ files (files or directories), and sets verbosity/streams.
	 * @param args Arguments including file paths and optional verbose flag */
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
					//n [RenameByHeader#002] LOW/dev: f.listFiles() can return null on an unreadable dir/I/O error (even past
					//n exists()&&isDirectory()) → NPE in this for-each. Unguarded (same shape as ConcatenateFiles#001). Dead tool → LOW.
					for(File f2 : f.listFiles()){
						String name=f2.getAbsolutePath();
						if(f2.isFile() && FileFormat.hasFastaOrFastqExtension(name)){
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
	
	/** Processes all collected files, applying header-based renaming.
	 * @param t Timer for tracking execution */
	void process(Timer t){
		for(String s : list){
			processFile(s);
		}
	}
	
	/**
	 * Processes a single FASTA/FASTQ file: reads first header, extracts taxonomy tokens, and renames file with genus/species.
	 * Preserves original extension and directory.
	 * @param path Absolute path to the file
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
			//NOTE [driver/RenameByHeader#001] LOW/dev: File.renameTo(...) returns a boolean success flag that is IGNORED here.
			//A failed rename (target already exists, cross-filesystem move, permission denied) is silently swallowed — the
			//file keeps its old name with no error/log. Should check the return. Dead dev tool (no .sh/callers) → LOW.
			//n GOOD (studied praise): the header-token parse above (split[1]/[2]/[3]) is wrapped in try/catch that logs the
			//n offending path and skips — deliberate defense against short/malformed headers, unlike most unguarded split[] tools.
			f.renameTo(new File(sb.toString()));
		}
	}
	
	/*--------------------------------------------------------------*/

	private ArrayList<String> list=new ArrayList<String>();
	private PrintStream outstream=System.err;
	private static boolean verbose=false;
	
}
