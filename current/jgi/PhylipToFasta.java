package jgi;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 3, 2014
 *
 */
public class PhylipToFasta {
	


	/** Program entry point for Phylip to FASTA conversion.
	 * @param args Command-line arguments specifying input and output files */
	public static void main(String[] args){
		Timer t=new Timer();
		PhylipToFasta x=new PhylipToFasta(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs PhylipToFasta converter and parses command-line arguments.
	 * Sets up input/output file formats and validates file paths.
	 * @param args Command-line arguments for input file, output file, and options
	 */
	public PhylipToFasta(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.PHYLIP, ".phylip", true, true);
	}
	
	/**
	 * Processes Phylip input file and converts to FASTA format.
	 * Reads sequence names from first occurrence of each sequence,
	 * then concatenates sequence data from subsequent interleaved blocks.
	 * Skips the header line containing sequence count and length.
	 *
	 * @param t Timer for tracking processing time and statistics
	 */
	void process(Timer t){
		
		ArrayList<StringBuilder> data=new ArrayList<StringBuilder>();
		long bases=0;
		
		{
			final TextFile tf=new TextFile(ffin1);
			String s=tf.nextLine(); //first line is some numbers

			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.startsWith("        ")){break;}
				StringBuilder sb=new StringBuilder();
				data.add(sb);
				sb.append('>');
				int pos=0;
				for(pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isWhitespace(c)){break;}
					sb.append(c);
				}
				sb.append('\n');
				while(pos<s.length() && Character.isWhitespace(s.charAt(pos))){pos++;}
				while(pos<s.length()){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
					pos++;
				}
			}

			final int mod=data.size();
			for(int i=0; s!=null; i++){
				StringBuilder sb=data.get(i%mod);
				for(int pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
					pos++;
				}
				s=tf.nextLine();
			}
			errorState|=tf.errorState;
		}
		final long reads=data.size();
		
		if(ffout1!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffout1);
			tsw.start();
			for(int i=0; i<data.size(); i++){
				StringBuilder sb=data.set(i, null);
				sb.append('\n');
				tsw.print(sb);
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;
	
	/*--------------------------------------------------------------*/
	
	/** Input file format handler for Phylip files */
	private final FileFormat ffin1;
	/** Output file format handler for FASTA files */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics */
	private PrintStream outstream=System.err;
	/** Controls verbosity of status output during processing */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
