package gff;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.PGMTools;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;

/**
 * Single-threaded tool for extracting genomic features from FASTA sequences
 * using GFF annotations with attribute-based filtering and strand-specific processing.
 * Supports both feature extraction and region masking modes.
 * @author Brian Bushnell
 */
public class CutGff_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CutGff_ST x=new CutGff_ST(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CutGff_ST(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=true;
		GffLine.parseAttributes=true;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, append, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff") || a.equals("gffin")){
				assert(b!=null);
				Tools.addFiles(b, gffList);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("invert")){
				invert=Parse.parseBoolean(b);
			}else if(a.equals("type") || a.equals("types")){
				types=b;
			}else if(a.equals("attributes") || a.equals("requiredattributes")){
				requiredAttributes=b.split(",");
			}else if(a.equals("banattributes") || a.equals("bannedattributes")){
				bannedAttributes=b.split(",");
			}else if(a.equals("banpartial")){
				banPartial=Parse.parseBoolean(b);
			}
			
			else if(a.equals("minlen")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxLen=Integer.parseInt(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		ArrayList<String> banned=new ArrayList<String>();
		if(banPartial){banned.add("partial=true");}
		if(bannedAttributes!=null){
			for(String s : bannedAttributes){banned.add(s);}
		}
		bannedAttributes=banned.isEmpty() ? null : banned.toArray(new String[0]);
		
		if(gffList.isEmpty()){
			for(String s : fnaList){
				String prefix=ReadWrite.stripExtension(s);
				String gff=prefix+".gff";
				File f=new File(gff);
				if(!f.exists()){
					String gz=gff+".gz";
					f=new File(gz);
					assert(f.exists() && f.canRead()) : "Can't read file "+gff;
					gff=gz;
				}
				gffList.add(gff);
			}
		}
		assert(gffList.size()==fnaList.size()) : "Number of fna and gff files do not match: "+fnaList.size()+", "+gffList.size();
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(gffList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(out);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Actual Code          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that handles all input file pairs.
	 * Processes each FASTA/GFF file pair sequentially using the configured parameters.
	 * @param t Timer for execution time tracking
	 */
	public void process(Timer t){
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		for(int i=0; i<fnaList.size(); i++){
			processFile(fnaList.get(i), gffList.get(i), types, bsw);
		}
		
		bsw.poisonAndWait();
	}
	
	/**
	 * Processes a single FASTA/GFF file pair for feature extraction.
	 * Loads GFF annotations and FASTA sequences, then processes both strands
	 * for feature matching and extraction based on configured criteria.
	 *
	 * @param fna Path to input FASTA file
	 * @param gff Path to input GFF file
	 * @param types Comma-separated list of feature types to process
	 * @param bsw Output writer for extracted sequences
	 */
	private void processFile(String fna, String gff, String types, ByteStreamWriter bsw){
		ArrayList<GffLine> lines=GffLine.loadGffFile(gff, types, false);
		
		ArrayList<Read> list=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		HashMap<String, Read> map=new HashMap<String, Read>();
		for(Read r : list){map.put(r.id, r);}
		
		processStrand(lines, map, 0, bsw, invert);
		for(Read r : list){r.reverseComplement();}
		processStrand(lines, map, 1, bsw, invert);
		
		if(invert){
			for(Read r : list){
				r.reverseComplement();
				bsw.println(r);
			}
		}
	}
	
	/**
	 * Determines if a GFF line passes all configured attribute and length filters.
	 * Checks feature length constraints, banned attributes, and required attributes.
	 * @param gline GFF line to evaluate
	 * @return true if the line passes all filters, false otherwise
	 */
	private boolean hasAttributes(GffLine gline){
		int len=gline.length();
		if(len<minLen || len>maxLen){return false;}
		if(hasAttributes(gline, bannedAttributes)){return false;}
		return requiredAttributes==null || hasAttributes(gline, requiredAttributes);
	}
	
	/**
	 * Checks if a GFF line contains any of the specified attributes.
	 * Used for both required and banned attribute filtering.
	 *
	 * @param gline GFF line to check
	 * @param attributes Array of attributes to search for
	 * @return true if any attribute is found in the line, false otherwise
	 */
	private boolean hasAttributes(GffLine gline, String[] attributes){
		if(attributes==null){return false;}
		for(String s : attributes){
			if(gline.attributes.contains(s)){
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Processes GFF features for a specific strand orientation.
	 * Extracts or masks genomic regions based on strand-specific coordinates
	 * and configured filtering criteria.
	 *
	 * @param lines List of GFF lines to process
	 * @param map Mapping from sequence IDs to Read objects
	 * @param strand Strand orientation (0=forward, 1=reverse)
	 * @param bsw Output writer for extracted sequences
	 * @param invert If true, mask regions with 'N' instead of extracting
	 */
	private void processStrand(ArrayList<GffLine> lines, HashMap<String, Read> map, int strand, ByteStreamWriter bsw, boolean invert){
		for(GffLine gline : lines){
			if(gline.strand==strand && hasAttributes(gline)){
				Read scaf=map.get(gline.seqid);
				assert(scaf!=null) : "Can't find "+gline.seqid+" in "+map.keySet();
				int start, stop;
				if(strand==0){
					start=gline.start-1;
					stop=gline.stop-1;
				}else{
					start=scaf.length()-gline.stop-1;
					stop=scaf.length()-gline.start-1;
				}
				if(invert){
					byte[] bases=scaf.bases;
					for(int i=start; i<stop; i++){
						if(i>=0 && i<bases.length){
							bases[i]='N';
						}
					}
				}else{
					if(start>=0 && stop<scaf.length()){
						Read r=new Read(Arrays.copyOfRange(scaf.bases, start, stop), null, gline.attributes, 1);
						bsw.println(r);
					}
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** List of input FASTA file paths */
	private ArrayList<String> fnaList=new ArrayList<String>();
	/** List of input GFF file paths corresponding to FASTA files */
	private ArrayList<String> gffList=new ArrayList<String>();
	/** Output file path for extracted sequences */
	private String out=null;
	/** Comma-separated list of GFF feature types to process (default: "CDS") */
	private String types="CDS";
	/** If true, mask matching regions with 'N' instead of extracting them */
	private boolean invert=false;
	/** If true, exclude features marked as partial=true */
	private boolean banPartial=true;
	/** Minimum feature length to include in processing */
	private int minLen=1;
	/** Maximum feature length to include in processing */
	private int maxLen=Integer.MAX_VALUE;

	/** Array of attributes that must be present in GFF features for inclusion */
	private String[] requiredAttributes;
	/** Array of attributes that cause GFF features to be excluded */
	private String[] bannedAttributes;
	
	/*--------------------------------------------------------------*/
	
	/** Number of bytes written to output (currently unused) */
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file format configuration */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Stream for status and error messages */
	private PrintStream outstream=System.err;
	/** Enable verbose logging output */
	public static boolean verbose=false;
	/** Indicates if an error occurred during processing */
	public boolean errorState=false;
	/** Allow overwriting existing output files */
	private boolean overwrite=true;
	/** Append output to existing files instead of overwriting */
	private boolean append=false;
	
}
