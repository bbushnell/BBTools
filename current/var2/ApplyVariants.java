package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Applies variants
 * 
 * @author Brian Bushnell
 * @date August 27, 2019
 *
 */
public class ApplyVariants {
	
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
		ApplyVariants x=new ApplyVariants(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public ApplyVariants(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.capBuffers(4); //Only for singlethreaded programs
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffvcf=FileFormat.testInput(inVcf, FileFormat.VCF, null, true, true);
		ffdepth=FileFormat.testInput(inDepth, FileFormat.TXT, null, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("vcf") || a.equals("vars")){
				inVcf=b;
			}else if(a.equals("name") || a.equals("sample") || a.equals("samplename")){
				sampleName=b;
			}else if(a.equals("addnumbers")){
				addContigNumbers=Parse.parseBoolean(b);
			}else if(a.equals("useprefix") || a.equals("prefix")){
				usePrefix=Parse.parseBoolean(b);
			}else if(a.equals("delimiter")){
				delimiter=Parse.parseSymbolToCharacter(b);
			}else if(a.equals("cov") || a.equals("depth") || a.equals("indepth") || a.equals("basecov")){
				inDepth=b;
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepth=Integer.parseInt(b);
			}else if(a.equals("maxindel")){
				maxIndel=Integer.parseInt(b);
				if(maxIndel<0){maxIndel=Integer.MAX_VALUE;}
			}else if(a.equals("noindels")){
				noIndels=Parse.parseBoolean(b);
			}else if(a.equals("noframeshifts") || a.equals("banframeshifts")){
				noFrameshifts=Parse.parseBoolean(b);
			}else if(a.equals("frameshifts")){
				noFrameshifts=!Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		
		//Ensure there is an input file
		if(in1==null || inVcf==null){throw new RuntimeException("Error - one sequence and one vcf file are required.");}
	
		if(minDepth>0 && inDepth==null){throw new RuntimeException("Error - mindepth requires a coverage file.");}
		
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		inVcf=Tools.fixExtension(inVcf);
		inDepth=Tools.fixExtension(inDepth);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, inVcf, inDepth)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, inVcf, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		FASTQ.FORCE_INTERLEAVED=false;
		FASTQ.TEST_INTERLEAVED=false;
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		VCFFile vfile=new VCFFile(ffvcf);
		
		//Load vcf
		if(ScafMap.defaultScafMap()==null){
			ScafMap.setDefaultScafMap(vfile.toScafMap(null), ffvcf.name());
		}
		
		if(ffdepth!=null && minDepth>0){
			depthMap=CoverageArray.loadDepth(ffdepth, CoverageArray2.class);
		}
		
		ArrayList<VCFLine> lines=vfile.lines(true);
		varMap=new HashMap<String, ArrayList<Var>>(ScafMap.defaultScafMap().size());
		for(VCFLine line : lines){
			ArrayList<Var> value=varMap.get(line.scaf);
			if(value==null){
				value=new ArrayList<Var>();
				varMap.put(line.scaf, value);
			}
			value.add(line.toVar());
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=4;
		
		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln, cris, ros);

				//Fetch a new list
				ln=cris.nextList();
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Process a list of Reads.
	 * @param ln The list.
	 * @param cris Read Input Stream
	 * @param ros Read Output Stream for reads that will be retained
	 */
	void processList(ListNum<Read> ln, final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){

		//Grab the actual read list from the ListNum
		final ArrayList<Read> reads=ln.list;
		
		//Loop through each read in the list
		for(int idx=0; idx<reads.size(); idx++){
			final Read r=reads.get(idx);
			
			//Validate reads in worker threads
			if(!r.validated()){r.validate(true);}

			//Track the initial length for statistics
			final int initialLength1=r.length();

			//Increment counters
			readsProcessed+=r.pairCount();
			basesProcessed+=initialLength1;

			Read mutant=processRead(r);
			reads.set(idx, mutant);
			readsOut+=mutant.pairCount();
			basesOut+=mutant.pairLength();
		}

		//Output reads to the output stream
		if(ros!=null){ros.add(reads, ln.id);}

		//Notify the input stream that the list was used
		cris.returnList(ln);
//		if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
	}
	
	@SuppressWarnings("unused")
	private void applyDepth(Read r){applyDepth(r, null);}
	
	private void applyDepth(Read r, CoverageArray ca){
		if(ca==null){ca=depthMap.get(r.id);}
		if(ca==null){
			String s=Tools.trimToWhitespace(r.id);
			ca=depthMap.get(s);
		}
		assert(ca!=null) : "Can't find "+r.id+" in depth map.";
		final byte[] bases=r.bases;
		for(int i=0; i<r.bases.length; i++){
			if(ca.get(i)<minDepth){bases[i]=noCovSymbol;}
		}
		
		ArrayList<Var> vars=varMap.get(r.id);
		if(vars==null){
			String sub=Tools.trimToWhitespace(r.id);
			vars=varMap.get(sub);//Handles truncated sequence names
		}
		if(vars==null){return;}
		int removed=0;
		for(int i=0; i<vars.size(); i++){
			Var v=vars.get(i);
			if(ca.get(v.start)<minDepth && v.indel()){
				vars.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(vars);}
	}
	
	private void filterIndels(final Read r){
		filterIndels(r, null);
	}
	
	private void filterIndels(final Read r, CoverageArray ca){
		if(minDepth<=0 && maxIndel==Integer.MAX_VALUE && !noIndels && !noFrameshifts){return;}
		
		if(minDepth>0 && ca==null){
			ca=depthMap.get(r.id);
			if(ca==null){
				String s=Tools.trimToWhitespace(r.id);
				ca=depthMap.get(s);
			}
			assert(ca!=null) : "Can't find "+r.id+" in depth map.";
		}
		
		ArrayList<Var> vars=varMap.get(r.id);
		if(vars==null){
			String sub=Tools.trimToWhitespace(r.id);
			vars=varMap.get(sub);//Handles truncated sequence names
		}
		if(vars==null){return;}
		int removed=0;
		for(int i=0; i<vars.size(); i++){
			Var v=vars.get(i);
			boolean remove=false;
			if(noFrameshifts && v.frameshift()){
				remove=true;
			}else if(v.indel()){
				if(noIndels || Tools.max(v.reflen(), v.readlen())>maxIndel 
						|| (ca!=null && ca.get(v.start)<minDepth))
				remove=true;
			}
			if(remove){
				vars.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(vars);}
	}
	
	/**
	 * Process a single read.
	 * @param r Read 1
	 */
	Read processRead(final Read r){
		
		CoverageArray ca=null;
		if(minDepth>0){
			assert(depthMap!=null) : "minDepth is "+minDepth+" but depthMap is null. "
					+ "You need a coverage file which includes this contig.";
			ca=depthMap.get(r.id);
			if(ca==null){
				String s=Tools.trimToWhitespace(r.id);
				ca=depthMap.get(s);
			}
			assert(ca!=null) : "Can't find "+r.id+" in depth map.";
			applyDepth(r, ca);
		}
		filterIndels(r, ca);
		
		String name=r.id;
		final byte[] bases=r.bases;
		ArrayList<Var> vars=varMap.get(name);
		if(vars==null){return r;}
		HashMap<Integer, Var> map=new HashMap<Integer, Var>();
		for(Var v : vars){
			Integer key=v.start;
			Var old=map.get(key);
			if(old==null || old.alleleCount()<v.alleleCount()){//TODO: allow retaining the higher quality one
				map.put(key, v);
			}
		}
		
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<bases.length; ){
			Var v=map.get(i);
			if(v==null){
				bb.append(bases[i]);
				i++;
			}else{
				
				final int reflen=v.reflen();
				final int readlen=v.readlen();
				
				applied++;
				
				//Old code, did not handle 'N' for low coverage
//				if(readlen>0){
//					bb.append(v.allele);
//				}
//				if(reflen==0){//insertion
//					bb.append(bases[i]);
//					i++;
//				}else{
//					i+=reflen;
//				}

				//low-depth indels should have already been removed
				if(reflen==0){//insertion
//					final int depth2=(ca==null || minDepth<1 ? 9999 : ca.get(i));
//					if(depth2>=minDepth){//Apply ins
//						bb.append(v.allele);
//					}
					bb.append(bases[i]);
					i++;
				}else if(readlen==0){//del
//					final int depth2=(ca==null || minDepth<1 ? 9999 : ca.get(i));
//					if(depth2>=minDepth){//Apply del
						i+=reflen;
//					}else{//Ignore del
//						bb.append(bases[i]);
//						i++;
//					}
				}else if(readlen!=reflen){//complex
					int lim=i+reflen;
					for(int j=0; i<lim; i++, j++){
						final int depth2=(ca==null || minDepth<1 ? 9999 : ca.get(i));
						byte b=(j<v.allele.length ? v.allele[j] : bases[i]);
						bb.append(depth2>=minDepth ? b : (byte)'N');
					}
				}else{//sub
					assert(v.reflen()==v.readlen()) : v;
					for(byte b : v.allele){
						final int depth2=(ca==null || minDepth<1 ? 9999 : ca.get(i));
						bb.append(depth2>=minDepth ? b : (byte)'N');
						i++;
					}
				}
			}
		}
		
		final String name2=rename(name, r.numericID);
		
		Read mutant=new Read(bb.toBytes(), null, name2, r.numericID);
		return mutant;
	}
	
	private String rename(String old, long number){
		if(sampleName==null){return old;}
		nameBuilder.clear();
		if(usePrefix){
			return nameBuilder.append(sampleName).append(delimiter).append(old).toString();
		}else if(addContigNumbers){
			return nameBuilder.append(sampleName).append(delimiter).append(number).toString();
		}
		return sampleName;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Variant input file */
	private String inVcf;
	/** Per-base coverage depth file (optional) */
	private String inDepth;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** If positive, change regions below this depth to N */
	private int minDepth=0;
	
	private byte noCovSymbol='N';
	
	HashMap<String, ArrayList<Var>> varMap;
	
	HashMap<String, CoverageArray> depthMap;
	
	/** Name of output sequences, if different from input */
	String sampleName=null;
	/** Add numbers as a suffix to the name of each contig, if they will be renamed */
	boolean addContigNumbers=true;
	/** Use samplename as prefix, rather than replacing the existing name. */
	boolean usePrefix=false;
	/** Delimits modified names. */
	char delimiter='_';
	
	/** Ignore non-multiple-of-3 indels */
	private boolean noFrameshifts=false;
	/** Ignore indels longer than this */
	private int maxIndel=Integer.MAX_VALUE;
	/** Ignore indels */
	private boolean noIndels=false;
	
	private ByteBuilder nameBuilder=new ByteBuilder();
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private long applied=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Variant input file */
	private final FileFormat ffvcf;
	/** Coverage input file */
	private final FileFormat ffdepth;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
