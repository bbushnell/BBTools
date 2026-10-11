package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import shared.Random;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;


/**
 * Generates chimeric reads by concatenating copied fragments of retained input reads.
 * Loads nonempty, unpaired reads into memory, then selects two parents with
 * replacement for each candidate; both parents may be the same read. Fragments
 * and the completed chimera may be reverse complemented independently.
 * Qualities are retained in the generated Read only when both fragments have them.
 * <p>Construction changes shared I/O settings and resolves formats. Processing
 * finishes input loading before opening the optional output writer. Generation
 * and counters are owned by the processing thread; output writing is asynchronous.
 *
 * @author Brian Bushnell
 * @date Oct 7, 2014
 */
public class MakeChimeras {

	/** Runs generation and closes a redirected diagnostic stream on normal completion. */
	public static void main(String[] args){
		Timer t=new Timer();
		MakeChimeras x=new MakeChimeras(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses options and resolves input/output formats without opening the output writer.
	 * Applies parsed overwrite/append flags. Interleaved input detection is initially
	 * disabled; processing asserts that source reads are unpaired.
	 * @param args Command-line options, including optional configuration files
	 */
	public MakeChimeras(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
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
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(a.equals("forcelength")){
				forceLength=Integer.parseInt(b);
			}else if(a.equals("readsout") || a.equals("chimeras")){
				readsOut=Parse.parseKMG(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			overwrite=parser.overwrite;
			append=parser.append;
			
			readsIn=parser.maxReads;
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	/**
	 * Loads input, retaining nonempty reads, then generates the requested chimera count.
	 * Input counters include reads discarded for zero length. Generation retries null
	 * candidates, and output counters advance only when a writer is present. A missing
	 * output therefore still runs generation but reports zero emitted reads and bases.
	 * <p>Requires a positive requested count and usable source reads. A positive forced
	 * length requires a longer source read, checked before opening the writer.
	 * Normal completion returns input batches,
	 * closes the reader, finishes the writer and reports accumulated I/O errors;
	 * exceptional exits have no general finally cleanup here.
	 * @param t Started timer; restarted to report loading and generation separately
	 */
	void process(Timer t){
		assert(readsOut>0) : "Please set the 'readsout' flag to a positive integer.";
		
		ArrayList<Read> source=new ArrayList<Read>();
		{
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(readsIn, false, ffin1, null, qfin1, null);
				if(verbose){outstream.println("Started cris");}
				cris.start(); //4567
			}
			assert(!cris.paired());
			
			long readsProcessed=0;
			long basesProcessed=0;
			
			{
				
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				//			outstream.println("Fetched "+reads);
				
				if(reads!=null && !reads.isEmpty()){
					Read r=reads.get(0);
					assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
				}
				
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

					for(int idx=0; idx<reads.size(); idx++){
						final Read r1=reads.get(idx);
						assert(r1.mate==null);

						final int initialLength1=r1.length();
						
						if(initialLength1>0){
							source.add(r1);
						}

						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				if(ln!=null){
					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				}
			}
			
			errorState|=ReadWrite.closeStream(cris);
			
			t.stop();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		
		
		if(readsOut>=0){
			t.start();
			//A positive second fragment requires a source longer than forceLength.
			if(readsOut>0 && forceLength>0){
				int maxLength=0;
				for(Read r : source){
					maxLength=Tools.max(maxLength, r.length());
					if(maxLength>forceLength){break;}
				}
				if(maxLength<=forceLength){
					throw new IllegalArgumentException("No source read is longer than forcelength="+
							forceLength+"; maximum source length="+maxLength);
				}
			}
			
			final TextStreamWriter tsw;
			if(ffout1==null){
				tsw=null;
			}else{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
			}
			
			final Random randy=Shared.threadLocalRandom();
			
			long readsProcessed=0;
			long basesProcessed=0;
			final int mod=source.size();
			for(long i=0; i<readsOut; i++){
				Read a=source.get(randy.nextInt(mod));
				Read b=source.get(randy.nextInt(mod));
				Read c=makeChimera(a, b, randy, i);
				if(c==null){
					i--;
				}else{
					if(tsw!=null && c!=null){
						tsw.println(c);
						readsProcessed++;
						basesProcessed+=c.length();
					}
				}
			}
			
			if(tsw!=null){errorState|=tsw.poisonAndWait();}
			
			t.stop();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Concatenates copied fragments, preserving the parent reads.
	 * With positive forceLength, requests that length from a (clipped to a's length)
	 * and b.length()-forceLength from b. Otherwise draws each fragment's length.
	 * The ID joins parent IDs with " ~ "; the result may be reverse complemented.
	 * @param a First parent, which may be the same object as b
	 * @param b Second parent
	 * @param randy Generator used for fragment selection and orientation
	 * @param numericID Numeric ID for the generated read
	 * @return New chimera, or null when a forced fragment has no bases
	 */
	private Read makeChimera(Read a, Read b, Random randy, long numericID) {
		final String id=a.id+" ~ "+b.id;

		final Read a2, b2;
		if(forceLength>0){
//			if(a.length()>b.length()){
//				Read c=a;
//				a=b;
//				b=c;
//			}
			a2=getPiece(a, randy, forceLength);
			b2=getPiece(b, randy, b.length()-forceLength);
			if(a2==null || b2==null){return null;}
		}else{
			a2=getPiece(a, randy);
			b2=getPiece(b, randy);
		}
		
		a=b=null;
		
		final byte[] abases=a2.bases, bbases=b2.bases, aquals=a2.quality, bquals=b2.quality;
		final int alen=a2.length(), blen=b2.length();
		final int len=a2.length()+b2.length();
		byte[] bases=new byte[len];
		byte[] quals=(aquals==null || bquals==null) ? null : new byte[len];

		for(int i=0; i<alen; i++){
			bases[i]=abases[i];
			if(quals!=null){quals[i]=aquals[i];}
		}
		for(int i=alen, j=0; j<blen; i++, j++){
			bases[i]=bbases[j];
			if(quals!=null){quals[i]=bquals[j];}
		}
		
		Read r=new Read(bases, quals, id, numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * Copies a positive-length fragment no longer than the source read.
	 * Chooses a prefix, suffix or random valid start, then may reverse complement
	 * the copy. Copies qualities when present and retains the source IDs.
	 * @param a Nonempty source read; its bases and qualities are not modified
	 * @param randy Generator used for length, start and orientation
	 * @return New fragment with independent base and optional quality arrays
	 */
	private static Read getPiece(Read a, Random randy) {
		int len=randy.nextInt(a.length())+1;
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, quals, a.id, a.numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * Copies a requested-length fragment without modifying the source read.
	 * Fragment length is constrained to the specified length or read length if shorter.
	 * Chooses a prefix, suffix or random valid start.
	 * The copy may be reverse complemented; source IDs and optional qualities are retained.
	 *
	 * @param a Source read to extract fragment from
	 * @param randy Random number generator for position selection
	 * @param len Desired fragment length in bases
	 * @return New fragment, or null if the length after clipping is less than one
	 */
	private Read getPiece(Read a, Random randy, int len) {
		len=Tools.min(len, a.length());
		if(len<1){return null;}
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, quals, a.id, a.numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Primary unpaired input filename. */
	private String in1=null;
	
	/** Optional separate quality input passed to the reader factory. */
	private String qfin1=null;

	/** Optional output filename; null suppresses writing, not generation. */
	private String out1=null;
	
	/** Optional input format override. */
	private String extin=null;
	/** Optional output format override. */
	private String extout=null;

	/** Positive values request this first-fragment length; nonpositive values draw lengths. */
	private int forceLength=0;
	
	/*--------------------------------------------------------------*/

	/** Limit passed to the source reader; negative means unlimited. */
	private long readsIn=-1;
	/** Requested successful chimera count; must be positive under assertions. */
	private long readsOut=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Input format descriptor used for loading the retained source pool. */
	private final FileFormat ffin1;

	/** Output descriptor, or null when no writer is requested. */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostic stream selected by PreParser. */
	private PrintStream outstream=System.err;
	/** Enables additional diagnostics in this tool and selected I/O components. */
	public static boolean verbose=false;
	/** Accumulated reader/writer completion errors. */
	public boolean errorState=false;
	/** Parsed permission to overwrite existing output. */
	private boolean overwrite=true;
	/** Parsed request to append output. */
	private boolean append=false;
	
}
