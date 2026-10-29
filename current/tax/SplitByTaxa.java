package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Splits FASTA/FASTQ/SAM/BAM reads by the taxonomy parsed from each read name.
 * Each primary read is assigned to the first ancestor at or above
 * {@code level=}; its mate, if present, is written to the same bucket. Output
 * paths replace the first {@code %} in the {@code out=} pattern with that
 * ancestor's sanitized name.
 *
 * @author Brian Bushnell
 * @date November 23, 2015
 */
public class SplitByTaxa{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses arguments and splits reads by taxonomy.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		SplitByTaxa x=new SplitByTaxa(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses read input, per-taxon output pattern, tree, optional GI table, and
	 * target taxonomic level settings.
	 * @param args BBTools-style {@code key=value} arguments
	 */
	public SplitByTaxa(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set some shared static variables
		Shared.capBuffers(4);
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_GZIP=false;
		ReadWrite.setZipThreads(Shared.threads());

		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		String giTableFile=null;
		String taxTreeFile=null;

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("taxlevel") || a.equals("tl") || a.equals("level") || a.equals("lv")){
				taxLevelE=TaxTree.parseLevelExtended(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			out2=parser.out2;

			extin=parser.extin;
			extout=parser.extout;
		}

		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		assert(out1==null || out1.contains("%")) : "Output filename must contain % symbol.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain % symbol.";
		assert(taxTreeFile!=null) : "This requires a taxtree file.  On NERSC systems, set tree=auto.";

		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		assert(FastaReadInputStream.settingsOK());

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		TaxFilter.loadGiTable(giTableFile);
		tree=TaxFilter.loadTree(taxTreeFile);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens the input stream, dispatches reads to per-taxon writers, closes all
	 * streams and reports read/base counts.
	 * @param t Timer started by {@link #main(String[])}
	 */
	public void process(Timer t){

		//Create a read input stream
		ConcurrentReadInputStream cris=null;
		final HashMap<String, ConcurrentReadOutputStream> map=
				new HashMap<String, ConcurrentReadOutputStream>();
		Throwable thrown=null;
		try{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
			boolean paired=cris.paired();
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

			//Reset counters
			readsProcessed=0;
			basesProcessed=0;

			//Process the read stream
			processInner(cris, map);
		}catch(Throwable e){
			thrown=e;
		}finally{
			if(verbose){outstream.println("Finished; closing streams.");}

			thrown=writeStats(thrown);
			thrown=closeInput(cris, thrown);
			for(ConcurrentReadOutputStream ros : map.values()){
				thrown=closeOutput(ros, thrown);
			}
		}
		if(thrown!=null){rethrow(thrown);}

		//Report timing and results
		t.stop();
		outstream.println("Reads In:               \t"+readsProcessed+" reads       \t"+basesProcessed+" bases");
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Reads batches from {@code cris} and sends each primary read and its mate, if
	 * present, to the primary read's taxon writer. Classification first parses
	 * supported GI/NCBI/accession-style primary headers, then falls back to tree
	 * name lookup, then to an {@code UNKNOWN} synthetic node. The ancestor climb
	 * stops at the requested level, a coarser ancestor, the root, or
	 * {@code UNKNOWN}'s synthetic self-parent. A missing writer is opened by
	 * replacing the first {@code %} in the output pattern with the taxon name after
	 * spaces become underscores and path separators are removed.
	 *
	 * @param cris Input stream containing reads to process
	 * @param map Writer map keyed by taxon name and mutated by this method
	 */
	void processInner(final ConcurrentReadInputStream cris, HashMap<String, ConcurrentReadOutputStream> map){

		final HashMap<String, OutputPathOwner> pathOwners=
				new HashMap<String, OutputPathOwner>();

		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;

					TaxNode tn=tree.parseNodeFromHeader(r1.id, true);
					if(tn==null){tn=tree.getNodeByName(r1.id);}
					if(tn==null){tn=unknown;}
					//climb toward taxLevelE; terminates: root has id==pid, and the `unknown` fallback node also has id==pid(-99) so it never spins.
					while(tn.levelExtended<taxLevelE && tn.id!=tn.pid){tn=tree.getNode(tn.pid);}

					if(out1!=null){
						ConcurrentReadOutputStream ros=map.get(tn.name);
						if(ros==null){
							final int buff=4;
							FileFormat ffout1=null, ffout2=null;
							final String cleanName=toOutputName(tn.name);
							ffout1=makeOutputFile(out1, cleanName);
							registerOutputPath(pathOwners, tn.name, ffout1);
							if(out2!=null){
								ffout2=makeOutputFile(out2, cleanName);
								registerOutputPath(pathOwners, tn.name, ffout2);
							}
							ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
							map.put(tn.name, ros);
							//TODO: If ros.start() fails after opening a writer backend,
							// close queues poison for a consumer that never started;
							// reliable handle release for that case belongs in CROS.
							ros.start(); //Start the stream
						}
						ArrayList<Read> temp=new ArrayList<Read>(1); //Kind of inefficient
						temp.add(r1);
						ros.add(temp, 0);
					}
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		//Do anything necessary after processing

	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes accumulated statistics and preserves unchecked write failures.
	 * @param thrown Primary failure to preserve, or null
	 * @return Primary or new failure with cleanup failures suppressed */
	private Throwable writeStats(Throwable thrown){
		try{
			errorState|=ReadStats.writeAll();
		}catch(RuntimeException e){
			thrown=recordFailure(thrown, e);
		}catch(Error e){
			thrown=recordFailure(thrown, e);
		}
		return thrown;
	}

	/** Closes the input stream and preserves unchecked close failures.
	 * @param cris Input stream to close, or null if it was never opened
	 * @param thrown Primary failure to preserve, or null
	 * @return Primary or new failure with cleanup failures suppressed */
	private Throwable closeInput(final ConcurrentReadInputStream cris, Throwable thrown){
		if(cris!=null){
			try{
				errorState|=ReadWrite.closeStream(cris);
			}catch(RuntimeException e){
				thrown=recordFailure(thrown, e);
			}catch(Error e){
				thrown=recordFailure(thrown, e);
			}
		}
		return thrown;
	}

	/** Closes one per-taxon output stream and preserves unchecked close failures.
	 * @param ros Output stream to close
	 * @param thrown Primary failure to preserve, or null
	 * @return Primary or new failure with cleanup failures suppressed */
	private Throwable closeOutput(final ConcurrentReadOutputStream ros, Throwable thrown){
		try{
			errorState|=ReadWrite.closeStream(ros);
		}catch(RuntimeException e){
			thrown=recordFailure(thrown, e);
		}catch(Error e){
			thrown=recordFailure(thrown, e);
		}
		return thrown;
	}

	/** Adds {@code failure} as suppressed to the current primary failure.
	 * @param primary Primary failure to preserve, or null
	 * @param failure Cleanup failure to record
	 * @return Primary or cleanup failure */
	private static Throwable recordFailure(final Throwable primary,
			final Throwable failure){
		if(primary==failure){return primary;}
		if(primary!=null){
			primary.addSuppressed(failure);
			return primary;
		}
		return failure;
	}

	/** Rethrows {@code failure} without wrapping unchecked failures.
	 * @param failure Throwable to rethrow */
	private static void rethrow(final Throwable failure){
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException(failure);
	}

	/** Builds an output descriptor after replacing the first percent pattern.
	 * @param pattern Output path pattern containing a percent placeholder
	 * @param cleanName Sanitized taxon name to substitute
	 * @return Output descriptor for the resolved path */
	private FileFormat makeOutputFile(final String pattern, final String cleanName){
		return FileFormat.testOutput(pattern.replaceFirst("%", Matcher.quoteReplacement(cleanName)),
				FileFormat.FASTQ, extout, false, overwrite, append, false);
	}

	/** Converts a raw taxon name to the path component used for split output.
	 * @param name Raw taxon name
	 * @return Path component for split output */
	private static String toOutputName(final String name){
		return name.replaceAll("\\s+", "_").replaceAll("[/\\\\]", "");
	}

	/** Registers one output path and rejects collisions before a writer opens.
	 * @param pathOwners Prior output path owners by normalized absolute path
	 * @param rawName Raw taxon name for the new writer
	 * @param ffout Output descriptor to register */
	private static void registerOutputPath(
			final HashMap<String, OutputPathOwner> pathOwners,
			final String rawName, final FileFormat ffout){
		final String path=ffout.name();
		final String normalized=new File(path).getAbsoluteFile().toPath()
				.normalize().toString();
		final OutputPathOwner owner=pathOwners.get(normalized);
		if(owner!=null){
			if(owner.rawName.equals(rawName)){
				throw new IllegalArgumentException(
						"SplitByTaxa mate output paths resolve to the same file for taxon '"+
						rawName+"': "+owner.path+" and "+path);
			}
			throw new IllegalArgumentException(
					"Distinct taxon names produce the same SplitByTaxa output path: "+
					path+" ("+owner.rawName+" -> "+owner.path+", "+
					rawName+" -> "+path+")");
		}
		pathOwners.put(normalized, new OutputPathOwner(rawName, path));
	}

	/** Owner of a normalized output path. */
	private static final class OutputPathOwner{
		/** Records one output path owner.
		 * @param rawName_ Raw taxon name that first used the path
		 * @param path_ Original output path */
		OutputPathOwner(final String rawName_, final String path_){
			rawName=rawName_;
			path=path_;
		}
		/** Raw taxon name that first used this path. */
		final String rawName;
		/** Original output path before normalization. */
		final String path;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary read input path. */
	private String in1=null;
	/** Secondary read input path for paired input, or null. */
	private String in2=null;

	/** Primary output path pattern; its first {@code %} is replaced per taxon. */
	private String out1=null;
	/** Secondary output path pattern for paired output, or null. */
	private String out2=null;

	/** Explicit input extension override. */
	private String extin=null;
	/** Explicit output extension override. */
	private String extout=null;

	/** Extended taxonomic level at which reads are split. */
	private int taxLevelE=TaxTree.stringToLevelExtended("phylum");

	/*--------------------------------------------------------------*/

	/** Reads consumed, counting pairs according to {@link Read#pairCount()}. */
	protected long readsProcessed=0;
	/** Bases consumed from reads and mates. */
	protected long basesProcessed=0;

	/** Maximum primary input fragments to process, or -1 for no limit. */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input descriptor. */
	private final FileFormat ffin1;
	/** Secondary input descriptor, or null. */
	private final FileFormat ffin2;

	/** Taxonomy tree used for header/name resolution and ancestor climbing. */
	private final TaxTree tree;

	/** Synthetic bucket for reads with unresolved names. */
	private final TaxNode unknown=new TaxNode(-99, -99, TaxTree.LIFE, TaxTree.LIFE_E, "UNKNOWN");

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for status messages and summary statistics. */
	private PrintStream outstream=System.err;
	/** Whether to print verbose stream lifecycle diagnostics. */
	public static boolean verbose=false;
	/** Set when input, output or statistics closing reports an error. */
	public boolean errorState=false;
	/** Whether output files may be overwritten. */
	private boolean overwrite=true;
	/** Whether output should append to existing files. */
	private boolean append=false;
	/** Legacy field; per-taxon writers are constructed with ordered=false. */
	private final boolean ordered=false;

}
