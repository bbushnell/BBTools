package stream;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.bam.BamIndexWriter;
import structures.ListNum;
import var2.SamFilter;
import var2.ScafMap;
import var2.Scaffold;

/**
 * Wrapper for streaming SAM/BAM files with optional filtering, conversion, and CIGAR normalization.
 * Supports emitting reads or SAM/BAM output, running SamFilter-based screening, and generating BAM
 * indexes when requested.
 *
 * @author Brian Bushnell, Isla
 * @date November 6, 2025
 */
public class SamStreamerWrapper{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Entry point for command-line execution.
	 * Initializes the wrapper and runs processing with a timer. */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		SamStreamerWrapper x=new SamStreamerWrapper(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs the wrapper from command-line arguments.
	 * Handles preparsing for config/help, parses options, initializes IO formats, and
	 * configures SAM parsing behavior based on requested output.
	 */
	SamStreamerWrapper(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			threadsIn=parser.threadsIn;
			threadsOut=parser.threadsOut;
			in1=parser.in1;
			out1=parser.out1;
		}

		//Do input/output setup
		fixExtensions();
		checkFileExistence();

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, null, true, true, false, true);
		ffout2=(out2==null ? null : FileFormat.testOutput(out2, FileFormat.SAM, null, true, true, false, true));

		//Determine if we need to parse SAM fields or can skip for performance.
		//BED filtering needs rname/pos/cigar, and the outu split needs intact records, so both force parsing.
		if(!forceParse && !fixCigar && !eqx && out2==null && bed==null && (ffout1==null || !ffout1.samOrBam())){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}

		//Enable optimizations for SAM->SAM conversion
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;
	}

	/**
	 * Parses command-line options and initializes the SamFilter and parser state.
	 */
	private Parser parse(String[] args){

		//Create filter
		filter=new SamFilter();
		filter.includeNonPrimary=true;
		filter.includeLengthZero=true;
		boolean doFilter=true;

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("reference")){
				ref=b;
			}else if(a.equals("rnameasbytes")){
				SamLine.RNAME_AS_BYTES=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				Parser.parseSam(arg, a, b);
				fixCigar=true;
				//SAM 1.4 introduced =/X; requesting it IS requesting eqx (convert M->=/X). SAM 1.3 stays M-form.
				if(SamLine.VERSION!=1.3f){eqx=true;}
			}else if(a.equals("eqx") || a.equals("eqxcigar") || a.equals("toeqx")){
				//Convert M->=/X (resolving via MD tag or ref=). Synonym for sam=1.4's cigar effect. Default false;
				//INDEPENDENT of normalization (eqx does NOT imply left-shift).
				eqx=Parse.parseBoolean(b);
			}else if(a.equals("normalize") || a.equals("canonicalize") || a.equals("canonicalise")
					|| a.equals("canonize") || a.equals("leftalign")){
				//Deterministic left-alignment of indels vs the reference (alignment canonicalization).
				//IMPLIES eqx (can't canonicalize an ambiguous M) and REQUIRES ref= (crash-loud if absent).
				normalize=Parse.parseBoolean(b);
				if(normalize){eqx=true;}
			}

			//Scaffold rename (old->new name TSV); rewrites @SQ SN: header lines + record RNAME/RNEXT
			else if(a.equals("rename") || a.equals("renamescaffolds") || a.equals("renamebylist")){
				renameFile=b;
			}

			//Output-split + BED filter parameters
			else if(a.equals("outu") || a.equals("outunmatched")){
				//Second output: reads that do NOT pass the filter(s) go here (raw, untransformed). Null disables the split.
				out2=b;
			}else if(a.equals("bed") || a.equals("bedfile")){
				//BED file for positional filtering (BedReadFilter). Reads are kept/routed by mof+include.
				bed=b;
			}else if(a.equals("minoverlapfraction") || a.equals("mof") || a.equals("overlap")){
				//Fraction of a read's reference span that must lie in the BED to match; 0=any overlap, 1=containment.
				mof=Double.parseDouble(b);
			}else if(a.equals("include") || a.equals("bedinclude") || a.equals("includebed")){
				//true=keep reads matching the BED; false=keep reads that do NOT (reverse/exclude).
				include=Parse.parseBoolean(b);
			}

			//Filter parameters
			else if(a.equals("filter")){
				doFilter=Parse.parseBoolean(b);
			}else if(filter.parse(arg, a, b)){
				//do nothing

			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && !arg.contains("=") && parser.in1==null &&
					FileFormat.isSamOrBamFile(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else if(i==1 && !arg.contains("=") && parser.out1==null && parser.in1!=null &&
					(FileFormat.isSequence(arg) || FileFormat.isBaiFile(arg))){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		//Disable filter if requested
		if(!doFilter){filter=null;}

		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Normalizes input/output filename extensions (adds/removes compression suffixes).
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}

	/**
	 * Validates that input exists, output paths are writable, and filenames are unique.
	 */
	private void checkFileExistence(){
		//Ensure input file exists
		if(in1==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}

		//Ensure output files can be written
		if(!Tools.testOutputFiles(true, false, false, out1, out2)){
			throw new RuntimeException("\nCan't write to output file "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read input file "+in1+"\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates Streamer/Writer components and processes all records according to parsed options.
	 */
	void process(Timer t){

		//Determine processing mode
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputBai=(ffout1!=null && ffout1.bai());
		final boolean useSharedHeader=inputSam && outputSam;
		final boolean makeReads=(outputReads || !inputSam);
		if(!inputSam) {
			System.err.println("Input is "+ffin1.formatString()+"; sam filter disabled.");
			filter=null;
			ref=null;
			bed=null; //BED filtering needs reference coordinates, unavailable for non-SAM input
			normalize=false; //no CIGARs to canonicalize for non-SAM input
			renameFile=null; //no scaffold names to rename for non-SAM input
		}

		//Alignment canonicalization needs reference context (rolling can reach beyond the read's footprint)
		if(normalize && ref==null){
			throw new RuntimeException("\nnormalize/canonicalize requires a reference: add ref=<fasta>.\n");
		}

		if(outputBai){
			assert(ffin1.bam()) : "bai output requires bam input.";
			try{
				BamIndexWriter.writeIndex(in1, out1);
			}catch(Throwable e){
				KillSwitch.exceptionKill(e);
			}
			t.stop();
			outstream.println("Time:                         \t"+t);
			return;
		}

		//Load reference if specified
		if(ref!=null){
			ScafMap.loadReference(ref, true);
			SamLine.RNAME_AS_BYTES=false;
		}

		//Build the BED filter if requested (reads failing it route to outu, or are dropped when outu is unset)
		if(bed!=null){bedFilter=new BedReadFilter(bed, mof, include);}

		//Build the scaffold renamer if requested
		if(renameFile!=null){renamer=new ScaffoldRenamer(renameFile);}

		//Create streamer and writers (fw=primary out; fw2=outu, the non-passing split)
		Streamer st=StreamerFactory.makeStreamer(ffin1, null, ordered, maxReads, useSharedHeader, makeReads, threadsIn);
		Writer fw=(ffout1==null ? null : WriterFactory.makeWriter(ffout1, null, threadsOut, null, useSharedHeader));
		Writer fw2=(ffout2==null ? null : WriterFactory.makeWriter(ffout2, null, threadsOut, null, useSharedHeader));

		//Process data
		st.start();
		//Rename @SQ SN: header lines AFTER the input has loaded the shared header, BEFORE the writer emits it
		if(renamer!=null && outputSam){renameSharedHeader();}
		if(fw!=null){fw.start();}
		if(fw2!=null){fw2.start();}

		if(outputReads || !inputSam){
			processAsReads(st, fw, fw2);
		}else{
			processAsSam(st, fw, fw2);
		}

		//Wait for writers to finish
		if(fw!=null){
			errorState|=fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		if(fw2!=null){
			errorState|=fw2.poisonAndWait();
			readsOutu=fw2.readsWritten();
			basesOutu=fw2.basesWritten();
		}

		//Check for errors
		errorState|=st.errorState();

		//Print statistics
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2f Mbp/sec", (basesProcessed/(double)(t.elapsed))*1000));
		if(ffout1!=null){
			outstream.println("Reads Out:          "+readsOut);
			outstream.println("Bases Out:          "+basesOut);
		}
		if(ffout2!=null){
			outstream.println("Reads Outu:         "+readsOutu);
			outstream.println("Bases Outu:         "+basesOutu);
		}

		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private void processAsReads(Streamer st, Writer fw, Writer fw2){
		final boolean splitting=(fw2!=null);
		for(ListNum<Read> ln=st.nextList(); ln!=null && ln.size()>0; ln=st.nextList()){
			ArrayList<Read> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			//Passthrough (out==list, no copy) only when nothing filters or splits the stream
			final boolean passthrough=(filter==null && bedFilter==null && !splitting);
			ArrayList<Read> out=(passthrough ? list : new ArrayList<Read>(list.size()));
			ArrayList<Read> outu=(splitting ? new ArrayList<Read>() : null);

			for(Read r : list){
				final int len=r.length();
				readsProcessed++;
				basesProcessed+=len;

				//A read passes only if it clears every active filter; failures route to outu (raw) or are dropped
				boolean passes=(filter==null || filter.passesFilter(r.samline)) &&
					(bedFilter==null || bedFilter.passes(r.samline));
				if(!passes){
					if(outu!=null){outu.add(r);}
					continue;
				}
				if(!passthrough){out.add(r);}
			}

			if(fw!=null){fw.addReads(new ListNum<Read>(out, ln.id));}
			if(fw2!=null){fw2.addReads(new ListNum<Read>(outu, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/**
	 * Processes records as SamLine objects for SAM/BAM-to-SAM/BAM workflows, applying filtering and CIGAR fixes.
	 */
	private void processAsSam(Streamer st, Writer fw, Writer fw2){
		final boolean splitting=(fw2!=null);
		final boolean transforming=(fixCigar || eqx);
		final ScafMap scafMap=(normalize ? ScafMap.defaultScafMap() : null); //reference for indel left-alignment
		for(ListNum<SamLine> ln=st.nextLines(); ln!=null && ln.size()>0; ln=st.nextLines()){
			ArrayList<SamLine> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			//Passthrough (out==list, no copy) only when nothing filters, transforms, renames, or splits the stream
			final boolean passthrough=(filter==null && bedFilter==null && renamer==null && !transforming && !splitting);
			ArrayList<SamLine> out=(passthrough ? list : new ArrayList<SamLine>(list.size()));
			ArrayList<SamLine> outu=(splitting ? new ArrayList<SamLine>() : null);

			for(SamLine sl : list){
				final int len=sl.lengthOrZero();
				readsProcessed++;
				basesProcessed+=len;

				//A read passes only if it clears every active filter (SamFilter AND BedReadFilter).
				//Reads that fail route to outu UNTRANSFORMED (the BED filter precedes cigar-fix), or are dropped if outu is unset.
				boolean passes=(filter==null || filter.passesFilter(sl)) &&
					(bedFilter==null || bedFilter.passes(sl)); //filters use ORIGINAL scaffold names
				if(!passes){
					if(outu!=null){
						if(renamer!=null){renamer.renameRecord(sl);} //outu records must match the renamed header too
						outu.add(sl);
					}
					continue;
				}

				//Length-zero guard for non-SAM output (kept reads only); unchanged behavior for SAM/BAM output
				if(!(len>0 || ffout1==null || ffout1.samOrBam())){continue;}

				//CIGAR transforms on KEPT reads. eqx (convert M->=/X) and sam-version conversion are INDEPENDENT.
				if(transforming && sl.cigar!=null){
					if(eqx){
						//eqx (= sam=1.4): emit =/X, resolving M via the MD tag or loaded reference (ScafMap from ref=).
						//Crash-loud under -ea on an M-only cigar with neither MD nor ref: do what was asked, or fail -- never
						//silently keep an unconverted cigar or drop it.
						rebuildCigar14(sl, false);
					}else{
						sl.setCigar(SamLine.toCigar13(sl.cigar)); //sam=1.3 only: version conversion (=/X->M)
					}
				}

				//Left-align indels vs the reference (runs AFTER eqx, so M is already resolved to =/X)
				if(normalize && sl.cigar!=null){
					final Scaffold scaf=scafMap.getScaffold(sl.rnameS());
					if(scaf!=null && scaf.bases!=null){CigarNormalizer.normalize(sl, scaf.bases);}
				}

				//Scaffold rename LAST, after BED filter + cigar transforms (which all use the original names)
				if(renamer!=null){renamer.renameRecord(sl);}

				if(!passthrough){out.add(sl);}
			}

			if(fw!=null){fw.addLines(new ListNum<SamLine>(out, ln.id));}
			if(fw2!=null){fw2.addLines(new ListNum<SamLine>(outu, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/** Rewrites the @SQ SN: scaffold names in the shared input header (in place) via the renamer, so the
	 * emitted output header matches the renamed records. Must run after the header is loaded and before the
	 * writer emits it. Both out and outu writers share this header. */
	private void renameSharedHeader(){
		final ArrayList<byte[]> hdr=SamReadInputStream.getSharedHeader(true);
		if(hdr==null){return;}
		for(int i=0; i<hdr.size(); i++){
			final String line=new String(hdr.get(i));
			final String renamed=renamer.renameHeaderLine(line);
			if(!renamed.equals(line)){hdr.set(i, renamed.getBytes());}
		}
	}

	/** Rebuilds a SamLine's CIGAR in SAM 1.4 (=/X) representation from its match string.
	 * @param allowM true keeps M ops (version-conversion only); false resolves M->=/X via MD tag or
	 * loaded reference (and crashes loud on an M-only cigar with neither -- never silently drops/keeps it). */
	private static void rebuildCigar14(SamLine sl, boolean allowM){
		byte[] shortMatch=sl.toShortMatch(allowM);
		byte[] longMatch=Read.toLongMatchString(shortMatch);
		int start=sl.pos-1;
		int stop=start+Read.calcMatchLength(longMatch)-1;
		sl.setCigar(SamLine.toCigar14(longMatch, start, stop, Integer.MAX_VALUE, sl.seq));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** SAM/BAM filter for quality/mapping criteria; null disables filtering. */
	private SamFilter filter;
	/** Positional BED-overlap filter; null disables BED filtering. Built in process() from bed/mof/include. */
	private BedReadFilter bedFilter;
	/** Scaffold renamer (old->new); null disables renaming. Built in process() from renameFile. */
	private ScaffoldRenamer renamer;

	/** Primary input SAM/BAM filename. */
	private String in1=null;
	/** Primary output filename (SAM/BAM/FASTQ/FASTA). */
	private String out1=null;
	/** Secondary output (outu): reads that do NOT pass the filter(s), emitted untransformed; null disables the split. */
	private String out2=null;
	/** Reference path for coordinate lookups or CIGAR reconstruction. */
	private String ref=null;
	/** BED file for positional filtering; null disables it. */
	private String bed=null;
	/** Scaffold-rename TSV (old<TAB>new); null disables renaming. */
	private String renameFile=null;
	/** Minimum overlap fraction for the BED filter (0=any overlap, 1=containment). */
	private double mof=0.0;
	/** BED filter sense: true keeps reads matching the BED, false keeps non-matching. */
	private boolean include=true;

	/** Input file format descriptor. */
	private FileFormat ffin1;
	/** Output file format descriptor. */
	private FileFormat ffout1;
	/** Secondary (outu) output file format descriptor; null when no split is requested. */
	private FileFormat ffout2;

	/** Thread count for input streaming (-1 to auto-detect). */
	private int threadsIn=-1;
	/** Thread count for output writing (-1 to auto-detect). */
	private int threadsOut=-1;

	/** Number of reads processed by the pipeline. */
	private long readsProcessed=0;
	/** Number of reads emitted to the primary writer. */
	private long readsOut=0;
	/** Number of reads emitted to the outu (non-passing) writer. */
	private long readsOutu=0;
	/** Number of bases processed by the pipeline. */
	private long basesProcessed=0;
	/** Number of bases emitted to the primary writer. */
	private long basesOut=0;
	/** Number of bases emitted to the outu (non-passing) writer. */
	private long basesOutu=0;

	/*--------------------------------------------------------------*/

	/** True if any stage of processing encountered an error. */
	public boolean errorState=false;
	/** Preserve input order in output when true. */
	public boolean ordered=true;
	/** Limit on reads to process; -1 means no limit. */
	private long maxReads=-1;
	/**
	 * Force parsing of all SAM fields even when not required by the output format.
	 */
	private boolean forceParse;
	/** Normalize CIGAR strings to the requested SAM version. */
	private boolean fixCigar;
	/** Convert CIGAR M ops to =/X form (via MD tag or ref). Default false; independent of sam= and of normalization. */
	private boolean eqx=false;
	/** Left-align (canonicalize) indels vs the reference. Default false; implies eqx and requires ref=. */
	private boolean normalize=false;

	/*--------------------------------------------------------------*/

	/** Output stream used for status messages and timing summaries. */
	private PrintStream outstream=System.err;
	/** Enables verbose debugging output. */
	public static boolean verbose=false;

}