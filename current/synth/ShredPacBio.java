package synth;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Converts long single-molecule reads (e.g. PacBio HiFi) into fake paired-end reads by walking
 * each molecule front-to-back and consuming every base exactly once.
 *
 * <p>Each consecutive {@code 2*readlen} window of a molecule becomes one pair: the first
 * {@code readlen} bases are read 1, the next {@code readlen} bases are read 2 (optionally
 * reverse-complemented, with quality reversed, so the pair is FR-oriented like real Illumina data).
 * Motivation (Brian, 2026-06-22): with the extra pairing constraint a ~{@code 2*readlen} span can be
 * placed using only two {@code readlen}-bp alignments, which is cheaper than one long alignment and
 * helps in moderately-repetitive regions, while letting CallVariants see a bbmap-style mapping that
 * (per the PacBio race probe) scores better than long-read alignments.</p>
 *
 * <p>Three output streams keep every input base accounted for in exactly one place:
 * <ul>
 *   <li>{@code out=}  proper pairs (interleaved fastq)</li>
 *   <li>{@code outs=} singletons: a leftover read with no mate of length &gt;= {@code minlen}</li>
 *   <li>{@code outd=} discards: leftover fragments shorter than {@code minlen}</li>
 * </ul>
 *
 * <p>Per molecule the walk is: take full {@code 2*readlen} pairs while they fit; then the residual
 * ({@code 0..2*readlen-1} bases) is classified once:
 * residual {@code < minlen} -&gt; discard; {@code [minlen, readlen)} -&gt; singleton (no room for a mate);
 * {@code [readlen, 2*readlen)} -&gt; r1 is a full {@code readlen} read and the rest is its mate, emitted as
 * a (possibly length-uneven) pair if the mate is {@code >= minlen}, else r1 becomes a singleton and the
 * sub-{@code minlen} tail is discarded. A {@code repOK}-style assertion checks that pair+singleton+discard
 * bases sum to the molecule length (every base used exactly once).</p>
 *
 * <p>First pass: sequential pairs only. Long-mate / unsequenced-insert mode is deferred (FakeReads
 * already covers that shape).</p>
 *
 * @author UMP45
 * @date June 22, 2026
 */
public class ShredPacBio {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		ShredPacBio x=new ShredPacBio(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses arguments and configures shredding parameters and input/output streams.
	 * @param args Command-line arguments
	 */
	public ShredPacBio(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("readlen") || a.equals("rl") || a.equals("length") || a.equals("len") || a.equals("readlength")){
				readLength=Parse.parseIntKMG(b);
			}else if(a.equals("minlen") || a.equals("ml") || a.equals("minlength")){
				minLength=Parse.parseIntKMG(b);
			}else if(a.equals("rcompmate") || a.equals("rcomp") || a.equals("rc") || a.equals("rcmate")){
				rcompMate=Parse.parseBoolean(b);
			}else if(a.equals("out") || a.equals("out1") || a.equals("outpairs") || a.equals("outp")){
				outPairs=b;
			}else if(a.equals("outs") || a.equals("outsingle") || a.equals("outsingles") || a.equals("outsingletons") || a.equals("singletons")){
				outSingles=b;
			}else if(a.equals("outd") || a.equals("outdiscard") || a.equals("outdiscards") || a.equals("discards")){
				outDiscards=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{//Process parser fields
			Parser.processQuality();

			if(maxReads<0){maxReads=parser.maxReads;}
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			if(outPairs==null){outPairs=parser.out1;}

			extin=parser.extin;
			extout=parser.extout;
		}

		assert(readLength>0) : "readlen must be positive: "+readLength;
		assert(minLength>0 && minLength<=readLength) : "minlen must be in 1.."+readLength+": "+minLength;
		fragLength=2*readLength;

		assert(FastaReadInputStream.settingsOK());

		if(in1==null){throw new RuntimeException("Error - an input file is required (in=).");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		//Treat "null" as a disabled stream
		if(outPairs!=null && outPairs.equalsIgnoreCase("null")){outPairs=null;}
		if(outSingles!=null && outSingles.equalsIgnoreCase("null")){outSingles=null;}
		if(outDiscards!=null && outDiscards.equalsIgnoreCase("null")){outDiscards=null;}

		if(!Tools.testOutputFiles(overwrite, append, false, outPairs, outSingles, outDiscards)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to one of the output files "
					+outPairs+", "+outSingles+", "+outDiscards+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, outPairs, outSingles, outDiscards)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		ffout1=FileFormat.testOutput(outPairs, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutS=FileFormat.testOutput(outSingles, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutD=FileFormat.testOutput(outDiscards, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates streams, drives processing, and reports throughput and per-stream counts.
	 * @param t Timer for throughput reporting
	 */
	void process(Timer t){

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		if(cris.paired()){KillSwitch.kill("ShredPacBio requires single-ended (unpaired) input.");}

		final ConcurrentReadOutputStream rosP=makeStream(ffout1);
		final ConcurrentReadOutputStream rosS=makeStream(ffoutS);
		final ConcurrentReadOutputStream rosD=makeStream(ffoutD);

		processInner(cris, rosP, rosS, rosD);

		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, rosP, rosS, rosD);
		if(verbose){outstream.println("Finished.");}

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println();
		outstream.println("Pairs out:       \t"+pairsOut+" pairs ("+(2*pairsOut)+" reads), "+basesPairs+" bases");
		outstream.println("Singletons out:  \t"+singlesOut+" reads, "+basesSingles+" bases");
		outstream.println("Discards out:    \t"+discardsOut+" reads, "+basesDiscards+" bases");

		assert(basesPairs+basesSingles+basesDiscards==basesProcessed) :
			"Base accounting mismatch: "+basesPairs+"+"+basesSingles+"+"+basesDiscards+" != "+basesProcessed;

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Opens and starts a single-file (interleaved-capable) output stream, or returns null if ff is null. */
	private static ConcurrentReadOutputStream makeStream(FileFormat ff){
		if(ff==null){return null;}
		final int buff=4;
		ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start();
		return ros;
	}

	/**
	 * Reads input lists and routes each molecule's shreds to the three output streams.
	 * One output list per input list per stream keeps each stream's ordering monotonic by list id.
	 */
	private void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream rosP,
			final ConcurrentReadOutputStream rosS, final ConcurrentReadOutputStream rosD){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			final ArrayList<Read> pairs=new ArrayList<Read>();
			final ArrayList<Read> singles=new ArrayList<Read>();
			final ArrayList<Read> discards=new ArrayList<Read>();

			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				readsProcessed++;
				basesProcessed+=r.length();
				processRead(r, pairs, singles, discards);
			}

			if(rosP!=null){rosP.add(pairs, ln.id);}
			if(rosS!=null){rosS.add(singles, ln.id);}
			if(rosD!=null){rosD.add(discards, ln.id);}

			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Shreds one molecule into pairs/singletons/discards, using every base exactly once.
	 * @param r Input molecule
	 * @param pairs Accumulator for interleaved proper pairs (r1 added with mate attached)
	 * @param singles Accumulator for singletons
	 * @param discards Accumulator for discards
	 */
	private void processRead(final Read r, final ArrayList<Read> pairs, final ArrayList<Read> singles,
			final ArrayList<Read> discards){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final int L=(bases==null ? 0 : bases.length);
		if(L==0){return;}
		final String origHeader=r.id;

		long basesAccounted=0;
		int pos=0;

		//Take full 2*readLength fragments as proper pairs
		while(L-pos>=fragLength){
			addPair(bases, quals, pos, pos+readLength, pos+fragLength, origHeader, pairs);
			basesAccounted+=fragLength;
			pos+=fragLength;
		}

		//Classify the residual (0..fragLength-1 bases) exactly once
		final int rem=L-pos;
		if(rem>0){
			if(rem<minLength){
				addSingle(bases, quals, pos, L, origHeader, "d", discards);
				discardsOut++; basesDiscards+=rem;
				basesAccounted+=rem;
			}else if(rem<readLength){
				addSingle(bases, quals, pos, L, origHeader, "s", singles);
				singlesOut++; basesSingles+=rem;
				basesAccounted+=rem;
			}else{//readLength <= rem < fragLength
				final int r1end=pos+readLength;
				final int mateLen=rem-readLength;
				if(mateLen>=minLength){
					addPair(bases, quals, pos, r1end, L, origHeader, pairs);
					basesAccounted+=rem;
				}else{
					addSingle(bases, quals, pos, r1end, origHeader, "s", singles);
					singlesOut++; basesSingles+=readLength;
					basesAccounted+=readLength;
					if(mateLen>0){
						addSingle(bases, quals, r1end, L, origHeader, "d", discards);
						discardsOut++; basesDiscards+=mateLen;
						basesAccounted+=mateLen;
					}
				}
			}
		}

		assert(basesAccounted==L) : "Per-molecule base accounting mismatch: "+basesAccounted+" != "+L+" for "+r.id;
	}

	/**
	 * Builds one proper pair (r1=[start,mid), r2=[mid,end)) and appends r1 (with mate attached) to the list.
	 * r2 is reverse-complemented (quality reversed) when {@code rcompMate}. Updates pair counters.
	 */
	private void addPair(final byte[] bases, final byte[] quals, final int start, final int mid, final int end,
			final String origHeader, final ArrayList<Read> pairs){
		final long nid=pairsOut;
		final byte[] b1=KillSwitch.copyOfRange(bases, start, mid);
		final byte[] b2=KillSwitch.copyOfRange(bases, mid, end);
		final byte[] q1=(quals==null ? null : KillSwitch.copyOfRange(quals, start, mid));
		final byte[] q2=(quals==null ? null : KillSwitch.copyOfRange(quals, mid, end));
		if(rcompMate){
			Vector.reverseComplementInPlaceFast(b2);
			if(q2!=null){Vector.reverseInPlace(q2);}
		}
		final Read r1=new Read(b1, q1, pairHeader(nid, '1', origHeader, start), nid, 0);
		final Read r2=new Read(b2, q2, pairHeader(nid, '2', origHeader, start), nid, Read.PAIRNUMMASK);
		r1.mate=r2;
		r2.mate=r1;
		pairs.add(r1);
		pairsOut++;
		basesPairs+=(end-start);
	}

	/**
	 * Builds one unpaired read ([start,end)) and appends it to the given list.
	 * @param tag Short label distinguishing the source ("s"=singleton, "d"=discard)
	 */
	private void addSingle(final byte[] bases, final byte[] quals, final int start, final int end,
			final String origHeader, final String tag, final ArrayList<Read> list){
		final long nid=("d".equals(tag) ? discardsOut : singlesOut);
		final byte[] b=KillSwitch.copyOfRange(bases, start, end);
		final byte[] q=(quals==null ? null : KillSwitch.copyOfRange(quals, start, end));
		final String name="shredpb_"+tag+nid+"\t"+(origHeader==null ? "." : origHeader)+"\t"+start;
		final Read r=new Read(b, q, name, nid, 0);
		list.add(r);
	}

	/**
	 * Builds a paired-read header: {@code shredpb_<nid> <pairDigit>:N:0\t<originalHeader>\t<originalPairStart>}.
	 * The {@code <pairDigit>:} after the first space is the Casava-style token BBTools uses to auto-detect
	 * pairs in an interleaved file; the tab-fields preserve full provenance (source molecule + start offset).
	 * @param nid Global pair index (shared by r1 and r2 so they are recognized as one pair)
	 * @param pairDigit '1' for read 1, '2' for read 2
	 * @param origHeader Full original molecule header (or null)
	 * @param start 0-based start of this pair's fragment within the original molecule
	 */
	private static String pairHeader(final long nid, final char pairDigit, final String origHeader, final int start){
		return "shredpb_"+nid+" "+pairDigit+":N:0\t"+(origHeader==null ? "." : origHeader)+"\t"+start;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String outPairs=null;
	private String outSingles=null;
	private String outDiscards=null;

	private String extin=null;
	private String extout=null;

	/** Length of each emitted read (r1 and full r2); a pair consumes 2*readLength bases. */
	private int readLength=500;
	/** Minimum length for any emitted read; shorter residuals go to discards. */
	private int minLength=200;
	/** 2*readLength; the span consumed per full pair. */
	private int fragLength;
	/** Reverse-complement r2 (and reverse its quality) so pairs are FR-oriented. */
	private boolean rcompMate=true;

	private long maxReads=-1;

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long pairsOut=0;
	protected long singlesOut=0;
	protected long discardsOut=0;
	protected long basesPairs=0;
	protected long basesSingles=0;
	protected long basesDiscards=0;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutS;
	private final FileFormat ffoutD;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

}
