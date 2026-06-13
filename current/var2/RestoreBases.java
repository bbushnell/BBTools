package var2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SamReadInputStream;
import stream.SamWriterST2;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Restores SEQ and QUAL fields onto secondary (0x100) and supplementary (0x800)
 * alignments by copying them from the primary (non-supplementary) alignment of the
 * same read, matched by read name.  Aligners such as minimap2 emit SEQ=* / QUAL=*
 * on secondary/supplementary records, which makes those records unusable for variant
 * calling (CallVariants evicts reads with null bases).  The only correct way to fill
 * them is to copy the true bases/quals from the primary record — restoring from the
 * MD tag would incur reference bias (a true het at AF 0.50 degrades to ~0.17).
 *
 * <p>Method: the input is partitioned into WAYS headerless sam.gz subfiles by
 * qname.hashCode()%WAYS so that all alignments of a read land in the same subfile.
 * Each subfile is then processed alone, single-threaded: load, name-sort, group by
 * (qname, pairnum), and for each group copy the primary's seq/qual onto its
 * secondaries/supplementaries.  Two modes: <b>fix</b> (write seq/qual directly;
 * supplementary hard-clips are converted H&rarr;S so the full read is valid) and
 * <b>tag</b> (attach OS:Z:/OQ:Z: tags in original sequencing orientation, leaving
 * SEQ=*).  The output header is rewritten with SO:unsorted.
 *
 * @author Brian Bushnell, UMP45
 * @date June 13, 2026
 */
public class RestoreBases {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point.
	 * @param args Command-line arguments specifying input/output files and options
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RestoreBases x=new RestoreBases(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor that parses command-line arguments and initializes processing parameters.
	 * @param args Command-line arguments containing file paths and processing options
	 */
	public RestoreBases(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.SET_FROM_OK=true;

		Parser parser=new Parser();

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ways")){
				ways=Integer.parseInt(b);
			}else if(a.equals("mode")){
				if(b.equalsIgnoreCase("fix")){mode=FIX;}
				else if(b.equalsIgnoreCase("tag")){mode=TAG;}
				else{throw new RuntimeException("Unknown mode "+b+"; valid modes are fix and tag.");}
			}else if(a.equals("tmp") || a.equals("tmpdir")){
				tmpDir=b;
			}else if(a.equals("in") || a.equals("in1")){
				parser.in1=b;
			}else if(a.equals("out") || a.equals("out1")){
				parser.out1=b;
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{//Process parser fields
			Parser.processQuality();
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			in1=parser.in1;
			out1=parser.out1;
		}

		if(in1==null){throw new RuntimeException("Error - an input file is required (in=).");}
		if(out1==null){throw new RuntimeException("Error - an output file is required (out=).");}
		assert(ways>0) : "ways must be positive: "+ways;

		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read input file "+in1+"\n");
		}
		if(!Tools.testOutputFiles(overwrite, false, false, out1)){
			throw new RuntimeException("\noverwrite="+overwrite+"; can't write to output file "+out1+"\n");
		}

		//Subfiles are headerless sam.gz; one in-memory subfile at a time bounds memory.
		if(tmpDir==null){tmpDir=out1+"_rbtmp";}

		ffin=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		assert(ffin.samOrBam()) : "Input must be sam or bam: "+in1;

		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method: partition by name-hash, then restore each subfile.
	 * @param t Timer for measuring execution time
	 */
	void process(Timer t){
		final File tdir=new File(tmpDir);
		tdir.mkdirs();
		assert(tdir.isDirectory()) : "Could not create temp directory "+tmpDir;

		final String[] subPaths=new String[ways];
		for(int i=0; i<ways; i++){subPaths[i]=tmpDir+"/part_"+i+".sam.gz";}

		//Phase 1: split input into WAYS headerless sam.gz subfiles by qname hash.
		partition(subPaths);

		//Phase 2: capture the input header and force SO:unsorted.
		final ArrayList<byte[]> header=patchHeaderSO(SamReadInputStream.getSharedHeader(true));

		//Phase 3: open the final output writer (header restored) and restore each subfile.
		final FileFormat ffout=FileFormat.testOutput(out1, FileFormat.SAM, null, true, overwrite, false, false);
		final SamWriterST2 out=new SamWriterST2(ffout, header, false);
		out.start();

		long outId=0;
		for(int i=0; i<ways; i++){
			outId=processSubfile(subPaths[i], out, outId);
		}
		errorState|=out.poisonAndWait();

		//Phase 4: clean up temp files.
		cleanup(subPaths, tdir);

		t.stop();
		printStats(t);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Streams the input once and routes every SamLine to subfile[qname.hashCode()%ways],
	 * buffering in chunks so all alignments of a read share a subfile.
	 * @param subPaths Paths of the headerless sam.gz subfiles to create
	 */
	private void partition(String[] subPaths){
		//Subfiles must be headerless; write with the Java codec (no per-subfile subprocess).
		final boolean oldNoHeader=ReadStreamWriter.NO_HEADER;
		ReadStreamWriter.NO_HEADER=true;
		final Writer[] subWriters=new Writer[ways];
		for(int i=0; i<ways; i++){
			FileFormat ff=FileFormat.testOutput(subPaths[i], FileFormat.SAM, null, false, overwrite, false, false);
			subWriters[i]=new SamWriterST2(ff, null, false);
			subWriters[i].start();
		}

		final Streamer cris=StreamerFactory.makeStreamer(ffin, 0, true, maxReads, true, false, -1);
		cris.start();

		@SuppressWarnings("unchecked")
		final ArrayList<SamLine>[] buffers=new ArrayList[ways];
		for(int i=0; i<ways; i++){buffers[i]=new ArrayList<SamLine>(CHUNK);}
		final long[] ids=new long[ways];

		for(ListNum<SamLine> ln=cris.nextLines(); ln!=null; ln=cris.nextLines()){
			if(ln.list==null){continue;}
			for(SamLine sl : ln.list){
				if(sl==null){continue;}
				readsIn++;
				basesIn+=(sl.seq==null ? 0 : sl.seq.length);
				final int part=(sl.qname==null ? 0 : (sl.qname.hashCode()&Integer.MAX_VALUE)%ways);
				ArrayList<SamLine> buf=buffers[part];
				buf.add(sl);
				if(buf.size()>=CHUNK){
					subWriters[part].addLines(new ListNum<SamLine>(buf, ids[part]++));
					buffers[part]=new ArrayList<SamLine>(CHUNK);
				}
			}
		}
		for(int i=0; i<ways; i++){
			if(!buffers[i].isEmpty()){
				subWriters[i].addLines(new ListNum<SamLine>(buffers[i], ids[i]++));
			}
		}
		errorState|=ReadWrite.closeStream(cris);
		for(Writer w : subWriters){errorState|=w.poisonAndWait();}
		ReadStreamWriter.NO_HEADER=oldNoHeader;
	}

	/**
	 * Loads one subfile, name-sorts it, restores bases on non-primary records, and writes
	 * the name-sorted records to the output stream in chunks of {@link #CHUNK}.
	 * @param path Subfile path
	 * @param out Final output writer
	 * @param outId Next list id to use for the output stream
	 * @return Updated output list id
	 */
	private long processSubfile(String path, SamWriterST2 out, long outId){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.SAM, null, true, true);
		final Streamer cris=StreamerFactory.makeStreamer(ff, 0, true, -1, false, false, -1);
		cris.start();
		final ArrayList<SamLine> all=new ArrayList<SamLine>();
		for(ListNum<SamLine> ln=cris.nextLines(); ln!=null; ln=cris.nextLines()){
			if(ln.list!=null){all.addAll(ln.list);}
		}
		errorState|=ReadWrite.closeStream(cris);

		if(all.isEmpty()){return outId;}

		Collections.sort(all, NAME_SORT);
		restoreGroups(all);

		for(int i=0; i<all.size(); i+=CHUNK){
			final int end=Math.min(i+CHUNK, all.size());
			final ArrayList<SamLine> chunk=new ArrayList<SamLine>(all.subList(i, end));
			out.addLines(new ListNum<SamLine>(chunk, outId++));
		}
		return outId;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Walks the name-sorted list, forming (qname, pairnum) groups, and restores each
	 * group's non-primary records from its primary record.
	 * @param all Name-sorted SamLines for one subfile
	 */
	private void restoreGroups(ArrayList<SamLine> all){
		final int n=all.size();
		int i=0;
		while(i<n){
			final SamLine primary=all.get(i);
			int j=i+1;
			while(j<n && sameGroup(primary, all.get(j))){j++;}
			//group is [i, j); primary is the rank-0 record (primary && !supplementary)
			final boolean usable=primary.primary() && !primary.supplementary() && primary.seq!=null;
			if(usable && hasHardClip(primary.cigar)){
				//primary itself lacks the full read (hard-clipped) - cannot restore this group
				skippedHardPrimary+=(j-i-1);
			}else if(usable){
				for(int k=i+1; k<j; k++){restoreOne(primary, all.get(k));}
			}else{
				groupsNoPrimary++;
			}
			i=j;
		}
	}

	/**
	 * Restores a single non-primary record from the primary, in either fix or tag mode.
	 * @param primary The primary (non-supplementary) record carrying the true bases
	 * @param np A non-primary (secondary or supplementary) record to restore
	 */
	private void restoreOne(SamLine primary, SamLine np){
		restoreAttempts++;
		if(np.cigar==null || !np.mapped()){skippedNoCigar++; return;}
		//Bases (soft+hard inclusive) the record's cigar accounts for must equal the full read.
		final int fullLen=SamLine.calcCigarBases(np.cigar, true, true);
		if(fullLen!=primary.seq.length){skippedLenMismatch++; return;}

		if(mode==TAG){
			np.addOptionalTag("OS:Z:"+new String(primary.seq));
			np.addOptionalTag("OQ:Z:"+qualToAscii(primary.qual));
			restored++;
			return;
		}

		//FIX mode: supplementary hard-clips become soft-clips so the full read is a valid SEQ.
		if(np.supplementary() && hasHardClip(np.cigar)){
			np.cigar=convertHardToSoft(np.cigar);
			hardConverted++;
		}
		assert(SamLine.calcCigarBases(np.cigar, true, false)==primary.seq.length) :
			"cigar query length "+SamLine.calcCigarBases(np.cigar, true, false)+
			" != read length "+primary.seq.length+" for cigar "+np.cigar;
		np.seq=primary.seq.clone();
		np.qual=(primary.qual==null ? null : primary.qual.clone());
		restored++;
	}

	/** True if both records share qname and pairnum.
	 * @param a First record
	 * @param b Second record
	 * @return True if same (qname, pairnum) group */
	private static boolean sameGroup(SamLine a, SamLine b){
		return a.pairnum()==b.pairnum() && compareQnames(a.qname, b.qname)==0;
	}

	/** True if the cigar contains a hard-clip operation.
	 * @param cigar CIGAR string
	 * @return True if an H operation is present */
	private static boolean hasHardClip(String cigar){
		return cigar!=null && cigar.indexOf('H')>=0;
	}

	/**
	 * Converts every hard-clip (H) operation to a soft-clip (S) and merges adjacent
	 * identical operations so the result is valid SAM.  After conversion the cigar's
	 * query-consuming length equals the full read length, matching the copied SEQ.
	 * @param cigar Original CIGAR string (must contain at least one H)
	 * @return Equivalent CIGAR with H replaced by S and adjacent ops merged
	 */
	static String convertHardToSoft(String cigar){
		if(cigar==null || cigar.indexOf('H')<0){return cigar;}
		final StringBuilder sb=new StringBuilder(cigar.length());
		final int n=cigar.length();
		int prevCount=0;
		char prevOp=0;
		int i=0;
		while(i<n){
			int count=0;
			while(i<n && cigar.charAt(i)>='0' && cigar.charAt(i)<='9'){
				count=count*10+(cigar.charAt(i)-'0');
				i++;
			}
			char op=cigar.charAt(i);
			i++;
			if(op=='H'){op='S';}
			if(op==prevOp){
				prevCount+=count;
			}else{
				if(prevOp!=0){sb.append(prevCount).append(prevOp);}
				prevOp=op;
				prevCount=count;
			}
		}
		if(prevOp!=0){sb.append(prevCount).append(prevOp);}
		return sb.toString();
	}

	/** Encodes phred qualities as ASCII (+33), or "*" if null.
	 * @param qual Phred quality scores (already decoded, 0-based)
	 * @return ASCII+33 quality string, or "*" */
	private static String qualToAscii(byte[] qual){
		if(qual==null){return "*";}
		final char[] c=new char[qual.length];
		for(int i=0; i<qual.length; i++){c[i]=(char)(qual[i]+33);}
		return new String(c);
	}

	/**
	 * Rewrites the SO: field of the @HD header line to SO:unsorted (adding @HD if absent),
	 * without mutating the shared header arrays.
	 * @param hdr Header lines from the input file
	 * @return A new header list with SO:unsorted
	 */
	private static ArrayList<byte[]> patchHeaderSO(ArrayList<byte[]> hdr){
		if(hdr==null){return null;}
		final ArrayList<byte[]> out=new ArrayList<byte[]>(hdr.size()+1);
		for(byte[] line : hdr){
			if(line!=null && Tools.startsWith(line, "@HD")){
				String s=new String(line);
				if(s.contains("SO:")){
					s=s.replaceAll("SO:[^\t]*", "SO:unsorted");
				}else{
					s=s+"\tSO:unsorted";
				}
				out.add(s.getBytes());
			}else{
				out.add(line);
			}
		}
		final boolean hasHD=!out.isEmpty() && Tools.startsWith(out.get(0), "@HD");
		if(!hasHD){out.add(0, "@HD\tVN:1.6\tSO:unsorted".getBytes());}
		return out;
	}

	/** Deletes the subfiles and the temp directory.
	 * @param subPaths Subfile paths
	 * @param tdir Temp directory */
	private void cleanup(String[] subPaths, File tdir){
		for(String p : subPaths){
			File f=new File(p);
			if(f.exists()){f.delete();}
		}
		tdir.delete();//only succeeds if empty
	}

	/** Prints a summary of the work performed, naming the tool that produced each number.
	 * @param t Stopped timer */
	private void printStats(Timer t){
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println();
		outstream.println("Mode:                  \t"+(mode==FIX ? "fix" : "tag"));
		outstream.println("Ways:                  \t"+ways);
		outstream.println("Alignments in:         \t"+readsIn);
		outstream.println("Restore attempts:      \t"+restoreAttempts);
		outstream.println("Records restored:      \t"+restored);
		outstream.println("Hardclips converted:   \t"+hardConverted);
		outstream.println("Skipped (hardclipped primary):\t"+skippedHardPrimary);
		outstream.println("Skipped (length mismatch):    \t"+skippedLenMismatch);
		outstream.println("Skipped (no cigar/unmapped):  \t"+skippedNoCigar);
		outstream.println("Groups with no primary:       \t"+groupsNoPrimary);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Comparator         ----------------*/
	/*--------------------------------------------------------------*/

	/** Null-safe qname comparison.
	 * @param a First qname
	 * @param b Second qname
	 * @return Comparison result */
	private static int compareQnames(String a, String b){
		if(a==null){return b==null ? 0 : -1;}
		if(b==null){return 1;}
		return a.compareTo(b);
	}

	/** Group-priority rank: primary-non-supplementary &lt; secondary &lt; supplementary.
	 * @param sl A SamLine
	 * @return 0 for primary-non-supplementary, 1 for secondary, 2 for supplementary */
	private static int rank(SamLine sl){
		if(sl.supplementary()){return 2;}
		if(!sl.primary()){return 1;}
		return 0;
	}

	/** Sort by qname, then pairnum, then group-priority rank (primary first). */
	private static final Comparator<SamLine> NAME_SORT=new Comparator<SamLine>(){
		@Override
		public int compare(SamLine a, SamLine b){
			int cmp=compareQnames(a.qname, b.qname);
			if(cmp!=0){return cmp;}
			cmp=a.pairnum()-b.pairnum();
			if(cmp!=0){return cmp;}
			return rank(a)-rank(b);
		}
	};

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String out1=null;
	private String tmpDir=null;

	/** Number of subfiles to partition into; larger values reduce per-subfile memory. */
	private int ways=31;
	/** Processing mode: FIX (write seq/qual) or TAG (attach OS:/OQ: tags). */
	private int mode=FIX;
	/** Maximum input alignments to process, or -1 for all. */
	private long maxReads=-1;
	private boolean overwrite=true;

	private final FileFormat ffin;

	/*--------------------------------------------------------------*/
	/*----------------           Counters           ----------------*/
	/*--------------------------------------------------------------*/

	private long readsIn=0;
	private long basesIn=0;
	private long restoreAttempts=0;
	private long restored=0;
	private long hardConverted=0;
	private long skippedHardPrimary=0;
	private long skippedLenMismatch=0;
	private long skippedNoCigar=0;
	private long groupsNoPrimary=0;

	/** True if an error was encountered during processing. */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Output/buffer chunk size for ListNum batches. */
	private static final int CHUNK=200;
	static final int FIX=0, TAG=1;

	/** Output stream for status messages. */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging. */
	public static boolean verbose=false;

}
