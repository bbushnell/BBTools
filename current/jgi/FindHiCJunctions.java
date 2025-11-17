package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongPair;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class FindHiCJunctions {
	
	/** Program entry point for HiC junction analysis.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		FindHiCJunctions x=new FindHiCJunctions(args);

		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs FindHiCJunctions with command-line argument parsing.
	 * Initializes parameters, file formats, and k-mer counting arrays.
	 * Sets up counts arrays for k-mer sizes 4, 6, 8, 10 and junction
	 * position analysis arrays for left/right flanking regions.
	 *
	 * @param args Command-line arguments for configuration
	 */
	public FindHiCJunctions(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBuffers(4); //Only for singlethreaded programs
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("minclip")){
				minClipLength=Integer.parseInt(b);
			}else if(a.equals("printkmers") || a.equals("printcounts")){
				printKmers=Parse.parseBoolean(b);
			}else if(a.equals("junctionfile") || a.equals("junctions") || a.equals("outk")){
				junctionFile=b;
				assert(b==null || b.indexOf('%')>0);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}
			
			else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
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
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		counts=new long[11][];
		counts[10]=new long[1024*1024];
		counts[8]=new long[256*256];
		counts[6]=new long[4096];
		counts[4]=new long[256];
		
		leftCounts=new long[6][];
		rightCounts=new long[6][];
		leftCounts[5]=new long[1024];
		leftCounts[4]=new long[256];
		leftCounts[3]=new long[64];
		leftCounts[2]=new long[16];
		rightCounts[5]=new long[1024];
		rightCounts[4]=new long[256];
		rightCounts[3]=new long[64];
		rightCounts[2]=new long[16];
	}
	
	/**
	 * Main processing pipeline for HiC junction detection.
	 * Processes mapped reads to identify junction sites, extracts k-mers
	 * around junctions, and optionally trims reads at junction positions.
	 * Outputs junction k-mer statistics and processed reads.
	 *
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, null, null, buff, null, false);
			ros.start();
		}else{ros=null;}

		long readsProcessed=0, readsOut=0;
		long basesProcessed=0, basesOut=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					final boolean keep=process(r1, r1.samline);
					if(keep || true){
						listOut.add(r1);
						readsOut++;
						basesOut+=r1.length();
					}
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);

		if(printKmers){
			printKmers(10, counts[10], "");
			printKmers(8, counts[8], "");
			printKmers(6, counts[6], "");
			printKmers(4, counts[4], "");

			printKmers(5, leftCounts[5], "L");
			printKmers(4, leftCounts[4], "L");
			printKmers(3, leftCounts[3], "L");

			printKmers(5, rightCounts[5], "R");
			printKmers(4, rightCounts[4], "R");
			printKmers(3, rightCounts[3], "R");
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Outputs k-mer counts to junction file in FASTA or TSV format.
	 * Filters k-mers by count threshold and frequency, then writes
	 * sorted results with sequence and abundance information.
	 *
	 * @param k K-mer length
	 * @param array Count array for this k-mer size
	 * @param direction Direction suffix for filename (L/R or empty)
	 */
	private void printKmers(int k, long[] array, String direction){
		boolean tsv=junctionFile.endsWith(".tsv") || junctionFile.endsWith(".tsv.gz");
		final String fname=junctionFile.replaceFirst("%", k+direction);
		final ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		final long sum=shared.Vector.sum(array);
		final double mult=1.0/(Tools.max(1, sum));
		final long thresh=Tools.max((2*sum)/array.length, (long)Math.ceil(Tools.max(minCount, minFraction*sum)));
//		if(k==4){System.err.println(Arrays.toString(array));}
//		
////		assert(k!=4) : thresh+", "+sum;
		
		ArrayList<LongPair> list=new ArrayList<LongPair>();
		for(int kmer=0; kmer<array.length; kmer++){
			final long count=array[kmer];
			if(count>=thresh){
				list.add(new LongPair(count, kmer));
			}
		}
		Collections.sort(list);
		Collections.reverse(list);
		

		ByteBuilder bb=new ByteBuilder(4200);
		for(LongPair pair : list){
			if(tsv){
				bb.append(AminoAcid.kmerToString(pair.b, k)).append('\t');
				bb.append(pair.a).append('\t').append(mult*pair.a, 5).append('\n');
			}else{
				bb.append('>').append(pair.a).append('\t').append(mult*pair.a, 5).append('\n');
				bb.append(AminoAcid.kmerToString(pair.b, k)).append('\n');
			}
			if(bb.length()>4096){
				bsw.print(bb);
				bb.clear();
			}
		}

		if(!bb.isEmpty()){
			bsw.print(bb);
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes individual read for junction detection and k-mer extraction.
	 * Identifies junction position from clipping patterns, extracts k-mers
	 * around the junction site, and optionally trims the read.
	 *
	 * @param r The read to process
	 * @param sl SAM line information for the read
	 * @return true if read should be kept in output, false otherwise
	 */
	private boolean process(Read r, SamLine sl){
		
		if(sl==null || !sl.mapped() || !sl.primary() || sl.supplementary()|| r.match==null || !r.containsNonNM()){
			return false;
		}
		
		boolean definiteJunction=(sl.nextMapped() && !sl.properPair());
		
//		final String cigar=sl.cigar;
		final byte[] bases=r.bases;
		final boolean rcomp=(r.strand()==Shared.MINUS);
		final int left, right;
		{
			if(r.shortmatch()){r.toLongMatchString(true);}
			int left0=SamLine.countLeadingClip(r.match);
			int right0=SamLine.countTrailingClip(r.match);
			if(left0==0 && right0==0){
				byte[] smatch=softClipMatch(r.match, minClipLength, true);
				left0=SamLine.countLeadingClip(smatch);
				right0=SamLine.countTrailingClip(smatch);
			}
			left=left0;
			right=right0;
		}
		
		if((left>1 && right>1) || (left<minClipLength && right<minClipLength)){return false;}

		final int pos=(left>right ? left-1 : bases.length-right-1);
		if(printKmers && definiteJunction && junctionFile!=null){
			if(rcomp){r.reverseComplement();}

			//Pos: position of base to left of junction
			//		System.err.println(left+", "+right+", "+pos+", "+rcomp);
			//		System.err.println(new String(bases).substring(0, pos+1)+"~"+new String(bases).substring(pos+1));
			for(int k=10, half=5, start=pos-4; k>=4; k-=2, half--, start++){
				int kmer=0;
				assert(start>=0) : left+", "+right+", "+pos+", "+start+", "+k+", "+bases.length+"\n"+new String(r.match);
				for(int i=start, len=0; len<k; i++, len++){
					final byte b=bases[i];
					final int num=AminoAcid.baseToNumber[b];
					kmer=(kmer<<2)|num;
				}
				if(kmer<0){return false;}
				//			System.err.println("Adding "+AminoAcid.kmerToString(kmer, k));
				counts[k][kmer]++;
				leftCounts[half][kmer>>k]++;
				rightCounts[half][kmer&~((-1)<<k)]++;
			}

			if(rcomp){r.reverseComplement();}
		}
		
		if(trim){
			int trimAmount=bases.length-pos-1;
			int remainingAmount=bases.length-trimAmount;
			if(remainingAmount>=minTrimLength){
				if(rcomp){
					TrimRead.trimByAmount(r, trimAmount, 0, 1, false);
				}else{
					TrimRead.trimByAmount(r, 0, trimAmount, 1, false);
				}
			}else if(trimAmount>=minTrimLength){
				if(rcomp){
					TrimRead.trimByAmount(r, 0, remainingAmount, 1, false);
				}else{
					TrimRead.trimByAmount(r, remainingAmount, 0, 1, false);
				}
			}else{
				//do nothing
			}
		}
		
		return true;
	}
	
	/**
	 * Converts alignment match string to soft-clipped version based on scoring.
	 * Uses dynamic programming approach to find optimal alignment region
	 * and converts flanking low-scoring regions to soft clips.
	 *
	 * @param match Original match string from alignment
	 * @param minClipLength Minimum bases required for clipping
	 * @param allowMutation Whether to modify the input array directly
	 * @return Modified match string with soft clips applied
	 */
	public static byte[] softClipMatch(byte[] match, int minClipLength, boolean allowMutation){

		final int matchScore=100;
		final int subScore=-200;
		final int subScore2=-100;
		final int insScore=-200;
		final int delScore=-200;
		final int delScore2=-10;
		final int clipScore=-1;
		final int nScore=1;

		int insCount=0;
		int delCount=0;
		
		long score=0;
		long maxScore=0;
		int maxPos=-1;
		int maxStart=-1;
		int currentStart=-1;
		byte current='?';
		
		for(int mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
//			long prevScore=score;
			
			if(m=='m' || m=='N' || m=='R'){
				if(score==0){currentStart=mpos;}
				
				score=score+(m=='m' ? matchScore : nScore);
				
				if(score>maxScore){
					maxScore=score;
					maxPos=mpos;
					maxStart=currentStart;
				}
			}else{
				if(m=='S' || m=='s'){
					score=score+(m==current ? subScore2 : subScore);
				}else if(m=='D'){
					score=score+(m==current ? delScore2 : delScore);
					delCount++;
				}else if(m=='I' || m=='X' || m=='Y'){
					score=score+insScore;
					insCount++;
				}else if(m=='C'){
					score=score+clipScore;
				}
				score=Tools.max(0, score);
			}
			current=m;
		}
		
		if(maxScore<1){return match;}
		final int leftClipM=maxStart;
		final int rightClipM=(match.length-maxPos-1);
		int leftClip=0, rightClip=0;
		for(int i=0; i<match.length; i++){
			byte m=match[i];
			if(i<maxStart){
				leftClip+=(m=='D' ? 0 : 1);
			}else if(i>maxPos){
				rightClip+=(m=='D' ? 0 : 1);
			}
		}
		if(leftClip<minClipLength && rightClip<minClipLength){return match;}
		if(delCount==0){
			final byte[] array=allowMutation ? match : match.clone();
			for(int i=0; i<leftClip; i++){array[i]='C';}
			for(int i=0, j=array.length-1; i<rightClip; i++, j--){array[j]='C';}
			return array;
		}
		
		ByteBuilder bb=new ByteBuilder(match.length);
		if(leftClip>=minClipLength){
			for(int mpos=0, processed=0; mpos<match.length; mpos++){
				byte m=match[mpos];
				if(m=='D'){
					if(mpos>=leftClipM){bb.append(m);}
				}else{
					bb.append(mpos<leftClipM ? (byte)'C' : m);
					processed++;
				}
			}
		}else{
			bb.append(match);
		}
		if(rightClip>=minClipLength){
			bb.reverseInPlace();
			byte[] temp=bb.toBytes();
			bb.clear();
			for(int mpos=0, processed=0; mpos<temp.length; mpos++){
				byte m=temp[mpos];
				if(m=='D'){
					if(mpos>=rightClipM){bb.append(m);}
				}else{
					bb.append(mpos<rightClipM ? (byte)'C' : m);
					processed++;
				}
			}
			bb.reverseInPlace();
		}
//		System.out.println(new String(match)+"\n"+bb.toString()+"\n");
//		System.out.println("length="+match.length);
//		System.out.println("maxScore="+maxScore);
//		System.out.println("maxPos="+maxPos);
//		System.out.println("maxStart="+maxStart);
//		System.out.println("leftClip="+leftClip);
//		System.out.println("rightClip="+rightClip);
//		System.out.println("leftClipM="+leftClipM);
//		System.out.println("rightClipM="+rightClipM);
//		System.out.println();
		
		return bb.toBytes();
	}
	
	/*--------------------------------------------------------------*/
	
	/** Input file path for reads to process */
	private String in1=null;

	/** Output file path for processed reads */
	private String out1=null;
	
	/** Input file extension override */
	private String extin=null;
	/** Output file extension override */
	private String extout=null;
	
	/**
	 * Output file pattern for junction k-mer counts (% replaced with k-mer size)
	 */
	private String junctionFile="junctions_k%.txt";
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Minimum clip length required to consider a junction */
	private int minClipLength=8;
	/** Minimum read length to retain after trimming */
	private int minTrimLength=25;
	/** Minimum count threshold for k-mer output */
	private int minCount=2;
	/** Minimum fraction of total k-mers for output threshold */
	private float minFraction=0.0005f;

	/** Whether to output k-mer statistics to junction files */
	boolean printKmers=true;
	/** Whether to trim reads at junction positions */
	boolean trim=true;
	/** K-mer count arrays for different k-mer lengths (4, 6, 8, 10) */
	private long[][] counts;
	/** Count arrays for left flanking k-mers at junction sites */
	private long[][] leftCounts;
	/** Count arrays for right flanking k-mers at junction sites */
	private long[][] rightCounts;
	
	/*--------------------------------------------------------------*/
	
	/** Input file format specification */
	private final FileFormat ffin1;

	/** Output file format specification */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and status information */
	public static boolean verbose=false;
	/** Tracks whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
