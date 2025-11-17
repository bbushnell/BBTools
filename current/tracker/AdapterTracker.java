package tracker;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Tools;
import stream.Read;
import structures.LongList;

/** Tracks counts of base types per position to make a consensus sequence.
 * Designed for use with BBMerge since adapters are inferred by insert size;
 * could be used with alignment too. */
public class AdapterTracker {
	
	/**
	 * Constructs an AdapterTracker with initialized count arrays.
	 * Creates LongList arrays for tracking A, C, G, T counts at each position
	 * for both read 1 and read 2 adapters.
	 */
	public AdapterTracker() {
		for(int i=0; i<counts.length; i++){
			for(int j=0; j<counts[i].length; j++){
				counts[i][j]=new LongList(151);
			}
		}
	}
	
	/**
	 * Stores adapter sequence information from a read with known insert size.
	 * Counts bases that extend beyond the insert size as potential adapter sequence.
	 * Filters out PhiX-like adapters if configured and increments position-specific
	 * base counts for consensus generation.
	 *
	 * @param r The read containing potential adapter sequence
	 * @param insert The insert size; bases beyond this position are considered adapters
	 */
	public void storeAdapterSequence(Read r, int insert){
		reads++;
		if(r.length()<=insert) {return;}
		shortInserts++;
		
		if(looksLikePhix(r, insert)){
			phixLike++;
			if(ignorePhixAdapters){return;}
		}
		
		LongList[] lists=counts[r.pairnum()];
		byte[] bases=r.bases;
		
		for(int i=insert, j=0; i<bases.length; i++, j++){
			byte b=bases[i];
			int num=AminoAcid.baseToNumber[b];
			if(num>=0){
				lists[num].increment(j);
			}
		}
	}
	
	/**
	 * Checks if a read pair contains PhiX-like adapter sequences.
	 * Tests both reads in the pair for PhiX adapter pattern.
	 *
	 * @param r The read to check (includes mate)
	 * @param insert The insert size
	 * @return true if either read contains PhiX-like adapter sequence
	 */
	private boolean looksLikePhix(Read r, int insert){
		return looksLikePhix(r.bases, insert) || looksLikePhix(r.mate.bases, insert);
	}
	
	/**
	 * Checks if a sequence contains PhiX adapter pattern starting at insert position.
	 * Compares sequence suffix against known PhiX prefix "AGATCGGAAGAGCG".
	 * Allows N bases to match any position in the pattern.
	 *
	 * @param bases The sequence bases to check
	 * @param insert The position where adapter sequence begins
	 * @return true if sequence matches PhiX adapter pattern
	 */
	private boolean looksLikePhix(byte[] bases, int insert){
		int len=bases.length-insert;
		if(len<phixPrefix.length){return false;}
		for(int i=insert, j=0; i<bases.length && j<phixPrefix.length; i++, j++){
			byte b=bases[i];
			if(b!='N' && b!=phixPrefix[j]){
				return false;
			}
		}
//		outstream.println(new String(bases).substring(insert));
//		outstream.println(new String(phixPrefix));
		return true;
	}
	
	/**
	 * Generates consensus adapter sequences from accumulated base counts.
	 * Creates adapter sequences for both read 1 and read 2 based on
	 * position-specific base count data with optional poly-A/G trimming.
	 * @return true if at least one valid adapter sequence was generated
	 */
	public boolean makeSequence() {
		seq1=seq2=null;
		seq1=toAdapterSequence(counts[0], trimPolyAorG);
		seq2=toAdapterSequence(counts[1], trimPolyAorG);
		return hasSequence();
	}
	
	/** Checks if valid adapter sequences have been generated.
	 * @return true if either seq1 or seq2 contains more than one base */
	public boolean hasSequence() {
		return (seq1!=null && seq1.length()>1) || (seq2!=null && seq2.length()>1);
	}
	
	/**
	 * Writes consensus adapter sequences to a FASTA file.
	 * Outputs Read1_adapter and Read2_adapter sequences based on accumulated counts.
	 * @param fname Output filename for the FASTA file
	 * @return Total count of bases used to generate the consensus sequences
	 */
	public long writeAdapterConsensus(String fname){
		StringBuilder sb=new StringBuilder();
		{
			sb.append(">Read1_adapter\n");
			String adapter=toAdapterSequence(counts[0], trimPolyAorG);
			sb.append(adapter).append('\n');
		}
		if(counts.length>1){
			sb.append(">Read2_adapter\n");
			String adapter=toAdapterSequence(counts[1], trimPolyAorG);
			sb.append(adapter).append('\n');
		}
		long count=counts[0][0].get(0)+counts[0][1].get(0)+
				counts[0][2].get(0)+counts[0][3].get(0);
//		outstream.println("Adapters counted: "+count);
		ReadWrite.writeString(sb, fname);
		return count;
	}
	
	/**
	 * Converts position-specific base count lists into consensus adapter sequence.
	 * Determines consensus base at each position using count thresholds and
	 * applies optional poly-A/G trimming and junk trimming.
	 *
	 * @param lists Arrays of base counts [A,C,G,T] for each position
	 * @param trimPolyAorG Whether to trim poly-A and poly-G sequences
	 * @return Consensus adapter sequence as a string
	 */
	private static String toAdapterSequence(LongList[] lists, boolean trimPolyAorG){
		StringBuilder adapter=new StringBuilder();
		long max=0;
		int lastBase=-1;
		for(int i=0; true; i++){
			long a=lists[0].get(i);
			long c=lists[1].get(i);
			long g=lists[2].get(i);
			long t=lists[3].get(i);
			long sum=(a+c+g+t);
			max=Tools.max(max, sum);
			if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
			long thresh=(max>100 ? 4+(sum*2)/3 : (sum*2)/3);
			if(a>thresh){
				adapter.append('A');
				lastBase=i;
			}else if(c>thresh){
				adapter.append('C');
				lastBase=i;
			}else if(g>thresh){
				adapter.append('G');
				lastBase=i;
			}else if(t>thresh){
				adapter.append('T');
				lastBase=i;
			}else{
				adapter.append('N');
			}
		}
		if(lastBase<0){return "N";}

		String trimmed=trimPoly2(adapter.toString(), 'N');
		if(trimPolyAorG){
			for(int len=-1; len!=trimmed.length(); ) {
				len=trimmed.length();
				trimmed=trimPoly2(trimmed, 'G');
				trimmed=trimPoly2(trimmed, 'A');
			}
		}
		if(trimJunk){
			trimmed=trimJunk(trimmed, 6);
		}
		
//		if(lastBase>=0){
//			char A=(trimPolyAorG ? 'A' : 'N');
//			while(lastBase>=0 && (adapter.charAt(lastBase)=='N' || adapter.charAt(lastBase)==A)){lastBase--;}
//		}
		
		if(trimmed.length()<1){return "N";}
		return trimmed;
	}
	
	/**
	 * Trims specified character and N's from the end of adapter sequence.
	 * Removes trailing occurrences of the specified base and ambiguous bases.
	 *
	 * @param adapter The adapter sequence to trim
	 * @param trim The character to trim from the end
	 * @return Trimmed adapter sequence
	 */
	private static String trimPoly(String adapter, char trim){
		int lastBase=-1;
		for(int i=0; i<adapter.length(); i++){
			char c=adapter.charAt(i);
			if(AminoAcid.isFullyDefined(c)){
				lastBase=i;
			}
		}
		
		int aCount=0;
		int nCount=0;
		int count=0;
		while(lastBase>=0){
			char c=adapter.charAt(lastBase);
			if(c=='N'){nCount++;}
			else if(c==trim){aCount++;}
			else{break;}
			count++;
			lastBase--;
		}
		
		if(lastBase<0){return "N";}
		if(count==nCount || (aCount>3)){
			return adapter.substring(0, lastBase+1);
		}
		return adapter;
	}
	
	/**
	 * Enhanced trimming method for homopolymer and N sequences.
	 * Removes trailing occurrences of specified character or N bases
	 * with improved logic for different trim thresholds.
	 *
	 * @param adapter The adapter sequence to trim
	 * @param poly The homopolymer character to trim
	 * @return Trimmed adapter sequence, or "N" if completely trimmed
	 */
	private static String trimPoly2(String adapter, char poly){
		int last=adapter.length()-1;
		int trim=0;
		while(last>=0) {
			char c=adapter.charAt(last);
//			System.err.println("c="+Character.toString(c)+", poly="+Character.toString(poly));
			if(c==poly || c=='N') {
				trim++;
				last--;
			}else{
				break;
			}
		}
		
		if(trim>3 || (trim>0 && poly=='N')) {
			adapter=adapter.substring(0, last+1);
		}
//		assert(poly=='N') : Character.toString(poly)+"\n"+adapter+"\n"+trim+"\n"+last;
		return adapter==null || adapter.length()<1 ? "N" : adapter;
	}
	
	/**
	 * Trims low-quality sequence from adapter end using scoring system.
	 * Uses scoring where N bases subtract 1 point and defined bases add 2 points.
	 * Trims until minimum score threshold is reached.
	 *
	 * @param s The sequence to trim
	 * @param minScore Minimum score required to retain sequence
	 * @return Trimmed sequence with junk removed
	 */
	private static String trimJunk(String s, int minScore) {
		int score=0, last=s.length()-1;
		for(; last>=0 && score<minScore; last--) {
			char c=s.charAt(last);
			if(c=='N') {
				score--;
			}else {
				score+=2;
			}
		}
		last++;
		while(last<s.length() && s.charAt(last)!='N' || (last<s.length()-1 && s.charAt(last+1)!='N')) {last++;}
		return (last<1 ? "N" : last>s.length() ? s : s.substring(0, last));
	}
	
	/**
	 * Merges another AdapterTracker's data into this one.
	 * Combines base count arrays and statistics from both trackers.
	 * @param b The AdapterTracker to merge into this one
	 */
	public void merge(AdapterTracker b){
		for(int x=0; x<counts.length; x++){
			for(int y=0; y<counts[x].length; y++){
				counts[x][y].incrementBy(b.counts[x][y]);
			}
		}
		reads+=b.reads;
		shortInserts+=b.shortInserts;
		phixLike+=b.phixLike;
	}
	
	/** Base count arrays [read1/2][A/C/G/T] for each position */
	final LongList[][] counts=new LongList[2][4];
	/** Consensus adapter sequence for read 1 */
	public String seq1=null;
	/** Consensus adapter sequence for read 2 */
	public String seq2=null;
	/** Total number of reads processed */
	public long reads=0;
	/** Number of reads with insert size shorter than read length */
	public long shortInserts=0;
	/** Number of reads containing PhiX-like adapter sequences */
	public long phixLike=0;
	
	/** PhiX adapter sequence prefix used for detection */
	private static final byte[] phixPrefix="AGATCGGAAGAGCG".getBytes();
	/** Whether to ignore PhiX-like adapters during processing */
	public static boolean ignorePhixAdapters=false;
	/** Whether to trim poly-A and poly-G sequences from adapter ends */
	public static boolean trimPolyAorG=true;
	/** Whether to apply junk trimming to remove low-quality adapter ends */
	public static boolean trimJunk=true;
	
}
