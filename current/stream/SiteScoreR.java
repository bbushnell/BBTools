package stream;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import parse.LineParserS2;
import shared.Shared;

/**
 * Represents a scored alignment site for a read with additional read-specific metadata.
 * Copies selected SiteScore values into a separate type, adding read identity and sorting metadata.
 * Natural equality uses score, paired score, chromosome, strand and start; positional
 * equality also checks stop and ignores scores. The overloads compare different fields.
 * Mutable scores/coordinates can change ordering and equality. Do not use this class as
 * a key in hash-based collections, or mutate comparison fields while sorting.
 * Constructor initialization makes perfect imply semiperfect; callers must maintain
 * that relationship when changing the public flags. Text output omits some working fields.
 * @author Brian Bushnell
 * @date Jul 16, 2012
 */
public final class SiteScoreR implements Comparable<SiteScoreR>{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies alignment fields and supplies read-specific identity.
	 * @param ss Nonnull source site; values are copied, not retained as a shared SiteScore
	 * @param readlen_ Read length
	 * @param numericID_ Numeric read identity
	 * @param pairnum_ Pair member identifier
	 */
	public SiteScoreR(SiteScore ss, int readlen_, long numericID_, byte pairnum_){
		this(ss.chrom, ss.strand, ss.start, ss.stop, readlen_, numericID_, pairnum_, ss.score, ss.pairedScore, ss.perfect, ss.semiperfect);
	}
	
	/** Initializes a scored inclusive alignment interval; perfect also sets semiperfect.
	 * @param chrom_ Reference chromosome
	 * @param strand_ Strand identifier
	 * @param start_ Inclusive start, no greater than stop_
	 * @param stop_ Inclusive stop
	 * @param readlen_ Read length
	 * @param numericID_ Numeric read identity
	 * @param pairnum_ Pair member identifier
	 * @param score_ Alignment score
	 * @param pscore_ Paired alignment score
	 * @param perfect_ Perfect-alignment flag
	 * @param semiperfect_ Semiperfect flag, promoted to true when perfect_ is true
	 */
	public SiteScoreR(int chrom_, byte strand_, int start_, int stop_, int readlen_, long numericID_, byte pairnum_, int score_, int pscore_, boolean perfect_, boolean semiperfect_){
		chrom=chrom_;
		strand=strand_;
		start=start_;
		stop=stop_;
		readlen=readlen_;
		numericID=numericID_;
		pairnum=pairnum_;
		score=score_;
		pairedScore=pscore_;
		perfect=perfect_;
		semiperfect=semiperfect_|perfect_;
		assert(start_<=stop_) : this.toText();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Orders by descending score/pairedScore, then ascending chrom/strand/start.
	 * Stop, read identity and other fields do not break ties. Integer subtraction
	 * assumes that the compared values have representable differences.
	 * @param other Nonnull site to compare
	 * @return Negative, zero or positive according to this ordering
	 */
	@Override
	public int compareTo(SiteScoreR other){
		//Descending scores then ascending coordinates; callers must keep integer differences representable.
		int x=other.score-score;
		if(x!=0){return x;}
		
		x=other.pairedScore-pairedScore;
		if(x!=0){return x;}
		
		x=chrom-other.chrom;
		if(x!=0){return x;}
		
		x=strand-other.strand;
		if(x!=0){return x;}
		
		x=start-other.start;
		return x;
	}
	
	/** Tests natural-order equality for another SiteScoreR; null and other types return false.
	 * This does not establish read identity or compare stop/perfection flags. */
	@Override
	public boolean equals(Object other){
		//[stream/SiteScoreR#001] guard null/wrong-type → honor equals contract (false, not NPE/CCE). Twin of SiteScore#001 (FIXED). hashCode asserts-false so no hash-collection use; reachability LOW.
		return other instanceof SiteScoreR && compareTo((SiteScoreR)other)==0;
	}
	
	/**
	 * Rejects hash use when assertions are enabled; do not use this type in hash-based collections.
	 * With assertions disabled the fallback is Object's identity hash, which is not
	 * a hash of the fields used by equals.
	 * @return Identity hash only when assertions are disabled
	 * @throws AssertionError When assertions are enabled
	 */
	@Override
	public int hashCode(){
		assert(false) : "This class should not be hashed.";
		return super.hashCode();
	}
	
	/** Compares chrom/strand/start/stop only, ignoring scores and read identity.
	 * @param other Nonnull base SiteScore
	 * @return true when all four positional fields match
	 */
	public boolean equals(SiteScore other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom){return false;}
		if(other.strand!=strand){return false;}
		return true;
	}
	
	/** Typed natural-order equality; unlike equals(Object), this overload requires a nonnull argument.
	 * @param other Nonnull SiteScoreR
	 * @return true when compareTo returns zero
	 */
	public boolean equals(SiteScoreR other){
		return compareTo(other)==0;
	}
	
	/** Returns string representation of this alignment site by delegating to toText().
	 * @return Comma-separated string of alignment data */
	@Override
	public String toString(){return toText().toString();}
	
	/** Serializes ten comma-separated fields, optionally prefixed with * when correct is true.
	 * Order: chrom, strand, start, stop, readlen, numericID, pairnum, two binary flag
	 * digits (semiperfect then perfect), pairedScore and score. normalizedScore and
	 * retainVotes are not included. The flags are written as currently stored.
	 * @return New builder containing this record without a line terminator
	 */
	public StringBuilder toText(){
		StringBuilder sb=new StringBuilder(50);
		if(correct){sb.append('*');}
		sb.append(chrom);
		sb.append(',');
		sb.append(strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append(readlen);
		sb.append(',');
		sb.append(numericID);
		sb.append(',');
		sb.append(pairnum);
		sb.append(',');
		//Two digits, no comma: [semiperfect][perfect]. Constructor-produced pairs are00/10/11.
		//Public flag changes can produce01; fromText reconstructs through the constructor and promotes it to11.
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
		return sb;
	}
	
	/** Tests inclusive overlap on the same chromosome and strand.
	 * @param ss Nonnull site with ordered endpoints
	 * @return true for an intersecting interval, including a shared endpoint
	 */
	public final boolean overlaps(SiteScoreR ss){
		return chrom==ss.chrom && strand==ss.strand && overlap(start, stop, ss.start, ss.stop);
	}
	/** Tests inclusive overlap of two ordered intervals. */
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/** Returns the ten labels in the order written by toText.
	 * @return Comma-separated labels without a line terminator
	 */
	public static String header(){
		//[stream/SiteScoreR#004 FIXED] Omit obsolete quickScore/slowScore labels; records contain ten fields.
		return "chrom,strand,start,stop,readlen,numericID,pairnum,semiperfect+perfect,pairedScore,score";
	}
	
	/** Parses the ten fields emitted by toText and restores an optional leading * marker.
	 * A legacy eleventh field is accepted but ignored. Constructor flag normalization
	 * applies; normalizedScore and retainVotes keep their new-object defaults.
	 * Uses a local sequential cursor and Java numeric conversions; trailing empty
	 * comma fields do not count toward the ten/eleven-field assertion.
	 * @param s Nonnull comma-separated record in this format
	 * @return New parsed site
	 */
	public static SiteScoreR fromText(String s){
		return fromText(s, new LineParserS2(','));
	}

	/** Parses one record after resetting the caller-owned comma cursor.
	 * @param s Nonnull comma-separated record
	 * @param lp Local comma cursor, reused only sequentially by this caller
	 * @return New parsed site
	 */
	private static SiteScoreR fromText(String s, LineParserS2 lp){
		lp.set(s);
		String first=lp.parseString();
		boolean correct=false;
		if(first.charAt(0)=='*'){
			correct=true;
			first=first.substring(1);
		}
		//Historical fix [stream/SiteScoreR#002]: parse chrom as int to match toText;
		//the former byte parser rejected values above127. Twin of SiteScore#002.
		int chrom=Integer.parseInt(first);
		byte strand=Byte.parseByte(lp.parseString());
		int start=Integer.parseInt(lp.parseString());
		int stop=Integer.parseInt(lp.parseString());
		int readlen=Integer.parseInt(lp.parseString());
		long numericID=Long.parseLong(lp.parseString());
		byte pairnum=Byte.parseByte(lp.parseString());
		int p=Integer.parseInt(lp.parseString(), 2);
		boolean perfect=(p&1)==1;
		boolean semiperfect=(p&2)==2;
		int pairedScore=Integer.parseInt(lp.parseString());
		int score=Integer.parseInt(lp.parseString());
		int terms=10;
		for(int field=10; lp.hasMore(); field++){
			if(lp.advance()>0){terms=field+1;}
		}
		assert(terms==10 || terms==11) : "SiteScoreR uses ten fields and accepts one ignored legacy field; found "+terms+": "+s;
		SiteScoreR ss=new SiteScoreR(chrom, strand, start, stop, readlen, numericID, pairnum, score, pairedScore, perfect, semiperfect);
		ss.correct=correct;
		
		return ss;
	}
	
	/** Parses tab-separated site records; trailing empty tab fields are omitted.
	 * Counts retained fields before allocating the result and reuses local cursors.
	 * @param s Nonnull text containing records accepted by fromText
	 * @return New array in input order
	 */
	public static SiteScoreR[] fromTextArray(String s){
		assert(s!=null) : "SiteScoreR.fromTextArray requires nonnull text; each retained tab field is passed to fromText";
		if(s.indexOf('\t')<0){return new SiteScoreR[] {fromText(s)};}
		final LineParserS2 records=new LineParserS2('\t').set(s);
		int terms=0;
		for(int field=0; records.hasMore(); field++){
			if(records.advance()>0){terms=field+1;}
		}
		final SiteScoreR[] out=new SiteScoreR[terms];
		if(terms==0){return out;}
		records.reset();
		final LineParserS2 fields=new LineParserS2(',');
		for(int i=0; i<terms; i++){out[i]=fromText(records.parseString(), fields);}
		return out;
	}
	
	/** Compares chrom/strand/start/stop, ignoring scores, flags and read identity.
	 * @param b Nonnull site
	 * @return true when all four positional fields match
	 */
	public boolean positionalMatch(SiteScoreR b){
		if(chrom!=b.chrom || strand!=b.strand || start!=b.start || stop!=b.stop){
			return false;
		}
		return true;
	}
	
	/** Ascending chrom/start/stop/strand, descending score, then perfect before nonperfect.
	 * Compared sites must be nonnull with representable integer differences. */
	public static class PositionComparator implements Comparator<SiteScoreR>{
		
		/** Used by the shared PCOMP instance. */
		private PositionComparator(){}
		
		/** Compares two nonnull sites in this comparator's documented order. */
		@Override
		public int compare(SiteScoreR a, SiteScoreR b){
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			return 0;
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	/** Descending int-cast normalizedScore, score, pairedScore and retainVotes, then
	 * perfect first and ascending chrom/start/stop/strand. Fractional normalized scores
	 * do not break ties. Compared values must have representable integer differences. */
	public static class NormalizedComparator implements Comparator<SiteScoreR>{
		
		/** Used by the shared NCOMP instance. */
		private NormalizedComparator(){}
		
		/** Compares two nonnull sites in this comparator's documented order. */
		@Override
		public int compare(SiteScoreR a, SiteScoreR b){
			if((int)a.normalizedScore!=(int)b.normalizedScore){return (int)b.normalizedScore-(int)a.normalizedScore;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.pairedScore!=b.pairedScore){return b.pairedScore-a.pairedScore;}
			if(a.retainVotes!=b.retainVotes){return b.retainVotes-a.retainVotes;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			return 0;
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	/** Ascending numericID/pairnum/chrom/start/stop/strand, descending score, then perfect first.
	 * Numeric IDs use direct relational comparison; integer differences must be representable. */
	public static class IDComparator implements Comparator<SiteScoreR>{
		
		/** Used by the shared IDCOMP instance. */
		private IDComparator(){}
		
		/**
		 * Compares two SiteScoreR objects by read identity first, then genomic position.
		 * [stream/SiteScoreR#003 DOC] Primary sort is numericID, then pairnum, then chromosome, start, stop, strand, score (descending), perfect — NOT position-primary (the old javadoc mis-described it as PositionComparator's order).
		 * @param a First SiteScoreR to compare
		 * @param b Second SiteScoreR to compare
		 * @return Negative, zero or positive according to the ordering above
		 */
		@Override
		public int compare(SiteScoreR a, SiteScoreR b){
			//Direct long comparison avoids subtraction overflow or narrowing an ID difference to int.
			if(a.numericID!=b.numericID){return a.numericID>b.numericID ? 1 : -1;}
			if(a.pairnum!=b.pairnum){return a.pairnum-b.pairnum;}
			
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			return 0;
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(ArrayList<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Shared.sort(list, this);
		}
		
		/** Sorts in place in this comparator's order; null or fewer than two entries is a no-op.
		 * Compared entries must be nonnull and must not be modified during sorting. */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}

	/** Shared positional comparator; it holds no per-sort scratch. */
	public static final PositionComparator PCOMP=new PositionComparator();
	/** Shared normalized-score comparator; fractional scores are truncated for ordering. */
	public static final NormalizedComparator NCOMP=new NormalizedComparator();
	/** Shared read-identity comparator. */
	public static final IDComparator IDCOMP=new IDComparator();
	
	/** Returns the inclusive interval length stop-start+1. */
	public int reflen(){return stop-start+1;}
	
	/** Inclusive reference start. */
	public int start;
	/** Inclusive reference stop. */
	public int stop;
	/** Read length, independent of the aligned reference span. */
	public int readlen;
	/** Alignment score; used by natural equality and ordering. */
	public int score;
	/** Paired alignment score; used by natural equality and ordering. */
	public int pairedScore;
	/** Reference chromosome. */
	public final int chrom;
	/** Strand identifier. */
	public final byte strand;
	/** Perfect flag; keep semiperfect true when setting this true. */
	public boolean perfect;
	/** Semiperfect flag, promoted by constructor when perfect is true. */
	public boolean semiperfect;
	/** Numeric read identity; ignored by natural equality. */
	public final long numericID;
	/** Read pair member; ignored by natural equality. */
	public final byte pairnum;
	/** Working score; NCOMP compares its int-cast value, and toText omits it. */
	public float normalizedScore;
	/** Correctness marker serialized as an optional leading *. */
	public boolean correct=false;
	/** Working retention vote count; omitted from toText. */
	public int retainVotes=0;
	
}
