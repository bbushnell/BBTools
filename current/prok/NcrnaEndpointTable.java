package prok;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import map.IntHashMap;
import prok.NcrnaTextIO.Output;
import prok.TrnaBoundaryFeatures.BoundaryType;
import prok.TrnaBoundaryFeatures.NinemerTable;
import structures.ByteBuilder;

/** Raw true-endpoint word counts with caller-compatible geometry and serialization.
 * Build worker-local tables, merge after joining, then share the result read-only.
 * Each vector worker owns a reusable leave-one-record-out view. That view removes
 * the record's contributed KEY at every matching position, including repeated
 * words away from the true coordinate. Validation uses the unmodified table.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointTable implements NinemerTable{

	/** Creates empty raw counts for one endpoint's sense-oriented window geometry. */
	public NcrnaEndpointTable(final BoundaryType type_, final int inside_, final int outside_){
		if(type_==null || inside_<0 || outside_<0 || (long)inside_+outside_<1 || (long)inside_+outside_>15){
			throw new IllegalArgumentException("Endpoint geometry needs a boundary type and one to fifteen total bases");
		}
		type=type_;
		inside=inside_;
		outside=outside_;
		k=inside+outside;
	}

	/** Adds one eligible record's true word, or counts its sequence limitation. */
	public void add(final byte[] bases, final int boundary){
		final int key=key(bases, windowStart(boundary));
		if(key==MISSING){missing++;}else if(key==AMBIGUOUS){ambiguous++;}else{
			increment(key, 1);
			counted++;
		}
	}

	/** Merges a completed worker table once; neither table may be mutated concurrently. */
	public void merge(final NcrnaEndpointTable other){
		if(other==null || other==this || other.type!=type || other.inside!=inside || other.outside!=outside){
			throw new IllegalArgumentException("Only distinct endpoint tables with identical geometry can merge");
		}
		final int[] keys=other.counts.keys(), values=other.counts.values();
		for(int i=0; i<keys.length; i++){
			if(keys[i]>=0){increment(keys[i], values[i]);}
		}
		counted+=other.counted;
		missing+=other.missing;
		ambiguous+=other.ambiguous;
	}

	/** Fails before the native map's saturating increment could lose count information. */
	private void increment(final int key, final int amount){
		assert(key>=0 && amount>0) : "Only packed valid words with positive counts belong in this table: key="+key+", amount="+amount;
		final int prior=counts.get(key);
		if((long)Math.max(0, prior)+amount>Integer.MAX_VALUE){
			throw new IllegalStateException("Endpoint word count exceeds the caller's integer range: "+key);
		}
		counts.increment(key, amount);
	}

	/** Returns a raw count; no smoothing or leave-one-out is applied here. */
	@Override
	public int rawCount(final byte[] bases, final int position){
		final int key=key(bases, position);
		return key<0 ? 0 : Math.max(0, counts.get(key));
	}

	/** Creates private mutable exclusion state over this completed table. */
	public LeaveOneOut leaveOneOut(){return new LeaveOneOut(this);}
	/** Number of eligible records that contributed a valid word. */
	public long counted(){return counted;}
	/** Number whose endpoint window lay outside the available context. */
	public long missing(){return missing;}
	/** Number whose available endpoint window contained an ambiguous base. */
	public long ambiguous(){return ambiguous;}

	/** Writes the existing inference loader's native format, honoring compression suffixes. */
	public void write(final String path){
		try(Output owner=new Output(path)){
			final ByteStreamWriter out=owner.writer;
			final ByteBuilder b=new ByteBuilder();
			b.append("#boundary\t").append(type==BoundaryType.START ? "start" : "stop").nl()
				.append("#insidecount\t").append(inside).nl().append("#outsidecount\t").append(outside).nl()
				.append("#k\t").append(k).append("\n#kmer\tcount\n");
			out.print(b);
			final int[] keys=counts.toKeyArray();
			Arrays.sort(keys);
			for(int key:keys){
				b.clear();
				for(int base=k-1; base>=0; base--){b.append(AminoAcid.numberToBase[(key>>>(2*base))&3]);}
				out.print(b.tab().append(counts.get(key)).nl());
			}
		}
	}

	/** Uses the same sense-oriented start/stop window calculation as inference. */
	private int windowStart(final int boundary){
		assert(boundary>=0) : "A true endpoint must index the sense-oriented context: boundary="+boundary;
		return TrnaBoundaryFeatures.kmerWindowStart(boundary, type, inside, outside);
	}

	/** Two-bit key; invalid context and ambiguous sequence are distinct counting outcomes. */
	private int key(final byte[] bases, final int position){
		if(bases==null || position<0 || position>bases.length-k){return MISSING;}
		int key=0;
		for(int i=position, end=position+k; i<end; i++){
			final int base=bases[i];
			final int value=base<0 || base>=AminoAcid.baseToNumber.length ? -1 : AminoAcid.baseToNumber[base];
			if(value<0){return AMBIGUOUS;}
			key=(key<<2)|value;
		}
		return key;
	}

	/** Mutable worker-local view over a completed, immutable parent table. */
	public static final class LeaveOneOut implements NinemerTable{
		private LeaveOneOut(final NcrnaEndpointTable table_){table=table_;}

		/** Binds one record that contributed to the parent, or an invalid word that contributed nothing. */
		public void set(final byte[] bases, final int trueBoundary){
			excludedKey=table.key(bases, table.windowStart(trueBoundary));
			if(excludedKey>=0 && table.counts.get(excludedKey)<1){
				throw new IllegalArgumentException("Leave-one-out record did not contribute to this endpoint table");
			}
		}

		@Override
		public int rawCount(final byte[] bases, final int position){
			final int key=table.key(bases, position);
			if(key<0){return 0;}
			final int raw=Math.max(0, table.counts.get(key));
			return raw-(key==excludedKey ? 1 : 0);
		}

		/** Borrowed completed counts; this view never changes the parent table. */
		private final NcrnaEndpointTable table;
		/** True endpoint word to subtract at every matching lookup, or a missing sentinel. */
		private int excludedKey=MISSING;
	}

	/** Transcriptional endpoint whose window geometry defines the counted word. */
	public final BoundaryType type;
	/** Number of bases inside and outside the true endpoint. */
	public final int inside, outside;
	/** Total word length, inside+outside, bounded to fifteen bases. */
	public final int k;
	/** Exact raw counts keyed by two-bit forward words; absent keys contribute zero. */
	private final IntHashMap counts=new IntHashMap();
	/** Mutually exclusive counts of valid, unavailable and ambiguous endpoint words. */
	private long counted, missing, ambiguous;
	/** Distinguishes unavailable geometry from sequence ambiguity during population accounting. */
	private static final int MISSING=-1, AMBIGUOUS=-2;
}
