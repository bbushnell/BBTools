package prok;

import dna.GeneticCode;

/** One outward, uninterrupted terminal peptide. Coordinates are genomic, zero-based,
 * half-open; this endpoint is never a real-stop key. No annotation selects the query.
 * @author Qiqi
 */
public final class ContigEdgeFragment{

	ContigEdgeFragment(final Scope scope_, final String contig_, final GeneticCode code_, final int length,
			final int strand, final int frame, final int sense, final int boundary_, final byte[] query_, final int[] starts_){
		assert(scope_!=null && code_!=null && sense>=MINIMUM && query_.length==Math.min(sense, MAXIMUM)) :
			"Only an eligible terminal sense interval may own a fragment; no stop or padding is added";
		scope=scope_; contig=contig_; code=code_; query=query_; starts=starts_;
		final int end=frame+3*((length-frame)/3), begin=end-3*sense, queryBegin=end-3*query.length;
		key=new Key(scope, strand==0 ? length : 0, strand, frame);
		intervalStart=strand==0 ? begin : length-end;
		intervalEnd=strand==0 ? end : length-begin;
		queryStart=strand==0 ? queryBegin : length-end;
		queryEnd=strand==0 ? end : length-queryBegin;
		trailingBases=length-end; intervalBoundary=boundary_;
		assert(intervalStart>=0 && intervalEnd<=length && queryStart>=intervalStart && queryEnd<=intervalEnd &&
			queryEnd-queryStart==3*query.length && trailingBases>=0 && trailingBases<=2) :
			"Fragment coordinates must retain complete codons in the oriented frame at the physical boundary";
	}
	public byte[] copyQuery(){return query.clone();}
	/** Recognized start codons as oriented, zero-based base offsets, in translation order.
	 * These describe the full uninterrupted interval, including starts before a capped query. */
	public int[] copyStartOffsets(){return starts.clone();}
	public int length(){return query.length;}
	public final Scope scope;
	public final Key key;
	public final String contig;
	public final GeneticCode code;
	public final int intervalStart, intervalEnd, queryStart, queryEnd, trailingBases, intervalBoundary;
	private final byte[] query;
	private final int[] starts;
	public static final int MINIMUM=40, MAXIMUM=100;
	public static final String POLICY="contig_edge_3p_40_100_v1";
	public enum Topology{UNKNOWN, LINEAR, CIRCULAR}

	/** Identity and declared topology supplied by the input manifest, never inferred from misses.
	 * inputRecord identifies a record within the frozen manifest, not just a file basename. */
	public static final class Scope{
		public Scope(final String genome_, final String inputRecord_, final long contigOrdinal_, final Topology topology_){
			if(genome_==null || genome_.isEmpty() || inputRecord_==null || inputRecord_.isEmpty() ||
					contigOrdinal_<0 || topology_==null){
				throw new IllegalArgumentException("Edge search requires genome, input-record identity, "
					+"contig ordinal and declared topology");
			}
			genome=genome_; inputRecord=inputRecord_; contigOrdinal=contigOrdinal_; topology=topology_;
		}
		public final String genome, inputRecord;
		public final long contigOrdinal;
		public final Topology topology;
	}

	/** Value identity for deduplicating repeated endpoint processing. Code/query/policy belong
	 * to the associated evidence and must also match before reusing a cached decision. */
	public static final class Key{
		public Key(final Scope scope_, final int boundary_, final int strand_, final int frame_){
			if(scope_==null || boundary_<0 || strand_<0 || strand_>1 || frame_<0 || frame_>2 ||
					strand_==1 && boundary_!=0){
				throw new IllegalArgumentException("An outward endpoint requires a scoped physical boundary, strand "
					+"and oriented frame");
			}
			scope=scope_; boundary=boundary_; strand=strand_; frame=frame_;
		}
		@Override public boolean equals(final Object other){
			if(this==other){return true;}
			if(!(other instanceof Key)){return false;}
			final Key k=(Key)other;
			return boundary==k.boundary && strand==k.strand && frame==k.frame
				&& scope.contigOrdinal==k.scope.contigOrdinal && scope.genome.equals(k.scope.genome)
				&& scope.inputRecord.equals(k.scope.inputRecord);
		}
		@Override public int hashCode(){
			int h=scope.genome.hashCode(); h=31*h+scope.inputRecord.hashCode();
			h=31*h+(int)(scope.contigOrdinal^(scope.contigOrdinal>>>32));
			return 31*(31*(31*h+boundary)+strand)+frame;
		}
		public final Scope scope;
		public final int boundary, strand, frame;
		public static final String KIND="CONTIG_EDGE_3P";
	}
}
