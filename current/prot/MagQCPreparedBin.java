package prot;

import java.util.Arrays;

import parse.LineParser1;

/**
 * Reusable whole-bin sufficient statistics for the frozen MAG-QC input layout.
 * Counts are raw integers, never recomputed normalization or truth labels.
 * @author Nilou
 */
final class MagQCPreparedBin {

	static final String SCHEMA="magqc_prepared_bin_v1";
	static final String COLUMNS="bin_id\tqc_status\tqc_domain\tqc_phylum\tlength_bp\tgc_bases\tacgt_bases"
		+"\tcoding_bp\tcds\tmapped_cds\tgene_length_sum\tgene_length_sq_sum\tr16\tr23\tr5\trother"
		+"\ttrna_total\tfamily_counts\tanticodon_counts\tdimers";

	/** Allocates scratch once at the resource-bound family width. */
	MagQCPreparedBin(int familyCount){
		if(familyCount<1){throw new IllegalArgumentException("Prepared bins require a positive bound family width");}
		families=new int[familyCount]; agg.setFam(families);
	}

	/** Parses strict integer fields and typed count lists without transient field arrays. */
	void parse(byte[] line){
		parser.set(line);
		if(parser.terms()!=20){throw new IllegalArgumentException("Prepared bin requires exactly 20 columns");}
		agg.reset(); Arrays.fill(families, 0);
		id=parser.parseString(0); status=parser.parseString(1);
		domain=parser.parseString(2); phylum=parser.parseString(3);
		length=number(4, Long.MAX_VALUE); agg.gc=number(5, Long.MAX_VALUE);
		agg.acgt=number(6, Long.MAX_VALUE); agg.coding=number(7, Long.MAX_VALUE);
		agg.cds=(int)number(8, Integer.MAX_VALUE); agg.mapped=(int)number(9, Integer.MAX_VALUE);
		agg.glenSum=number(10, Long.MAX_VALUE); agg.glenSq=number(11, Long.MAX_VALUE);
		agg.r16=(int)number(12, Integer.MAX_VALUE); agg.r23=(int)number(13, Integer.MAX_VALUE);
		agg.r5=(int)number(14, Integer.MAX_VALUE); agg.rother=(int)number(15, Integer.MAX_VALUE);
		agg.trna=(int)number(16, Integer.MAX_VALUE);
		sparse(17, families); sparse(18, agg.anti); denseDimers(19);
		validate();
	}

	/** Checks only invariants supported by the cache/aggregation contract, not assumed mapped-count identities. */
	void validate(){
		if(id==null || id.isEmpty() || domain==null || domain.isEmpty() || phylum==null || phylum.isEmpty()){
			throw new IllegalArgumentException("Prepared bin identity and QuickClade fields must be explicit");
		}
		if(!"classified".equals(status) && !"partial".equals(status) && !"unknown".equals(status)){
			throw new IllegalArgumentException("Unknown QuickClade status for "+id);
		}
		// CacheBuilder sequence scan counts G/C as ACGT and each ACGT base as one base of length.
		if(length<=0 || agg.gc<0 || agg.gc>agg.acgt || agg.acgt>length){
			throw new IllegalArgumentException("Prepared bin violates 0<=GC<=ACGT<=positive length: "+id);
		}
		// MVM.copyAnticodonSnapshot derives unknown=trna-classified; a negative residual is invalid.
		long classified=0;
		for(int count:agg.anti){classified+=count;}
		if(classified>agg.trna){throw new IllegalArgumentException("Structural anticodons exceed total tRNA for "+id);}
		// Existing MVM.computeContext/encodeCount add one in int arithmetic. Refuse that endpoint
		// until the legacy arithmetic is changed separately; this is not a biological count bound.
		if(agg.cds==Integer.MAX_VALUE || agg.r16==Integer.MAX_VALUE || agg.r23==Integer.MAX_VALUE ||
			agg.r5==Integer.MAX_VALUE || agg.rother==Integer.MAX_VALUE || agg.trna==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Prepared CDS/RNA counts must leave room for the existing count+1 encoding: "+id);
		}
		// Coding is summed CDS length (CacheBuilder.Acc.addCds), so overlap may exceed assembly length.
		// No family-sum=mapped or mapped<=CDS assertion is implied by the serving contract.
	}

	/** Parses a nonnegative decimal integer with an explicit overflow bound. */
	private long number(int field, long max){
		parser.setBounds(field);
		return number(parser.line(), parser.a(), parser.b(), max);
	}

	/** Shared byte-range integer parser; malformed fields never wrap or silently become zero. */
	private static long number(byte[] line, int from, int to, long max){
		assert(from>=0 && to<=line.length) : "LineParser bounds and list delimiters must remain inside the source row";
		if(from>=to){throw new IllegalArgumentException("Missing prepared integer");}
		long value=0;
		for(int i=from; i<to; i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || value>(max-digit)/10){throw new IllegalArgumentException("Malformed or overflowing prepared integer");}
			value=value*10+digit;
		}
		return value;
	}

	/** Rank lists are strictly increasing and unique; a dash explicitly denotes all zero counts. */
	private void sparse(int field, int[] counts){
		parser.setBounds(field);
		final byte[] line=parser.line(); final int end=parser.b();
		int start=parser.a(), previous=-1;
		if(end-start==1 && line[start]=='-'){return;}
		if(start==end){throw new IllegalArgumentException("Use - for an empty prepared count list");}
		while(start<end){
			int colon=start; while(colon<end && line[colon]!=':' && line[colon]!=';'){colon++;}
			if(colon==end || line[colon]!=':'){throw new IllegalArgumentException("Prepared count requires rank:count");}
			int stop=colon+1; while(stop<end && line[stop]!=';'){stop++;}
			final int rank=(int)number(line, start, colon, Integer.MAX_VALUE);
			final int count=(int)number(line, colon+1, stop, Integer.MAX_VALUE);
			if(rank<=previous || rank>=counts.length){throw new IllegalArgumentException("Prepared count ranks must be unique, ordered and in range");}
			counts[rank]=count; previous=rank;
			if(stop==end){return;}
			start=stop+1;
			if(start==end){throw new IllegalArgumentException("Trailing prepared count delimiter");}
		}
	}

	/** KmerTracker uses the cache's native sixteen-dimer order; no precomputed ratios are accepted. */
	private void denseDimers(int field){
		parser.setBounds(field);
		final byte[] line=parser.line(); final int end=parser.b();
		int start=parser.a();
		for(int i=0; i<16; i++){
			int stop=start; while(stop<end && line[stop]!=','){stop++;}
			agg.dimer[i]=(int)number(line, start, stop, Integer.MAX_VALUE);
			if((i==15)!=(stop==end)){throw new IllegalArgumentException("Prepared dimer field requires exactly 16 counts");}
			start=stop+1;
		}
	}

	String id,status,domain,phylum;
	long length;
	final int[] families;
	final MagQCVectorMaker.Agg agg=new MagQCVectorMaker.Agg();
	private final LineParser1 parser=new LineParser1((byte)'\t');
}
