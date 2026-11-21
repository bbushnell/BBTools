package jgi;

import java.util.Arrays;
import java.util.HashMap;

import shared.Tools;

/**
 * Mutable coverage statistics for one sequence, parsed with a shared column map.
 * Call initializeHeader before constructing rows for a table. The map is process-wide,
 * unsynchronized and replaced in place; callers must coordinate schema changes with all parsing.
 * Constructed rows retain their values independently of later header changes.
 *
 * Constructors reject nonfinite scalar input, but do not validate finite ranges or count
 * relationships. Aggregation checks count sums and weighted-result finiteness before mutation;
 * public mutation remains outside that finite-input guarantee.
 * Missing numeric columns retain zero; Covered_percent is derived rather than parsed.
 * Median_fold and Std_Dev also accept literal NA. Merged medians and standard deviations
 * are unavailable and serialize as NA; known singleton values remain numeric.
 * Unavailable values use zero numeric storage, never floating-point NaN.
 * Public length/coveredBases/underMin totals are long; clients compiled against their
 * older int field descriptors must be recompiled, even though method signatures are unchanged.
 *
 * @author Brian Bushnell, Shinobu, Sucrose
 * @date Oct 10, 2014
 */
public class CovStatsLine{

	/**
	 * Parses a tab-delimited data row using the current shared header map.
	 * String.split discards trailing empty terms; this constructor does not preserve them.
	 * @param s Data row, with an identifier in the first column
	 * @throws NumberFormatException If an accessed number is malformed, outside long range or nonfinite
	 * @throws ArrayIndexOutOfBoundsException If a configured column is absent from the row
	 * @throws NullPointerException If s is null
	 */
	public CovStatsLine(String s){
		this(s.split("\t"));
	}

	/**
	 * Parses selected columns from a caller-owned array without retaining the array.
	 * Unmapped terms are ignored; absent numeric columns leave their fields at zero.
	 * Assertions require at least eight terms, a data row rather than a header, and ID at
	 * column zero in an initialized map. They do not verify every mapped index or value.
	 * @param split Data terms under the current shared schema; accessed terms must be nonnull
	 * @throws NumberFormatException If an accessed number is malformed, outside long range or nonfinite
	 * @throws ArrayIndexOutOfBoundsException If a configured column is absent from the array
	 * @throws NullPointerException If split is null
	 */
	public CovStatsLine(String[] split){
		assert(split.length>=8) : Arrays.toString(split);
		assert(!split[0].startsWith("#")) : Arrays.toString(split);
		assert(id_FNUM==0) : "Not initialized with header: "+id_FNUM;

		if(id_FNUM>=0){id=split[id_FNUM];}
		if(length_FNUM>=0){length=Long.parseLong(split[length_FNUM]);}
		if(coveredBases_FNUM>=0){coveredBases=Long.parseLong(split[coveredBases_FNUM]);}
		if(plusReads_FNUM>=0){plusReads=Long.parseLong(split[plusReads_FNUM]);}
		if(minusReads_FNUM>=0){minusReads=Long.parseLong(split[minusReads_FNUM]);}
		if(avgFold_FNUM>=0){avgFold=parseFinite(split[avgFold_FNUM], "Avg_fold", id);}
		if(refGC_FNUM>=0){refGC=parseFinite(split[refGC_FNUM], "Ref_GC", id);}
		if(median_FNUM>=0){
			final String token=split[median_FNUM];
			medianAvailable=!"NA".equals(token);
			if(medianAvailable){median=Long.parseLong(token);}
		}
		if(underMin_FNUM>=0){underMin=Long.parseLong(split[underMin_FNUM]);}
		if(readGC_FNUM>=0){readGC=parseFinite(split[readGC_FNUM], "Read_GC", id);}
		if(stdDev_FNUM>=0){
			final String token=split[stdDev_FNUM];
			stdDevAvailable=!"NA".equals(token);
			if(stdDevAvailable){stdDev=parseFinite(token, "Std_Dev", id);}
		}
	}

	/**
	 * Retains lexical parse failures and every finite value, including signed zero and sentinels.
	 * @param token Original scalar term
	 * @param column Coverage column name for a nonfinite-value diagnostic
	 * @param id Record identifier for a nonfinite-value diagnostic
	 * @return Parsed finite value without range normalization
	 * @throws NumberFormatException If parsing fails or the result is NaN or infinity
	 * @throws NullPointerException If token is null
	 */
	private static double parseFinite(final String token, final String column, final String id){
		final double value=Double.parseDouble(token);
		if(Double.isNaN(value) || Double.isInfinite(value)){
			throw new NumberFormatException("Nonfinite coverage value for record '"+id+"', column '"+column+"': "+token);
		}
		return value;
	}

	/** Returns 100*coveredBases/max(1,length), without clamping the result to a percentage range. */
	public final double coveredPercent(){
		return (100.0*coveredBases)/Tools.max(1, length);
	}

	/** Returns plusReads+minusReads, throwing ArithmeticException if the total cannot fit in long. */
	public final long reads(){return Math.addExact(plusReads, minusReads);}

	/**
	 * Merges the other row's current values into this receiver, retaining the receiver's ID.
	 * Fold coverage/reference GC use length weights; read GC uses read-count weights, not
	 * observed base counts. Weight denominators are clamped to at least one, not the inputs.
	 * Checks base totals, both rows' read sums, strand totals, combined reads and weighted
	 * result finiteness before mutating any receiver state. Finite operation order is retained.
	 *
	 * Per-row summaries cannot reconstruct a pooled median or standard deviation:
	 * successful merges mark both unavailable and clear their numeric storage to zero,
	 * including zero-length merges. This method does not generally validate
	 * public mutations. Finite ranges and count relationships are not normalized.
	 * Neither row may be concurrently mutated during this operation.
	 * @param csl Row to merge; its fields are read directly, without a defensive copy
	 * @throws ArithmeticException If a count sum overflows long or a weighted result is nonfinite; receiver remains unchanged
	 * @throws NullPointerException If csl is null
	 */
	public void add(CovStatsLine csl){
		final long combinedLength=Math.addExact(length, csl.length);
		final long combinedCovered=Math.addExact(coveredBases, csl.coveredBases);
		final long combinedUnder=Math.addExact(underMin, csl.underMin);
		final long thisReads=reads(), otherReads=csl.reads();
		final long combinedPlus=Math.addExact(plusReads, csl.plusReads);
		final long combinedMinus=Math.addExact(minusReads, csl.minusReads);
		final long combinedReads=Math.addExact(thisReads, otherReads);
		final double invlen2=1d/Tools.max(1L, combinedLength);
		final double combinedFold=((avgFold*length)+(csl.avgFold*csl.length))*invlen2;
		final double combinedRefGC=((refGC*length)+(csl.refGC*csl.length))*invlen2;
		final double combinedReadGC=((readGC*thisReads)+(csl.readGC*otherReads))*1.0/(Tools.max(1, combinedReads));
		checkFiniteAggregate(combinedFold, "Avg_fold", csl);
		checkFiniteAggregate(combinedRefGC, "Ref_GC", csl);
		checkFiniteAggregate(combinedReadGC, "Read_GC", csl);

		avgFold=combinedFold;
		refGC=combinedRefGC;
		readGC=combinedReadGC;
		length=combinedLength;
		coveredBases=combinedCovered;
		plusReads=combinedPlus;
		minusReads=combinedMinus;
		median=0;
		medianAvailable=false;
		stdDev=0;
		stdDevAvailable=false;
		underMin=combinedUnder;
	}

	/** Rejects a nonfinite candidate before add assigns fields, naming both records and the weighted column. */
	private void checkFiniteAggregate(final double value, final String column, final CovStatsLine other){
		if(Double.isNaN(value) || Double.isInfinite(value)){
			throw new ArithmeticException("Nonfinite aggregated coverage value for records '"+id+"' and '"+other.id
					+"', column '"+column+"': "+value);
		}
	}

	/** Returns the known median or literal NA when parsing/aggregation marked it unavailable. */
	public final String medianText(){return medianAvailable ? Long.toString(median) : "NA";}

	/** Returns the known deviation to four decimal places, or literal NA when unavailable. */
	public final String stdDevText(){return stdDevAvailable ? Tools.format("%.4f", stdDev) : "NA";}

	/**
	 * Serializes twelve fields in a fixed order, with four decimal places for finite doubles.
	 * Order: ID, Avg_fold, Length, Ref_GC, derived Covered_percent, Covered_bases, Plus_reads,
	 * Minus_reads, Median_fold, Under_min, Read_GC, Std_Dev. The active input schema and its
	 * original Under threshold/window label do not determine this output layout.
	 * Median_fold and Std_Dev use literal NA when unavailable, retaining their column positions.
	 * @return Tab-delimited row without a header or trailing newline
	 */
	@Override
	public String toString(){
		return Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%s\t%d\t%.4f\t%s", id, avgFold, length,
				refGC, coveredPercent(), coveredBases, plusReads, minusReads, medianText(), underMin, readGC, stdDevText());
	}

	/**
	 * Replaces the process-wide parsing indices for subsequent constructors.
	 * Removes leading # characters, splits on tabs, and lowercases names using the default
	 * locale without trimming. Any Under_* name maps to Under_min; duplicate normalized
	 * names use their last position and unknown columns are ignored. Missing names map to -1.
	 * An assertion requires ID at column zero, but runs after replacing the shared indices;
	 * failure does not restore the previous schema. No schema is retained in constructed rows.
	 * @param header Column labels for the next serially parsed table
	 * @throws NullPointerException If header is null
	 */
	public static void initializeHeader(String header){
		while(header.startsWith("#")){header=header.substring(1);}
		String[] split=header.split("\t");
		HashMap<String, Integer> map=new HashMap<String, Integer>(23);
		for(int i=0; i<split.length; i++){
			String s=split[i].toLowerCase();
			if(s.startsWith("under_")){s="under_min";}
			map.put(s, i);
		}
		id_FNUM=map.containsKey("id") ? map.get("id") : -1;
		avgFold_FNUM=map.containsKey("avg_fold") ? map.get("avg_fold") : -1;
		length_FNUM=map.containsKey("length") ? map.get("length") : -1;
		refGC_FNUM=map.containsKey("ref_gc") ? map.get("ref_gc") : -1;
		coveredBases_FNUM=map.containsKey("covered_bases") ? map.get("covered_bases") : -1;
		plusReads_FNUM=map.containsKey("plus_reads") ? map.get("plus_reads") : -1;
		minusReads_FNUM=map.containsKey("minus_reads") ? map.get("minus_reads") : -1;
		median_FNUM=map.containsKey("median_fold") ? map.get("median_fold") : -1;
		underMin_FNUM=map.containsKey("under_min") ? map.get("under_min") : -1;
		readGC_FNUM=map.containsKey("read_gc") ? map.get("read_gc") : -1;
		stdDev_FNUM=map.containsKey("std_dev") ? map.get("std_dev") : -1;

		assert(id_FNUM==0) : "Bad header: "+id_FNUM+"\n"+header+"\n"+map;
	}

	/** Sequence identifier copied from the mapped ID term. */
	public String id;
	/** Reported sequence length in bases; count ranges/relationships are not validated. */
	public long length;
	/** Reported covered-base count, used to derive coveredPercent. */
	public long coveredBases;
	/** Reported plus-strand read count. */
	public long plusReads;
	/** Reported minus-strand read count. */
	public long minusReads;
	/** Reported mean fold coverage; construction requires a finite value. */
	public double avgFold;
	/** Reported reference GC fraction; construction requires a finite value, without a range check. */
	public double refGC;
	/** Reported median fold coverage, meaningful only when medianAvailable is true; otherwise zero. */
	public long median;
	/** Whether median is available; absent columns retain the legacy known-zero default. */
	public boolean medianAvailable=true;
	/** Count from an Under_* column; its threshold/window metadata is not retained here. */
	public long underMin;
	/** Reported read GC fraction; construction requires a finite value, without a range check. */
	public double readGC;
	/** Reported coverage standard deviation, meaningful only when stdDevAvailable is true; otherwise zero. */
	public double stdDev;
	/** Whether standard deviation is available; absent columns retain the legacy known-zero default. */
	public boolean stdDevAvailable=true;

	/** Shared column indices, initially -1; initializeHeader replaces them without synchronization. */
	private static int id_FNUM=-1;
	private static int length_FNUM=-1;
	private static int coveredBases_FNUM=-1;
	private static int plusReads_FNUM=-1;
	private static int minusReads_FNUM=-1;
	private static int avgFold_FNUM=-1;
	private static int refGC_FNUM=-1;
	private static int median_FNUM=-1;
	private static int underMin_FNUM=-1;
	private static int readGC_FNUM=-1;
	private static int stdDev_FNUM=-1;
}
