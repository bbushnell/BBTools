package var2;

import java.util.HashSet;

import fileIO.ReadWrite;
import parse.Parse;
import shared.Tools;
import stream.SamLine;

/**
 * Mutable SAM, VCF, variant and contig-name filter configuration.
 * Public passesFilter methods reject null before applying inversion. SAM records
 * use length, mapped status and selected alignment criteria; VCF/Var use region
 * and name criteria; strings use names only. These paths are not interchangeable.
 * Region bounds are zero-based and inclusive. Nonempty Var spans are converted
 * from exclusive stops; empty spans retain point-overlap behavior.
 *
 * Configure before concurrent use. Fields and the contig set are mutable and
 * unsynchronized. Identity calculation delegates to SamLine and may modify its
 * sequence while resolving an ambiguous CIGAR. External BAM filtering is a separate
 * global setting installed by setSamtoolsFilter, not by parse or passesFilter.
 * @author Brian Bushnell
 * @contributor Isla
 */
public class SamFilter{

	/**
	 * Applies one recognized option; false means this parser did not handle it.
	 * Numeric position bounds are taken literally without a one-based conversion.
	 * Identity values above one are divided by 100 once. Assertions check position
	 * ordering and identity upper bounds only, after assignment; failure can leave
	 * the new value stored. mappedonly=f and unmappedonly=f do not change settings.
	 * @param arg Original argument, unused
	 * @param a Non-null option name, already lowercased by the caller
	 * @param b Option value; contigs accepts null, other branches use their parser contracts
	 * @return Whether the option name was recognized
	 */
	public boolean parse(String arg, String a, String b){

		if(a.equals("min") || a.equals("minpos")){
			minPos=Parse.parseIntKMG(b);
			assert(minPos<=maxPos) : "minPos>maxPos";
		}else if(a.equals("max") || a.equals("maxpos")){
			maxPos=Parse.parseIntKMG(b);
			assert(minPos<=maxPos) : "minPos>maxPos";
		}else if(a.equals("minreadmapq") || a.equals("minsammapq") || a.equals("minmapq")){
			minMapq=Parse.parseIntKMG(b);
		}else if(a.equals("maxreadmapq") || a.equals("maxsammapq") || a.equals("maxmapq")){
			maxMapq=Parse.parseIntKMG(b);
		}else if(a.equals("mappedonly")){
			if(Parse.parseBoolean(b)){includeMapped=true; includeUnmapped=false;}
		}else if(a.equals("unmappedonly")){
			if(Parse.parseBoolean(b)){includeMapped=false; includeUnmapped=true;}
		}else if(a.equals("mapped")){
			includeMapped=Parse.parseBoolean(b);
		}else if(a.equals("unmapped")){
			includeUnmapped=Parse.parseBoolean(b);
		}else if(a.equals("secondary") || a.equals("nonprimary")){
			includeNonPrimary=Parse.parseBoolean(b);
		}else if(a.equals("supplementary") || a.equals("supplimentary")){
			includeSupplementary=Parse.parseBoolean(b);
		}else if(a.equals("duplicate") || a.equals("duplicates")){
			includeDuplicate=Parse.parseBoolean(b);
		}else if(a.equals("qfail") || a.equals("samqfail")){
			includeQfail=Parse.parseBoolean(b);
		}else if(a.equals("lengthzero")){
			includeLengthZero=Parse.parseBoolean(b);
		}else if(a.equals("invert")){
			invert=Parse.parseBoolean(b);
		}else if(a.equals("minid")){
			minId=Float.parseFloat(b);
			if(minId>1f){minId/=100;}
			assert(minId<=1f) : "minid should be between 0 and 1.";
		}else if(a.equals("maxid")){
			maxId=Float.parseFloat(b);
			if(maxId>1f){maxId/=100;}
			assert(maxId<=1f) : "maxid should be between 0 and 1.";
		}else if(a.equals("contigs")){
			addContig(b);
		}else{
			return false;
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Filters            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds exact, case-sensitive names and selected aliases to the shared set.
	 * Comma-separated elements are added recursively, and the original comma-bearing
	 * string is also retained. Underscores produce a space-replaced alias, which can
	 * also contribute its first whitespace-separated token. Spaces are not converted
	 * to underscores. Repeated calls accumulate; null leaves the current set unchanged.
	 * @param s Name or comma-separated names; may be null
	 */
	void addContig(String s){
		if(s==null){return;}
		if(s.indexOf(',')>=0){
			for(String s2 : s.split(",")){
				addContig(s2);
			}
		}
		if(contigs==null){contigs=new HashSet<String>();}
		contigs.add(s);
		if(s.indexOf('_')>=0){addContig(s.replace('_', ' '));}
		String[] split=s.split("\\s+");
		if(split.length>0 && !split[0].equals(s)){contigs.add(split[0]);}
	}

	/** Rejects null even under inversion; otherwise inverts the SAM match if requested.
	 * Identity filtering inherits SamLine.calcIdentity's possible sequence mutation. */
	public boolean passesFilter(SamLine sl){
		if(sl==null){return false;}
		return invert^matchesFilter(sl);
	}

	/**
	 * Tests a SAM record before inversion. QC failure and length are checked first;
	 * an unmapped record then returns includeUnmapped without checking other flags, region, mapQ,
	 * identity or contig. Mapped records use the remaining checks in source order.
	 * The region includes soft clips but excludes hard clips via SamLine.start/stop.
	 * Identity limits are skipped for null CIGAR; otherwise calcIdentity's CIGAR/MD/
	 * reference policy applies and may modify sequence or terminate on unresolved M.
	 * @param sl SAM record, possibly null
	 * @return Whether the record matches the implemented criteria, before inversion
	 */
	boolean matchesFilter(SamLine sl){
		if(sl==null){return false;}
		//VAR2-017: QC failure applies to mapped and unmapped records before inversion,
		//matching the existing SamStreamerWrapper STR-043 policy.
		if(!includeQfail && sl.discarded()){return false;}
		if(!includeLengthZero && sl.length()<1){return false;}

		if(!sl.mapped()){return includeUnmapped;}
		else if(!includeMapped){return false;}

		if(!includeNonPrimary && !sl.nonSecondary()){return false;}
		if(!includeSupplementary && sl.supplementary()){return false;}
		if(!includeDuplicate && sl.duplicate()){return false;}

		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=sl.start(true, false);
			final int stop=sl.stop(start, true, false);
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}

		if(minMapq>Integer.MIN_VALUE || maxMapq<Integer.MAX_VALUE){
			if(sl.mapq>maxMapq || sl.mapq<minMapq){return false;}
		}

		if(sl.cigar!=null && (minId>0 || maxId<1)){
			float identity=sl.calcIdentity();
			if(identity<minId || identity>maxId){return false;}
		}

		if(contigs!=null){
			String rname=sl.rnameS();
			if(rname==null){return false;}
			return contigs.contains(rname);
		}

		return true;
	}

	/** Rejects null even under inversion; otherwise applies region/name matching and inversion. */
	public boolean passesFilter(VCFLine vl){
		if(vl==null){return false;}
		return invert^matchesFilter(vl);
	}

	/**
	 * Tests only region and contig criteria, before inversion.
	 * Converts one-based VCF pos to zero-based start and uses an inclusive stop of
	 * start+max(0,reflen-1). Thus a missing/zero-length REF is treated as one position;
	 * this is a reference-span test, not an ALT/event-span interpretation.
	 * @param vl VCF record, possibly null
	 * @return Whether region and exact contig name match
	 */
	boolean matchesFilter(VCFLine vl){
		if(vl==null){return false;}

		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=vl.pos-1;
			final int stop=start+(Tools.max(0, vl.reflen-1));
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}

		if(contigs!=null){
			String rname=vl.scaf;
			if(rname==null){return false;}
			return contigs.contains(rname);
		}

		return true;
	}

	/** Rejects null even under inversion; otherwise applies the Var region/name path.
	 * map is required only when contig names are configured. */
	public boolean passesFilter(Var v, ScafMap map){
		if(v==null){return false;}
		return invert^matchesFilter(v, map);
	}

	/**
	 * Tests only region and contig criteria, before inversion.
	 * Converts a nonempty Var span's exclusive stop to its last included base.
	 * Empty spans, including insertions, are tested as points at start. Scaffold
	 * lookup is needed only for name filtering.
	 * @param v Variant, possibly null
	 * @param map Required for contig lookup when contigs is non-null
	 * @return Whether the implemented region and name checks pass
	 */
	boolean matchesFilter(Var v, ScafMap map){
		if(v==null){return false;}

		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=v.start;
			//VAR2-018: Tools.overlap includes both endpoints; retain empty spans as points.
			final int stop=(v.stop>start ? v.stop-1 : start);
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}

		if(contigs!=null){
			String rname=map.getName(v.scafnum);
			if(rname==null){return false;}
			return contigs.contains(rname);
		}

		return true;
	}

	/** Rejects null even under inversion; otherwise applies exact-name matching and inversion. */
	public boolean passesFilter(String name){
		if(name==null){return false;}
		return invert^matchesFilter(name);
	}

	/** Accepts any non-null name if contigs is null, otherwise requires exact set membership. */
	boolean matchesFilter(String name){
		if(name==null){return false;}
		if(contigs!=null){
			return contigs.contains(name);
		}
		return true;
	}

	/** Removes only mapping-quality bounds. All other settings, including invert
	 * and the contig set, remain unchanged; this does not reset the global BAM mask. */
	public void clear(){
		minMapq=Integer.MIN_VALUE;
		maxMapq=Integer.MAX_VALUE;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Inclusive zero-based lower region bound; MIN_VALUE disables the lower limit. */
	public int minPos=Integer.MIN_VALUE;
	/** Inclusive zero-based upper region bound; MAX_VALUE disables the upper limit. */
	public int maxPos=Integer.MAX_VALUE;
	/** Inclusive minimum SAM mapping quality; MIN_VALUE disables it. */
	public int minMapq=Integer.MIN_VALUE;
	/** Inclusive maximum SAM mapping quality; MAX_VALUE disables it. */
	public int maxMapq=Integer.MAX_VALUE;
	/** Minimum SAM identity fraction; filtering runs when {@code minId>0 || maxId<1}. */
	public float minId=Integer.MIN_VALUE;
	/** Maximum SAM identity fraction; initial unbounded value disables this limit. */
	public float maxId=Integer.MAX_VALUE;
	/** Unmapped-record result after QC/length checks, before inversion. */
	public boolean includeUnmapped=true;
	/** Whether mapped records may proceed to the remaining SAM checks. */
	public boolean includeMapped=true;
	/** Allow supplementary mapped records; also controls external BAM filtering. */
	public boolean includeSupplementary=true;
	/** Allow QC-failed records in Java filtering and the external BAM mask. */
	public boolean includeQfail=false;
	/** Allow duplicate mapped records; also controls external BAM filtering. */
	public boolean includeDuplicate=true;
	/** Allow secondary mapped records; supplementary status is separate. */
	public boolean includeNonPrimary=false;
	/** Skip the SamLine.length check when true. */
	public boolean includeLengthZero=false;
	/** Mutable exact-name set with configured aliases; null disables name filtering. */
	public HashSet<String> contigs=null;
	/** Invert non-null match results; null inputs always fail. */
	public boolean invert=false;

	/**
	 * Replaces ReadWrite.SAMTOOLS_IGNORE_FLAG from five include settings.
	 * The external samtools/sambamba BAM paths use this mask to discard records
	 * before Java filtering. This is global process state, not a per-filter value.
	 * Inverted filters clear the mask so input survives for Java-side inversion.
	 * Otherwise it ignores includeMapped, length, region, identity and mapQ. The
	 * secondary/supplementary/duplicate mask bits also apply to unmapped records,
	 * unlike matchesFilter's early return; the mask is not a complete Java filter.
	 * Direct SAM filtering does not consult this mask. Native BAM-reader behavior
	 * must be checked separately; selecting the native branch bypasses the external
	 * command construction in ReadWrite.getBamInputStream.
	 */
	public void setSamtoolsFilter(){
		ReadWrite.SAMTOOLS_IGNORE_FLAG=0;
		//VAR2-019: input discarded here cannot be recovered by later Java inversion.
		if(invert){return;}
		if(!includeUnmapped){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;}
		if(!includeNonPrimary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!includeSupplementary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SUPPLEMENTARY;}
		if(!includeQfail){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_QFAIL;}
		if(!includeDuplicate){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_DUPLICATE;}
	}
}
