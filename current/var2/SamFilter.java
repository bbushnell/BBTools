package var2;

import java.util.HashSet;

import fileIO.ReadWrite;
import shared.Parse;
import shared.Tools;
import stream.SamLine;

/**
 * Filters SAM/BAM alignments, VCF variants, and other genomic data based on 
 * configurable criteria including mapping quality, position ranges, alignment 
 * identity, and SAM flags. Supports both inclusive and exclusive filtering modes.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class SamFilter {
	
	/**
	 * Parses command-line arguments to configure filter parameters.
	 * 
	 * @param arg Original argument string
	 * @param a Argument key (lowercase)
	 * @param b Argument value
	 * @return true if argument was recognized and parsed
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
		}else if(a.equals("mapped")){
			includeMapped=Parse.parseBoolean(b);
		}else if(a.equals("unmapped")){
			includeUnmapped=Parse.parseBoolean(b);
		}else if(a.equals("secondary") || a.equals("nonprimary")){
			includeNonPrimary=Parse.parseBoolean(b);
		}else if(a.equals("supplimentary")){
			includeSupplimentary=Parse.parseBoolean(b);
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
	 * Adds contig names to the filter whitelist. Handles comma-separated lists
	 * and automatically adds common name variants (underscore/space conversion).
	 * 
	 * @param s Contig name or comma-separated list of names
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
	
	/**
	 * Tests whether a SAM line passes all filter criteria.
	 * 
	 * @param sl SAM line to test
	 * @return true if line passes filters (accounting for invert flag)
	 */
	public boolean passesFilter(SamLine sl){
		if(sl==null){return false;}
		return invert^matchesFilter(sl);
	}
	
	/**
	 * Internal filter logic for SAM lines. Tests all configured criteria
	 * including position, mapping quality, SAM flags, and alignment identity.
	 * 
	 * @param sl SAM line to test
	 * @return true if line matches filter criteria (before inversion)
	 */
	boolean matchesFilter(SamLine sl){
		if(sl==null){return false;}
		if(!includeLengthZero && sl.length()<1){return false;}
		
		if(!sl.mapped()){return includeUnmapped;}
		else if(!includeMapped){return false;}

		if(!includeNonPrimary && !sl.primary()){return false;}
		if(!includeSupplimentary && sl.supplementary()){return false;}
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
	
	/**
	 * Tests whether a VCF line passes position and contig filters.
	 * 
	 * @param vl VCF line to test
	 * @return true if line passes filters (accounting for invert flag)
	 */
	public boolean passesFilter(VCFLine vl){
		if(vl==null){return false;}
		return invert^matchesFilter(vl);
	}
	
	/**
	 * Internal filter logic for VCF lines. Only position and contig
	 * filters are applicable to VCF data.
	 * 
	 * @param vl VCF line to test
	 * @return true if line matches filter criteria (before inversion)
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
	
	/**
	 * Tests whether a Var object passes position and contig filters.
	 * 
	 * @param v Var object to test
	 * @param map ScafMap for contig name resolution
	 * @return true if variant passes filters (accounting for invert flag)
	 */
	public boolean passesFilter(Var v, ScafMap map){
		if(v==null){return false;}
		return invert^matchesFilter(v, map);
	}
	
	/**
	 * Internal filter logic for Var objects. Only position and contig
	 * filters are applicable to variant data.
	 * 
	 * @param v Var object to test
	 * @param map ScafMap for contig name resolution
	 * @return true if variant matches filter criteria (before inversion)
	 */
	boolean matchesFilter(Var v, ScafMap map){
		if(v==null){return false;}
		
		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=v.start;
			final int stop=v.stop;
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}
		
		if(contigs!=null){
			String rname=map.getName(v.scafnum);
			if(rname==null){return false;}
			return contigs.contains(rname);
		}
		
		return true;
	}
	
	/**
	 * Tests whether a contig name passes the contig filter.
	 * 
	 * @param name Contig name to test
	 * @return true if name passes filter (accounting for invert flag)
	 */
	public boolean passesFilter(String name){
		if(name==null){return false;}
		return invert^matchesFilter(name);
	}
	
	/**
	 * Internal filter logic for contig names.
	 * 
	 * @param name Contig name to test
	 * @return true if name matches filter criteria (before inversion)
	 */
	boolean matchesFilter(String name){
		if(name==null){return false;}
		if(contigs!=null){
			return contigs.contains(name);
		}
		return true;
	}
	
	/**
	 * Resets mapping quality filters to default (no filtering).
	 */
	public void clear() {
		minMapq=Integer.MIN_VALUE;
		maxMapq=Integer.MAX_VALUE;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum position for coordinate filtering */
	public int minPos=Integer.MIN_VALUE;
	/** Maximum position for coordinate filtering */
	public int maxPos=Integer.MAX_VALUE;
	/** Minimum mapping quality threshold */
	public int minMapq=Integer.MIN_VALUE;
	/** Maximum mapping quality threshold */
	public int maxMapq=Integer.MAX_VALUE;
	/** Minimum alignment identity (0.0-1.0) */
	public float minId=Integer.MIN_VALUE;
	/** Maximum alignment identity (0.0-1.0) */
	public float maxId=Integer.MAX_VALUE;
	/** Whether to include unmapped reads */
	public boolean includeUnmapped=true;
	/** Whether to include mapped reads */
	public boolean includeMapped=true;
	/** Whether to include supplementary alignments */
	public boolean includeSupplimentary=true;
	/** Whether to include reads that failed quality checks */
	public boolean includeQfail=false;
	/** Whether to include duplicate reads */
	public boolean includeDuplicate=true;
	/** Whether to include non-primary (secondary) alignments */
	public boolean includeNonPrimary=false;
	/** Whether to include zero-length alignments */
	public boolean includeLengthZero=false;
	/** Set of allowed contig names (null means all allowed) */
	public HashSet<String> contigs=null;
	/** Whether to invert all filter results */
	public boolean invert=false;

	/**
	 * Configures samtools filtering flags based on current filter settings.
	 * Sets ReadWrite.SAMTOOLS_IGNORE_FLAG to exclude unwanted alignment types.
	 */
	public void setSamtoolsFilter(){
		ReadWrite.SAMTOOLS_IGNORE_FLAG=0;
		if(!includeUnmapped){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;}
		if(!includeNonPrimary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!includeSupplimentary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SUPPLIMENTARY;}
		if(!includeQfail){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_QFAIL;}
		if(!includeDuplicate){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_DUPLICATE;}
	}
}