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
			if(Parse.parseBoolean(b)) {includeMapped=true; includeUnmapped=false;}
		}else if(a.equals("unmappedonly")){
			if(Parse.parseBoolean(b)) {includeMapped=false; includeUnmapped=true;}
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
	 * Adds contig names to the filter whitelist. Handles comma-separated lists
	 * and automatically adds common name variants (underscore/space conversion).
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
	
	public boolean passesFilter(SamLine sl){
		if(sl==null){return false;}
		return invert^matchesFilter(sl);
	}
	
	/**
	 * Internal filter logic for SAM lines. Tests all configured criteria
	 * including position, mapping quality, SAM flags, and alignment identity.
	 * @param sl SAM line to test
	 * @return true if line matches filter criteria (before inversion)
	 */
	boolean matchesFilter(SamLine sl){
		if(sl==null){return false;}
		if(!includeLengthZero && sl.length()<1){return false;}
		
		if(!sl.mapped()){return includeUnmapped;}
		else if(!includeMapped){return false;}

		if(!includeNonPrimary && !sl.primary()){return false;}
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
	
	public boolean passesFilter(VCFLine vl){
		if(vl==null){return false;}
		return invert^matchesFilter(vl);
	}
	
	/**
	 * Internal filter logic for VCF lines. Only position and contig
	 * filters are applicable to VCF data.
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
	
	public boolean passesFilter(String name){
		if(name==null){return false;}
		return invert^matchesFilter(name);
	}
	
	boolean matchesFilter(String name){
		if(name==null){return false;}
		if(contigs!=null){
			return contigs.contains(name);
		}
		return true;
	}
	
	public void clear() {
		minMapq=Integer.MIN_VALUE;
		maxMapq=Integer.MAX_VALUE;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	public int minPos=Integer.MIN_VALUE;
	public int maxPos=Integer.MAX_VALUE;
	public int minMapq=Integer.MIN_VALUE;
	public int maxMapq=Integer.MAX_VALUE;
	public float minId=Integer.MIN_VALUE;
	public float maxId=Integer.MAX_VALUE;
	public boolean includeUnmapped=true;
	public boolean includeMapped=true;
	public boolean includeSupplementary=true;
	public boolean includeQfail=false;
	public boolean includeDuplicate=true;
	public boolean includeNonPrimary=false;
	public boolean includeLengthZero=false;
	public HashSet<String> contigs=null;
	public boolean invert=false;

	/** Configures samtools filtering flags based on current filter settings.
	 * Sets ReadWrite.SAMTOOLS_IGNORE_FLAG to exclude unwanted alignment types. */
	public void setSamtoolsFilter(){
		ReadWrite.SAMTOOLS_IGNORE_FLAG=0;
		if(!includeUnmapped){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;}
		if(!includeNonPrimary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!includeSupplementary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SUPPLEMENTARY;}
		if(!includeQfail){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_QFAIL;}
		if(!includeDuplicate){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_DUPLICATE;}
	}
}