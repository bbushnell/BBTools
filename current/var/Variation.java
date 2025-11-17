package var;
import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import dna.Data;
import dna.Gene;
import dna.GeneSet;
import driver.Search;
import shared.Shared;
import shared.Tools;
import structures.Range;



/**
 * Represents a genomic variation such as SNPs, insertions, deletions, or complex variations.
 * Provides coordinate-based comparisons, overlap detection, and gene interaction analysis.
 * Supports various variation types from reference matches to complex substitutions.
 * @author Brian Bushnell
 */
public class Variation implements Comparable<Variation>, Serializable, Cloneable {
	
//	>locus  ploidy  haplotype       chromosome      begin   end     varType reference       alleleSeq       totalScore      hapLink xRef

	/**
	 * 
	 */
	private static final long serialVersionUID = -3847258470952802740l;

	/**
	 * Creates a variation from a parsed VarLine object.
	 * Validates consistency between begin/end locations and variation type.
	 * @param line The VarLine containing variation data to convert
	 */
	public Variation(VarLine line){
//		this(line.chromosome, line.beginLoc, line.endLoc, line.xRef, line.varType, line.ref, line.call);
		this(line.chromosome, line.beginLoc, line.endLoc, line.varType, line.ref, line.call);
		
		assert(!((varType==INS || varType==DELINS || varType==SNP) && call==null)) : "\n"+line+"\n"+this+
		"\nline.ref="+line.ref+"\nline.call="+line.call+"\nref="+ref+"\ncall="+call;
		
		assert(beginLoc<=endLoc) : line.toString();

		assert(this.equals(line)) : "\n\n"+this+"\n!=\n"+line;
		assert(line.equals(this)) : "\n\n"+this+"\n!=\n"+line;
		

//		if(xRef==11429487){
//			System.out.println("\n"+this.toString());
//		}
	}
	
//	public Variation(GeneVarLine line){
////		this(line.chromosome, line.beginLoc, line.endLoc, line.xRef, line.varType, line.ref, line.call);
//		this(line.chromosome, line.beginLoc, line.endLoc, line.xRef, line.xRefArray, line.varType, line.ref, line.call);
//
//		assert(beginLoc<=endLoc) : line.toString();
//
//		assert(this.equals(line)) : "\n\n"+this+"\n!=\n"+line.toSuperString()+"\n\n"+line;
//		assert(line.equals(this)) : "\n\n"+this+"\n!=\n"+line.toSuperString()+"\n\n"+line;
//
//	}
	
	/**
	 * Copy constructor that creates a new variation from an existing one.
	 * Performs deep validation to ensure consistency of copied data.
	 * @param line The variation to copy
	 */
	public Variation(Variation line){
		this(line.chromosome, line.beginLoc, line.endLoc, line.varType, line.ref, line.call);
		
		assert(beginLoc<=endLoc) : line.toString();

		assert(this.equals(line)) : "\n\n"+this+"\n!=\n"+line.toSuperString()+"\n\n"+line;
		assert(line.equals(this)) : "\n\n"+this+"\n!=\n"+line.toSuperString()+"\n\n"+line;
		
	}
	
	/**
	 * Creates a variation with specified genomic coordinates and sequence details.
	 *
	 * @param chr Chromosome number (0-based indexing)
	 * @param bLoc Begin location on chromosome (0-based)
	 * @param eLoc End location on chromosome (0-based, inclusive)
	 * @param vType Variation type constant (REF, SNP, INS, DEL, etc.)
	 * @param rf Reference sequence at this location
	 * @param ca Called/alternate sequence at this location
	 */
	public Variation(int chr, int bLoc, int eLoc, byte vType, String rf, String ca){
		chromosome=chr;
		beginLoc=bLoc;
		endLoc=eLoc;
		varType=vType;
		
		setDetails(vType, rf, ca);
		
		assert(beginLoc<=endLoc) : toString();
		
	}
	
	/** Default constructor for creating empty variation objects */
	public Variation(){}
	
	
	@Override
	public Variation clone(){
		Variation v=null;
		try {
			v=(Variation) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return v;
	}
	
	
	/**
	 * Converts an array of VarLine objects to a set of unique Variation objects.
	 * Optionally filters out reference-equal variations to focus on actual changes.
	 *
	 * @param array Array of VarLine objects to convert
	 * @param retainEqual Whether to include variations that match reference
	 * @return Set containing unique variations from the input array
	 */
	public static final HashSet<Variation> toVariations(VarLine[] array, boolean retainEqual){
		HashSet<Variation> set=new HashSet<Variation>(array.length);
		for(VarLine line : array){
			Variation var=new Variation(line);
			if(retainEqual || var.varType!=Variation.REF){
				if(!set.contains(var)){
					set.add(var);
				}
			}
		}
		return set;
	}
	
	/**
	 * Converts multiple VarLine arrays into a single sorted array of unique variations.
	 * Combines all input arrays before deduplication and sorting.
	 *
	 * @param array Two-dimensional array of VarLine objects
	 * @param retainEqual Whether to include reference-equal variations
	 * @return Sorted array of unique variations from all input arrays
	 */
	public static final Variation[] toVariationArray(VarLine[][] array, boolean retainEqual){
		HashSet<Variation> set=toVariations(array[0], retainEqual);
		for(int i=1; i<array.length; i++){
			set.addAll(toVariations(array[i], retainEqual));
		}
		Variation[] vars=set.toArray(new Variation[set.size()]);
		Arrays.sort(vars);
		return vars;
	}
	
	/**
	 * Converts a single VarLine array into a sorted array of unique variations.
	 * @param array Array of VarLine objects to convert
	 * @param retainEqual Whether to include reference-equal variations
	 * @return Sorted array of unique variations from the input
	 */
	public static final Variation[] toVariationArray(VarLine[] array, boolean retainEqual){
		HashSet<Variation> set=toVariations(array, retainEqual);
		Variation[] vars=set.toArray(new Variation[set.size()]);
		Arrays.sort(vars);
		return vars;
	}
	
	/**
	 * Converts a Set of comparable objects to a sorted array of the specified type.
	 * Uses reflection to create appropriately typed arrays.
	 *
	 * @param c The class type for the returned array
	 * @param set Set of objects to convert to array
	 * @return Sorted array containing all elements from the set
	 */
	@SuppressWarnings("unchecked")
	public static final <X extends Comparable<? super X>> X[] toArray(Class<X> c, Set<X> set){
		
		set.getClass().getTypeParameters();
		X[] array=(X[])Array.newInstance(c,set.size());
		
		array=set.toArray(array);
		int i=0;
		for(X x : set){
			array[i]=x;
			i++;
		}
		Arrays.sort(array);
		return array;
	}
	
	
	/**
	 * Filters variations to include only those intersecting coding regions or nearby areas.
	 * Uses gene annotation data to determine which variations affect coding sequences.
	 *
	 * @param variances Array of variations to filter
	 * @param chrom Chromosome number to analyze
	 * @param nearby If true, include variations near coding regions; if false, only overlapping
	 * @return Array containing only coding-relevant variations
	 */
	public static VarLine[] filterCodingVariances(VarLine[] variances, int chrom, boolean nearby){
		Range[] ranges=(nearby ? Data.geneNearbyRangeMatrix(chrom) : Data.geneCodeAndExonRangeMatrix(chrom));
		
		ArrayList<VarLine> list=new ArrayList<VarLine>(8+variances.length/8);
		
		for(VarLine var : variances){

			if(var.varType!=VarLine.REF && var.varType!=VarLine.NOREF){
				int loc=var.beginLoc;
				int rnum=Search.findPointBinary(loc, ranges);
				
				if(ranges[rnum].intersects(var.beginLoc, var.endLoc)){
					list.add(var);
				}
				
				for(int i=rnum; i<ranges.length; i++){
					Range r=ranges[i];
					if(r.a>var.endLoc){break;} //Out of range

					if(r.intersects(var.beginLoc, var.endLoc)){
						list.add(var);
						break;
					}
				}
			}
		}
		
		return list.toArray(new VarLine[list.size()]);
	}
	
	


	/**
	 * Generates an array of non-overlapping Ranges, sorted by index, ascending.
	 * To each is attached a list of all overlapping Variations from the input array.
	 * @param va
	 * @return The array of ranges
	 */
	public static Range[] makeVarRanges(Variation[] va){
		//			System.out.println("va.length="+va.length);
		
		if(va==null || va.length==0){
			return new Range[0];
		}
		
		ArrayList<Range> ra=new ArrayList<Range>(va.length);
		for(Variation v : va){
			Range r=new Range(v.beginLoc, v.endLoc);
			r.obj1=new ArrayList<Variation>();
			((ArrayList<Variation>)r.obj1).add(v);
			ra.add(r);
		}
		Shared.sort(ra);
		ArrayList<Range> ra2=new ArrayList<Range>(va.length);
		Range current=null;
		//			System.out.println("ra.size="+ra.size());
		for(Range r : ra){
			//				System.out.println("\ncurrent="+current+", r="+r);
			if(current==null){current=r;}
			else if(current.intersects(r)){
				//					System.out.println("merged");
				Range temp=current.merge(r);
				temp.obj1=current.obj1;
				((ArrayList<Variation>)temp.obj1).addAll((ArrayList<Variation>)r.obj1);
				current=temp;
			}else{
				//					System.out.println("added");
				ra2.add(current);
				current=r;
			}
			//				System.out.println("current="+current+", r="+r);
		}
		//			System.out.println("\ncurrent="+current);
		//			System.out.println("ra2.size="+ra2.size());
		assert(current!=null); //Note: this could be null if input was empty, I guess...
		assert(ra2.size()==0 || ra2.get(ra2.size()-1)!=current);
		ra2.add(current);
		return ra2.toArray(new Range[ra2.size()]);
	}
	
	/**
	 * Extracts rsID number from a reference identifier string.
	 * Handles various dbSNP reference formats including prefixes.
	 * @param s Reference string containing rsID information
	 * @return Numeric rsID, or -1 if not found or invalid
	 */
	public static final int toRsid(String s){return xRefToId(s);}
	/**
	 * Converts external reference strings to numeric identifiers.
	 * Parses strings containing colons, prefixes, and other formatting to extract IDs.
	 * @param s External reference string to parse
	 * @return Numeric identifier, or -1 if parsing fails
	 */
	public static final int xRefToId(String s){
//		System.out.println(s);
		if(s==null || s.length()==0){return -1;}
//		assert(s.startsWith("dbsnp:rs")) : s;
		
		if(s.contains(":")){
			s=s.substring(s.indexOf(':')+1);
		}
		
		int i=0, max=s.length();
//		System.err.println(s);
		while(i<max && !Tools.isDigit(s.charAt(i))){i++;}
		if(i>=max){assert(s.equals(".")) : s; return -1;}
		s=s.substring(i);
		
		return Integer.parseInt(s);
	}
	
	/**
	 * Extracts multiple rsID numbers from a comma or semicolon separated string.
	 * @param s String containing multiple reference identifiers
	 * @return Array of numeric rsIDs, or null if no valid IDs found
	 */
	public static final int[] toRsidArray(String s){return xRefToIdArray(s);}
	/**
	 * Converts multiple external reference strings to an array of numeric identifiers.
	 * Splits on commas and semicolons to handle multi-ID strings.
	 * @param s String containing multiple external references
	 * @return Array of numeric identifiers, or null if no valid IDs found
	 */
	public static final int[] xRefToIdArray(String s){
		if(s==null || s.length()<1){return null;}
		String[] array=s.split("[,;]");
		int[] r=new int[array.length];
		for(int i=0; i<array.length; i++){
			r[i]=xRefToId(array[i]);
			if(r[i]==-1){
				if(r.length==1){return null;}
				
				//This can be safely disabled.  But it is best to fix this case by making the array smaller.
				assert(false) : "Not a real rsID: "+s;
			}
		}
		return r;
	}
	
	/**
	 * Determines if this variation exactly matches another variation.
	 * Compares chromosome, coordinates, type, and sequence content for equality.
	 * @param line Variation to compare against
	 * @return true if variations are identical, false otherwise
	 */
	public boolean matches(Variation line){
		if(line==null || chromosome!=line.chromosome || beginLoc!=line.beginLoc || endLoc!=line.endLoc || varType!=line.varType){
			return false;
		}
		return matches(line.varType, line.ref, line.call);
	}
	
//	public boolean matchesLoose(VarLine line){
//		if(line==null || chromosome!=line.chromosome || !intersects(line)){
//			return false;
//		}
//		if(isEqual() && line.isEqual()){return true;}
//		if(varType!=line.varType){return false;}
//		return matches(line.varType, line.ref, line.call);
//	}
	
	/** Overlap and don't contradict each other */
	public boolean matchesLoose(VarLine line){
		if(line==null || chromosome!=line.chromosome || !intersects(line)){
			return false;
		}
		
		if(isTrueVariation()){
			if(varType!=line.varType){return false;}
			return matches(line.varType, line.ref, line.call);
		}else if(isRef()){
			return line.isRef();
		}else{
			assert(isUnsureVariation()) : this;
			return line.isUnsureVariation();
		}
	}
	
	/**
	 * Internal method to check sequence compatibility for specific variation types.
	 * Handles type-specific matching logic for different variation classes.
	 *
	 * @param type Variation type to check
	 * @param ref2 Reference sequence to compare
	 * @param call2 Called sequence to compare
	 * @return true if sequences are compatible for the given type
	 */
	private boolean matches(int type, String ref2, String call2){
		if(type==REF || type==REFPOINT || type==DEL || type==NOCALL){
			return true;
		}
		return call.equals(call2);
		
	}
	
	/**
	 * Internal method to set reference and call sequences based on variation type.
	 * Applies type-specific rules for which sequences should be stored.
	 * Automatically interns strings to reduce memory usage.
	 *
	 * @param vt Variation type constant
	 * @param rf Reference sequence string
	 * @param ca Called/alternate sequence string
	 */
	private void setDetails(byte vt, String rf, String ca){
		
		ref=null;
		call=null;
		
		switch(vt){
			
			case REF: {
			}break;
			case SNP: {
				ref=rf; call=ca;
			}break;
			case INS: {
				call=ca;
			}break;
			case DEL: {
				ref=rf;
			}break;
			case DELINS: {
				ref=rf; call=ca;
			}break;
			case REFCON: {
				ref=rf; call=ca;
			}break;
			case REFINCON: {
				ref=rf; call=ca;
			}break;
			case NOCALL: {
				 //I can't remember if nocalls need N or null calls
//				ref=rf;
//				call=ca;
//				assert(ref!=null && call!=null && ref.length()==call.length()) : ref+", "+call;
			}break;
			case NOREF: {
				 //I can't remember if norefs need N or null refs
//				ref=rf;
				call=ca;
//				assert(ref!=null && call!=null && ref.length()==call.length()) : ref+", "+call;
			}break;
			case PAR: {
				ref=rf; call=ca;
			}break;
			case NULL: {
				ref=rf; call=ca;
			}break;
			case REFPOINT: {
				ref=call=null;
			}break;

			default: {assert(false);}
		}
		intern();
	}
	
	/**
	 * Returns a formatted string representing the genomic location.
	 * Uses 0-based coordinates by default.
	 * @return Location string in format "(start)" or "(start - end)"
	 */
	public String locationString(){
		return locationString(0);
	}
	
	/**
	 * Returns a formatted string representing the genomic location with specified base offset.
	 * @param base Base offset (0 for 0-based, 1 for 1-based coordinates)
	 * @return Location string adjusted by the base offset
	 */
	public String locationString(int base){
		assert(base==0 || base==1);
		
		if(beginLoc==endLoc){
			return "("+(beginLoc+base)+")";
		}
		return "("+(beginLoc+base)+" - "+(endLoc+base)+")";
		
//		if(beginLoc==endLoc){
//			return (beginLoc+base)+"";
//		}
//		return (beginLoc+base)+"-"+(endLoc+base);
	}
	
	/** Returns the same output as toString() for compatibility */
	public String toSuperString(){return toString();}
	
	@Override
	public String toString(){
		return toString(0);
	}
	
	/**
	 * Returns a formatted string representation with specified coordinate base.
	 * @param base Base offset for coordinates (0 for 0-based, 1 for 1-based)
	 * @return Tab-separated string with adjusted coordinates
	 */
	public String toString(int base){
		StringBuilder sb=new StringBuilder();
		
		sb.append("chr"+Gene.chromCodes[chromosome]);
		while(sb.length()<5){sb.append(' ');}
		sb.append('\t');
		sb.append(locationString(base)+"\t");
		
		sb.append(varTypeMap[varType]);
		
		sb.append("\t"+(ref==null ? "" : ref));
		sb.append("\t"+(call==null ? "" : call));
		
		sb.append('\t');
		
		return sb.toString();
	}
	
	/**
	 * Returns a string formatted for source file output compatibility.
	 * Uses specific coordinate adjustments for different variation types.
	 * @return Tab-separated string in source file format
	 */
	public String toSourceString(){
		StringBuilder sb=new StringBuilder(64);
		
		sb.append("chr"+Gene.chromCodes[chromosome]+"\t");
		sb.append(beginLoc+"\t");
		
		if(varType==INS){
			sb.append(beginLoc+"\t");
		}else{
			sb.append((endLoc+1)+"\t");
		}
		
		sb.append(varTypeMap[varType]+"\t");
		sb.append((ref==null ? "" : ref)+"\t");
		sb.append((call==null ? "" : call)+"\t");
		
		return sb.toString();
	}
	
	/** Returns the standard header line for variation file output.
	 * @return Tab-separated header string */
	public static String header(){
		
		return "chrom\tstart\tstop\ttype\tref\tcall\trsID";
	}
	
	/** Returns a compact string representation containing only essential information.
	 * @return Shortened tab-separated string with location, type, and sequences */
	public String toShortString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append(locationString()+"\t");
		
		sb.append(varTypeMap[varType]);

		if(ref!=null){sb.append("\t"+ref);}
		if(call!=null){sb.append("\t"+call);}
		
		return sb.toString();
	}
	
	
	/**
	 * Finds the index of a string within an array using exact matching.
	 * @param a String to search for
	 * @param array Array to search within
	 * @return Index of the string, or -1 if not found
	 */
	public static final int find(String a, String[] array){
		for(int i=0; i<array.length; i++){
			if(a.equals(array[i])){return i;}
		}
		assert(false) : "Can't find "+a+" in "+Arrays.toString(array);
		return -1;
	}
	
	/**
	 * Returns the length of the reference sequence affected by this variation.
	 * Handles type-specific length calculations for different variation classes.
	 * @return Length of reference sequence in bases
	 */
	public final int lengthRef(){
		switch(varType){
			case SNP: {
				assert(endLoc-beginLoc+1==1) : "\n"+endLoc+"-"+beginLoc+"+1 = "+(endLoc-beginLoc+1)+" != "+1+"\n"+this.toString()+"\n";
				assert(call!=null && call.length()==1) : "\ncall= '"+call+"'\n"+this.toString();
				return 1;
			}
			case INS: {
				return 0;
			}
			case REFPOINT: {
				return 0;
			}
//			case NOREF: {
//				throw new RuntimeException();
//			}
//			case NULL: {
//				throw new RuntimeException();
//			}
			default: {
				break;
			}
		}
		return endLoc-beginLoc+1;
	}

	/** Returns the maximum of reference length and variation length */
	public final int lengthMax(){return max(lengthRef(), lengthVar());}
	/** Returns the minimum of reference length and variation length */
	public final int lengthMin(){return min(lengthRef(), lengthVar());}
	/** Returns the difference between variation length and reference length */
	public final int lengthDif(){return isNR_or_NC() ? 0 : lengthVar()-lengthRef();}
	
	/**
	 * Returns the length of the variant sequence for this variation.
	 * Handles type-specific calculations including insertions and complex variations.
	 * @return Length of variant sequence in bases
	 */
	public final int lengthVar(){
		switch(varType){
			
			case REF: {
				return endLoc-beginLoc+1;
			}
			case SNP: {
				assert(endLoc-beginLoc+1==1);
				assert(call!=null && call.length()==1);
				return 1;
			}
			case INS: {
				assert(call!=null);
				return call.length();
			}
			case REFPOINT: {
				return 0;
			}
			case DEL: {
				assert(call==null);
				return 0;
			}
			case DELINS: {
				assert(call!=null);
				return call.length();
			}
			case REFCON: {
				return endLoc-beginLoc+1;
			}
			case REFINCON: {
				assert(false) : "Warning - Length cannot be known for certain.";
				return endLoc-beginLoc+1;
			}
			case NOCALL: {
				assert(false) : "Warning - Length cannot be known for certain.";
				return endLoc-beginLoc+1;
			}
			case NOREF: {
				assert(false) : "Warning - Length cannot be known for certain.";
				return endLoc-beginLoc+1;
			}
			case PAR: {
				assert(call!=null);
				return call.length();
			}
			case NULL: {
				assert(false);
				throw new RuntimeException();
			}

			default: {throw new RuntimeException();}
		}
	}
	
	
//	//TODO Note that this may be wrong for e.g. insertions, deletions, and if/when changed to half-open numbering.
//	public final int length(){
//		if(varType==INS){return 0;}
//		return endLoc-beginLoc+1;
//	}
//	public final int length2(){
//		if(varType==INS){return call==null ? 0 : call.length();}
//		if(varType==DELINS){return call==null ? (endLoc-beginLoc+1) : max(call.length(), endLoc-beginLoc+1);}
//		return endLoc-beginLoc+1;
//	}
	
	/** Returns true if this is a point variation (insertion or reference point) */
	public final boolean isPoint(){
		return varType==INS || varType==REFPOINT;
	}
	
	/** Returns true if this variation matches the reference sequence */
	public final boolean isRef(){
		return varType==REF || varType==REFPOINT;
	}
	
	/** Returns true if this represents an actual sequence change from reference */
	public final boolean isTrueVariation(){
		return varType==SNP || varType==INS || varType==DEL || varType==DELINS;
	}
	
	/** Returns true if this variation represents a no-call or uncertain region */
	public final boolean isNoCall(){
		return varType==NOCALL || varType==REFCON || varType==REFINCON;
	}
	
	/** Returns true if this is a no-ref or no-call variation type */
	public final boolean isNR_or_NC(){
		return varType==NOCALL || varType==NOREF || varType==REFCON || varType==REFINCON;
	}
	
	/** Returns true if this variation has uncertain or ambiguous status */
	public final boolean isUnsureVariation(){
		return varType==NOCALL || varType==NOREF || varType==REFINCON || varType==REFCON;
	}
	
	
//	/** TODO May be slow.  Perhaps add a boolean field. */
//	public boolean isCoding(){
//		int middle=((beginLoc+endLoc)/2);
//		GeneSet[] sets=Data.getNearestGeneSets(chromosome, middle);
//		for(GeneSet gs : sets){
//			for(Gene g : gs.genes){
//				if(g.intersectsCodeAndExon(beginLoc, endLoc)){
//					return true;
//				}
//			}
//		}
//		return false;
//	}

	
	/** Does this variation intersect within (range) of a coding region or splice site? */
	public boolean isNearCodingOrSplice(int range, boolean includeExonsForUntranslatedGenes, boolean includeSplice){
		assert(beginLoc<=endLoc);
		int a=beginLoc-range, b=endLoc+range;
		return isNearCodingOrSplice(range, includeExonsForUntranslatedGenes, Data.getNearestGeneSets(chromosome, a, b), includeSplice);
	}

	
	/** Does this variation intersect within (range) of a coding region or splice site? */
	public boolean isNearCodingOrSplice(int range, boolean includeExonsForUntranslatedGenes){
		return isNearCodingOrSplice(range, includeExonsForUntranslatedGenes, true);
	}
	
	
	/** Does this variation lie at least partially within an intron? */
	public boolean intersectsIntron(){
		assert(beginLoc<=endLoc);
		int a=beginLoc, b=endLoc;
		return intersectsIntron(Data.getNearestGeneSets(chromosome, a, b));
	}
	
	
	/** Does this variation intersect within (range) of a coding region or splice site? */
	public boolean isNearCodingOrSplice(int range, boolean includeExonsForUntranslatedGenes, GeneSet[] sets, boolean includeSplice){
		assert(beginLoc<=endLoc);
		int a=beginLoc-range, b=endLoc+range;
		
//		int middle=((beginLoc+endLoc)/2);
//		GeneSet[] sets=Data.getNearestGeneSets(chromosome, middle);

//		boolean flag=(chromosome==21 && intersects(9929078));//TODO UNDO
		
//		if(flag){System.out.println("Found: "+Arrays.toString(sets));}
//		assert(false);
		
		if(sets==null){
			assert(chromosome>=25);
			return true;
		}
		
		
		for(GeneSet gs : sets){
//			if(flag){System.out.println("### "+gs);}//TODO UNDO
			for(Gene g : gs.genes){
				
				if(!g.untranslated){

//					if(flag){System.out.println("*** "+g);}//TODO UNDO

//					if(flag){
//						System.out.println("intersectsCodeAndExon: "+g.intersectsCodeAndExon(a, b));
//						System.out.println("intersectsCode: "+g.intersectsCode(a, b));
//						System.out.println("intersectsExon: "+g.intersectsExon(a, b));
//					}

					if(g.intersectsCodeAndExon(a, b)){
						return true;
					}
					
				}else if(includeExonsForUntranslatedGenes){
//					if(flag){System.out.println("*** "+g);}//TODO UNDO
//
//					if(flag){
//						System.out.println("intersectsExon: "+g.intersectsExon(a, b));
//					}

					if(g.intersectsExon(a, b)){
						return true;
					}
					
				}
				
				if(includeSplice){
					int[] array=g.nearestSpliceSite(beginLoc, endLoc);
					if(array[0]<=range){return true;}
				}
				
			}
		}
		return false;
	}
	
	
	/** Does this variation lie at least partially within an intron? */
	public boolean intersectsIntron(GeneSet[] sets){
		assert(beginLoc<=endLoc);
		int a=beginLoc, b=endLoc;
		
//		int middle=((beginLoc+endLoc)/2);
//		GeneSet[] sets=Data.getNearestGeneSets(chromosome, middle);
		
		if(sets==null){
			assert(chromosome>=25);
			return true;
		}
		
		
		for(GeneSet gs : sets){
			for(Gene g : gs.genes){
				if(g.intersectsIntron(a, b)){return true;}
			}
		}
		return false;
	}
	
	/** Start coordinate of the variation (0-based, inclusive) */
	public int beginLoc=-2;
	/** End coordinate of the variation (0-based, inclusive) */
	public int endLoc=-2;
	
	/** Chromosome number (0-based indexing) */
	public int chromosome=-1;
	/** Variation type constant (REF, SNP, INS, DEL, etc.) */
	public byte varType=-1;
	
	/** Reference sequence at this genomic location */
	public String ref=null;
	/** Called/alternate sequence at this genomic location */
	public String call=null;
	
	
	/** Static mapping between ploidy values and their string representations */
	public static final HashMap<Object, Object> ploidyMap=makePloidyMap();
	/** Array mapping haplotype indices to string representations */
	public static final String[] haploMap={"0","1","2","all"};
	
	/** Array mapping variation type constants to string names */
	public static final String[] varTypeMap={"ref","snp","ins","del","sub",
		"no-call-rc","no-call-ri","no-call","no-ref","PAR-called-in-X","null","refpoint"};
	
	/** Reverse mapping from variation type strings to byte constants */
	public static final HashMap<String, Byte> varTypeMap2=makeVarTypeMap();
	
	/**
	 * Creates the reverse mapping from variation type strings to byte constants.
	 * Includes aliases and alternative names for variation types.
	 * @return HashMap mapping string names to variation type constants
	 */
	private static final HashMap<String, Byte> makeVarTypeMap(){
		HashMap<String, Byte> r=new HashMap<String, Byte>(32);
		
		for(byte i=0; i<varTypeMap.length; i++){r.put(varTypeMap[i], i);}
		r.put("=", REF);
		r.put("ref-consistent", REFCON);
		r.put("ref-inconsistent", REFINCON);
//		r.put("no-call-rc", REFCON);
//		r.put("no-call-ri", REFINCON);
		r.put("delins", DELINS);
		
		return r;
	}

	/** Constant for reference-matching variations */
	public static final byte REF=0;
	/** Constant for single nucleotide polymorphisms */
	public static final byte SNP=1;
	/** Constant for insertion variations */
	public static final byte INS=2;
	/** Constant for deletion variations */
	public static final byte DEL=3;
	/** Constant for deletion-insertion (complex substitution) variations */
	public static final byte DELINS=4;
	/** Constant for reference-consistent no-call variations */
	public static final byte REFCON=5;
	/** Constant for reference-inconsistent no-call variations */
	public static final byte REFINCON=6;
	/** Constant for general no-call variations */
	public static final byte NOCALL=7;
	/** Constant for no-reference variations */
	public static final byte NOREF=8;
	/** Constant for pseudoautosomal region variations */
	public static final byte PAR=9;
	/** Constant for null/undefined variations */
	public static final byte NULL=10;
	/** Constant for reference point variations (zero-length) */
	public static final byte REFPOINT=11;
	
	/** Interns reference and call strings to reduce memory usage.
	 * Uses Data.intern() to share common sequence strings across variations. */
	public void intern(){
//		assert(false) : ref+", "+call+", "+this;
		if(ref!=null){ref=Data.intern(ref);}
		if(call!=null){call=Data.intern(call);}
	}
	
	/**
	 * Creates bidirectional mapping between ploidy values and string representations.
	 * Supports both byte and integer keys with string values and reverse mapping.
	 * @return HashMap containing ploidy mappings
	 */
	private static HashMap<Object, Object> makePloidyMap(){
		HashMap<Object, Object> hashy=new HashMap<Object, Object>(64);
		for(int i=0; i<10; i++){
			hashy.put((Byte)(byte)i, i+"");
			hashy.put((Integer)i, i+"");
			hashy.put(i+"", (Byte)(byte)i);
		}
		hashy.put((Byte)(byte)-1, "?");
		hashy.put((Integer)(-1), "?");
		hashy.put("?",(Byte)(byte)-1);
		return hashy;
	}
	
	/** Returns the smaller of two integers */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integers */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	@Override
	public final int hashCode(){
		long x=chromosome;
		x=x<<4;
		x^=varType;
		x=x<<28;
		x^=beginLoc;
		x=x<<16;
		x^=(endLoc-beginLoc+1);
		return (int)(x^(x>>>32));
	}
	
	@Override
	public int compareTo(Variation other) {
		if(chromosome!=other.chromosome){return other.chromosome>chromosome ? -1 : 1;}
		if(beginLoc!=other.beginLoc){return other.beginLoc>beginLoc ? -1 : 1;}
		if(endLoc!=other.endLoc){return other.endLoc>endLoc ? -1 : 1;}
		if(varType!=other.varType){return other.varType>varType ? -1 : 1;}
		if(varType==REF || varType==NOCALL){return 0;}
		
		if(call==null){
			return other.call==null ? 0 : -1;
		}
		return other.call==null ? 1 : call.compareTo(other.call);
	}
	
	@Override
	public boolean equals(Object other){
		return equals((Variation)other);
	}
	
	/**
	 * Tests equality with another variation using compareTo for consistency.
	 * @param other Variation to compare against
	 * @return true if variations are equivalent
	 */
	public boolean equals(Variation other){
		return compareTo(other)==0;
	}
	
	/**
	 * Tests if a genomic coordinate intersects this variation's span.
	 * @param point Genomic coordinate to test
	 * @return true if point falls within this variation's coordinates
	 */
	public boolean intersects(int point){
		return point>=beginLoc && point<=endLoc;
	}
	
	/**
	 * Tests if a genomic coordinate is adjacent to this variation.
	 * Includes positions one base before or after the variation.
	 * @param point Genomic coordinate to test
	 * @return true if point is within one base of this variation
	 */
	public boolean touches(int point){
		return point>=beginLoc-1 && point<=endLoc+1;
	}
	
	/** This is quite clever. */
	public static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	/**
	 * Tests if two genomic intervals overlap or are adjacent.
	 * Extends overlap detection to include touching intervals.
	 *
	 * @param a1 Start of first interval
	 * @param b1 End of first interval
	 * @param a2 Start of second interval
	 * @param b2 End of second interval
	 * @return true if intervals overlap or touch
	 */
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=(b1+1) && b2>=(a1-1);
	}
	
	/** Is (a1, b1) within (a2, b2) ? */
	public static boolean isWithin(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2;
	}
	
	/**
	 * Tests if the first interval is contained within the second without touching edges.
	 * Requires strict containment with gaps at both ends.
	 *
	 * @param a1 Start of first interval
	 * @param b1 End of first interval
	 * @param a2 Start of second interval
	 * @param b2 End of second interval
	 * @return true if first interval is strictly within second interval
	 */
	public static boolean isWithinNotTouching(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>a2 && b1<b2;
	}
	
	//Slow if not inlined
	/**
	 * Tests if this variation overlaps with the specified genomic interval.
	 * @param a2 Start coordinate of interval
	 * @param b2 End coordinate of interval
	 * @return true if this variation overlaps the specified interval
	 */
	public boolean intersects(int a2, int b2){return overlap(beginLoc, endLoc, a2, b2);}

	/**
	 * Tests if this variation is completely contained within the specified interval.
	 * @param a2 Start coordinate of containing interval
	 * @param b2 End coordinate of containing interval
	 * @return true if this variation is within the specified interval
	 */
	public boolean isWithin(int a2, int b2){return isWithin(beginLoc, endLoc, a2, b2);}
	
	/**
	 * Tests if this variation is strictly contained within the specified interval.
	 * @param a2 Start coordinate of containing interval
	 * @param b2 End coordinate of containing interval
	 * @return true if this variation is strictly within the specified interval
	 */
	public boolean isWithinNotTouching(int a2, int b2){return isWithinNotTouching(beginLoc, endLoc, a2, b2);}
	
	/**
	 * Complex method to determine if two variations intersect or overlap.
	 * Handles special cases for different variation types including point mutations,
	 * insertions, deletions, and no-call regions with sophisticated overlap logic.
	 *
	 * @param v Variation to test intersection against
	 * @return true if variations intersect or overlap
	 */
	public boolean intersects(Variation v){
		
		if(v.chromosome!=chromosome){
			return false;
		}
		
		int len1=lengthRef();
		int len2=v.lengthRef();
		
		if(len1<len2){
			return v.intersects(this);
		}
		
//		if(v.beginLoc==46397336 || v.beginLoc==46397348){
//			System.err.println(len1+": "+this+"\n"+len2+": "+v+"\n");
//		}
		//Now, this is at least as long (ref-wise) as v.
		
//		if(varType==EQUAL && v.varType==INS){
//			assert(false);
//		}

		if(!touch(beginLoc, endLoc, v.beginLoc, v.endLoc)){return false;}
		
//		if(v.beginLoc==46397336 || v.beginLoc==46397348){
//			System.err.println("Touch("+beginLoc+", "+endLoc+", "+v.beginLoc+", "+v.endLoc+")");
//		}
		
		if(v.isPoint()){
//			if(v.beginLoc==46397336 || v.beginLoc==46397348){System.out.println("v");}
			if(isPoint()){
//				assert(beginLoc==v.beginLoc && endLoc==v.endLoc) : this+"\n"+v;
//				return true;

				return beginLoc==v.beginLoc;
			}
//			if(v.beginLoc==46397336 || v.beginLoc==46397348){System.out.println("w");}
			if(this.isNoCall()){
				if(len1>0){return overlap(beginLoc, endLoc, v.beginLoc, v.endLoc);} //Normal case
				else{
					//TODO: Bad news! Original MAY have been a length 0 no-call in half-open coordinates.
					return overlap(beginLoc, endLoc+1, v.beginLoc, v.endLoc);
				}
			}
//			if(v.beginLoc==46397336 || v.beginLoc==46397348){System.out.println("x");}
			if(v.beginLoc<=beginLoc){return false;}
//			if(v.beginLoc==46397336 || v.beginLoc==46397348){System.out.println("y");}
		}
//		if(v.beginLoc==46397336 || v.beginLoc==46397348){System.out.println("z");}
		
		return overlap(beginLoc, endLoc, v.beginLoc, v.endLoc);
	}
}
