package prok;

import stream.Read;

/** Reusable worker-local view of one sense-oriented flanked learning record.
 * Every record is joined to the original population, including excluded loci.
 * Missing geometry or GC is an explicit learning limitation; contradictory
 * identifiers or declared geometry fail rather than silently dropping a locus.
 * The Read and its bases are borrowed only until the next set call.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningContext{

	/** Binds one record; metadata keys may be separated by whitespace or semicolons. */
	public void set(final Read read_, final NcrnaLearningPopulation population){
		if(read_==null || read_.id==null || read_.mate!=null || population==null){
			throw new IllegalArgumentException("Learning contexts need unpaired named records and the original population");
		}
		read=read_;
		final String header=read.id;
		final String id=NcrnaSupport.recordName(header);
		gene=population.index(id);
		if(gene<0){throw bad("Unknown original locus");}
		start=stop=-1;
		gc=Float.NaN;
		status=NO_FLANK;
		int left=-1, right=-1, seen=0;
		for(int from=id.length(); from<header.length();){
			while(from<header.length() && separator(header.charAt(from))){from++;}
			int to=from;
			while(to<header.length() && !separator(header.charAt(to))){to++;}
			int flag=0;
			if(header.startsWith("source=", from)){
				flag=1;
				final String expected=population.sourceName(gene);
				if(to-from-7!=expected.length() || !header.regionMatches(from+7, expected, 0, expected.length())){
					throw bad("Context source disagrees with original locus");
				}
			}else if(header.startsWith("lflank=", from)){
				flag=2;
				left=integer(header, from+7, to);
			}else if(header.startsWith("rflank=", from)){
				flag=4;
				right=integer(header, from+7, to);
			}else if(header.startsWith("contig_gc=", from)){
				flag=8;
				gc=fraction(header, from+10, to);
			}else{
				final int equals=header.indexOf('=', from);
				if(equals>=from && equals<to && !legacyMetadata(header, from)){
					throw bad("Unknown context metadata key "+header.substring(from, equals));
				}
			}
			if((seen&flag)!=0){throw bad("Duplicate context metadata");}
			seen|=flag;
			from=to;
		}
		if((seen&1)==0){throw bad("Missing original source metadata");}
		if((left>=0 && left!=population.leftFlank[gene]) || (right>=0 && right!=population.rightFlank[gene])){
			throw bad("Supplied flank length disagrees with original locus");
		}
		if(left<0 || right<0){return;}
		if(read.bases==null || (long)read.length()-left-right<16){
			status=MALFORMED_SPAN;
			return;
		}
		if(read.length()-left-right!=population.length[gene]){
			throw bad("Context flank/span geometry disagrees with original locus");
		}
		start=left;
		stop=read.length()-right-1;
		status=Float.isNaN(gc) ? NO_GC : AVAILABLE;
		assert(start>=0 && stop<read.length() && stop-start+1==population.length[gene])
			: "Sense coordinates are L through N-R-1 and must retain the original gene length: gene="+gene+
			", start="+start+", stop="+stop+", contextLength="+read.length()+", geneLength="+population.length[gene];
	}

	/** True when the flanked sequence has a valid gene span, even if GC is unavailable. */
	public boolean hasGeometry(){return status==AVAILABLE || status==NO_GC;}
	/** True only when both span geometry and whole-contig GC are available. */
	public boolean hasFeatures(){return status==AVAILABLE;}

	/** Adds the original locus to strict decimal flank diagnostics. */
	private int integer(final String text, final int from, final int to){
		try{
			return NcrnaNumbers.unsignedInt(text, from, to);
		}catch(IllegalArgumentException failure){
			throw bad("Invalid flank length: "+failure.getMessage());
		}
	}

	/** Whole-token input includes legacy exponent notation; absent/NA/NaN is unavailable. */
	private float fraction(final String text, final int from, final int to){
		if(from==to || (to-from==2 && text.regionMatches(from, "NA", 0, 2))
				|| (to-from==3 && text.regionMatches(from, "NaN", 0, 3))){return Float.NaN;}
		final float value;
		try{
			// Parse's substring overload needs a process-wide switch for exponents.
			// One GC token per context is parsed here, never per candidate vector.
			value=Float.parseFloat(text.substring(from, to));
		}catch(NumberFormatException failure){throw bad("Malformed contig GC fraction");}
		if(!(value>=0 && value<=1)){throw bad("Contig GC is outside [0,1]");}
		return value;
	}

	private IllegalArgumentException bad(final String message){return new IllegalArgumentException(message+": "+read.id);}
	private static boolean separator(final char c){return c==';' || Character.isWhitespace(c);}

	/** Extractor metadata retained for provenance, without affecting learning features.
	 * Unkeyed description words remain legal; other key=value tokens are errors.
	 */
	private static boolean legacyMetadata(final String header, final int from){
		return header.startsWith("contig=", from) || header.startsWith("strand=", from)
			|| header.startsWith("start=", from) || header.startsWith("stop=", from)
			|| header.startsWith("family=", from) || header.startsWith("status=", from);
	}

	/** Borrowed current context, valid until the next set call. */
	public Read read;
	/** Original-population locus index for this context. */
	public int gene=-1;
	/** Zero-based inclusive sense-oriented gene start in read.bases. */
	public int start=-1;
	/** Zero-based inclusive sense-oriented gene stop in read.bases. */
	public int stop=-1;
	/** Availability code describing span geometry and GC. */
	public int status=NO_FLANK;
	/** Whole-contig GC fraction, or NaN when unavailable. */
	public float gc=Float.NaN;
	/** Valid span geometry and a finite whole-contig GC fraction. */
	public static final int AVAILABLE=0;
	/** At least one flank length is missing; gene coordinates remain unavailable. */
	public static final int NO_FLANK=1;
	/** No bases, or fewer than sixteen bases remain between declared flanks. */
	public static final int MALFORMED_SPAN=2;
	/** Valid span geometry but no usable whole-contig GC value. */
	public static final int NO_GC=3;
}
