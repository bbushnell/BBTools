package var2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import fileIO.FileFormat;
import ml.CellNet;
import shared.KillSwitch;
import shared.LineParser2;
import shared.Parse;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/**
* Represents a single genomic variant with associated quality metrics and statistics.
* Stores variant position, type, allele information, and comprehensive read-level data
* for variant calling quality assessment and filtering.
* 
* Core variant calling class that accumulates evidence from multiple reads,
* calculates statistical scores for variant confidence, and outputs results
* in standard formats (VAR, VCF). Handles all major variant types including
* substitutions, insertions, deletions, and junction variants.
* 
* Key features:
* - Multi-threaded variant accumulation from read alignments
* - Sophisticated statistical scoring algorithms for quality assessment
* - Bias detection (strand bias, read bias, positional bias)
* - Homopolymer and repeat region analysis
* - Support for multiple output formats with comprehensive metadata
* 
* @author Brian Bushnell
* @contributor Isla
* @date November 4, 2016
*/
public class Var implements Comparable<Var>, Serializable, Cloneable {

	private static final long serialVersionUID=3328626403863586829L;

	/**
	 * Main method for testing and loading variant files.
	 * Supports both VAR and VCF format input files for validation and analysis.
	 * 
	 * @param args Command line arguments: [0] = input file path, [1+] = optional parameters
	 */
	public static void main(String[] args){
		if(args.length>1){
			for(int i=1; i<args.length; i++){			
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				Parser.parseCommonStatic(arg, a, b);
			}
		}
		Timer t=new Timer();
		VarMap vmap;
		FileFormat ff=FileFormat.testInput(args[0], FileFormat.VAR, null, true, true);
		if(ff.vcf()){
			ScafMap smap=ScafMap.loadVcfHeader(ff);
			vmap=VcfLoader.loadVcfFile(args[0], smap, false, false);
		}else{
			vmap=VcfLoader.loadVarFile(args[0], null);
		}
		t.stop("Loaded "+vmap.size()+" variants.\nTime: \t");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Var clone(){
		try{
			return (Var)super.clone();
		}catch(CloneNotSupportedException e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	/**
	 * Creates variant from basic parameters using allele as integer ASCII code.
	 * Convenience constructor that converts single-character alleles from ASCII codes
	 * to byte arrays using the pre-allocated AL_MAP lookup table. Used primarily
	 * for single-base substitutions where the allele is known as a character code.
	 * 
	 * @param scafnum_ Scaffold/chromosome number (0-based indexing)
	 * @param start_ Start position in reference coordinates (0-based, inclusive)
	 * @param stop_ Stop position in reference coordinates (0-based, exclusive)
	 * @param allele_ Single-character allele as ASCII integer code
	 * @param type_ Variant type constant (SUB, INS, DEL, etc.)
	 */
	public Var(int scafnum_, int start_, int stop_, int allele_, int type_){
		this(scafnum_, start_, stop_, AL_MAP[allele_], type_);
	}

	/**
	 * Copy constructor - creates new variant with identical properties.
	 * @param v Source variant to copy
	 */
	public Var(Var v){
		this(v.scafnum, v.start, v.stop, v.allele, v.type);
	}
	
	/**
	* Primary constructor for creating variants.
	* Validates coordinates and allele data, computes hash code for efficient storage.
	* 
	* @param scafnum_ Scaffold/chromosome number
	* @param start_ Start position (0-based, inclusive)
	* @param stop_ Stop position (0-based, exclusive) 
	* @param allele_ Allele sequence as byte array
	* @param type_ Variant type (SUB, INS, DEL, etc.)
	*/
	public Var(int scafnum_, int start_, int stop_, byte[] allele_, int type_){
		scafnum=scafnum_;
		start=start_;
		stop=stop_;
		allele=allele_;
		type=type_;
		hashcode=hash();
		
		// Validate allele is either empty, single base, or multi-base sequence
		assert(allele.length>1 || allele==AL_0 ||
				allele==AL_A || allele==AL_C || allele==AL_G || allele==AL_T || allele==AL_N) : 
				new String(allele_)+", "+allele_.length;
		assert(start<=stop) : "\n"+VarHelper.toBasicHeader()+"\n"+this+"\n";
	}
	
	/**
	 * Constructs variant by parsing delimited text line from VAR format file.
	 * Parses all 24+ fields including position, type, counts, and quality statistics.
	 * Used for loading existing variant data from disk storage.
	 * 
	 * @param line Byte array containing tab-delimited variant data
	 * @param delimiter Field separator character (typically tab)
	 */
	public Var(final byte[] line, final byte delimiter){
	    LineParser2 lp=new LineParser2(delimiter);
	    lp.set(line);
	    
	    // Parse core variant information
	    scafnum=lp.parseInt(0); // Field 0: Scaffold number
	    start=lp.parseInt(1); // Field 1: Start position (0-based)
	    stop=lp.parseInt(2); // Field 2: Stop position (exclusive)
	    type=typeInitialArray[lp.parseByte(3, 0)]; // Field 3: Variant type
	    
	    // Field 4: Allele sequence with special handling for empty/single bases
	    int alleleLen=lp.length(4);
	    if(alleleLen==0){
	        allele=AL_0; // Empty allele for deletions
	    }else if(alleleLen==1){
	        allele=AL_MAP[lp.parseByte(4, 0)]; // Single base - use pre-allocated array
	    }else{
	        allele=lp.parseByteArray(4); // Multi-base sequence
	    }
	    
	    // Parse read count statistics by strand and read pair
	    r1plus=lp.parseInt(5); // Field 5: Read 1 plus strand count
	    r1minus=lp.parseInt(6); // Field 6: Read 1 minus strand count  
	    r2plus=lp.parseInt(7); // Field 7: Read 2 plus strand count
	    r2minus=lp.parseInt(8); // Field 8: Read 2 minus strand count
	    properPairCount=lp.parseInt(9); // Field 9: Proper pair count
	    
	    // Parse quality and mapping statistics  
	    lengthSum=lp.parseLong(10); // Field 10: Sum of supporting read lengths
	    mapQSum=lp.parseLong(11); // Field 11: Sum of mapping qualities
	    mapQMax=lp.parseInt(12); // Field 12: Maximum mapping quality
	    baseQSum=lp.parseLong(13); // Field 13: Sum of base qualities
	    baseQMax=lp.parseInt(14); // Field 14: Maximum base quality
	    
	    // Parse positional and identity statistics
	    endDistSum=lp.parseLong(15); // Field 15: Sum of distances from read ends
	    endDistMax=lp.parseInt(16); // Field 16: Maximum distance from read ends
	    idSum=lp.parseLong(17); // Field 17: Sum of alignment identities  
	    idMax=lp.parseInt(18); // Field 18: Maximum alignment identity
	    
	    // Parse coverage and annotation data
	    coverage=lp.parseInt(19); // Field 19: Total coverage at position
	    minusCoverage=lp.parseInt(20); // Field 20: Minus strand coverage
	    nearbyVarCount=lp.parseInt(21); // Field 21: Count of nearby variants
	    flagged=lp.parseInt(22)>0; // Field 22: Flagged status (boolean)
	    
	    // Fields 23-24 (contig end distance, phred score) are parsed but not stored
	    // These are calculated dynamically when needed
	    
	    // Compute hash and validate object state
	    hashcode=hash();
	    assert(allele.length>1 || allele==AL_0 || 
	           allele==AL_A || allele==AL_C || allele==AL_G || 
	           allele==AL_T || allele==AL_N);
	    assert(start<=stop) : this.toString();
	    assert(type>=LJUNCT || type==type_old()) : type+", "+type_old()+", "+this.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resets all accumulated statistics to initial state while preserving variant identity.
	 * Clears read counts, quality metrics, and flags but maintains position and allele data.
	 * Used when reprocessing variants or clearing cached statistics.
	 * 
	 * @return This Var object for method chaining
	 */
	public Var clear(){
		// Reset coverage tracking
		coverage=-1;
		minusCoverage=-1;
		
		// Clear all read count statistics
		r1plus=0;
		r1minus=0;
		r2plus=0;
		r2minus=0;
		properPairCount=0;
		
		// Reset mapping quality statistics
		mapQSum=0;
		mapQMax=0;
		
		// Reset base quality statistics
		baseQSum=0;
		baseQMax=0;
		
		// Reset positional statistics
		endDistSum=0;
		endDistMax=0;
		
		// Reset identity statistics
		idSum=0;
		idMax=0;
		
		lengthSum=0;
		
		// Reset derived statistics and flags
		revisedAlleleFraction=-1;
		nearbyVarCount=-1;
		forced=false;
		flagged=false;
		return this;
	}

	/**
	 * Adds coverage data from another variant (used in multi-sample processing).
	 * Only merges coverage statistics, not read-level data or quality metrics.
	 * Variants must be equivalent (same position and allele).
	 * 
	 * @param b Source variant to merge coverage from
	 */
	public void addCoverage(Var b){
		assert(this.equals(b));
		coverage+=b.coverage; // Add total coverage
		minusCoverage+=b.minusCoverage; // Add minus strand coverage
	}

	/**
	 * Merges complete statistical data from another equivalent variant.
	 * Combines read counts, quality metrics, and positional statistics.
	 * Used when consolidating variants from multiple processing threads.
	 * 
	 * @param b Source variant with data to merge
	 */
	public void add(Var b){
		final int oldReads=alleleCount();
		
		// Validate current state before merging
		assert(oldReads==0 || baseQSum/oldReads<=60) : this;
		assert(this.equals(b)); // Must be same variant
		
		// Merge read count statistics by strand and pair
		r1plus+=b.r1plus;
		r1minus+=b.r1minus;
		r2plus+=b.r2plus;
		r2minus+=b.r2minus;
		properPairCount+=b.properPairCount;
		lengthSum+=b.lengthSum;
		
		// Merge mapping quality statistics (sum and max)
		mapQSum+=b.mapQSum;
		mapQMax=Tools.max(mapQMax, b.mapQMax);
		
		// Merge base quality statistics (sum and max)
		baseQSum+=b.baseQSum;
		baseQMax=Tools.max(baseQMax, b.baseQMax);

		// Merge positional statistics (distance from read ends)
		endDistSum+=b.endDistSum;
		endDistMax=Tools.max(endDistMax, b.endDistMax);

		// Merge identity statistics
		idSum+=b.idSum;
		idMax=Tools.max(idMax, b.idMax);

		// Validate merged state
		assert(alleleCount()>=oldReads) : "\n"+this+"\n"+b;
		assert(alleleCount()==oldReads+b.alleleCount()) : "\n"+this+"\n"+b;
		assert(alleleCount()==0 || baseQSum/alleleCount()<=60) : "\n"+this+"\n"+b;
	}

	/**
	 * Adds evidence from a single read supporting this variant.
	 * Calculates read-specific statistics and updates accumulated metrics.
	 * 
	 * @param r Read containing this variant
	 */
	public void add(Read r){
		final SamLine sl=r.samline;
		final int bstart=calcBstart(r, sl); // Find variant start in read
		final int bstop=calcBstop(bstart, r); // Find variant end in read
		add(r, bstart, bstop);
	}

	/**
	 * Adds evidence from a read with pre-calculated variant boundaries.
	 * Updates strand counts, quality metrics, and positional statistics.
	 * Core method for accumulating variant evidence from aligned reads.
	 * 
	 * @param r Read supporting this variant
	 * @param bstart Start position of variant within read bases
	 * @param bstop Stop position of variant within read bases  
	 */
	public void add(Read r, final int bstart, final int bstop){
		final int oldReads=alleleCount();
		final SamLine sl=r.samline;
		
		// Update strand-specific read counts
		if(sl.strand()==0){ // Plus strand
			if(sl.pairnum()==0){
				r1plus++; // Read 1 plus
			}else{
				r2plus++; // Read 2 plus
			}
		}else{ // Minus strand
			if(sl.pairnum()==0){
				r1minus++; // Read 1 minus
			}else{
				r2minus++; // Read 2 minus
			}
		}
		
		// Update read-level statistics
		lengthSum+=r.length();
		properPairCount+=(sl.properPair() ? 1 : 0);
		
		// Update mapping quality statistics
		mapQSum+=sl.mapq;
		mapQMax=Tools.max(mapQMax, sl.mapq);
		
		// Calculate and update base quality for this variant
		int baseQ=calcBaseQ(bstart, bstop, r, sl);
		baseQSum+=baseQ;
		baseQMax=Tools.max(baseQMax, baseQ);
		
		// Calculate and update distance from read ends
		int endDist=calcEndDist(bstart, bstop, r);
		endDistSum+=endDist;
		endDistMax=Tools.max(endDistMax, endDist);
		
		// Calculate and update alignment identity
		int id=(int)(1000*Read.identitySkewed(r.match, false, false, false, true));
		idSum+=id;
		idMax=Tools.max(idMax, id);
		
		// Validate updated state
		assert(alleleCount()>0) : this;
		assert(alleleCount()==oldReads+1) : this;
		assert(baseQSum/alleleCount()<=60) : this;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts reads to variant list using scaffold mapping for coordinate resolution.
	 * Convenience wrapper that resolves scaffold name to number.
	 * 
	 * @param r Read to analyze for variants
	 * @param sl SAM alignment line
	 * @param callNs Whether to call variants at N positions
	 * @param scafMap Scaffold mapping for name resolution
	 * @return List of variants found in the read, or null if none
	 */
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, ScafMap scafMap){
		final int scafnum=scafMap.getNumber(sl.rnameS());
		return toVars(r, sl, callNs, scafnum);
	}

	/**
	 * Extracts all variants from a read alignment.
	 * Processes match string to identify substitutions, insertions, deletions, and junctions.
	 * Handles read orientation and creates appropriate Var objects for each variant type.
	 * 
	 * @param r Read with alignment information
	 * @param sl SAM line with mapping details
	 * @param callNs Whether to report variants at N positions in reference
	 * @param scafnum Scaffold number for variant coordinates
	 * @return List of Var objects representing all variants in the read
	 */
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, final int scafnum){
		final boolean hasV=r.containsVariants(); // Check for substitutions/indels
		final boolean callSID=(CALL_DEL || CALL_INS || CALL_SUB) && hasV;
		final boolean callJ=(CALL_JUNCTION) && (hasV || r.containsClipping());
		if(!callSID && !callJ){return null;}
		
		// Prepare read for variant analysis
		r.toLongMatchString(false);
		if(sl.strand()==1 && !r.swapped()){
			r.reverseComplement(); // Orient read to match reference
			r.setSwapped(true);
		}
		
		// Extract different variant types
		ArrayList<Var> sidList=null, jList=null;
		if(callSID){sidList=toSubsAndIndels(r, sl, callNs, scafnum);}
		if(callJ){jList=VarHelper.toJunctions(r, sl, scafnum, hasV, 8);}
		
		// Combine results
		if(sidList!=null){
			if(jList!=null){sidList.addAll(jList);}
			return sidList;
		}else{
			return jList;
		}
	}

	/**
	 * Extracts substitutions, insertions, and deletions from read alignment.
	 * Parses match string character by character to identify variant positions.
	 * Creates Var objects for each variant with appropriate coordinates and alleles.
	 * 
	 * @param r Read with variant information
	 * @param sl SAM alignment line
	 * @param callNs Whether to call variants at N positions
	 * @param scafnum Scaffold number for coordinates
	 * @return List of substitution and indel variants
	 */
	private static ArrayList<Var> toSubsAndIndels(Read r, SamLine sl, boolean callNs, final int scafnum){
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final int rpos0=sl.pos-1; // Reference start position (0-based)
		ArrayList<Var> list=new ArrayList<Var>();

		// Track current variant state
		int bstart=-1, bstop=-1; // Base positions in read
		int rstart=-1, rstop=-1; // Reference positions
		int mode=-1; // Current match state
		
		// Parse match string to find variants
		int mpos=0, bpos=0, rpos=rpos0;
		for(; mpos<match.length; mpos++){
			byte m=match[mpos];
			
			// Handle end of current variant
			if(m!=mode){
				if(mode=='D'){ // End of deletion
					bstop=bpos;
					rstop=rpos;
					if(CALL_DEL){
						Var v=new Var(scafnum, rstart, rstop, 0, DEL);
						v.add(r, bstart, bstop);
						list.add(v);
					}
					bstart=bstop=rstart=rstop=-1;
				}else if(mode=='I'){ // End of insertion
					bstop=bpos;
					rstop=rpos;
					int blen=bstop-bstart;
					if(CALL_INS){
						Var v;
						if(blen==1){
							v=new Var(scafnum, rstart, rstop, bases[bstart], INS);
						}else{
							v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop), INS);
						}
						v.add(r, bstart, bstop);
						list.add(v);
					}
					bstart=bstop=rstart=rstop=-1;
				}
			}
			
			// Process current match character
			if(m=='C'){ // Clipping (soft/hard)
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){ // Match, substitution, or N
				if(m=='S' || (m=='N' && callNs)){
					if(CALL_SUB){
						Var v=new Var(scafnum, rpos, rpos+1, bases[bpos], SUB);
						v.add(r, bpos, bpos+1);
						list.add(v);

						// Validate substitution against reference if available
						if(TEST_REF_VARIANTS && v.type()==SUB && v.allele.length==1){
							final byte call=v.allele[0];
							final Scaffold scaf=ScafMap.defaultScafMap().getScaffold(scafnum);
							final byte ref=scaf.bases[v.start];
							assert(ref!=call) : (char)call+"="+(char)ref+" at scaf "+scafnum+" pos "+v.start+"\n"
							+sl+"\n"+ScafMap.defaultScafMap().getScaffold(scafnum).getSequence(sl)+"\n";
						}
					}
				}
				bpos++;
				rpos++;
			}else if(m=='D'){ // Deletion
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				rpos++;
			}else if(m=='I'){ // Insertion
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
			mode=m;
		}
		
		// Handle variant at end of match string
		if(mode=='D'){
			bstop=bpos;
			rstop=rpos;
			if(CALL_DEL){
				Var v=new Var(scafnum, rstart, rstop, 0, DEL);
				v.add(r, bstart, bstop);
				list.add(v);
			}
		}else if(mode=='I'){
			bstop=bpos;
			rstop=rpos-1;
			int blen=bstop-bstart;
			if(CALL_INS){
				Var v;
				assert(rstart<=rstop) : "\n"+rstart+", "+rstop+", "+rpos+
				"\n"+bstart+", "+bstop+", "+bpos+
				"\n"+r+"\n"+sl;
				if(blen==1){
					v=new Var(scafnum, rstart, rstop, bases[bstart], INS);
				}else{
					v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop), INS);
				}
				v.add(r, bstart, bstop);
				list.add(v);
			}
		}
		
		return list;
	}

	/**
	 * Calculates the starting position of this variant within the read bases.
	 * Parses match string to find where the variant begins in the read sequence.
	 * 
	 * @param r Read containing the variant
	 * @param sl SAM line with alignment information
	 * @return Starting base position of variant in read (0-based)
	 */
	public int calcBstart(Read r, SamLine sl){
		r.toLongMatchString(false);
		byte[] match=r.match;
		final int rstart=sl.pos-1; // Reference start position
		final int type=type();
		
		int bstart=-1;
		
		// Parse match string to find variant position
		for(int mpos=0, rpos=rstart, bpos=0; mpos<match.length; mpos++){
			byte m=match[mpos];
			if(m=='C'){ // Clipping
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){ // Match/substitution/N
				if(rpos==rstart){ // Found variant position
					assert(type==SUB || type==NOCALL) : type+", "+bpos+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				bpos++;
				rpos++;
			}else if(m=='D'){ // Deletion
				if(rpos==rstart){
					assert(type==DEL) : type+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				rpos++;
			}else if(m=='I'){ // Insertion
				if(rpos==rstart && type==INS){
					bstart=bpos;
					break;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
		}
		assert(bstart>=0);
		return bstart;
	}

	/**
	 * Calculates the ending position of this variant within the read bases.
	 * 
	 * @param bstart Starting position of variant in read
	 * @param r Read containing the variant
	 * @return Ending base position of variant in read (exclusive)
	 */
	public int calcBstop(int bstart, Read r){
		assert(bstart>=0);
		int bstop=bstart+readlen(); // Add variant length
		assert(bstop<=r.length());
		return bstop;
	}

	/**
	 * Calculates distance from variant to nearest read end.
	 * Used for detecting variants near read ends which may be less reliable.
	 * 
	 * @param bstart Starting position of variant in read
	 * @param bstop Ending position of variant in read
	 * @param r Read containing the variant
	 * @return Distance to nearest read end
	 */
	public int calcEndDist(int bstart, int bstop, Read r){
		int dist=Tools.min(bstart, r.length()-bstop); // Distance to nearest end
		assert(dist>=0 && dist<=r.length()/2) : dist+", "+r.length()+", "+bstart+", "+bstop+"\n"+this+"\n"+
			"\n"+new String(r.match)+"\n"+r.samline+"\n";
		assert(dist<=(r.length()-readlen())/2) : "\ndist="+dist+", r.len="+r.length()+", readlen="+readlen()+", allele='"+new String(allele)+
			"', allele.length="+allele.length+"\n"+"allele2="+toString(allele)+"\n(r.length()-readlen())/2="+((r.length()-readlen())/2)+"bstart="+bstart+", bstop="+bstop+
			"\nthis="+this+"\nmatch="+new String(r.match)+"\nsamline="+r.samline+"\n";
		return dist;
	}

	/**
	 * Debug method for converting byte array to readable string.
	 * @param array Byte array to convert
	 * @return Comma-separated string of byte values
	 */
	String toString(byte[] array){
		StringBuilder sb=new StringBuilder();
		for(byte b : array){
			sb.append((int)b).append(',');
		}
		return sb.toString();
	}

	/**
	 * Calculates average base quality for this variant position.
	 * Handles strand orientation and special cases for deletions.
	 * 
	 * @param bstart0 Starting base position (forward orientation)
	 * @param bstop0 Ending base position (forward orientation)
	 * @param r Read with quality information
	 * @param sl SAM line with strand information
	 * @return Average base quality score for the variant
	 */
	public int calcBaseQ(final int bstart0, final int bstop0, Read r, SamLine sl){
		final byte[] quals=r.quality;
		if(quals==null){return Shared.FAKE_QUAL;} // No quality data available
		final int type=type();
		final int bstart, bstop;
		final int len=r.length();
		
		// Adjust coordinates for read orientation
		if(sl.strand()==0 || (sl.strand()==1 && r.swapped())){
			bstart=bstart0; // Forward orientation
			bstop=bstop0;
		}else{
			bstart=len-bstop0-1; // Reverse complement coordinates
			bstop=len-bstart0-1;
			assert(bstop-bstart==bstop0-bstart0);
		}
		
		int sum=0, avg=0;
		if(type==DEL){ // Special handling for deletions
			if(bstart==0){
				sum=avg=quals[0]; // Use first base quality
			}else if(bstop>=len-1){
				sum=avg=quals[len-1]; // Use last base quality
			}else{
				assert(bstop==bstart) : bstart0+", "+bstop0+", "+bstart+", "+bstop+"\n"+
						r.length()+", "+r.swapped()+", "+type()+", "+readlen()+", "+reflen()+
						"\n"+this+"\n"+new String(r.match)+"\n"+r.samline+"\n";
				sum=quals[bstart]+quals[bstop+1]; // Average flanking bases
				avg=sum/2;
			}
		}else{ // Sum qualities for variant bases
			for(int i=bstart; i<bstop; i++){
				sum+=quals[i];
			}
			avg=Math.round(sum/(bstop-bstart));
		}
		return avg;
	}

	/**
	 * Calculates reference sequence length affected by this variant.
	 * @return Number of reference bases spanned (stop - start)
	 */
	public int reflen(){
		return stop-start;
	}

	/**
	 * Calculates read sequence length of this variant's allele.
	 * @return Length of alternative allele sequence
	 */
	int readlen(){
		return (allele==null || allele.length==0 || allele[0]=='.' ? 0 : allele.length);
	}

	/**
	 * Gets the variant type constant for this variant.
	 * @return Type constant (SUB, INS, DEL, etc.)
	 */
	public int type(){return type;}

	/**
	 * Legacy method for calculating variant type from coordinates and allele.
	 * Used for validation and backward compatibility.
	 * @return Calculated variant type
	 */
	int type_old(){
		int reflen=reflen(), readlen=readlen();
		return typeReadlenReflen(readlen, reflen, allele);
	}

	/**
	 * Determines variant type from start/stop coordinates and allele.
	 * @param start Variant start position
	 * @param stop Variant stop position  
	 * @param allele Alternative allele sequence
	 * @return Variant type constant
	 */
	static int typeStartStop(int start, int stop, byte[] allele){
		final int readlen=(allele.length==0 || allele[0]=='.' ? 0 : allele.length);
		final int reflen=stop-start;
		return typeReadlenReflen(readlen, reflen, allele);
	}

	/**
	 * Determines variant type from read length, reference length, and allele content.
	 * Core typing logic used by multiple variant type calculation methods.
	 * 
	 * @param readlen Length of alternative allele
	 * @param reflen Length of reference sequence affected
	 * @param allele Alternative allele sequence
	 * @return Variant type constant (INS, DEL, SUB, or NOCALL)
	 */
	static int typeReadlenReflen(int readlen, int reflen, byte[] allele){
		if(reflen<readlen){return INS;} // Insertion: read longer than reference
		if(reflen>readlen){return DEL;} // Deletion: reference longer than read
		for(byte b : allele){ // Same length: check for substitutions
			if(b!='N'){return SUB;} // Non-N substitution
		}
		return NOCALL; // All N bases = no-call
	}

	/**
	 * Gets human-readable string representation of variant type.
	 * @return Type name (e.g., "SUB", "INS", "DEL")
	 */
	String typeString(){
		return typeArray[type()];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	* Tests equality with another object (must be a Var).
	* Uses optimized comparison with hash code pre-check for performance.
	* 
	* @param b Object to compare against
	* @return True if objects represent the same variant
	*/
	@Override
	/**
	* Tests equality with another object (must be a Var).
	* Uses optimized comparison with hash code pre-check for performance.
	* 
	* @param b Object to compare against
	* @return True if objects represent the same variant
	*/
	public boolean equals(Object b){
		return equals((Var)b);
	}

	/**
	 * Tests equality with another Var object.
	 * Two variants are equal if they have the same position, type, and allele.
	 * Uses hash code comparison for fast inequality detection.
	 * 
	 * @param b Var to compare against
	 * @return True if variants are equivalent
	 */
	public boolean equals(Var b){
		return hashcode==b.hashcode && compareTo(b)==0; // Hash pre-check + full comparison
	}

	/**
	* Returns hash code for this variant.
	* Used for efficient storage in hash-based collections.
	* 
	* @return Pre-computed hash code value
	*/
	@Override
	/**
	* Returns hash code for this variant.
	* Used for efficient storage in hash-based collections.
	* 
	* @return Pre-computed hash code value
	*/
	public int hashCode(){
		return hashcode;
	}

	/**
	* Generates unique key for this variant for use in bloom filters and hash maps.
	* Combines type, allele hash, length, and position into compact 64-bit key.
	* Uses bit shifting to pack multiple fields efficiently.
	* 
	* @return 64-bit key with high bit cleared (always positive)
	*/
	public long toKey(){
		final int len=(type==DEL ? reflen() : readlen()); // Use appropriate length for type
		final long key=type^((hash(allele)&0x3F)>>alleleShift)^(len>>lenShift)^(Long.rotateRight(start, startShift));
		return key&0x7FFFFFFFFFFFFFFFL; // Clear high bit to ensure positive
	}

	/**
	 * Compares this variant to another for sorting.
	 * Primary sort: scaffold number
	 * Secondary sort: adjusted start position (deletions sort slightly earlier)
	 * Tertiary sort: variant type  
	 * Quaternary sort: stop position
	 * Final sort: allele sequence
	 * 
	 * @param v Variant to compare against
	 * @return Negative, zero, or positive for less than, equal, or greater than
	 */
	@Override
	/**
	* Compares this variant to another for sorting.
	* Primary sort: scaffold number
	* Secondary sort: adjusted start position (deletions sort slightly earlier)
	* Tertiary sort: variant type  
	* Quaternary sort: stop position
	* Final sort: allele sequence
	* 
	* @param v Variant to compare against
	* @return Negative, zero, or positive for less than, equal, or greater than
	*/
	public int compareTo(Var v){
		if(scafnum!=v.scafnum){return scafnum-v.scafnum;} // Sort by scaffold first
		
		final int typeA=type(), typeB=v.type();
		int stA=start+(typeA==DEL ? -1 : 0); // Adjust deletion positions slightly earlier
		int stB=v.start+(typeB==DEL ? -1 : 0);
		if(stA!=stB){return stA-stB;} // Sort by adjusted start position
		if(typeA!=typeB){return typeA-typeB;} // Sort by variant type
		if(stop!=v.stop){return stop-v.stop;} // Sort by stop position
		return compare(allele, v.allele); // Sort by allele sequence
	}

	/**
	 * Compares two byte arrays lexicographically.
	 * Used for allele comparison in variant sorting.
	 * 
	 * @param a First byte array
	 * @param b Second byte array  
	 * @return Comparison result (negative/zero/positive)
	 */
	public int compare(byte[] a, byte[] b){
		if(a==b){return 0;} // Same reference
		if(a.length!=b.length){return b.length-a.length;} // Sort by length (longer first)
		for(int i=0; i<a.length; i++){
			byte ca=a[i], cb=b[i];
			if(ca!=cb){return ca-cb;} // Lexicographic comparison
		}
		return 0;
	}

	/**
	 * Returns string representation of this variant.
	 * Uses quick formatting with default parameters for debugging.
	 * 
	 * @return Formatted variant string
	 */
	@Override
	/**
	* Returns string representation of this variant.
	* Uses quick formatting with default parameters for debugging.
	* 
	* @return Formatted variant string
	*/
	public String toString(){
		return toTextQuick(new ByteBuilder()).toString();
	}

	/**
	 * Generates formatted text representation with default parameters.
	 * Quick version for debugging and logging purposes.
	 * 
	 * @param bb ByteBuilder to append formatted text to
	 * @return ByteBuilder with formatted variant data
	 */
	public ByteBuilder toTextQuick(ByteBuilder bb){
		return toText(bb, 0.99f, 30, 30, 150, 1, 2, null, null); // Default parameters for quick output
	}

	/**
	* Generates comprehensive formatted text representation in VAR format.
	* Includes all statistical data, quality metrics, and calculated scores.
	* 
	* @param bb ByteBuilder to append formatted output
	* @param properPairRate Overall proper pair rate for dataset
	* @param totalQualityAvg Average base quality across dataset
	* @param totalMapqAvg Average mapping quality across dataset  
	* @param readLengthAvg Average read length across dataset
	* @param rarity Minimum variant frequency threshold
	* @param ploidy Expected ploidy level
	* @param map Scaffold mapping for reference information
	* @return ByteBuilder with complete formatted variant data
	*/
	public ByteBuilder toText(ByteBuilder bb, double properPairRate, double totalQualityAvg, double totalMapqAvg, 
			double readLengthAvg, double rarity, int ploidy, ScafMap map, CellNet net){
		useIdentity=true; // Enable identity scoring for output
		
		// Core variant identification
		bb.append(scafnum).tab();
		bb.append(start).tab();
		bb.append(stop).tab();
		bb.append(typeArray[type()]).tab();
		for(byte b : allele){bb.append(b);} // Output allele sequence
		bb.tab();
		
		// Read count statistics by strand and pair
		bb.append(r1plus).tab();
		bb.append(r1minus).tab();
		bb.append(r2plus).tab();
		bb.append(r2minus).tab();
		bb.append(properPairCount).tab();
		bb.append(lengthSum).tab();

		// Quality statistics (sums and maximums)
		bb.append(mapQSum).tab();
		bb.append(mapQMax).tab();
		bb.append(baseQSum).tab();
		bb.append(baseQMax).tab();
		bb.append(endDistSum).tab();
		bb.append(endDistMax).tab();
		bb.append(idSum).tab();
		bb.append(idMax).tab();

		// Coverage and annotation data
		bb.append(coverage).tab();
		bb.append(minusCoverage).tab();
		bb.append(nearbyVarCount).tab();
		bb.append(flagged ? 1 : 0).tab();
		
		// Calculate distance from contig/scaffold ends
		final int scafEndDist=!doNscan ? nScan : (map==null ? start : contigEndDist(map));
		bb.append(scafEndDist).tab();

		// Calculate and output Phred-scaled quality score
		final double score=score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, net);
		bb.append(VarHelper.toPhredScore(score), 2).tab();
		
		// Extended statistics if enabled
		if(extendedText){
			bb.append(alleleCount()).tab();
			final double af=alleleFraction();
			bb.append(af, 4).tab();
			bb.append(revisedAlleleFraction(af, readLengthAvg), 4).tab();
			bb.append(strandRatio(), 4).tab();
			bb.append(baseQAvg(), 2).tab();
			bb.append(mapQAvg(), 2).tab();
			bb.append(edistAvg(), 2).tab();
			bb.append(identityAvg(), 2).tab();
			
			// Individual scoring components
			bb.append(edistScore(), 4).tab();
			bb.append(identityScore(), 4).tab();
			bb.append(qualityScore(totalQualityAvg, totalMapqAvg), 4).tab();
			bb.append(pairedScore(properPairRate, scafEndDist), 4).tab();
			bb.append(biasScore(properPairRate, scafEndDist), 4).tab();
			bb.append(coverageScore(ploidy, rarity, readLengthAvg), 4).tab();
			bb.append(homopolymerScore(map), 4).tab();
			bb.append(score, 4).tab();
		}
		
		bb.length--; // Remove final tab
		return bb;
	}
	
	/**
	* Generates VCF format output for this variant with comprehensive INFO fields.
	* Includes all statistical data, quality metrics, and sample-specific information.
	* Handles coordinate conversion, allele normalization, and proper VCF formatting.
	* 
	* @param bb ByteBuilder to append VCF line to
	* @param properPairRate Dataset proper pair rate for scoring
	* @param totalQualityAvg Dataset average base quality
	* @param mapqAvg Dataset average mapping quality
	* @param readLengthAvg Dataset average read length
	* @param ploidy Expected organism ploidy level
	* @param map Scaffold mapping for reference sequence access
	* @param filter Variant filter for pass/fail determination
	* @param trimWhitespace Whether to trim scaffold names
	* @return ByteBuilder with complete VCF line
	*/
	public ByteBuilder toVCF(ByteBuilder bb, double properPairRate, double totalQualityAvg, double mapqAvg, double readLengthAvg,
			int ploidy, ScafMap map, VarFilter filter, CellNet net, boolean trimWhitespace){
		
		final Scaffold scaf=map.getScaffold(scafnum);
		final byte[] bases=scaf.bases;
		final int reflen=reflen(), readlen=readlen(), type=type();
		//TODO: Scoring is done twice here
		final double score=phredScore(properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, filter.rarity, ploidy, map, net);
		final boolean pass=(filter==null ? true :
			filter.passesFilter(this, properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, ploidy, map, net, true));
		
		// CHROM field
		bb.append(trimWhitespace ? Tools.trimWhitespace(scaf.name) : scaf.name).tab();
		
		// POS field - convert to 1-based VCF coordinates
		boolean indel=(type==INS || type==DEL);
		boolean addPrevBase=true;
		final int vcfStart=start+(indel && addPrevBase ? 0 : 1);
		bb.append(vcfStart).tab();
		
		// ID field
		bb.append('.').tab();
		
		// REF and ALT fields with proper indel handling
		byte prevBase=(bases==null ? (byte)'N' : bases[Tools.mid(start-1, 0, bases.length-1)]);
		if(UPPER_CASE_ALLELES){prevBase=(byte) Tools.toUpperCase(prevBase);}
		
		if(addPrevBase){
			// REF field with leading base for indel normalization
			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			for(int i=0, rpos=start; i<reflen; i++, rpos++){
				bb.append(bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos]);
			}
			bb.tab();

			// ALT field with leading base
			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			bb.append(allele).tab();
		}else{
			// REF field without leading base
			if(reflen==0){
				bb.append('.');
			}else{
				for(int i=0, rpos=start; i<reflen; i++, rpos++){
					char refBase=bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos];
					if(UPPER_CASE_ALLELES){refBase=Tools.toUpperCase(refBase);}
					bb.append(refBase);
				}
			}
			bb.tab();

			// ALT field
			if(allele.length<1){
				bb.append('.').tab();
			}else{
				bb.append(allele).tab();
			}
		}
		
		// QUAL field
		bb.append(score, 2).tab();
		
		// FILTER field
		bb.append(pass ? "PASS\t" : "FAIL\t");
		
		// Calculate derived statistics
		final int scafEndDist=!doNscan ? nScan : (map==null ? start : contigEndDist(map));
		final int count=alleleCount();
		final double af=alleleFraction();
		final double raf=revisedAlleleFraction(af, readLengthAvg);
		final double strandBias=strandBiasScore(scafEndDist);

		// INFO field with comprehensive variant statistics
		{
			assert(Scaffold.trackStrand()==(minusCoverage>=0)) : Scaffold.trackStrand()+", "+minusCoverage;
			final int covMinus=(Scaffold.trackStrand() ? minusCoverage : coverage/2);
			final int covPlus=Tools.max(0, coverage-covMinus);
			final int refMinus=Tools.max(0, covMinus-alleleMinusCount());
			final int refPlus=Tools.max(0, covPlus-allelePlusCount());
			
			bb.append("SN=").append(scafnum).append(';');
			bb.append("STA=").append(start).append(';');
			bb.append("STO=").append(stop).append(';');
			bb.append("TYP=").append(typeArray[type()]).append(';');
			
			bb.append("R1P=").append(r1plus).append(';');
			bb.append("R1M=").append(r1minus).append(';');
			bb.append("R2P=").append(r2plus).append(';');
			bb.append("R2M=").append(r2minus).append(';');
			
			bb.append("AD=").append(count).append(';');
			bb.append("DP=").append(Tools.max(coverage, count)).append(';');
			bb.append("MCOV=").append(minusCoverage).append(';');
			bb.append("PPC=").append(properPairCount).append(';');
			
			bb.append("AF=").append(af,4).append(';');
			bb.append("RAF=").append(raf,4).append(';');
			bb.append("LS=").append(lengthSum).append(';');
			
			bb.append("MQS=").append(mapQSum).append(';');
			bb.append("MQM=").append(mapQMax).append(';');
			bb.append("BQS=").append(baseQSum).append(';');
			bb.append("BQM=").append(baseQMax).append(';');
			
			bb.append("EDS=").append(endDistSum).append(';');
			bb.append("EDM=").append(endDistMax).append(';');
			bb.append("IDS=").append(idSum).append(';');
			bb.append("IDM=").append(idMax).append(';');

			bb.append("NVC=").append(nearbyVarCount).append(';');
			bb.append("FLG=").append(flagged ? 1 : 0).append(';');
			bb.append("CED=").append(scafEndDist).append(';');
			bb.append("HMP=").append(homopolymerCount(map)).append(';');
			bb.append("SB=").append(strandBias,4).append(';');
			
			if(Scaffold.trackStrand()){
				bb.append("DP4=").append(refPlus).append(',').append(refMinus).append(',');
				bb.append(allelePlusCount()).append(',').append(alleleMinusCount()).append(';');
			}
			
			bb.length--; // Remove final semicolon
		}
		
		// FORMAT field
		{
			bb.tab();
			bb.append("GT:DP:AD:AF:RAF:NVC:FLG:SB:SC:PF");
			bb.tab();

			// Sample data
			bb.append(genotype(ploidy, pass));
			bb.append(':');
			bb.append(Tools.max(coverage, count));
			bb.append(':');
			bb.append(count);
			bb.append(':');
			bb.append(af,4);
			bb.append(':');

			bb.append(raf,4);
			bb.append(':');
			bb.append(nearbyVarCount);
			bb.append(':');
			bb.append(flagged ? 1 : 0);
			bb.append(':');
			bb.append(strandBias,4);
			bb.append(':');
			
			bb.append(score,2);
			bb.append(':');
			bb.append(pass ? "PASS" : "FAIL");
		}
		
		return bb;
	}

	/**
	 * Calculates number of variant copies based on allele frequency and ploidy.
	 * Uses allele frequency thresholds to determine heterozygous vs homozygous calls.
	 * Ensures minimum variant copy requirements are met.
	 * 
	 * @param ploidy Expected organism ploidy level
	 * @return Number of variant copies (0 to ploidy)
	 */
	public int calcCopies(int ploidy){
		final double af=alleleFraction();
		if(ploidy==1){
			return af<0.4 ? 0 : 1; // Haploid: threshold at 40%
		}else if(ploidy==2){
			if(af<0.2){return 0;} // Homozygous reference
			if(af<0.8){return 1;} // Heterozygous
			return 2; // Homozygous variant
		}
		
		// General ploidy handling
		int copies=(int)Math.round(ploidy*af); // Round to nearest integer
		if(af>=0.5){copies=Tools.max(copies, 1);} // At least 1 copy if AF >= 50%
		copies=Tools.mid(MIN_VAR_COPIES, copies, ploidy); // Clamp to valid range
		return copies;
	}/**
	 * Generates VCF genotype string based on variant copies and filter status.
	 * Uses allele frequency and ploidy to determine genotype call.
	 * Can output dot genotype for failed variants if configured.
	 * 
	 * @param ploidy Expected organism ploidy level
	 * @param pass Whether variant passed quality filters
	 * @return VCF genotype string (e.g., "0/1", "1/1", "./.")
	 */
	private String genotype(int ploidy, boolean pass){
		// Handle failed variants with dot genotype if enabled
		if(!pass && noPassDotGenotype){
			if(ploidy==1){return ".";} // Haploid no-call
			else if(ploidy==2){return "./.";} // Diploid no-call
			StringBuilder sb=new StringBuilder(ploidy*2-1);
			sb.append('.');
			for(int i=1; i<ploidy; i++){
				sb.append('/').append('.'); // Multi-ploid no-call
			}
			return sb.toString();
		}
		
		int copies=calcCopies(ploidy); // Calculate variant copies
		
		// Handle common ploidy cases
		if(ploidy==1){return copies==0 ? "0" : "1";} // Haploid: 0 or 1
		if(ploidy==2){
			if(copies==0){return "0/0";} // Homozygous reference
			if(copies==1){return "0/1";} // Heterozygous
			return "1/1"; // Homozygous variant
		}
		
		// General ploidy handling
		StringBuilder sb=new StringBuilder(ploidy*2);
		int refCopies=ploidy-copies; // Reference allele copies
		
		// Add reference alleles
		for(int i=0; i<refCopies; i++){
			sb.append(0).append('/');
		}
		// Add variant alleles
		for(int i=0; i<copies; i++){
			sb.append(1).append('/');
		}
		sb.setLength(sb.length()-1); // Remove final slash
		return sb.toString();
	}

	/**
	 * Computes hash code for this variant based on position and allele.
	 * Uses position rotation and allele hash combination for good distribution.
	 * 
	 * @return Hash code value for this variant
	 */
	private int hash(){
		return scafnum^Integer.rotateLeft(start, 9)^Integer.rotateRight(stop, 9)^hash(allele);
	}

	/**
	 * Computes hash code for byte array using position-dependent mixing.
	 * Provides good hash distribution for allele sequences of varying lengths.
	 * 
	 * @param a Byte array to hash (typically allele sequence)
	 * @return Hash code for the byte array
	 */
	public static final int hash(byte[] a){
		int code=123456789; // Prime seed
		for(byte b : a){
			code=Integer.rotateLeft(code, 3)^codes[b]; // Mix with pre-computed random codes
		}
		return code&Integer.MAX_VALUE; // Ensure positive result
	}

	/**
	 * Calculates or retrieves coverage at this variant position.
	 * Uses cached value if available, otherwise queries scaffold coverage data.
	 * Also calculates strand-specific coverage if strand tracking is enabled.
	 * 
	 * @param map Scaffold mapping containing coverage information
	 * @return Total coverage depth at this position
	 */
	public int calcCoverage(ScafMap map){
		if(coverage>=0){return coverage;} // Return cached value if available
		
		Scaffold scaf=map.getScaffold(scafnum);
		coverage=scaf.calcCoverage(this); // Calculate total coverage
		if(Scaffold.trackStrand()){
			minusCoverage=scaf.minusCoverage(this); // Calculate minus strand coverage
		}
		return coverage;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Scoring Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Scoring Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Calculates Phred-scaled quality score for this variant.
	 * Converts internal scoring (0.0-1.0) to standard Phred scale for output.
	 * 
	 * @param properPairRate Overall proper pair rate for dataset
	 * @param totalQualityAvg Average base quality across all reads
	 * @param totalMapqAvg Average mapping quality across all reads
	 * @param readLengthAvg Average read length across dataset
	 * @param rarity Minimum variant frequency threshold
	 * @param ploidy Expected organism ploidy level
	 * @param map Scaffold mapping for reference sequence access
	 * @param net Optional CellNet for prediction
	 * @return Phred-scaled quality score (higher = more confident)
	 */
	public double phredScore(double properPairRate, double totalQualityAvg, double totalMapqAvg, 
			double readLengthAvg, double rarity, int ploidy, ScafMap map, CellNet net){
		double score=score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, net);
		return VarHelper.toPhredScore(score);
	}

	/**
	 * Calculates composite variant quality score from multiple evidence types.
	 * Combines coverage, quality, bias, positional, and sequence context scores.
	 * Uses geometric mean (power 0.2) to balance different evidence types.
	 * 
	 * @param properPairRate Dataset-wide proper pair rate for normalization
	 * @param totalQualityAvg Dataset-wide average base quality
	 * @param totalMapqAvg Dataset-wide average mapping quality
	 * @param readLengthAvg Dataset-wide average read length
	 * @param rarity Minimum expected variant frequency
	 * @param ploidy Expected number of chromosome copies
	 * @param map Scaffold map for sequence context analysis
	 * @param net Optional CellNet for prediction
	 * @return Composite quality score (0.0 to 1.0, higher = better)
	 */
	public double score(double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map, CellNet net){
		int scafEndDist=(map==null ? start : contigEndDist(map)); // Distance from contig ends
		
		// Calculate individual scoring components
		double cs=coverageScore(ploidy, rarity, readLengthAvg); // Coverage adequacy
		if(cs==0){return 0;} // No coverage = no call
		double es=(useEdist ? edistScore() : 1); // Distance from read ends
		double qs=qualityScore(totalQualityAvg, totalMapqAvg); // Base and mapping quality
		double ps=(usePairing ? pairedScore(properPairRate, scafEndDist) : 1); // Proper pairing rate
		double bs=(useBias ? biasScore(properPairRate, scafEndDist) : 1); // Strand/read bias
		double is=(useIdentity ? identityScore() : 1); // Alignment identity
		double hs=(useHomopolymer ? homopolymerScore(map) : 1); // Homopolymer context
		
		// Geometric mean of all components (power 0.2 = 5th root)
		double gMean=Math.pow(es*qs*ps*bs*cs*is*hs, 0.2);
		double score=gMean;
		if(net!=null) {
			float[] vec=FeatureVectorMaker.toVector(this, properPairRate, totalQualityAvg, 
					totalMapqAvg, readLengthAvg, ploidy, map);
			float output=net.applyInput(vec).feedForward();
			
		}
		return score;
	}

	/**
	 * Scores variant based on distance from read ends.
	 * Variants near read ends are less reliable due to sequencing quality decline.
	 * Uses read length and position to calculate confidence penalty.
	 * 
	 * @return Score from 0.05 to 1.0 (higher = farther from ends)
	 */
	public double edistScore(){
		double lengthAvg=lengthAvg(); // Average supporting read length
		double edistAvg=((edistAvg()*2+endDistMax))*0.333333333333; // Weighted distance average
		double constant=5+Tools.min(20, lengthAvg*0.1)+lengthAvg*0.01; // Length-dependent threshold
		double weighted=Tools.max(0.05, edistAvg-Tools.min(constant, edistAvg*0.95)); // Apply threshold
		weighted=weighted*weighted; // Square for non-linear penalty
		return weighted/(weighted+4); // Normalize to 0-1 range
	}

	/**
	 * Scores variant based on alignment identity of supporting reads.
	 * Higher identity reads provide more reliable variant evidence.
	 * Adjusts for the variant's own contribution to identity calculation.
	 * 
	 * @return Score from 0.75 to 1.0 (higher = better alignment identity)
	 */
	public double identityScore(){
		double lengthAvg=lengthAvg(); // Average read length
		double idAvg=0.001f*(((identityAvg()+idMax))*0.5f); // Average identity (0-1 scale)
		// Diminish impact of this variant on overall identity
		double weighted=Tools.min(1, (idAvg*lengthAvg+(0.65f*Tools.max(1, readlen())))/lengthAvg);
		weighted=0.75f+0.25f*weighted; // Compress range to 0.75-1.0
		return weighted;
	}

	/**
	 * Combines base quality and mapping quality scores.
	 * Both quality types contribute to variant reliability assessment.
	 * 
	 * @param totalBaseqAvg Dataset average base quality for normalization
	 * @param totalMapqAvg Dataset average mapping quality for normalization
	 * @return Combined quality score (product of individual scores)
	 */
	public double qualityScore(double totalBaseqAvg, double totalMapqAvg){
		return baseQualityScore(totalBaseqAvg)*mapQualityScore(totalMapqAvg);
	}

	/**
	 * Scores variant based on base quality of supporting reads.
	 * Compares variant's base quality against dataset average.
	 * Includes correction for recalibrated quality scores.
	 * 
	 * @param totalBaseqAvg Dataset-wide average base quality
	 * @return Base quality score (0.0 to 1.0)
	 */
	public double baseQualityScore(double totalBaseqAvg){
		double bqAvg=baseQAvg(); // This variant's average base quality
		
		// Fudge factor for recalibrated quality scores (not well-tested)
		if(totalBaseqAvg<32 && bqAvg<32){
			double fudgeFactor1=0.75*(32-totalBaseqAvg); // Dataset adjustment
			double fudgeFactor2=0.75*(32-bqAvg); // Variant adjustment
			totalBaseqAvg+=fudgeFactor1;
			bqAvg+=Tools.min(fudgeFactor1, fudgeFactor2);
		}
		
		// Apply penalty if variant quality is below dataset average
		final double delta=totalBaseqAvg-bqAvg; // Quality deficit
		if(delta>0){
			bqAvg=Tools.max(bqAvg*0.5, bqAvg-0.5*delta); // Reduce effective quality
		}
		
		// Transform quality with threshold and multiplier
		double mult=0.25;
		double thresh=12;
		if(bqAvg>thresh){
			bqAvg=bqAvg-thresh+(thresh*mult); // Linear above threshold
		}else{
			bqAvg=bqAvg*mult; // Scaled below threshold
		}
		
		// Convert to probability and square for emphasis
		double baseProbAvg=1-Math.pow(10, 0-.1*bqAvg); // Phred to probability
		double d=baseProbAvg*baseProbAvg; // Square the probability
		return d;
	}

	/**
	* Scores variant based on distance from read ends.
	* Variants near read ends are less reliable due to sequencing quality decline.
	* Uses read length and position to calculate confidence penalty.
	* 
	* @return Score from 0.05 to 1.0 (higher = farther from ends)
	*/
	public double mapQualityScore(double totalMapqAvg){
		double mqAvg=0.5f*(mapQAvg()+mapQMax); // Average of mean and max
		double mapProbAvg=1-Math.pow(10, 0-.1*(mqAvg+2)); // Phred to probability (+2 bonus)
		double d=mapProbAvg; // Direct probability score
		return d;
	}

	/**
	 * Scores variant based on proper pairing rate of supporting reads.
	 * Compares variant's pairing rate against dataset baseline.
	 * Applies positional correction for variants near contig ends.
	 * 
	 * @param properPairRate Dataset-wide proper pair rate
	 * @param scafEndDist Distance from nearest contig end
	 * @return Pairing score (0.1 to 1.0)
	 */
	public double pairedScore(double properPairRate, int scafEndDist){
		if(properPairRate<0.5){return 0.98;} // Skip if dataset has poor pairing
		final double count=alleleCount();
		if(count==0){return 0;} // No reads = no score
		double rate=properPairCount/count; // Variant's pairing rate
		rate=rate*(count/(0.1+count)); // Weight by read count
		if(rate*1.05>=properPairRate){ // Good pairing rate
			return Tools.max(rate, 1-0.001*properPairRate);
		}
		double score=((rate*1.05)/properPairRate)*0.5+0.5; // Scale poor pairing
		score=Tools.max(0.1, score); // Minimum score
		return modifyByEndDist(score, scafEndDist); // Adjust for position
	}

	/**
	 * Adjusts scores based on distance from contig ends.
	 * Variants near contig ends may have artifacts due to assembly issues.
	 * Provides score bonus for variants in problematic regions.
	 * 
	 * @param x Original score to modify
	 * @param scafEndDist Distance from nearest contig end
	 * @return Modified score (may be increased near contig ends)
	 */
	public double modifyByEndDist(double x, int scafEndDist){
		if(x>=0.99 || !doNscan || scafEndDist>=nScan){return x;} // No adjustment needed
		if(scafEndDist<minEndDistForBias){ // Very close to end
			return Tools.max(x, 0.98+0.02*x); // Boost score significantly
		}
		double delta=1-x; // Score deficit
		delta=delta*(scafEndDist*scafEndDist)/(nScan*nScan); // Scale by distance squared
		return 1-delta; // Apply adjusted penalty
	}

	/**
	 * Scores variant based on coverage depth and allele fraction.
	 * Considers expected ploidy and minimum variant frequency thresholds.
	 * Includes adjustment for insertion length effects on coverage.
	 * 
	 * @param ploidy Expected organism ploidy level
	 * @param rarity Minimum expected variant frequency
	 * @param readLengthAvg Average read length for insertion adjustments
	 * @return Coverage adequacy score (0.0 to 1.0)
	 */
	public double coverageScore(int ploidy, double rarity, double readLengthAvg){
		int count=alleleCount();
		if(count==0){return 0;} // No supporting reads
		double rawScore=count/(lowCoveragePenalty+count); // Coverage adequacy (may be severe)
		
		double ratio=0.98; // Default allele fraction
		if(coverage>0){
			double dif=coverage-count; // Reference read count
			if(dif>0){
				// Adjust for expected sequencing errors and biases
				dif=dif-coverage*.01f-Tools.min(0.5f, coverage*.1f);
				dif=Tools.max(0.1f, dif);
			}
			ratio=(coverage-dif)/coverage; // Adjusted allele fraction
			
			// Use revised allele fraction for substitutions if available
			if(type()==SUB && revisedAlleleFraction!=-1 && revisedAlleleFraction<ratio){
				ratio=revisedAlleleFraction;
			}else{
				ratio=adjustForInsertionLength(ratio, readLengthAvg); // Adjust for insertion artifacts
			}
			
			// Handle rare variants with ploidy considerations
			if(rarity<1 && ratio>rarity){
				double minExpected=1f/ploidy; // Minimum expected for heterozygote
				if(ratio<minExpected){
					ratio=minExpected-((minExpected-ratio)*0.1); // Modest boost for low-frequency variants
				}
			}
		}
		
		double ratio2=Tools.min(1, ploidy*ratio); // Scale by ploidy
		return rawScore*ratio2; // Combine coverage and fraction scores
	}
	
	/**
	 * Revises allele fraction for insertions by adjusting nearby substitutions.
	 * Insertions can create false substitutions when reads span the insertion boundary.
	 * This method reduces allele fractions of nearby substitutions that may be artifacts.
	 * Only applies to insertion variants with sufficient length and valid positions.
	 * 
	 * @param readLengthAvg Average read length for calculating adjustment factors
	 * @param scaffold Reference scaffold containing this insertion
	 * @param map VarMap containing other variants that may need adjustment
	 */
	public void reviseAlleleFraction(double readLengthAvg, Scaffold scaffold, VarMap map){
		assert(type()==INS);
		final int ilen=readlen();
		if(ilen<3 || start<1 || start>=scaffold.length-2){return;} // Skip short insertions or edge cases
		final byte[] bases=scaffold.bases;
		
		// Calculate adjustment based on insertion's revised allele fraction
		final double afIns=alleleFraction(); // Original insertion AF
		final double rafIns=revisedAlleleFraction(afIns, readLengthAvg); // Corrected insertion AF
		final double revisedDif=0.55*(rafIns-afIns); // Half left, half right on average
		final double mult=revisedDif/allele.length; // Per-base adjustment factor
		
		// Adjust substitutions to the right of insertion
		for(int i=0, j=start; i<allele.length && j<scaffold.bases.length; i++, j++){
			final byte b=allele[i];
			if(b!=bases[j]){ // Insertion differs from reference
				Var key=new Var(scaffold.number, j, j+1, b, SUB);
				Var affectedSub=map.get(key);
				if(affectedSub!=null){
					assert(key.type()==SUB);
					final double subModifier=revisedDif-mult*i; // Distance-dependent adjustment
					synchronized(affectedSub){
						double afSub=affectedSub.alleleFraction();
						double rafSub=affectedSub.revisedAlleleFraction;
						double modified=afSub-subModifier; // Reduce substitution AF
						if(rafSub==-1){
							affectedSub.revisedAlleleFraction=Tools.max(afSub*0.05, modified);
						}else{
							affectedSub.revisedAlleleFraction=Tools.min(rafSub, Tools.max(afSub*0.05, modified));
						}
					}
				}
			}
		}
		
		// Adjust substitutions to the left of insertion (reverse order)
		for(int i=0, j=start-1; i<allele.length && j>=0; i++, j--){
			final byte b=allele[allele.length-1-i]; // Process insertion sequence backwards
			if(b!=bases[j]){
				Var key=new Var(scaffold.number, j, j+1, b, SUB);
				Var affectedSub=map.get(key);
				if(affectedSub!=null){
					assert(key.type()==SUB);
					final double subModifier=revisedDif-mult*i;
					synchronized(affectedSub){
						double afSub=affectedSub.alleleFraction();
						double rafSub=affectedSub.revisedAlleleFraction;
						double modified=afSub-subModifier;
						if(rafSub==-1){
							affectedSub.revisedAlleleFraction=Tools.max(afSub*0.05, modified);
						}else{
							affectedSub.revisedAlleleFraction=Tools.min(rafSub, Tools.max(afSub*0.05, modified));
						}
					}
				}
			}
		}
	}

	/**
	 * Gets revised allele fraction, calculating it if not already computed.
	 * For insertions, adjusts for read length bias where long insertions
	 * extending beyond read boundaries appear at artificially low frequencies.
	 * 
	 * @param af Original allele fraction
	 * @param readLengthAvg Average read length for insertion length adjustment
	 * @return Revised allele fraction accounting for technical biases
	 */
	public double revisedAlleleFraction(double af, double readLengthAvg){
		if(revisedAlleleFraction!=-1){
			return revisedAlleleFraction; // Return cached value
		}else if(type()==INS){
			return revisedAlleleFraction=adjustForInsertionLength(af, readLengthAvg);
		}
		return af; // No adjustment needed for non-insertions
	}

	/**
	 * Adjusts allele fraction for insertion length bias.
	 * Long insertions near read ends won't be fully captured, causing
	 * systematic underestimation of insertion allele frequency.
	 * 
	 * @param ratio Original allele fraction
	 * @param rlen0 Average read length from dataset
	 * @return Adjusted allele fraction accounting for insertion length bias
	 */
	public double adjustForInsertionLength(final double ratio, final double rlen0){
		if(type()!=INS){return ratio;} // Only applies to insertions
		final int ilen=readlen();
		if(ilen<2){return ratio;} // Skip very short insertions
		
		final double rlen=Tools.max(ilen*1.2+6, rlen0); // Effective read length
		final double sites=rlen+ilen-1; // Total possible observation sites
		final double goodSites=rlen-ilen*1.1-6; // Sites where insertion fully observable
		
		final double expectedFraction=goodSites/sites; // Expected observable fraction
		final double revisedRatio=Tools.min(ratio/expectedFraction, 1-(1-ratio)*0.1); // Upward adjustment
		return revisedRatio;
	}

	/**
	 * Calculates homopolymer penalty score for this variant.
	 * Homopolymer regions are prone to sequencing errors, especially indels.
	 * Longer homopolymer runs receive progressively lower confidence scores.
	 * 
	 * @param map Scaffold mapping for reference sequence access
	 * @return Homopolymer score (1.0 = no homopolymer, lower = longer homopolymer)
	 */
	public double homopolymerScore(ScafMap map){
		if(map==null){return 1;}
		
		int count=homopolymerCount(map); // Get homopolymer length
		if(count<2){return 1;} // No penalty for short runs
		return 1f-(count*0.1f/9); // Linear penalty up to count=9
	}

	/**
	 * Counts homopolymer run length around this variant position.
	 * Different calculation methods for substitutions, insertions, and deletions.
	 * Used to assess likelihood of sequencing errors in repetitive sequence.
	 * 
	 * @param map Scaffold mapping for reference sequence access
	 * @return Length of homopolymer run (0 = no homopolymer)
	 */
	public int homopolymerCount(ScafMap map){
		if(map==null){return 0;}
		final byte[] bases=map.getScaffold(scafnum).bases;
		if(bases==null){return 0;}
		
		final int type=type();
		if(type==SUB){
			assert(start==stop-1) : start+", "+stop;
			final byte base=allele[0]; // Substituted base
			int x=VarHelper.homopolymerCountSub(bases, start, base);
			return x;
		}else if(type==INS){
			final byte base1=allele[0], base2=allele[allele.length-1]; // First and last inserted bases
			int i=0;
			// Check if entire insertion is homopolymer
			while(i<allele.length && allele[i]==base1){i++;}
			while(i<allele.length && allele[i]==base2){i++;}
			if(i<bases.length){return 0;} // Mixed sequence insertion
			// Count flanking homopolymer
			int left=VarHelper.homopolymerCountLeft(bases, start, base1);
			int right=VarHelper.homopolymerCountRight(bases, stop+1, base2);
			return left+right+1;
		}else if(type==DEL){
			if(start<0 || start+1>=bases.length || stop<=0 || stop>=bases.length){return 0;}
			final byte base1=bases[start+1], base2=bases[stop-1]; // First and last deleted bases
			int pos=start+1;
			// Check if entire deletion is homopolymer
			while(pos<=stop && bases[pos]==base1){pos++;}
			while(pos<=stop && bases[pos]==base2){pos++;}
			if(pos<=stop){return 0;} // Mixed sequence deletion
			// Count flanking homopolymer
			int left=VarHelper.homopolymerCountLeft(bases, start, base1);
			int right=VarHelper.homopolymerCountRight(bases, stop, base2);
			return left+right+1;
		}else{
			return 0; // No homopolymer analysis for other types
		}
	}

	/**
	 * Calculates combined bias score from strand and read biases.
	 * Uses geometric mean (square root of product) to balance both bias types.
	 * 
	 * @param properPairRate Dataset proper pair rate for read bias calculation
	 * @param scafEndDist Distance from contig ends for positional adjustment
	 * @return Combined bias score (0.0 to 1.0, higher = less biased)
	 */
	public double biasScore(double properPairRate, int scafEndDist){
		double strandBias=strandBiasScore(scafEndDist); // Plus vs minus strand bias
		double readBias=readBiasScore(properPairRate); // Read 1 vs read 2 bias
		return Math.sqrt(strandBias*readBias); // Geometric mean
	}

	/**
	 * Calculates strand bias score using statistical significance testing.
	 * Tests whether plus/minus strand distribution is significantly skewed.
	 * Includes special handling for high-coverage variants with mild bias.
	 * 
	 * @param scafEndDist Distance from contig ends for positional correction
	 * @return Strand bias score (0.0 to 1.0, higher = less biased)
	 */
	public double strandBiasScore(int scafEndDist){
		int plus=allelePlusCount(); // Plus strand supporting reads
		int minus=alleleMinusCount(); // Minus strand supporting reads
		final double x=VarProb.eventProb(plus, minus); // Statistical significance of bias
		final double x2=modifyByEndDist(x, scafEndDist); // Adjust for position
		
		double result=x2;
		// Relaxed stringency for high-coverage variants seen on both strands
		if(plus+minus>=20 && x2<0.9){
			int min=Tools.min(plus, minus);
			int max=Tools.max(plus, minus);
			if(min>1 && min>0.06f*max){ // Present on both strands
				double y=0.15+(0.2*min)/max; // Relaxation factor
				result=y+(1-y)*x2; // Blend with original score
			}
		}
		return result;
	}

	/**
	 * Calculates read bias score (Read 1 vs Read 2 distribution).
	 * Tests whether variant appears preferentially in one read of pair.
	 * Includes relaxed scoring for high-coverage variants.
	 * 
	 * @param properPairRate Dataset proper pair rate (affects scoring)
	 * @return Read bias score (0.0 to 1.0, higher = less biased)
	 */
	public double readBiasScore(double properPairRate){
		if(properPairRate<0.5){return 0.95f;} // Skip for unpaired data
		final int r1=r1AlleleCount(), r2=r2AlleleCount(); // Read 1 vs Read 2 counts
		final double x=VarProb.eventProb(r1, r2); // Statistical significance
		
		final double x2=0.10+0.90*x; // Compress range
		double result=x2;
		// Relaxed stringency for high-coverage variants
		if(r1+r2>=20 && x2<0.9){
			int min=Tools.min(r1, r2);
			int max=Tools.max(r1, r2);
			if(min>1 && min>0.07f*max){ // Present in both reads
				double y=0.15+(0.2*min)/max;
				result=y+(1-y)*x2;
			}
		}
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns count of supporting reads on plus strand */
	public int allelePlusCount(){return r1plus+r2plus;}
	/** Returns count of supporting reads on minus strand */
	public int alleleMinusCount(){return r1minus+r2minus;}
	/** Returns count of supporting reads from Read 1 */
	public int r1AlleleCount(){return r1plus+r1minus;}
	/** Returns count of supporting reads from Read 2 */
	public int r2AlleleCount(){return r2plus+r2minus;}
	/** Returns total count of reads supporting this variant */
	public int alleleCount(){return r1plus+r1minus+r2plus+r2minus;}

	/**
	* Calculates allele fraction (variant reads / total coverage).
	* Uses maximum of variant count and total coverage for robustness.
	* 
	* @return Allele fraction (0.0 to 1.0)
	*/
	public double alleleFraction(){
		int count=alleleCount();
		int cov=Tools.max(count, coverage, 1); // Avoid division by zero
		return count/(double)cov;
	}

	/**
	 * Calculates strand ratio balance (minority strand / majority strand).
	 * Perfect balance returns 1.0, complete bias approaches 0.0.
	 * 
	 * @return Strand balance ratio (0.0 to 1.0)
	 */
	public double strandRatio(){
		int plus=allelePlusCount();
		int minus=alleleMinusCount();
		if(plus==minus){return 1;} // Perfect balance
		return (Tools.min(plus,  minus)+1)/(double)Tools.max(plus, minus);
	}

	/** Returns average base quality of supporting reads */
	public double baseQAvg(){return baseQSum/(double)alleleCount();}
	/** Returns average mapping quality of supporting reads */
	public double mapQAvg(){return mapQSum/(double)alleleCount();}
	/** Returns average distance from read ends */
	public double edistAvg(){return endDistSum/(double)alleleCount();}
	/** Returns average alignment identity of supporting reads */
	public double identityAvg(){return idSum/(double)alleleCount();}
	/** Returns average length of supporting reads */
	public double lengthAvg(){return lengthSum/(double)alleleCount();}
	/** Returns fraction of supporting reads that are properly paired */
	public double properPairRate(){return properPairCount/(double)alleleCount();}

	/**
	 * Sets coverage values for this variant position.
	 * @param coverage_ Total depth of coverage
	 * @param minusCoverage_ Coverage on minus strand
	 */
	public void setCoverage(int coverage_, int minusCoverage_){
		coverage=coverage_;
		minusCoverage=minusCoverage_;
	}

	/**
	 * Returns total coverage at this position.
	 * @return Coverage depth (must be calculated first)
	 */
	public int coverage(){
		assert(coverage>-1) : coverage+", "+this;
		return coverage;
	}

	/**
	 * Checks if coverage has been calculated for this variant.
	 * @return True if coverage data is available
	 */
	public boolean hasCoverage(){
		return coverage>-1;
	}

	/**
	 * Calculates distance from nearest contig end.
	 * Variants near contig ends may be less reliable due to assembly issues.
	 * Scans for runs of N bases to identify true contig boundaries.
	 * 
	 * @param map Scaffold mapping for sequence access
	 * @return Distance to nearest contig end (in bases)
	 */
	public int contigEndDist(ScafMap map){
		Scaffold scaf=map.getScaffold(scafnum);
		int len=scaf.length;
		byte[] bases=scaf.bases;
		
		int scafEndDist=Tools.max(0, Tools.min(start, len-stop)); // Distance to scaffold end
		if(bases==null || nScan<1){return scafEndDist;}
		int limit=Tools.min(nScan, scafEndDist);
		int contigEndDist=leftContigEndDist(bases, limit); // Check left side for N runs
		limit=Tools.min(limit, contigEndDist);
		contigEndDist=rightContigEndDist(bases, limit); // Check right side for N runs
		return Tools.min(scafEndDist, contigEndDist);
	}

	/**
	 * Calculates distance to nearest contig end on the left side.
	 * Searches for runs of 10+ N bases indicating contig boundaries.
	 * 
	 * @param bases Reference sequence
	 * @param maxDist Maximum distance to search
	 * @return Distance to left contig end
	 */
	public int leftContigEndDist(byte[] bases, int maxDist){
		if(start>=bases.length){return Tools.min(bases.length, maxDist+1);}
		int ns=0; // Count of consecutive Ns
		for(int i=start, lim=Tools.max(0, start-maxDist); i>=lim; i--){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0; // Reset N count
			}else{
				ns++;
				if(ns>=10){ // Found contig boundary
					int x=start-i-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1; // No boundary found
	}

	/**
	 * Calculates distance to nearest contig end on the right side.
	 * Searches for runs of 10+ N bases indicating contig boundaries.
	 * 
	 * @param bases Reference sequence
	 * @param maxDist Maximum distance to search
	 * @return Distance to right contig end
	 */
	public int rightContigEndDist(byte[] bases, int maxDist){
		if(stop<0){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=stop, lim=Tools.min(bases.length-1, stop+maxDist); i<=lim; i++){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=i-stop-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}

	/**
	 * Gets scaffold name using default scaffold mapping.
	 * @return Scaffold/chromosome name
	 */
	public String scafName(){
		return scafName(ScafMap.defaultScafMap());
	}

	/**
	 * Gets scaffold name using specified mapping.
	 * @param map Scaffold mapping to use
	 * @return Scaffold/chromosome name
	 */
	public String scafName(ScafMap map){
		return map.getScaffold(scafnum).name;
	}

	/** Sets forced variant status (from input VCF) */
	public Var setForced(boolean b){forced=b; return this;}
	/** Returns whether this variant was forced from input */
	public boolean forced(){return forced;}

	/** Sets flagged status for this variant */
	public Var setFlagged(boolean b){flagged=b; return this;}
	/** Returns whether this variant has been flagged */
	public boolean flagged(){return flagged;}

	/** Returns true if this is an insertion */
	public final boolean ins(){return type==INS;}
	/** Returns true if this is a deletion */
	public final boolean del(){return type==DEL;}
	/** Returns true if this is an insertion or deletion */
	public final boolean indel(){return type==INS || type==DEL;}
	/** Returns true if this is a substitution */
	public final boolean sub(){return type==SUB;}
	/** Returns true if this is a complex variant */
	public final boolean complex(){return type==COMPLEX;}
	/** 
	 * Returns true if this variant causes a frameshift.
	 * Frameshifts occur when indel length is not divisible by 3.
	 */
	public final boolean frameshift(){
		int delta=Tools.absdif(reflen(), readlen()); // Length difference
		return delta%3!=0; // Not divisible by 3
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Scaffold/chromosome number for this variant */
	public final int scafnum;
	/** Start position of variant (0-based, inclusive) */
	public final int start;
	/** Stop position of variant (0-based, exclusive) */
	public final int stop;
	/** Alternative allele sequence */
	public final byte[] allele;
	/** Pre-computed hash code for efficient storage and comparison */
	public final int hashcode;
	/** Variant type constant (SUB, INS, DEL, etc.) */
	public final int type;

	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Total coverage depth at this position (-1 = not calculated) */
	int coverage=-1;
	/** Coverage on minus strand (-1 = not calculated) */
	int minusCoverage=-1;

	/** Count of Read 1 supporting reads on plus strand */
	int r1plus;
	/** Count of Read 1 supporting reads on minus strand */
	int r1minus;
	/** Count of Read 2 supporting reads on plus strand */
	int r2plus;
	/** Count of Read 2 supporting reads on minus strand */
	int r2minus;
	/** Count of properly paired supporting reads */
	int properPairCount;

	/** Sum of mapping qualities from supporting reads */
	long mapQSum;
	/** Maximum mapping quality among supporting reads */
	public int mapQMax;

	/** Sum of base qualities from supporting reads */
	long baseQSum;
	/** Maximum base quality among supporting reads */
	public int baseQMax;

	/** Sum of distances from read ends for supporting reads */
	long endDistSum;
	/** Maximum distance from read ends among supporting reads */
	public int endDistMax;

	/** Sum of alignment identity scores from supporting reads */
	long idSum;
	/** Maximum alignment identity among supporting reads */
	int idMax;

	/** Sum of read lengths from supporting reads */
	long lengthSum;

	/** Count of nearby variants within scanning distance (-1 = not calculated) */
	int nearbyVarCount=-1;

	/** Revised allele fraction accounting for technical biases (-1 = not calculated) */
	double revisedAlleleFraction=-1;
	/** Whether this variant was forced from input VCF file */
	private boolean forced=false;
	/** Whether this variant has been flagged for special attention */
	boolean flagged=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Enable calling of insertion variants */
	public static boolean CALL_INS=true;
	/** Enable calling of deletion variants */
	public static boolean CALL_DEL=true;
	/** Enable calling of substitution variants */
	public static boolean CALL_SUB=true;
	/** Enable calling of no-call variants (all N bases) */
	public static boolean CALL_NOCALL=false;
	/** Enable calling of junction variants from clipped reads */
	public static boolean CALL_JUNCTION=false;

	/** Include extended statistics in text output */
	public static boolean extendedText=true;

	/** Use dot genotype for variants that fail filters */
	public static boolean noPassDotGenotype=false;

	/** Enable homopolymer context scoring */
	public static boolean useHomopolymer=true;
	/** Enable alignment identity scoring */
	public static boolean useIdentity=true;
	/** Enable proper pairing rate scoring */
	public static boolean usePairing=true;
	/** Enable strand/read bias scoring */
	public static boolean useBias=true;
	/** Enable distance from read ends scoring */
	public static boolean useEdist=true;
	/** Enable scanning for nearby N bases (contig ends) */
	public static boolean doNscan=true;
	/** Penalty factor for low coverage variants */
	public static double lowCoveragePenalty=0.8;
	/** Maximum distance to scan for contig ends */
	public static int nScan=600;
	/** Minimum distance from contig ends before applying bias penalties */
	public static int minEndDistForBias=200;

	/** Minimum number of variant copies required for calling */
	public static int MIN_VAR_COPIES=0;

	/** Convert allele sequences to uppercase in output */
	public static final boolean UPPER_CASE_ALLELES=true;
	/** Verify that variants don't match reference (debugging) */
	private static final boolean TEST_REF_VARIANTS=false;

	/** Semicolon character for VCF INFO field separation */
	private static final byte colon=';';
	/** Tab character for field separation */
	private static final byte tab='\t';

	/** Human-readable names for variant types */
	public static final String[] typeArray=new String[] {"INS","NOCALL","SUB","DEL","LJUNCT","RJUNCT","BJUNCT","MULTI","COMPLEX"};

	/** Insertion variant type constant */
	public static final int INS=0;
	/** No-call variant type (length-neutral, reference N) */
	public static final int NOCALL=1;
	/** Substitution variant type (length-neutral) */
	public static final int SUB=2;
	/** Deletion variant type */
	public static final int DEL=3;
	/** Left-junction variant (left side clipped, right side normal) */
	public static final int LJUNCT=4;
	/** Right-junction variant (right side clipped, left side normal) */
	public static final int RJUNCT=5;
	/** Bidirectional junction (both sides clipped) */
	public static final int BJUNCT=6;
	/** Multiallelic variant (dominates all other types) */
	public static final int MULTI=7;
	/** Complex variant type */
	public static final int COMPLEX=8;

	/** Total number of variant types */
	public static final int VAR_TYPES=COMPLEX+1;

	/** Maps type initial letters to type constants for parsing */
	static final byte[] typeInitialArray=new byte[128];

	/** Pre-allocated empty allele array for deletions */
	static final byte[] AL_0=new byte[0];
	/** Pre-allocated single A allele */
	static final byte[] AL_A=new byte[] {(byte)'A'};
	/** Pre-allocated single C allele */
	static final byte[] AL_C=new byte[] {(byte)'C'};
	/** Pre-allocated single G allele */
	static final byte[] AL_G=new byte[] {(byte)'G'};
	/** Pre-allocated single T allele */
	static final byte[] AL_T=new byte[] {(byte)'T'};
	/** Pre-allocated single N allele */
	static final byte[] AL_N=new byte[] {(byte)'N'};
	/** Maps ASCII codes to pre-allocated allele arrays */
	static final byte[][] AL_MAP=makeMap();
	/** Random codes for hash function */
	static final int[] codes=makeCodes();

	/*--------------------------------------------------------------*/
	/*----------------     Bit Packing Constants    ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of bits for variant type in hash keys */
	private static final int typeBits=2;
	/** Number of bits for allele hash in hash keys */
	private static final int alleleBits=6;
	/** Number of bits for scaffold number in hash keys */
	private static final int scafBits=16;
	/** Number of bits for length in hash keys */
	private static final int lenBits=8;

	/** Bit shift for allele hash in hash keys */
	private static final int alleleShift=typeBits;
	/** Bit shift for scaffold number in hash keys */
	private static final int scafShift=alleleShift+alleleBits;
	/** Bit shift for length in hash keys */
	private static final int lenShift=scafShift+scafBits;
	/** Bit shift for start position in hash keys */
	private static final int startShift=lenShift+lenBits;

	/*--------------------------------------------------------------*/
	/*----------------    Static Initialization     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Generates array of random integers for hash function.
	 * Uses fixed seed for reproducible hash codes across runs.
	 * 
	 * @return Array of 256 random integers for hash mixing
	 */
	static final int[] makeCodes(){
		Random randy=new Random(1); // Fixed seed for reproducibility
		int[] array=new int[256];
		for(int i=0; i<array.length; i++){
			array[i]=randy.nextInt(); // Generate random hash codes
		}
		return array;
	}

	/**
	 * Creates mapping from ASCII codes to pre-allocated allele arrays.
	 * Avoids repeated allocation of single-base allele sequences.
	 * Maps both uppercase and lowercase bases to same arrays.
	 * 
	 * @return Mapping array where AL_MAP[ascii_code] = allele_array
	 */
	static final byte[][] makeMap(){
		byte[][] map=new byte[128][];
		
		// Map special characters and empty alleles
		map[0]=map['.']=map['\t']=AL_0; // Empty allele for deletions
		
		// Map DNA bases (both upper and lowercase)
		map['A']=map['a']=AL_A;
		map['C']=map['c']=AL_C;
		map['G']=map['g']=AL_G;
		map['T']=map['t']=AL_T;
		map['N']=map['n']=AL_N;
		
		// Fill remaining entries with single-byte arrays
		for(int i=0; i<map.length; i++){
			if(map[i]==null){
				map[i]=new byte[]{(byte)i}; // Single-byte array for unmapped characters
			}
		}
		return map;
	}

	/**
	 * Static initialization block for type parsing arrays.
	 * Maps variant type initial letters to type constants.
	 */
	static {
		Arrays.fill(typeInitialArray, (byte)-1); // Initialize to invalid
		
		// Map type initial letters to constants
		/**
		* Static initialization block for type parsing arrays.
		* Maps variant type initial letters to type constants.
		*/
		typeInitialArray['I']=INS; // Insertion
		typeInitialArray['N']=NOCALL; // No-call
		typeInitialArray['S']=SUB; // Substitution
		typeInitialArray['D']=DEL; // Deletion
		typeInitialArray['L']=LJUNCT; // Left junction
		typeInitialArray['R']=RJUNCT; // Right junction
		typeInitialArray['B']=BJUNCT; // Bidirectional junction
		typeInitialArray['M']=MULTI; // Multiallelic
		typeInitialArray['C']=COMPLEX; // Complex
	}

	/** Version string for VAR format output */
	public static final String varFormat="1.3";

}

