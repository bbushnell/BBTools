package var2;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import parse.Parse;
import shared.Tools;

/**
 * Mutable parsed-variant record with identity and observation fields that can be refilled
 * by {@link #set(byte[], ScafMap)}. Parsing allocates a scaffold-name String and may copy
 * allele arrays; reuse avoids allocating another VarM, not all temporary allocation.
 * This class inherits Object's identity equality/hashCode rather than Var's value-key
 * semantics. A caller must coordinate mutation and avoid sharing partially refilled state.
 *
 * <p>The parser follows the extended-statistics path of {@link VcfToVar#fromVCF}, including
 * its deletion-coordinate convention. The accessors were ported from Var; later changes
 * to Var do not automatically update this separate implementation. The non-vector
 * score()/phredScore path is not mirrored.
 *
 * <p><b>Historical performance note — do NOT pool VarM for speed.</b> The original implementation
 * reported that reusing one VarM per thread was measured
 * SLOWER than allocating a fresh Var per variant on the giab net9b set (~204M candidate variants, 16
 * workers): about <b>+5.7% CPU</b> in the NN-scoring path and about <b>+29% CPU</b> in the no-NN
 * filter path (CROSSOVER results bit-identical either way). Those runs favored allocating Var, so the
 * gradevcf reusevar flag and the pipeline overloads VectorUMP45.makeVector(VarM) and
 * VarFilter.passesFilter(VarM) were removed. VarM is kept only as a mutable parsed-variant building
 * block; re-add those overloads if a genuine consumer appears.
 *
 * @author Ady
 * @contributor Sucrose
 */
public class VarM{

	/**
	 * Refills the record using VcfToVar's coverage/extended parsing rules.
	 * Decodes CHROM as US-ASCII, parses POS and copies or canonicalizes ALT.
	 * When ALT and REF lengths differ and ALT is nonempty, ALT loses its first byte
	 * and start becomes POS. Only a one-base REF has its length decremented, retaining
	 * the historical deletion span. Otherwise start is POS-1; stop is start plus REF length.
	 * This is not general allele normalization or multiallelic splitting.
	 *
	 * Searches for custom statistics after ALT using exact keys and semicolon
	 * delimiters; the R1P and PPC block markers must be found before their helpers
	 * are called. Missing individual statistics use the helpers' sentinel values,
	 * with most count fields floored at zero. DP must be positive under assertions.
	 * Successful parsing assigns every stored field, including resetting forced.
	 * Exceptions can leave a partially updated object; no rollback is performed.
	 *
	 * @param line Nonnull record in the expected single-ALT/custom-INFO layout
	 * @param scafMap Nonnull registry containing CHROM
	 */
	public void set(byte[] line, ScafMap scafMap){
		int a=0, b=0;

		// Field 0: CHROM - scaffold name
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		String scaf=new String(line, a, b-a, StandardCharsets.US_ASCII);
		b++;
		a=b;

		// Field 1: POS - position (1-based in VCF)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		int pos=Parse.parseInt(line, a, b);
		b++;
		a=b;

		// Field 2: ID - skip
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		b++;
		a=b;

		// Field 3: REF - reference allele (length only)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		int reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;

		// Field 4: ALT - alternate allele
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		byte[] alt;
		if(b<=a+1){
			alt=Var.AL_MAP[line[a]];
		}else{
			alt=Arrays.copyOfRange(line, a, b);
		}
		b++;
		a=b;

		// Search begins at QUAL, immediately after ALT; INFO is not separately bounded.
		int infoStart=b;

		// Convert VCF coordinates to internal format (identical to VcfToVar.fromVCF).
		final int start;
		if(alt.length!=reflen && alt.length>0){
			alt=Arrays.copyOfRange(alt, 1, alt.length);
			start=pos;
			if(alt.length==0){alt=Var.AL_0;}
			else if(alt.length==1 && Var.AL_MAP[alt[0]]!=null){alt=Var.AL_MAP[alt[0]];}
			if(reflen==1){reflen--;}
		}else{
			start=pos-1;
		}
		final int stop=start+reflen;

		final int scafNum=scafMap.getNumber(scaf);
		assert(scafNum>=0) : scaf+"\n"+scafMap.keySet()+"\n"+scafMap.altKeySet()+"\n";

		// Type: prefer the explicit TYP= field, else derive from coordinates+allele.
		int type=parseVcfType(line, infoStart);
		if(type<0){type=Var.typeStartStop(start, stop, alt);}

		// Identity (mutable here, unlike Var).
		this.scafnum=scafNum;
		this.start=start;
		this.stop=stop;
		this.allele=alt;
		this.type=type;

		// Coverage block: R1P=2;R1M=0;R2P=0;R2M=0;...DP=24  (MCOV is not read by the vector)
		infoStart=Tools.indexOfDelimited(line, "R1P=", infoStart, (byte)';');
		r1plus=Tools.max(0, parseVcfIntDelimited(line, "R1P=", infoStart));
		r1minus=Tools.max(0, parseVcfIntDelimited(line, "R1M=", infoStart));
		r2plus=Tools.max(0, parseVcfIntDelimited(line, "R2P=", infoStart));
		r2minus=Tools.max(0, parseVcfIntDelimited(line, "R2M=", infoStart));
		coverage=parseVcfIntDelimited(line, "DP=", infoStart);
		assert(coverage>0) : new String(line, infoStart, line.length-infoStart);
		minusCoverage=parseVcfIntDelimited(line, "MCOV=", infoStart);

		// Extended block: PPC=0;RAF=...;LS=...;MQS=...;FLG=...
		infoStart=Tools.indexOfDelimited(line, "PPC=", infoStart, (byte)';');
		properPairCount=Tools.max(0, parseVcfIntDelimited(line, "PPC=", infoStart));
		revisedAlleleFraction=parseVcfDoubleDelimited(line, "RAF=", infoStart);
		lengthSum=Tools.max(0, parseVcfLongDelimited(line, "LS=", infoStart));
		mapQSum=Tools.max(0, parseVcfLongDelimited(line, "MQS=", infoStart));
		mapQMax=Tools.max(0, parseVcfIntDelimited(line, "MQM=", infoStart));
		baseQSum=Tools.max(0, parseVcfLongDelimited(line, "BQS=", infoStart));
		baseQMax=Tools.max(0, parseVcfIntDelimited(line, "BQM=", infoStart));
		endDistSum=Tools.max(0, parseVcfLongDelimited(line, "EDS=", infoStart));
		endDistMax=Tools.max(0, parseVcfIntDelimited(line, "EDM=", infoStart));
		idSum=Tools.max(0, parseVcfLongDelimited(line, "IDS=", infoStart));
		idMax=Tools.max(0, parseVcfIntDelimited(line, "IDM=", infoStart));
		nearbyVarCount=Tools.max(0, parseVcfIntDelimited(line, "NVC=", infoStart));
		flagged=(Tools.max(0, parseVcfIntDelimited(line, "FLG=", infoStart))>0);
		forced=false; //fromVCF never sets forced; reset each call for reuse safety
	}

	/*--------------------------------------------------------------*/
	/*----------------    Accessors (ported from Var)   -----------*/
	/*--------------------------------------------------------------*/

	/** Variant type constant (Var.SUB/INS/DEL/...). */
	public int type(){return type;}

	/** Total allele observation count across strands and read pairs. */
	public int alleleCount(){return r1plus+r1minus+r2plus+r2minus;}

	/** Allele fraction = allele count / coverage, guarded against divide-by-zero. */
	public double alleleFraction(){
		int count=alleleCount();
		int cov=Tools.max(count, coverage, 1);
		return count/(double)cov;
	}

	/** Total read depth at this site. */
	public int coverage(){
		assert(coverage>-1) : coverage+", "+this;
		return coverage;
	}

	/** Returns cached RAF unless it is -1; otherwise caches an insertion adjustment.
	 * Non-insertions return af without caching it. */
	public double revisedAlleleFraction(double af, double readLengthAvg){
		if(revisedAlleleFraction!=-1){
			return revisedAlleleFraction;
		}else if(type()==Var.INS){
			return revisedAlleleFraction=adjustForInsertionLength(af, readLengthAvg);
		}
		return af;
	}

	/** Mapping-quality sum divided by allele count; zero-count division is not guarded. */
	public double mapQAvg(){return mapQSum/(double)alleleCount();}

	/** Base-quality sum divided by allele count; zero-count division is not guarded. */
	public double baseQAvg(){return baseQSum/(double)alleleCount();}

	/** Identity sum (per-mille units) divided by allele count, without a zero-count guard. */
	public double identityAvg(){return idSum/(double)alleleCount();}

	/** Read-end-distance sum divided by allele count, without a zero-count guard. */
	public double edistAvg(){return endDistSum/(double)alleleCount();}

	/** Read-length sum divided by allele count, without a zero-count guard. */
	public double lengthAvg(){return lengthSum/(double)alleleCount();}

	/** Strand balance ratio in (0,1]; 1.0 when plus and minus counts are equal. */
	public double strandRatio(){
		int plus=allelePlusCount();
		int minus=alleleMinusCount();
		if(plus==minus){return 1;}
		return (Tools.min(plus, minus)+1)/(double)Tools.max(plus, minus);
	}

	/** Allele observations on the plus strand (R1+R2). */
	public int allelePlusCount(){return r1plus+r2plus;}

	/** Allele observations on the minus strand (R1+R2). */
	public int alleleMinusCount(){return r1minus+r2minus;}

	/** Allele observations from read 1 (both strands). */
	public int r1AlleleCount(){return r1plus+r1minus;}

	/** Allele observations from read 2 (both strands). */
	public int r2AlleleCount(){return r2plus+r2minus;}

	/** Proper-pair count divided by allele count, without a zero-count guard. */
	public double properPairRate(){return properPairCount/(double)alleleCount();}

	/** Reference-span length (stop-start). */
	public int reflen(){return stop-start;}

	/** Alternate-allele length (0 for empty/'.'). */
	public int readlen(){
		return (allele==null || allele.length==0 || allele[0]=='.' ? 0 : allele.length);
	}

	/** Insertion-length bias correction for allele fraction; copied from Var. */
	public double adjustForInsertionLength(final double ratio, final double rlen0){
		if(type()!=Var.INS){return ratio;}
		final int ilen=readlen();
		if(ilen<2){return ratio;}

		final double rlen=Tools.max(ilen*1.2+6, rlen0);
		final double sites=rlen+ilen-1;
		final double goodSites=rlen-ilen*1.1-6;

		final double expectedFraction=goodSites/sites;
		final double revisedRatio=Tools.min(ratio/expectedFraction, 1-(1-ratio)*0.1);
		return revisedRatio;
	}

	/** Calculates the type-specific, bounded homopolymer statistic through VarHelper.
	 * Returns zero for absent map/bases or unsupported types. Deletions use this
	 * parser's historical coordinate convention; this is not an unrestricted run length. */
	public int homopolymerCount(ScafMap map){
		if(map==null){return 0;}
		final byte[] bases=map.getScaffold(scafnum).bases;
		if(bases==null){return 0;}

		final int type=type();
		if(type==Var.SUB){
			assert(start==stop-1) : start+", "+stop;
			final byte base=allele[0];
			int x=VarHelper.homopolymerCountSub(bases, start, base);
			return x;
		}else if(type==Var.INS){
			final byte base1=allele[0], base2=allele[allele.length-1];
			int i=0;
			while(i<allele.length && allele[i]==base1){i++;}
			while(i<allele.length && allele[i]==base2){i++;}
			if(i<allele.length){return 0;}
			int left=VarHelper.homopolymerCountLeft(bases, start, base1);
			int right=VarHelper.homopolymerCountRight(bases, stop+1, base2);
			return left+right+1;
		}else if(type==Var.DEL){
			if(start<0 || start+1>=bases.length || stop<=0 || stop>=bases.length){return 0;}
			final byte base1=bases[start+1], base2=bases[stop-1];
			int pos=start+1;
			while(pos<=stop && bases[pos]==base1){pos++;}
			while(pos<=stop && bases[pos]==base2){pos++;}
			if(pos<=stop){return 0;}
			int left=VarHelper.homopolymerCountLeft(bases, start, base1);
			int right=VarHelper.homopolymerCountRight(bases, stop, base2);
			return left+right+1;
		}else{
			return 0;
		}
	}

	/** Starts with distance to the scaffold edge, then searches for runs of ten bases
	 * that are not fully defined. Var.nScan bounds the search, not necessarily the
	 * returned value: an unsuccessful scan can return its limit plus one. Without
	 * bases or with scanning disabled, returns the scaffold-edge distance. */
	public int contigEndDist(ScafMap map){
		Scaffold scaf=map.getScaffold(scafnum);
		int len=scaf.length;
		byte[] bases=scaf.bases;

		int scafEndDist=Tools.max(0, Tools.min(start, len-stop));
		if(bases==null || Var.nScan<1){return scafEndDist;}
		int limit=Tools.min(Var.nScan, scafEndDist);
		int contigEndDist=leftContigEndDist(bases, limit);
		limit=Tools.min(limit, contigEndDist);
		contigEndDist=rightContigEndDist(bases, limit);
		return Tools.min(scafEndDist, contigEndDist);
	}

	/** Scans left from start for ten consecutive not-fully-defined bases.
	 * Returns distance to the near edge of the first complete run found, or maxDist+1.
	 * A start beyond the sequence returns min(bases.length, maxDist+1). */
	public int leftContigEndDist(byte[] bases, int maxDist){
		if(start>=bases.length){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=start, lim=Tools.max(0, start-maxDist); i>=lim; i--){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=start-i-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}

	/** Scans right from stop for ten consecutive not-fully-defined bases.
	 * Returns distance to the near edge of the first complete run found, or maxDist+1.
	 * A negative stop returns min(bases.length, maxDist+1). */
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

	/** Describes identity, allele, coverage and strand/read counts for diagnostics. */
	@Override
	public String toString(){
		return "VarM[scafnum="+scafnum+" start="+start+" stop="+stop+" type="+type
				+" allele="+(allele==null ? "null" : new String(allele))+" cov="+coverage
				+" r="+r1plus+","+r1minus+","+r2plus+","+r2minus+"]";
	}

	/*--------------------------------------------------------------*/
	/*----------------   INFO parse helpers (from VcfToVar)  ------*/
	/*--------------------------------------------------------------*/

	/** Parses a long, caps it at Integer.MAX_VALUE, then casts to int.
	 * No lower-bound clamp or long-overflow detection is performed; absent keys yield -1. */
	private static int parseVcfIntDelimited(byte[] line, String query, int start){
		return (int)Tools.min(Integer.MAX_VALUE, parseVcfLongDelimited(line, query, start));
	}

	/** Searches at start or after a semicolon, then parses digits with optional '-'.
	 * Stops at tab/semicolon (asserting on other terminators); returns -1 if absent.
	 * Accumulation has no overflow check. A match at byte zero returns zero unparsed.
	 * start must be a valid nonnegative search offset. */
	private static long parseVcfLongDelimited(byte[] line, String query, final int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
		if(loc<0){return -1;}
		long current=0;
		long mult=1;
		if(loc>0){
			if(line[loc+query.length()]=='-'){mult=-1; loc++;}
			for(int i=loc+query.length(); i<line.length; i++){
				final byte x=line[i];
				if(Tools.isDigit(x)){
					current=current*10+(x-'0');
				}else{
					assert(x==tab || x==colon) : x+", "+query+", "+new String(line, loc, i-loc+1);
					break;
				}
			}
		}
		return mult*current;
	}

	/** Parses the span up to tab/semicolon using Parse.parseDouble.
	 * Returns -1 for an absent key or a match at byte zero; start must be nonnegative. */
	private static double parseVcfDoubleDelimited(byte[] line, String query, int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
		if(loc<0){return -1;}
		if(loc>0){
			loc=loc+query.length();
			int loc2=loc+1;
			while(loc2<line.length){
				byte b=line[loc2];
				if(b==tab || b==colon){break;}
				loc2++;
			}
			return Parse.parseDouble(line, loc, loc2);
		}
		return -1;
	}

	/** Parses numeric TYP or maps only its first byte through Var.typeInitialArray.
	 * Numeric codes are not range-checked here; absent/unrecognized initials yield -1. */
	private static int parseVcfType(byte[] line, int start){
		int idx=Tools.indexOfDelimited(line, "TYP=", start, (byte)';');
		if(idx<0){return -1;}
		idx+=4;
		if(idx>=line.length){return -1;}
		byte b=line[idx];
		if(b>='0' && b<='9'){
			return parseVcfIntDelimited(line, "TYP=", start);
		}
		if(b<Var.typeInitialArray.length){
			int t=Var.typeInitialArray[b];
			if(t>=0){return t;}
		}
		return -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Fields (all mutable)  ------------*/
	/*--------------------------------------------------------------*/

	/** Scaffold/chromosome number. */
	public int scafnum;
	/** Start from the parser's POS/REF/ALT conversion; see set for the deletion convention. */
	public int start;
	/** Parser start plus adjusted REF length; see set for the deletion convention. */
	public int stop;
	/** Alternate allele; may reference a shared Var canonical array, so do not mutate its contents. */
	public byte[] allele;
	/** Variant type constant (Var.SUB/INS/DEL/...). */
	public int type;

	/** Total read depth (DP). */
	int coverage=-1;
	/** Minus-strand coverage (MCOV); parity with Var, not read by the vector. */
	int minusCoverage=-1;
	/** Strand/read-pair allele counts (R1P/R1M/R2P/R2M). */
	int r1plus, r1minus, r2plus, r2minus;
	/** Proper-pair allele count (PPC). */
	int properPairCount;
	/** Mapping-quality sum (MQS) and max (MQM). */
	long mapQSum;
	public int mapQMax;
	/** Base-quality sum (BQS) and max (BQM). */
	long baseQSum;
	public int baseQMax;
	/** Read-end-distance sum (EDS) and max (EDM). */
	long endDistSum;
	public int endDistMax;
	/** Alignment-identity sum (IDS) and max (IDM). */
	long idSum;
	int idMax;
	/** Read-length sum (LS). */
	long lengthSum;
	/** Nearby variant count (NVC). */
	int nearbyVarCount=-1;
	/** Cached revised allele fraction (RAF); -1 = not yet computed. */
	double revisedAlleleFraction=-1;
	/** Forced-call flag (parity with Var; not parsed from VCF, reset each set()). */
	boolean forced=false;
	/** Flagged status (FLG>0); parity with Var, not read by the vector. */
	boolean flagged=false;

	/** Tab and semicolon delimiters for the INFO parse helpers. */
	private static final byte tab='\t', colon=';';
}
