package var2;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import parse.Parse;
import shared.Tools;

/**
 * Mutable, reusable parallel of Var. Var's identity fields (scafnum/start/stop/allele/type/
 * hashcode) are final so Var can safely serve as a hash key; that immutability blocks in-place
 * reuse. VarM instead has all-mutable fields and a {@link #set(byte[], ScafMap)} method that
 * refills it from a VCF line with no allocation. VarM is NEVER a map key and has no hashcode.
 *
 * <p>VarM mirrors the FULL mutable field set of Var (not just the vector inputs) so it can also
 * serve other per-variant consumers such as VarFilter — primitive fields are cheap. Only Var's
 * final hashcode is intentionally omitted (a mutable object has no stable hash identity; carrying
 * one would be misleading, not an optimization).
 *
 * <p>The parse mirrors {@link VcfToVar#fromVCF} (parseCoverage+parseExtended) and the accessors are
 * copied verbatim from Var, so a VarM yields values bit-identical to VcfToVar.fromVCF(...). If a Var
 * accessor changes, change the copy here too. The non-vector score()/phredScore path is not mirrored.
 *
 * <p><b>Performance note — do NOT pool VarM for speed.</b> Reusing one VarM per thread was measured
 * SLOWER than allocating a fresh Var per variant on the giab net9b set (~204M candidate variants, 16
 * workers): about <b>+5.7% CPU</b> in the NN-scoring path and about <b>+29% CPU</b> in the no-NN
 * filter path (CROSSOVER results bit-identical either way). The JVM handles short-lived Var
 * allocation better than a hand-pooled mutable object (escape analysis / cheap young-gen GC), so the
 * gradevcf reusevar flag and the pipeline overloads VectorUMP45.makeVector(VarM) and
 * VarFilter.passesFilter(VarM) were removed. VarM is kept only as a mutable parsed-variant building
 * block; re-add those overloads if a genuine consumer appears.
 *
 * @author Ady
 */
public class VarM {

	/**
	 * Refills this VarM from a VCF line in place (no allocation), mirroring
	 * VcfToVar.fromVCF with parseCoverage=true, parseExtended=true. EVERY vector-relevant
	 * field is re-set each call so nothing is stale from a previously-parsed variant.
	 *
	 * @param line VCF record as a byte array
	 * @param scafMap Scaffold map for CHROM-name resolution (and later base lookups)
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

		// INFO begins after QUAL; record the position for the field searches below.
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
	/*----------------   Accessors (verbatim from Var)  -----------*/
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

	/** Insertion-length-adjusted allele fraction; caches into revisedAlleleFraction. */
	public double revisedAlleleFraction(double af, double readLengthAvg){
		if(revisedAlleleFraction!=-1){
			return revisedAlleleFraction;
		}else if(type()==Var.INS){
			return revisedAlleleFraction=adjustForInsertionLength(af, readLengthAvg);
		}
		return af;
	}

	/** Mean mapping quality over allele observations. */
	public double mapQAvg(){return mapQSum/(double)alleleCount();}

	/** Mean base quality over allele observations. */
	public double baseQAvg(){return baseQSum/(double)alleleCount();}

	/** Mean alignment identity (per-mille) over allele observations. */
	public double identityAvg(){return idSum/(double)alleleCount();}

	/** Mean read-end distance over allele observations. */
	public double edistAvg(){return endDistSum/(double)alleleCount();}

	/** Mean read length over allele observations. */
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

	/** Fraction of allele observations from proper pairs. */
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

	/** Homopolymer run length around this variant, or 0 if none/unavailable; copied from Var. */
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
			if(i<bases.length){return 0;}
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

	/** Distance to the nearest contig end (N-run boundary), capped at Var.nScan; copied from Var. */
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

	/** Distance to the nearest contig end on the left (10+ Ns = boundary); copied from Var. */
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

	/** Distance to the nearest contig end on the right (10+ Ns = boundary); copied from Var. */
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

	@Override
	public String toString(){
		return "VarM[scafnum="+scafnum+" start="+start+" stop="+stop+" type="+type
				+" allele="+(allele==null ? "null" : new String(allele))+" cov="+coverage
				+" r="+r1plus+","+r1minus+","+r2plus+","+r2minus+"]";
	}

	/*--------------------------------------------------------------*/
	/*----------------   INFO parse helpers (from VcfToVar)  ------*/
	/*--------------------------------------------------------------*/

	/** Parses an int from an INFO key=value, clamped to int range; -1 if absent. */
	private static int parseVcfIntDelimited(byte[] line, String query, int start){
		return (int)Tools.min(Integer.MAX_VALUE, parseVcfLongDelimited(line, query, start));
	}

	/** Parses a long from an INFO key=value (handles a leading '-'); -1 if absent. */
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

	/** Parses a double from an INFO key=value; -1.0 if absent. */
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

	/** Parses TYP= as an integer code or a string type name; -1 if absent/unrecognized. */
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
	/** Start position (0-based, inclusive). */
	public int start;
	/** Stop position (0-based, exclusive). */
	public int stop;
	/** Alternate allele sequence. */
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
