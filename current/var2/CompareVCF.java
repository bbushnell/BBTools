package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;

/**
 * Performs set operations on multiple VCF (Variant Call Format) files.
 * Supports difference, union, and intersection operations to compare genetic
 * variants across multiple samples. Can filter variants by quality score, restrict
 * to the intervals of a BED file, left-align indels (normalize), and optionally
 * split complex variants into components.
 *
 * @author Brian Bushnell
 * @contributor UMP45
 * @date January 14, 2017
 */
public class CompareVCF {
	
	public static void main(String[] args){
		Timer t=new Timer();
		CompareVCF x=new CompareVCF(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs a CompareVCF instance from command-line arguments.
	 * Parses operation mode (difference/union/intersection), input files,
	 * quality filters, and variant splitting options.
	 * @param args Command-line arguments containing file paths and options
	 */
	public CompareVCF(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		int mode_=DIFFERENCE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("difference") || a.equals("minus") || a.equals("dif") || a.equals("diff") || a.equals("subtraction") || a.equals("subtract")){
				mode_=DIFFERENCE;
			}else if(a.equals("union") || a.equals("plus")){
				mode_=UNION;
			}else if(a.equals("intersection") || a.equals("shared")){
				mode_=INTERSECTION;
			}else if(a.equals("addsamples")){
				addSamples=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("splitalleles")) {
				splitAlleles=Parse.parseBoolean(b);
			}else if(a.equals("splitsubs") || a.equals("splitsnps")) {
				splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitcomplex")) {
				splitComplex=Parse.parseBoolean(b);
			}else if(a.equals("sass") || a.equals("split")) {
				splitAlleles=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitall") || a.equals("sascsss")) {
				splitAlleles=splitComplex=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("ploidyout") || a.equals("outploidy")) {
				ploidyOut=Integer.parseInt(b);
			}else if(a.equals("minscore") || a.equals("minqual") || a.equals("minq")) {
				minScore=Double.parseDouble(b);
			}else if(a.equals("trimtocanonical") || a.equals("canonicalize") || a.equals("canonicize") || a.equals("canonize")){
				VCFLine.TRIM_TO_CANONICAL=Parse.parseBoolean(b);
			}else if(a.equals("normalize") || a.equals("leftalign") || a.equals("norm")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("bed") || a.equals("bedfile")){
				bedFile=b;
			}else if(a.equals("invertbed") || a.equals("excludebed") || a.equals("bedexclude") || a.equals("bedinvert")){
				invertBed=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		{//Process parser fields
			in1=(parser.in1==null ? null : parser.in1.split(","));
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null || in1.length<2){throw new RuntimeException("Error - at least two input files are required.");}

		if(ploidyOut>0){
			if(mode!=UNION && mode!=INTERSECTION){
				throw new RuntimeException("Error - ploidyout requires union or intersection; subtraction is undefined for ploidy merging.");
			}
			//Do NOT force splitalleles here: upconvert decomposes multiallelic lines itself (per-alt
			//dosage from the original GT), which requires that GT intact.  Pre-splitting would discard it.
		}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		ffin1=new FileFormat[in1.length];
		for(int i=0; i<in1.length; i++){
			ffin1[i]=FileFormat.testInput(in1[i], FileFormat.VCF, null, true, true);
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.VCF, null, true, overwrite, append, false);
		if(ffout1!=null && ffout1.type()==FileFormat.VAR){outputVar=true;}
		if(normalize && ref==null){
			throw new RuntimeException("Error - normalize/leftalign requires ref= (the reference fasta) for left-alignment.");
		}
		if(ref!=null){ScafMap.loadReference(ref, null, null, true);}
		if(bedFile!=null){
			bedMask=new BedMask(bedFile);
			outstream.println("Loaded "+bedMask.intervalsLoaded()+" BED intervals across "+bedMask.scaffolds()+" scaffolds.");
		}
	}
	
	/**
	 * Loads VCF variants from a file into a HashSet.
	 * Optionally splits complex variants and filters by minimum quality score.
	 * Updates processing statistics and header information.
	 *
	 * @param ff FileFormat for the input VCF file
	 * @param set Existing HashSet to add variants to, or null to create new set
	 * @return HashSet containing variants from the file
	 */
	public HashSet<VCFLine> getSet(FileFormat ff, HashSet<VCFLine> set){
		if(set==null){set=new HashSet<VCFLine>();}
		VCFFile vfile=new VCFFile(ff);
		samples.addAll(vfile.sampleNames);
		if(header==null){
			header=vfile.header;
			if(ScafMap.defaultScafMap()==null){
				ScafMap.setDefaultScafMap(vfile.toScafMap(null), ff.name());
			}
		}
		for(Entry<VCFLine, VCFLine> e : vfile.map.entrySet()){
			VCFLine v=e.getValue();
			ArrayList<VCFLine> list=null;
			if(splitAlleles || splitComplex || splitSubs){list=v.split(splitAlleles, splitComplex, splitSubs);}
			if(list==null || list.isEmpty()){
				if(normalize){v.leftAlign(refBases(v.scaf));}
				if(inRegion(v) && !set.contains(v) && v.qual>=minScore){set.add(v);}
			}else{
				for(VCFLine line : list){
					if(normalize){line.leftAlign(refBases(line.scaf));}
					if(inRegion(line) && !set.contains(line) && v.qual>=minScore){set.add(line);}
				}
			}
		}
		
		linesProcessed+=vfile.linesProcessed();
		headerLinesProcessed+=vfile.header.size();
		variantLinesProcessed+=vfile.map.size();
		bytesProcessed+=vfile.bytesProcessed();
		
		errorState|=vfile.errorState;
		return set;
	}

	/** Returns reference bases for a scaffold (for left-alignment), cached by name. */
	private byte[] refBases(String scaf){
		if(scaf!=lastScaf && !scaf.equals(lastScaf)){
			Scaffold sc=ScafMap.defaultScafMap().getScaffold(scaf);
			lastBases=(sc==null ? null : sc.bases);
			lastScaf=scaf;
		}
		return lastBases;
	}

	/** Whether a variant passes the BED region filter (always true when no bed= was given).
	 * @param v Variant to test
	 * @return true if the variant should be included in the comparison */
	private boolean inRegion(VCFLine v){
		return bedMask==null || (bedMask.contains(v.scaf, v.pos)^invertBed);
	}

	/**
	 * Creates the union of all variants from input VCF files.
	 * Combines all unique variants from all input files into a single set.
	 * @return HashSet containing all unique variants from all input files
	 */
	public HashSet<VCFLine> union(){
		if(ploidyOut>0){return upconvert(UNION);}
		final HashSet<VCFLine> set=new HashSet<VCFLine>();
		for(FileFormat ff : ffin1){
			getSet(ff, set);
		}
		return set;
	}
	
	/**
	 * Creates the intersection of variants across all input VCF files.
	 * Returns only variants that are present in every input file.
	 * @return HashSet containing variants common to all input files
	 */
	public HashSet<VCFLine> intersection(){
		if(ploidyOut>0){return upconvert(INTERSECTION);}
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.retainAll(set);}
		}
		return set0;
	}
	
	/**
	 * Creates the difference between the first VCF file and all subsequent files.
	 * Returns variants present in the first file but absent from all others.
	 * @return HashSet containing variants unique to the first input file
	 */
	public HashSet<VCFLine> difference(){
		assert(ploidyOut<=0) : "ploidyout is undefined for subtraction (should have been rejected at parse).";
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.removeAll(set);}
		}
		return set0;
	}

	/**
	 * Merges variants across all inputs (union or intersection) and rewrites each
	 * output line's genotype to the summed ALT dosage, producing a polyploid call set.
	 * Each diploid input contributes 0, 1, or 2 ALT copies per variant; the output GT is
	 * (ploidyOut-sum) reference alleles followed by (sum) alternate alleles.  Only UNION
	 * and INTERSECTION are valid here; subtraction is rejected at parse time.
	 * @param opMode UNION (variant present in any input) or INTERSECTION (present in all)
	 * @return Representative VCFLines with rewritten polyploid genotypes
	 */
	private HashSet<VCFLine> upconvert(int opMode){
		assert(opMode==UNION || opMode==INTERSECTION) : opMode;
		final int nInputs=ffin1.length;
		final HashMap<VCFLine, VCFLine> rep=new HashMap<VCFLine, VCFLine>();
		final HashMap<VCFLine, int[]> acc=new HashMap<VCFLine, int[]>();//value = {ALT dosage sum, input count}
		final int[] inputPloidy=new int[nInputs];

		for(int fi=0; fi<nInputs; fi++){
			final HashSet<VCFLine> set=getSet(ffin1[fi], null);
			int ploidyThis=-1;
			for(VCFLine v : set){
				final int alleles=gtAlleleCount(v);
				if(ploidyThis<0){ploidyThis=alleles;}
				if(v.isMulti()){
					decomposeMulti(v, rep, acc);//Multiallelic: split per-ALT, dosage = per-allele copy count
				}else{
					accumulate(rep, acc, v, altDosage(v));//Biallelic: dosage = ALT copies in the GT
				}
			}
			assert(ploidyThis>=1) : "Input "+fi+" had no readable genotypes.";
			inputPloidy[fi]=ploidyThis;
		}

		//Required invariant: the output ploidy equals the sum of the input ploidies.
		int ploidySum=0;
		for(int p : inputPloidy){ploidySum+=p;}
		assert(ploidySum==ploidyOut) : "ploidyout ("+ploidyOut+") must equal the sum of input ploidies ("+ploidySum+").";

		final HashSet<VCFLine> out=new HashSet<VCFLine>();
		for(VCFLine v : rep.values()){
			final int[] a=acc.get(v);
			final int sum=a[0], count=a[1];
			if(opMode==INTERSECTION && count<nInputs){continue;}//Not present in every input
			assert(sum>=1 && sum<=ploidyOut) : "ALT dosage "+sum+" out of [1,"+ploidyOut+"] at "+v;
			setPolyploidGT(v, sum, ploidyOut);
			out.add(v);
		}
		return out;
	}

	/**
	 * Decomposes a multiallelic line into one biallelic line per ALT allele actually carried by the
	 * genotype, accumulating each into the polyploid merge keyed by (scaf,pos,ref,alt).  The per-allele
	 * dosage is how many copies of that allele's index the GT contains, so 1/2 yields each of the two
	 * ALTs at dosage 1 while 2/2 yields the second ALT at dosage 2.  Keying by allele SEQUENCE (not
	 * index) makes this agree across inputs that list the same alleles in a different order, and lets a
	 * decomposed allele merge with the same allele called biallelically in another input.
	 * @param v Multiallelic source line (ALT contains a comma)
	 * @param rep Representative-line map being built
	 * @param acc Per-key {dosage sum, input count} accumulator
	 */
	private void decomposeMulti(VCFLine v, HashMap<VCFLine, VCFLine> rep, HashMap<VCFLine, int[]> acc){
		final int[] gt=parseGtIndices(v);
		final String[] alts=new String(v.alt).split(",");
		for(int k=1; k<=alts.length; k++){//Allele indices are 1-based (0 is the reference)
			int count=0;
			for(int idx : gt){if(idx==k){count++;}}
			if(count==0){continue;}//Listed in ALT but not present in this genotype
			final VCFLine lk=v.isolateAllele(alts[k-1].getBytes());
			if(lk==null){continue;}//Allele equals reference (defensive; real ALTs never do)
			if(normalize){lk.leftAlign(refBases(lk.scaf));}//Native lines were left-aligned in getSet
			accumulate(rep, acc, lk, count);
		}
	}

	/** Adds one variant's ALT dosage into the polyploid accumulator: creates the entry on first sight
	 * (recording it as the representative line) or sums into it and bumps the input count.
	 * @param rep Representative-line map
	 * @param acc Per-key {dosage sum, input count} accumulator
	 * @param line Variant key (and representative on first sight)
	 * @param dosage ALT copies to add */
	private static void accumulate(HashMap<VCFLine, VCFLine> rep, HashMap<VCFLine, int[]> acc, VCFLine line, int dosage){
		final int[] a=acc.get(line);
		if(a==null){rep.put(line, line); acc.put(line, new int[]{dosage, 1});}
		else{a[0]+=dosage; a[1]++;}
	}

	/**
	 * Parses the first sample's GT into its allele indices (0 = reference, 1+ = ALT positions).
	 * Reads only the leading GT subfield, accepts phased ('|') or unphased ('/') separators, and
	 * handles multi-digit indices; a missing allele ('.') becomes -1, which matches no ALT index.
	 * @param v Variant whose first sample GT is read
	 * @return Allele indices in genotype order
	 */
	private static int[] parseGtIndices(VCFLine v){
		assert(!v.samples.isEmpty()) : "No genotype to read: "+v;
		final byte[] s=v.samples.get(0);
		int end=0;
		while(end<s.length && s[end]!=':'){end++;}//GT is the first FORMAT subfield
		int alleles=1;
		for(int i=0; i<end; i++){if(s[i]=='/' || s[i]=='|'){alleles++;}}
		final int[] out=new int[alleles];
		int slot=0, num=-1;
		for(int i=0; i<end; i++){
			final byte c=s[i];
			if(c=='/' || c=='|'){out[slot++]=num; num=-1;}
			else if(c>='0' && c<='9'){num=(num<0 ? 0 : num)*10+(c-'0');}
			//'.' (or any other byte) leaves num at -1: a no-call allele, which matches no ALT index
		}
		out[slot]=num;
		return out;
	}

	/**
	 * Counts ALT-allele copies in the first sample's biallelic genotype.  Only called on biallelic
	 * lines (multiallelic lines are routed to {@link #decomposeMulti}); a residual multiallelic allele
	 * index (&gt;1) trips an assertion rather than being guessed at, guarding that invariant.
	 * @param v Variant whose first sample GT is read
	 * @return ALT copies: 1 for a het (0/1), 2 for a hom-alt (1/1)
	 */
	private static int altDosage(VCFLine v){
		assert(!v.samples.isEmpty()) : "No genotype to read: "+v;
		final byte[] s=v.samples.get(0);
		int count=0;
		for(int i=0; i<s.length && s[i]!=':'; i++){
			final byte c=s[i];
			if(c=='/' || c=='|' || c=='.'){continue;}//Separator or missing allele
			assert(c=='0' || c=='1') : "Residual multiallelic GT (allele index >1); run splitalleles: "+new String(s);
			if(c=='1'){count++;}
		}
		return count;
	}

	/**
	 * Returns the number of alleles in the first sample's GT (separators + 1), i.e. the input ploidy.
	 * @param v Variant whose first sample GT is read
	 * @return Allele count of the genotype
	 */
	private static int gtAlleleCount(VCFLine v){
		assert(!v.samples.isEmpty()) : "No genotype to read: "+v;
		final byte[] s=v.samples.get(0);
		int alleles=1;
		for(int i=0; i<s.length && s[i]!=':'; i++){
			if(s[i]=='/' || s[i]=='|'){alleles++;}
		}
		return alleles;
	}

	/**
	 * Rewrites the variant's first sample to the polyploid genotype: (ploidy-sum) reference
	 * alleles then (sum) alternate alleles, unphased, preserving any trailing FORMAT subfields
	 * (e.g. :DP).  Collapses to a single merged sample column.
	 * @param v Representative line to mutate
	 * @param sum Total ALT dosage across inputs
	 * @param ploidy Output ploidy
	 */
	private static void setPolyploidGT(VCFLine v, int sum, int ploidy){
		final String sample=new String(v.samples.get(0));
		final int colon=sample.indexOf(':');
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<ploidy; i++){
			if(i>0){bb.append('/');}
			bb.append(i<ploidy-sum ? '0' : '1');//Reference alleles first, then ALT
		}
		if(colon>=0){bb.append(sample.substring(colon));}//Preserve ":<trailing subfields>"
		v.samples.clear();
		v.samples.add(bb.toBytes());
	}

	/**
	 * Rewrites the #CHROM header line's sample columns into one merged column named by
	 * joining all input sample names with '+', recording the polyploid provenance
	 * (e.g. HG001+HG002).  Only invoked when ploidyout is active.
	 */
	private void renameMergedSample(){
		if(header==null || samples.isEmpty()){return;}
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<samples.size(); i++){
			if(i>0){sb.append('+');}
			sb.append(samples.get(i));
		}
		final String merged=sb.toString();
		for(int i=0; i<header.size(); i++){
			final byte[] line=header.get(i);
			if(Tools.startsWith(line, "#CHROM")){
				final String[] cols=new String(line).split("\t");
				assert(cols.length>=10) : "Malformed #CHROM line (need 9 fixed columns + sample): "+new String(line);
				final ByteBuilder bb=new ByteBuilder();
				for(int c=0; c<9; c++){//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
					if(c>0){bb.append('\t');}
					bb.append(cols[c]);
				}
				bb.append('\t').append(merged);
				header.set(i, bb.toBytes());
				break;
			}
		}
	}
	
	/**
	 * Performs the specified set operation and returns results as a sorted list.
	 * Executes difference, union, or intersection based on configured mode.
	 * @return Sorted ArrayList of variants after applying the set operation
	 */
	ArrayList<VCFLine> toList(){
		final HashSet<VCFLine> set;
		if(mode==DIFFERENCE){
			set=difference();
		}else if(mode==UNION){
			set=union();
		}else if(mode==INTERSECTION){
			set=intersection();
		}else{
			throw new RuntimeException("Unknown mode "+mode);
		}
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(set.size());
		list.addAll(set);
		Shared.sort(list);
		return list;
	}
	
	/**
	 * Main processing method that executes the VCF comparison and writes output.
	 * Performs the configured set operation and writes results to output file.
	 * Prints processing statistics and timing information.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		ArrayList<VCFLine> list=toList();
		
		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		{//Processing block
			if(ploidyOut>0){renameMergedSample();}
			for(byte[] line : header){
				headerLinesOut++;
				if(bsw!=null){bsw.println(line);}
			}
			ByteBuilder bb=new ByteBuilder(33000);
			for(VCFLine line : list){
				variantLinesOut++;
				if(outputVar){
					assert(false) : "TODO";
//					Var v=line.toVar();
//					v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, ScafMap.defaultScafMap());
				}else{
					line.toText(bb);
				}
				bb.nl();
				if(bb.length>=32000){
					if(bsw!=null){bsw.print(bb);}
					bb.clear();
				}
			}
			if(bb.length>0){
				if(bsw!=null){bsw.print(bb);}
				bb.clear();
			}
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Header Lines In:   \t"+headerLinesProcessed);
		outstream.println("Variant Lines In:  \t"+variantLinesProcessed);
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/

	private long linesProcessed=0;
	private long headerLinesProcessed=0;
	private long variantLinesProcessed=0;
	private long headerLinesOut=0;
	private long variantLinesOut=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

	public ArrayList<byte[]> header=null;
	public ArrayList<String> samples=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	
	private String in1[]=null;
	private String out1=null;
	private String ref=null;

	private final FileFormat ffin1[];
	private final FileFormat ffout1;
	
	public final int mode;
	
	public boolean addSamples=true;
	private boolean outputVar=false;

	boolean splitAlleles=false;
	boolean splitSubs=false;
	boolean splitComplex=false;
	boolean normalize=false;
	double minScore=-99999;
	int ploidyOut=0;

	/** BED file restricting which variants are compared; null disables region filtering */
	private String bedFile=null;
	/** When true, keep variants OUTSIDE the BED intervals instead of inside */
	private boolean invertBed=false;
	/** Loaded region mask, or null when no bed= was given */
	private BedMask bedMask=null;

	/** Cached scaffold name and reference bases for the most recent leftAlign lookup */
	private String lastScaf=null;
	private byte[] lastBases=null;
	
	/*--------------------------------------------------------------*/
	
	/** Constant for intersection operation mode (value 2) */
	/** Constant for union operation mode (value 1) */
	public static int DIFFERENCE=0, UNION=1, INTERSECTION=2;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
