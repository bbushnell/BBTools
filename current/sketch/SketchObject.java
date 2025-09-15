package sketch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import aligner.IDAligner;
import aligner.SingleStateAlignerFlat2;
import aligner.SingleStateAlignerFlatFloat;
import dna.AminoAcid;
import dna.Data;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GeneModelParser;
import prok.ProkObject;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.IntList3;
import tax.TaxTree;
import tracker.EntropyTracker;

public class SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		for(int i=0; i<1000000; i++){
			wkidToAniExact(0.5, 32, 23, 0.000005);
		}
		t.stop();
		System.out.println(t);
	}
	
	public static boolean parseSketchFlags(String arg, String a, String b){
		
		if(parseCoding(arg, a, b)){
			//
		}
		
		else if(a.equalsIgnoreCase("k")){
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "Bad argument "+arg;
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
				if(k==k2){k2=0;}
			}else{
				k=Integer.parseInt(b);
				k2=0;
			}
			setK=true;
		}else if(a.equalsIgnoreCase("rcomp")){
			rcomp=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("minfakeid")){
			minFakeID=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("hashNames")){
			hashNames=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("preferSSUMap")){
			preferSSUMap=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("preferSSUMapForEuks") || a.equalsIgnoreCase("preferSSUMapEuks")){
			preferSSUMapForEuks=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useSSUMapOnly")){
			useSSUMapOnly=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useSSUMapOnlyForEuks") || a.equalsIgnoreCase("SSUMapOnlyEuks")){
			useSSUMapOnlyForEuks=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("skipNonCanonical")){
			skipNonCanonical=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useSizeEstimateInScore") || a.equalsIgnoreCase("useSizeEstimate")){
			useSizeEstimate=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("blacklist") || a.equalsIgnoreCase("bl")){
			blacklist=b;
		}else if(a.equalsIgnoreCase("whitelist") || a.equalsIgnoreCase("usewhitelist")){
			useWhitelist=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sampleseed")){
			sampleseed=Long.parseLong(b);
		}
		
		else if(a.equalsIgnoreCase("size") || a.equalsIgnoreCase("sketchsize")){
			if("auto".equalsIgnoreCase(b)){
				AUTOSIZE=true;
				AUTOSIZE_LINEAR=false;
				SET_AUTOSIZE=true;
			}else if("linear".equalsIgnoreCase(b)){
				AUTOSIZE=false;
				AUTOSIZE_LINEAR=true;
				SET_AUTOSIZE=true;
			}else{
				AUTOSIZE=false;
				AUTOSIZE_LINEAR=false;
				targetSketchSize=Parse.parseIntKMG(b);
				SET_TARGET_SIZE=true;
			}
		}else if(a.equalsIgnoreCase("maxfraction") || a.equalsIgnoreCase("maxgenomefraction") || a.equalsIgnoreCase("mgf")){
			maxGenomeFraction=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("smallsketchsize")){
			smallSketchSize=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minSketchSize")){
			minSketchSize=Tools.max(1, Integer.parseInt(b));
		}else if(a.equalsIgnoreCase("autosize")){
			AUTOSIZE=Parse.parseBoolean(b);
			if(AUTOSIZE){AUTOSIZE_LINEAR=false;}
			SET_AUTOSIZE=true;
		}else if(a.equalsIgnoreCase("autosizefactor") || a.equalsIgnoreCase("autosizefraction") || a.equalsIgnoreCase("autosizemult") || a.equalsIgnoreCase("sizemult")){
			AUTOSIZE_FACTOR=Float.parseFloat(b);
			SET_AUTOSIZE_FACTOR=true;
			AUTOSIZE_LINEAR=false;
			SET_AUTOSIZE=true;
		}else if(a.equalsIgnoreCase("linear") || a.equalsIgnoreCase("lineardensity") || a.equalsIgnoreCase("density")){
			if(Tools.isNumeric(b)){
				AUTOSIZE_LINEAR_DENSITY=Double.parseDouble(b);
				AUTOSIZE_LINEAR=true;
				AUTOSIZE=false;
				SET_AUTOSIZE=true;
			}else{
				AUTOSIZE_LINEAR=Parse.parseBoolean(b);
				if(AUTOSIZE_LINEAR){AUTOSIZE=false;}
				SET_AUTOSIZE=true;
			}
		}else if(a.equalsIgnoreCase("maxGenomeFractionSmall")){
			maxGenomeFractionSmall=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("keyFraction")){
			double d=Double.parseDouble(b);
			setKeyFraction(d);
		}else if(a.equalsIgnoreCase("prealloc")){
			if(b==null || Character.isLetter(b.charAt(0))){
				if(Parse.parseBoolean(b)){
					prealloc=1.0f;
				}else{
					prealloc=0;
				}
			}else{
				prealloc=Float.parseFloat(b);
			}
		}
		
		else if(a.equalsIgnoreCase("intmap")){
			SketchIndex.useIntMap=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("intmapsize")){
			SketchIndex.intMapSize=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("bitsetbits")){
			assert(false) : "bitsetbits should be 2.";
//			bitSetBits=Integer.parseInt(b);
		}
		
//		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount")){
//			minKeyOccuranceCount=Integer.parseInt(b);
//		}
		else if(a.equalsIgnoreCase("sketchHeapFactor")){
			sketchHeapFactor=Tools.max(1.0, Double.parseDouble(b));
		}
		
		else if(a.equalsIgnoreCase("killok")){
			KILL_OK=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("ssa") || a.equalsIgnoreCase("usessa")){
			useSSA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("ssa3") || a.equalsIgnoreCase("usessa3")){
			useSSA3=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("exactani")){
			EXACT_ANI=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("translate") || a.equalsIgnoreCase("callgenes")){
			translate=Parse.parseBoolean(b);
			defaultParams.translate=translate;
		}else if(a.equalsIgnoreCase("sixframes") || a.equalsIgnoreCase("6frames")){
			sixframes=Parse.parseBoolean(b);
			defaultParams.sixframes=sixframes;
			if(sixframes){
				translate=defaultParams.translate=true;
			}
		}else if(a.equalsIgnoreCase("processSSU") || a.equalsIgnoreCase("process16S") || a.equalsIgnoreCase("16S") || a.equalsIgnoreCase("SSU")){
			processSSU=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("hashSeed")){
			hashSeed=Long.parseLong(b);
			//remakeCodes(hashSeed); //Happens in postParse
		}
		
		else if(a.equalsIgnoreCase("forceDisableMultithreadedFastq")){
			forceDisableMultithreadedFastq=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("verifyentropy")){
			EntropyTracker.verify=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("entropyK")){
			entropyK=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("fastentropy") || a.equalsIgnoreCase("fentropy")){
			if(Parse.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.FAST;}
		}else if(a.equalsIgnoreCase("mediumentropy") || a.equalsIgnoreCase("mentropy")){
			if(Parse.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.MEDIUM;}
		}else if(a.equalsIgnoreCase("slowentropy") || a.equalsIgnoreCase("sentropy")){
			if(Parse.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.SLOW;}
		}else if(a.equalsIgnoreCase("superslowentropy") || a.equalsIgnoreCase("ssentropy")){
			if(Parse.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.SUPERSLOW;}
		}else if(a.equalsIgnoreCase("verbose2")){
			verbose2=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("loadSketchesFromSketchFile2")){
			LOADER2=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("useToValue2") || a.equalsIgnoreCase("ToValue2")){
			useToValue2=Parse.parseBoolean(b);
		}else if(a.equals("ssu") || a.equals("ssufile")){
			assert(false) : "ssufile is deprecated; please specify 16Sfile and/or 18Sfile independently";
		}else if(a.equalsIgnoreCase("16Sfile")/* || a.equalsIgnoreCase("ssufile") || a.equalsIgnoreCase("silvafile")*/){
			SSUMap.r16SFile=b;
			if("auto".equalsIgnoreCase(SSUMap.r16SFile)){SSUMap.r16SFile=TaxTree.default16SFile();}
			assert(SSUMap.r16SMap==null);
		}else if(a.equalsIgnoreCase("18Sfile")){
			SSUMap.r18SFile=b;
			if("auto".equalsIgnoreCase(SSUMap.r18SFile)){SSUMap.r18SFile=TaxTree.default18SFile();}
			assert(SSUMap.r18SMap==null);
		}
		
//		else if(a.equalsIgnoreCase("minHashValue")){
//			minHashValue=Tools.max(1, Long.parseLong(b));
//		}
		
		else{
			return false;
		}
		return true;
	}
	
	private static boolean parseCoding(String arg, String a, String b){
		if(a.equalsIgnoreCase("deltaout") || a.equalsIgnoreCase("delta")){
			deltaOut=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("a33") || a.equalsIgnoreCase("a48")){
			boolean x=Parse.parseBoolean(b);
			if(x){CODING=A48;}
			else if(CODING==A48){CODING=HEX;}
		}else if(a.equalsIgnoreCase("hex")){
			boolean x=Parse.parseBoolean(b);
			if(x){CODING=HEX;}
			else if(CODING==HEX){CODING=A48;}
		}else if(a.equalsIgnoreCase("raw")){
			boolean x=Parse.parseBoolean(b);
			if(x){CODING=RAW;}
			else if(CODING==RAW){CODING=A48;}
		}else{
			return false;
		}
		return true;
	}
	
	static int parseMode(String[] args){
		int mode=defaultParams.mode;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			int x=parseMode(arg, a, b);
			if(x>-1){mode=x;}
		}
		return mode;
	}
	
	static int parseMode(String arg, String a, String b){
		int mode_=-1;
		if(a.equalsIgnoreCase("mode")){
			assert(b!=null) : "Bad parameter: "+arg;
			if(Tools.isDigit(b.charAt(0))){
				mode_=Integer.parseInt(b);
			}else if(b.equalsIgnoreCase("single") || b.equalsIgnoreCase("onesketch")){
				mode_=ONE_SKETCH;
			}else if(b.equalsIgnoreCase("taxa") || b.equalsIgnoreCase("pertaxa")){
				mode_=PER_TAXA;
			}else if(b.equalsIgnoreCase("sequence") || b.equalsIgnoreCase("persequence") || b.equalsIgnoreCase("header") || b.equalsIgnoreCase("perheader")){
				mode_=PER_SEQUENCE;
//			}else if(b.equalsIgnoreCase("header") || b.equalsIgnoreCase("perheader")){
//				mode_=PER_HEADER;
			}else if(b.equalsIgnoreCase("img")){
				mode_=PER_IMG;
			}else if(b.equalsIgnoreCase("perfile") || b.equalsIgnoreCase("file")){
				mode_=PER_FILE;
			}else{
				assert(false) : "Bad parameter: "+arg;
			}
		}else if(a.equalsIgnoreCase("single") || a.equalsIgnoreCase("onesketch")){
			mode_=ONE_SKETCH;
		}else if(a.equalsIgnoreCase("taxa") || a.equalsIgnoreCase("pertaxa")){
			mode_=PER_TAXA;
		}else if(a.equalsIgnoreCase("sequence") || a.equalsIgnoreCase("persequence") || a.equalsIgnoreCase("header") || a.equalsIgnoreCase("perheader")){
			mode_=PER_SEQUENCE;
//		}else if(a.equalsIgnoreCase("header") || a.equalsIgnoreCase("perheader")){
//			mode_=PER_HEADER;
		}else if(a.equalsIgnoreCase("perfile") || a.equalsIgnoreCase("file")){
			mode_=PER_FILE;
		}
		return mode_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	static synchronized void setTaxtree(String taxTreeFile, PrintStream outstream){
		if(taxTreeFile==null){
			taxtree=null;
			return;
		}
		if(treefile!=null){
			assert(!treefile.equals(taxTreeFile));
			if(treefile.equals(taxTreeFile)){return;}
			treefile=taxTreeFile;
		}
		taxtree=TaxTree.loadTaxTree(taxTreeFile, outstream, hashNames, false);
	}
	
	public static void reset(){
		postparsed=false;
		blacklist=null;
		useWhitelist=false;
		defaultParams=new DisplayParams();
		amino=false;
		translate=false;
		sixframes=false;
	}
	
	public static void postParse(){
		if(postparsed){return;}
		postparsed=true;
		IntList3.defaultMode=IntList3.UNIQUE; //Not really safe, if Seal uses Sketch...
		
		if(defaultParams.amino()){amino=true;}
		if(amino){Shared.AMINO_IN=true;}
		
		if(amino8){
			AminoAcid.AMINO_SHIFT=3;
			System.err.println("Set AMINO_SHIFT to "+AminoAcid.AMINO_SHIFT);
		}
		
		if(amino || translate){
			rcomp=false;
			if(k>12){
				int oldk=k;
				int oldk2=k2;
//				k=Tools.min(k, 63/AminoAcid.AMINO_SHIFT);
//				k2=Tools.min(k2, k);
				k=defaultKAmino;
				k2=defaultK2Amino;
				if(k==k2){k2=0;}
				if(k!=oldk || k2!=oldk2){System.err.println("Set k to "+k+(k2>0 ? ","+k2 : ""));}
			}
//			AUTOSIZE_FACTOR=(AUTOSIZE_FACTOR*defaultK)/k;
//			System.err.println("Set AUTOSIZE_FACTOR to "+AUTOSIZE_FACTOR);
		}
		
		if(aminoOrTranslate()){
			bitsPerBase=5;//Maybe it is safe to keep these at 4 and 8 or 5 and 8; need to check.
			bitsPerCycle=10;
		}else{
			bitsPerBase=2;
			bitsPerCycle=8;
		}
		basesPerCycle=bitsPerCycle/bitsPerBase;
		hashCycles=((k*bitsPerBase)+bitsPerCycle-1)/bitsPerCycle;
		
		cycleMask=~((-1L)<<bitsPerCycle);
		maxCycles=(64+bitsPerCycle-1)/bitsPerCycle;
		codeIncrement=(int)(cycleMask+1);
		codes=makeCodes(maxCycles, codeIncrement, hashSeed, false);
		codes1D=makeCodes1D(codes);
		
		if(k2>0){
			assert(k2<k) : "k2 must be less than k.";
//			assert(k2%basesPerCycle==0) : "k2 must be a multiple of "+basesPerCycle; //No longer required!  Still recommended for speed though.
			
			k2mask=~((-1L)<<(bitsPerBase*k2));
			k2submask=~((-1L)<<(bitsPerBase*(k2%basesPerCycle)));
			k2shift=(k-k2); //for toValue2
			k2midmask=(k2mask<<((k2shift*bitsPerBase)/2)); //for toValue2
			hashCycles2=k2/basesPerCycle;
		}else{
			useToValue2=false;
		}
		codeMax=hashCycles*codeIncrement;
		codeMax2=hashCycles2*codeIncrement;
		
//		hasher=k2<1 ? new Hasher1() : k2submask==0 ? new Hasher2() : new Hasher3();
		
		if(translate || processSSU){
			if(pgmFile==null){
				pgmFile=Data.findPath("?model.pgm");
			}
			pgm=GeneModelParser.loadModel(pgmFile);
			GeneCaller.call16S=processSSU;
			GeneCaller.call18S=false;//processSSU;
			GeneCaller.call23S=false;
			GeneCaller.call5S=false;
			GeneCaller.calltRNA=false;
			GeneCaller.callCDS=translate;
			
			if(GeneCaller.call16S || GeneCaller.call18S || GeneCaller.call23S){
				ProkObject.loadLongKmers();
				ProkObject.loadConsensusSequenceFromFile(true, true);
			}
		}
		
//		if(HASH_VERSION<2 && useToValue2){HASH_VERSION=2;}
//		else if(HASH_VERSION==2 && !useToValue2){HASH_VERSION=1;}
//		assert(blacklist!=null) : blacklist+", "+(blacklist==null);
		if(blacklist!=null){Blacklist.addFiles(blacklist);}
		
//		System.err.println("amino="+amino+", k="+k+", k2="+k2+", bitsPerCycle="+bitsPerCycle+", bitsPerBase="+bitsPerBase+
//				", basesPerCycle="+basesPerCycle+", hashCycles="+hashCycles+", k2mask="+k2mask+", k2submask="+k2submask+", hashCycles2="+hashCycles2+
//				", codeMax="+codeMax+", codeMax2="+codeMax2);
//		structures.IntList3.ascending=false;//123 for testing
		
		alignerPool=new AlignmentThreadPool(Shared.threads());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void remakeCodes(long hashSeed){
		long[][] codes2=makeCodes(maxCycles, codeIncrement, hashSeed, false);
		long[] codes1D2=makeCodes1D(codes2);
		for(int i=0; i<codes2.length; i++){
			for(int j=0; j<codes2[i].length; j++){
				codes[i][j]=codes2[i][j];
			}
		}
		for(int i=0; i<codes1D2.length; i++){
				codes1D[i]=codes1D2[i];
		}
	}
	
	public static long[][] makeCodes(int symbols, int modes, long seed, boolean positive){
		Random randy=new Random(seed);
		long mask=positive ? Long.MAX_VALUE : -1L;
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong()&mask;
			}
		}
		for(int i=0; i<3; i++) {
			antialias(r, randy);
		}
		return r;
	}
	
	private static void antialias(long[][] matrix, Random randy){
		for(long[] array : matrix){
			antialias(array, randy);
		}
	}
	
	private static void antialias(long[] array, Random randy){
		for(int i=0; i<64; i++){
			antialiasNumbers(array, randy);
			antialiasBit(array, randy, i);
		}
	}
	
	private static void antialiasBit(long[] array, Random randy, int bit){
		int half=array.length/2;
		long ones=0;
		for(int i=0; i<array.length; i++){
			ones+=(array[i]>>bit)&1L;
		}
		final long orMask=1L<<bit;
		final long andMask=~orMask;
		while(ones<half-1){
			int loc=randy.nextInt(array.length);
			while((array[loc]&orMask)!=0){
				loc=randy.nextInt(array.length);
			}
			array[loc]|=orMask;
			ones++;
		}
		while(ones>half+1){
			int loc=randy.nextInt(array.length);
			while((array[loc]&orMask)==0){
				loc=randy.nextInt(array.length);
			}
			array[loc]&=andMask;
			ones--;
		}
	}
	
	private static void antialiasNumbers(long[] array, Random randy){
		for(int i=0; i<array.length; i++){
			array[i]=antialiasNumber(array[i], randy);
		}
	}
	
	private static long antialiasNumber(long number, Random randy){
		int ones=Long.bitCount(number);
		while(Long.bitCount(number)<31){
			number=number|(1L<<randy.nextInt(64));
		}
		while(Long.bitCount(number)>33){
			number=number&~(1L<<randy.nextInt(64));
		}
		return number;
	}
	
//	public static long[] makeCodes1D(int symbols, int modes, long seed, boolean positive){
//		Random randy=Shared.threadLocalRandom(seed);
//		long mask=positive ? Long.MAX_VALUE : -1L;
//		long[] r=new long[symbols*modes];
//		for(int i=0; i<r.length; i++){
//			r[i]=randy.nextLong()&mask;
//		}
//		return r;
//	}
	
	public static long[] makeCodes1D(long[][] codes2D){
		final int len=codes2D.length*codes2D[0].length;
		long[] codes1D=new long[len];
		int k=0;
		for(long[] array : codes2D){
			for(long x : array){
				codes1D[k]=x;
				k++;
			}
		}
		return codes1D;
	}
	
	public static final long hash(long kmer){//Avoid using this.
		return rcomp ? hash(kmer, AminoAcid.reverseComplementBinaryFast(kmer,  k), rcomp) : 
			hash(kmer, kmer, rcomp);
	}
	
	public static final long hash(long kmer, long rkmer){
		assert(postparsed);
		if(useToValue2){return hashToValue2(kmer, rkmer, rcomp);}
		final long key=toValue(kmer, rkmer);
		return k2<1 ? hash1(key) : k2submask==0 ? hash2(key) : hash3(key);
	}
	
	public static final long hash(long kmer, long rkmer, boolean rcomp){
		assert(postparsed);
		if(useToValue2){return hashToValue2(kmer, rkmer, rcomp);}
		final long key=rcomp ? toValue(kmer, rkmer) : kmer;
		return k2<1 ? hash1(key) : k2submask==0 ? hash2(key) : hash3(key);
	}
	
	/**
	 * New version.
	 * Generates a hash code from a kmer.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hashToValue2(final long kmer0, final long rkmer0, final boolean rcomp){
//		return Tools.max(Tools.hash64shift(kmer0), Tools.hash64shift(kmer0&0xFFFFFFFFFFFFL));//Something like this would be around 10% faster; not worth it.
		long kmer=kmer0, rkmer=rkmer0;
		final long key;
		final boolean useK1;
		if(rcomp){
			assert(!amino && !translate);
			final long kmer2=((kmer&k2midmask)>>>k2shift);
			final long rkmer2=((rkmer&k2midmask)>>>k2shift);
			final long max2=Tools.max(kmer2, rkmer2);
//			long cmasked=max2&kmerChoiceMask;
//			useK1=kmerChoiceBitset.get((int)cmasked);
			useK1=((max2%4999L)&1L)==0L;
			key=useK1 ? Tools.max(kmer, rkmer) : max2;
//			System.err.println("\n"+longer+", "+max2+", "+(max2%4999L)+", "+((max2%4999L)&1L)+
//					", "+Long.toHexString(kmer)+", "+Long.toHexString(max2)+", "+Long.toHexString(value));
//			long value=longer ? (kmer2>rkmer2 ? kmer : rkmer) : max2;
//			System.err.println("\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+
//					"\n"+AminoAcid.kmerToString(kmer2, k2)+"\n"+AminoAcid.kmerToString(rkmer2, k2)+"\n"+AminoAcid.kmerToString(value, k)+"\n");
//			assert(kmer2==AminoAcid.reverseComplementBinaryFast(rkmer2, k2));)
		}else if(amino || translate){
			assert(amino || translate);
			final long kmer2=(kmer&k2mask);
			useK1=((kmer2%4999)&1L)==0L;
			key=useK1 ? kmer : kmer2;
		}else{//if for whatever reason you want forward kmers
			final long kmer2=((kmer&k2midmask)>>>k2shift);
			useK1=((kmer2%4999L)&1L)==0L;
			key=useK1 ? kmer : kmer2;
		}
//		if(key%5!=0){return -1;}//This makes it a little faster.
		final long bit=useK1 ? 0 : 1; //Note that this gets reversed later during the subtraction process
		
//		System.err.println(bit+", "+Long.toHexString(kmer)+", "+Long.toHexString(k2mask));
//		assert(bit==0) : bit+", "+Long.toHexString(kmer);
		
		long hashcode=key, data=key;
//		for(int i=0; i<codeMax; i+=codeIncrement){
//			int x=(int)(data&cycleMask);
//			data>>>=bitsPerCycle;
//			hashcode^=codes1D[i+x];
//		}
//		for(int i=0; data!=0; i+=codeIncrement){
//			int x=(int)(data&cycleMask);
//			data>>>=bitsPerCycle;
//			hashcode^=codes1D[i+x];
//		}
		{//5% faster than other loops.
			int i=0;
			do{
				int x=(int)(data&cycleMask);
				data>>>=bitsPerCycle;
				hashcode^=codes1D[i+x];
				i+=codeIncrement;
			}while(data!=0);
		}
		hashcode=(hashcode&(-2L))|bit;
		
		return hashcode;
	}
	
	/**
	 * Fastest version!
	 * Generates a hash code from a kmer.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash1(long kmer){
//		if(kmer%5!=0){return -1;}
		long code=kmer;
		for(int i=0; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		
		return code;
	}
	
	
	/**
	 * Generates a hash code from a kmer, using dual kmer lengths.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash2(final long kmer0){
//		return Tools.hash64shift(kmer0);
		long kmer=kmer0;
		long code=0;
		
		for(int i=0; i<codeMax2; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		final long code2=code^(k2mask&kmer0);
		for(int i=codeMax2; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
	}
	
	
	/**
	 * Generates a hash code from a kmer, using dual kmer lengths, allowing K2 to be a non-multiple-of-4.
	 * Uses one additional and, xor, and lookup.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash3(final long kmer0){
//		return Tools.hash64shift(kmer0);
		
		long kmer=kmer0;
		long code=0;
		assert(k2submask!=0);
		
		int i=0;
		for(; i<codeMax2; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		final int residual=(int)(kmer&k2submask);
		final long code2=code^(k2mask&kmer0)^codes1D[i+residual];
		for(; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------      Convenience Methods     ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final float align(byte[] query, byte[] ref){
		IDAligner ssa=GeneCaller.getSSA();
		return ssa.align(query, ref);
	}
	
	public static final float alignAndMakeMatch(Read r, byte[] ref){
		byte[] query=r.bases;
		SingleStateAlignerFlat2 ssa=GeneCaller.getSSA();
		int a=0, b=ref.length-1;
		int[] max=ssa.fillUnlimited(query, ref, a, b, -9999);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		byte[] match=ssa.traceback(query, ref, a, b, rows, maxCol, 0);
		int[] score=ssa.score(query, ref, a, b, rows, maxCol, 0);
		final float id=Read.identity(match);
		r.match=match;
		r.start=score[1];
		r.stop=score[2];
		r.setMapped(true);
		return id;
	}
	
	public static final float alignAndMakeMatch(Read r, byte[] ref, float[] refWeights/*, float[] insWeights, float[] delWeights*/){
		byte[] query=r.bases;
		SingleStateAlignerFlatFloat ssa=GeneCaller.getSSAF();
//		ssa.setWeights(refWeights, insWeights, delWeights);
		ssa.setWeights(refWeights);
		int a=0, b=ref.length-1;
		int[] max=ssa.fillUnlimited(query, ref, a, b, -9999);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		byte[] match=ssa.traceback(query, ref, a, b, rows, maxCol, 0);
		int[] score=ssa.score(query, ref, a, b, rows, maxCol, 0);
		final float id=Read.identity(match);
		r.match=match;
		r.start=score[1];
		r.stop=score[2];
		r.setMapped(true);
		return id;
	}
	
	static String fixMeta(String s){
		if(s==null){return null;}
		int colon=s.indexOf(':');
		assert(colon>=0);
		if(s.length()==colon+5 && s.endsWith(":null")){return null;}
		return s.replace('\t', ' ');
	}
	
	static ArrayList<String> fixMeta(ArrayList<String> list){
		if(list==null || list.isEmpty()){return null;}
		for(int i=0; i<list.size(); i++){
			String s=list.get(i);
			s=fixMeta(s);
			if(s==null){
				list.remove(i);
				i--;
			}else{
				list.set(i, s);
			}
		}
		if(list.isEmpty()){return null;}
		list.trimToSize();
		Collections.sort(list);
		return list;
	}
	
	public static final float aniToWkid(final double ani){
		final float wkid;
		if(ani<=0){
			wkid=0;
		}else if(k2<1){
			wkid=(float)Math.pow(ani, k);
		}else{
			wkid=0.5f*(float)(Math.pow(ani, k)+Math.pow(ani, k2));
		}
		return wkid;
	}
	
	public static final float wkidToAniExact(final double wkid, final int A, final int B, final double maxDeviation){
		assert(A>B);
		assert(maxDeviation>0);
		final double logWkid=Math.log(wkid);
		final double invA=1.0/A;
		final double invB=1.0/B;
		
		double aniUpper=Math.exp(logWkid*invA);
		double aniLower=Math.exp(logWkid*invB);
		assert(aniLower<=aniUpper);
		double aniEst=(aniLower+aniUpper)*0.5f;
//		System.err.println(aniEst);
//		if(B>16){//Fast mode
//			aniUpper=(aniUpper+3*aniEst)*0.25;
//			aniLower=(aniLower+3*aniEst)*0.25;
//		}
		double wkidEst=aniToWkid(aniEst, A, B);
		int iters=1;
//		System.err.println(aniLower+"\t"+aniUpper+"\t"+aniEst+"\t"+aniToWkid(aniEst, A, B)+"\t"+(wkidEst-wkid));
		while(Math.abs(wkidEst-wkid)>maxDeviation && iters<40){
			if(wkidEst<wkid){aniLower=aniEst;}
			else{aniUpper=aniEst;}
			aniEst=(aniLower+aniUpper)*0.5f;
			wkidEst=aniToWkid(aniEst, A, B);
//			System.err.println(aniLower+"\t"+aniUpper+"\t"+aniEst+"\t"+wkidEst+"\t"+(wkidEst-wkid));
			iters++;
		}
//		System.err.println("Iterations: "+iters);
		return (float)aniEst;
	}
	
	public static double aniToWkid(double ani, int A, int B){
		if(A<B){int C=A; A=B; B=C;}//Swap
		double aPow=Math.pow(ani, A);
		double bPow=Math.pow(ani, B);
//		return 0.5*(aPow+bPow);
		return useToValue2 ? 0.5*(aPow+bPow) : 0.5*(aPow+(bPow*0.4+aPow*0.6));
	}
	
	//From Kayla
//	public static double aniToWkid(double ANI, int K2, int K1){
//		if(K2<K1){int C=K2; K2=K1; K1=C;}//Swap
//		
//		return 0.5*((1 - (1 - Math.pow(ANI, K2-K1))*Math.pow(ANI, K1))*Math.pow(ANI, K2) +
//                (1 + (1 - Math.pow(ANI, K2-K1))*Math.pow(ANI, K1))*Math.pow(ANI, K1));
//	}
	
	public static double aniToWkid(double ani, int A){
		return Math.pow(ani, A);
	}
	
	public static final float wkidToAniExact(final double wkid, final int k){
		return (float)Math.exp(Math.log(wkid)/(k));
	}
	
	public static final float wkidToAni(final double wkid){
		final float ani;
		if(wkid<=0){
			ani=0;
		}else if(k2<1){
			ani=(float)Math.exp(Math.log(wkid)/k);
		}else{
			if(EXACT_ANI){return wkidToAniExact(wkid, k, k2, 0.00005);}
			
			//This is linear interpolation which is not correct.  But for some reason it works really well.
			final double log=Math.log(wkid);
			
//			double ani1=Math.exp(log/k);
//			double ani2=Math.exp(log/k2);
//			ani=(float)(0.5*(ani1+ani2));
			
			//alternatively...  this one seems to work better.
//			ani=(float)Math.exp(2*log/(k*1.1+k2*0.9));
			ani=(float)Math.exp(2*log/(1.2*k+.8*k2));
		}
		return ani;
	}

//	public static void kmerArrayToSketchArray(long[] kmers){
//		for(int i=0; i<kmers.length; i++){
//			long kmer=kmers[i];
//			long hash=SketchTool.hash(kmer);
//			assert(hash>=minHashValue);
//			hash=Long.MAX_VALUE-hash;
//			kmers[i]=hash;
//		}
//		Shared.sort(kmers);
//	}
	
	public static void hashArrayToSketchArray(long[] keys){
		for(int i=0; i<keys.length; i++){
			long hash=keys[i];
			assert(hash>=minHashValue);
			hash=Long.MAX_VALUE-hash;
			keys[i]=hash;
		}
		Shared.sort(keys);
	}
	
	public static final long genomeSizeEstimate(long min, int length) {
		if(length==0){return 0;}
		double est=((double)Long.MAX_VALUE)*2*length/(Tools.max(min, 1));
//		if(k2>0){est*=0.5;} //This is necessary if the hash function uses max, but not random.
//		System.err.println("max="+Long.MAX_VALUE+", min="+min+", len="+length+", est="+(long)est);
//		new Exception().printStackTrace(System.err);
		return (long)Math.ceil(est);
	}
	
	public static final int toSketchSize(long genomeSizeBases, long genomeSizeKmers, long genomeSizeEstimate, int maxSketchSize){
//		assert(false) : genomeSizeBases+", "+genomeSizeKmers+", "+genomeSizeEstimate+", "+maxSketchSize+", "+useSizeEstimate;
		if(genomeSizeEstimate>0 && useSizeEstimate){
			return toSketchSizeKmers(genomeSizeEstimate, maxSketchSize);
		}
		if(genomeSizeKmers>0){
			return toSketchSizeKmers(genomeSizeKmers, maxSketchSize);
		}
		assert(genomeSizeBases>0) : "BBSketch does not currently support empty files.\n"
				+genomeSizeBases+", "+genomeSizeKmers+", "+genomeSizeEstimate+", "+maxSketchSize+", "+useSizeEstimate;
		return toSketchSizeBases(genomeSizeBases, maxSketchSize);
	}
	
	private static final int toSketchSizeBases(long genomeSizeBases, int maxSketchSize) {
		return toSketchSizeKmers(Tools.max(0, genomeSizeBases-k+1), maxSketchSize);
	}
	
	private static final int toSketchSizeKmers(long genomeSizeKmers, int maxSketchSize) {
//		System.err.println(genomeSizeKmers+", "+maxSketchSize);
//		new Exception().printStackTrace();
		if(AUTOSIZE){
			if(aminoOrTranslate()){
				float linear1=Tools.min(60+0.5f*(float)Math.sqrt(genomeSizeKmers),
						maxGenomeFractionSmall*genomeSizeKmers*2);
				float linear2=genomeSizeKmers*maxGenomeFraction*2;
				float poly=0+1f*(float)Math.sqrt(genomeSizeKmers)+90f*(float)Math.pow(genomeSizeKmers, 0.3);
				float log=Tools.max(1000, -4000+3500*(float)Math.log(Tools.max(1, genomeSizeKmers)));
				float min=Tools.min(Tools.max(linear1, linear2), poly, log);
				assert(min>=Integer.MIN_VALUE && min<=Integer.MAX_VALUE) : min;

				//			System.err.println(genomeSizeKmers+" -> "+linear1+", "+linear2+", "+poly+", "+log);

				return (int)Tools.min(genomeSizeKmers*keyFraction2, min*AUTOSIZE_FACTOR);
				//			return (int)Tools.max(1, min*AUTOSIZE_FACTOR);
			}else{
				float linear1=Tools.min(smallSketchSize+0.5f*(float)Math.sqrt(genomeSizeKmers),
						maxGenomeFractionSmall*genomeSizeKmers-10);
				float linear2=genomeSizeKmers*maxGenomeFraction;
				float poly=0+1f*(float)Math.sqrt(genomeSizeKmers)+90f*(float)Math.pow(genomeSizeKmers, 0.3);
				float log=Tools.max(1000, -4000+3500*(float)Math.log(Tools.max(1, genomeSizeKmers))+8f*(float)Math.pow(genomeSizeKmers, 0.3));
				float min=Tools.min(Tools.max(linear1, linear2), poly, log);
				assert(min>=Integer.MIN_VALUE && min<=Integer.MAX_VALUE) : min;

				//			System.err.println(genomeSizeKmers+" -> "+linear1+", "+linear2+", "+poly+", "+log);
				
				int result=(int)Tools.min(genomeSizeKmers*keyFraction2, min*AUTOSIZE_FACTOR);
//				System.err.println(result);
				return result;
				//			return (int)Tools.max(1, min*AUTOSIZE_FACTOR);
			}
		}else if(AUTOSIZE_LINEAR){
			if(aminoOrTranslate()){
				return (int)(Tools.min(10000000, Tools.min(maxGenomeFractionSmall, AUTOSIZE_LINEAR_DENSITY)*genomeSizeKmers));
			}else{
				return (int)(Tools.min(10000000, Tools.min(maxGenomeFractionSmall, AUTOSIZE_LINEAR_DENSITY)*genomeSizeKmers));
			}
		}else{
			if(aminoOrTranslate()){//kmers are roughly 3x smaller...
				return Tools.min((int)(2+3*maxGenomeFraction*genomeSizeKmers), maxSketchSize);
			}else{
				return Tools.min((int)(2+maxGenomeFraction*genomeSizeKmers), maxSketchSize);
			}
		}
	}
	
//	static final long toValue2(long kmer, long rkmer){
//		assert(k2>0 && k2<k) : "k="+k+","+k2;
//		
//		final long value;
//		if(rcomp){
//			assert(!amino && !translate);
////			if(!rcomp){return kmer;}
//			long kmer2=((kmer&k2midmask)>>>k2shift);
//			long rkmer2=((rkmer&k2midmask)>>>k2shift);
//			
////			assert(kmer2>=0);
////			assert(kmer<0 || kmer2<=kmer);
////			assert(rkmer<0 || rkmer2<=kmer);
//			
////			assert(false) : "\n"+Long.toBinaryString(kmer)+
////				"\n"+Long.toBinaryString(rkmer)+
////				"\n"+Long.toBinaryString(kmer2)+
////				"\n"+Long.toBinaryString(rkmer2);
//			
//			//TODO: Slow
////			assert(kmer==AminoAcid.reverseComplementBinaryFast(rkmer, k)) : 
////				"\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+"\n"+AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(rkmer, k), k)+"\n";
////			assert(kmer2==AminoAcid.reverseComplementBinaryFast(rkmer2, k2));
//			long max2=Tools.max(kmer2, rkmer2);
////			long cmasked=max2&kmerChoiceMask;
////			boolean longer=kmerChoiceBitset.get((int)cmasked);
//			boolean longer=((max2%4999L)&1L)==0L;
//			value=longer ? Tools.max(kmer, rkmer) : max2;
////			System.err.println("\n"+longer+", "+max2+", "+(max2%4999L)+", "+((max2%4999L)&1L)+
////					", "+Long.toHexString(kmer)+", "+Long.toHexString(max2)+", "+Long.toHexString(value));
////			long value=longer ? (kmer2>rkmer2 ? kmer : rkmer) : max2;
////			System.err.println("\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+
////					"\n"+AminoAcid.kmerToString(kmer2, k2)+"\n"+AminoAcid.kmerToString(rkmer2, k2)+"\n"+AminoAcid.kmerToString(value, k)+"\n");
////			assert(kmer2==AminoAcid.reverseComplementBinaryFast(rkmer2, k2));)
//		}else{
//			assert(amino || translate);
////			long kmer2=(kmer&k2mask);
//			long kmer2=((kmer&k2midmask)>>>k2shift);
//			boolean longer=((kmer2%4999)&1)==0;
//			value=longer ? kmer : kmer2;
//		}
//		return value;
//	}
	
	private static final long toValue(long kmer, long rkmer){
//		assert(toValue2(kmer, rkmer)==toValue2(rkmer, kmer));
		assert(!useToValue2);
//		if(useToValue2){return toValue2(kmer, rkmer);}
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Plaintext decimal */
	public static final int RAW=0;
	/** Hexadecimal (4-bit) */
	public static final int HEX=1;
	/** ASCII-48 6-bit encoding */
	public static final int A48=2;
	
	public static final char[] codingArray={'R', 'H', 'A'};
	
	static final byte[] hexTable=new byte[128];
	static {
		Arrays.fill(hexTable, (byte)-1);
		for(int i='0'; i<='9'; i++){
			hexTable[i]=(byte)(i-'0');
		}
		for(int i='A'; i<='F'; i++){
			hexTable[i]=hexTable[i+'a'-'A']=(byte)(i-'A'+10);
		}
		hexTable['x']=hexTable['X']=hexTable['-']=hexTable['+']=0;
	}

	public static final int ONE_SKETCH=1, PER_SEQUENCE=2, PER_TAXA=3, PER_IMG=4/*, PER_HEADER=5*/, PER_FILE=6;
	
	/*--------------------------------------------------------------*/
	/*----------------      Default Locations       ----------------*/
	/*--------------------------------------------------------------*/

//	public static ArrayList<String> IMG_PATH(){return makePaths(IMG_PATH, 31);}
//	public static ArrayList<String> NT_PATH(){return makePaths(NT_PATH, 31);}
//	public static ArrayList<String> NR_PATH(){throw new RuntimeException("NR is not currently available.");}
//	public static ArrayList<String> REFSEQ_PATH(){return makePaths(REFSEQ_PATH, 31);}
//	public static ArrayList<String> SILVA_PATH(){return makePaths(SILVA_PATH, 31);}
	
	private static ArrayList<String> makePaths(String pattern, int files){
		ArrayList<String> list=new ArrayList<String>(files);
		for(int i=0; i<files; i++){
			list.add(pattern.replace("#", ""+i));
		}
		return list;
	}
	
	private static final String IMG_PATH="/global/cfs/cdirs/bbtools/img/current/img#.sketch";
	private static final String NT_PATH="/global/cfs/cdirs/bbtools/nt/current/taxa#.sketch";
	private static final String NR_PATH="/global/cfs/cdirs/bbtools/nr/current/taxa#.sketch";
	private static final String REFSEQ_PATH="/global/cfs/cdirs/bbtools/refseq/current/taxa#.sketch";
	private static final String REFSEQ_PATH_BIG="/global/cfs/cdirs/bbtools/refseq/current/big#.sketch";
	private static final String SILVA_PATH="/global/cfs/cdirs/bbtools/silva/latest/both_taxa#.sketch";
	private static final String PROKPROT_PATH="/global/cfs/cdirs/bbtools/refseq/current/prot/taxa#.sketch";
	private static final String PROKPROT_PATH_BIG="/global/cfs/cdirs/bbtools/refseq/current/prot/big#.sketch";
	private static final String MITO_PATH="/global/cfs/cdirs/bbtools/mito2/taxa#.sketch";
	private static final String FUNGI_PATH="/global/cfs/cdirs/bbtools/mito2/fungi#.sketch";

	private static final String IMG_PATH_IGBVM="/data/sketch/img/current/img#.sketch";
	private static final String NT_PATH_IGBVM="/data/sketch/nt/current/taxa#.sketch";
	private static final String NR_PATH_IGBVM="/data/sketch/nr/current/taxa#.sketch";
	private static final String REFSEQ_PATH_IGBVM="/data/sketch/refseq/current/taxa#.sketch";
	private static final String REFSEQ_PATH_BIG_IGBVM="/data/sketch/refseq/current/big#.sketch";
	private static final String SILVA_PATH_IGBVM="/data/sketch/silva/current/both_taxa#.sketch";
	private static final String PROKPROT_PATH_IGBVM="/data/sketch/refseq/current/prot/taxa#.sketch";
	private static final String PROKPROT_PATH_BIG_IGBVM="/data/sketch/refseq/current/prot/big#.sketch";
	private static final String MITO_PATH_IGBVM="/data/sketch/mito2/taxa#.sketch";
	private static final String FUNGI_PATH_IGBVM="/data/sketch/mito2/fungi#.sketch";

	private static final String IMG_PATH_AWS=null;
	private static final String NT_PATH_AWS="/test1/sketch/latest/nt/taxa#.sketch";
	private static final String NR_PATH_AWS=null;
	private static final String REFSEQ_PATH_AWS="/test1/sketch/latest/refseq/taxa#.sketch";
	private static final String SILVA_PATH_AWS="/test1/sketch/latest/ribo/both_taxa#.sketch";
	private static final String PROKPROT_PATH_AWS="/test1/sketch/latest/protein/taxa#.sketch";
	private static final String MITO_PATH_AWS=null;
	private static final String FUNGI_PATH_AWS=null;

	public static final String IMG_PATH(){return Shared.IGBVM ? IMG_PATH_IGBVM : Shared.AWS ? IMG_PATH_AWS : IMG_PATH;}
	public static final String NT_PATH(){return Shared.IGBVM ? NT_PATH_IGBVM : Shared.AWS ? NT_PATH_AWS : NT_PATH;}
	public static final String NR_PATH(){return Shared.IGBVM ? NR_PATH_IGBVM : Shared.AWS ? NR_PATH_AWS : NR_PATH;}
	public static final String REFSEQ_PATH(){return Shared.IGBVM ? REFSEQ_PATH_IGBVM : Shared.AWS ? REFSEQ_PATH_AWS : REFSEQ_PATH;}
	public static final String REFSEQ_PATH_BIG(){return Shared.IGBVM ? REFSEQ_PATH_BIG_IGBVM : REFSEQ_PATH_BIG;}
	public static final String SILVA_PATH(){return Shared.IGBVM ? SILVA_PATH_IGBVM : Shared.AWS ? SILVA_PATH_AWS : SILVA_PATH;}
	public static final String PROKPROT_PATH(){return Shared.IGBVM ? PROKPROT_PATH_IGBVM : Shared.AWS ? PROKPROT_PATH_AWS : PROKPROT_PATH;}
	public static final String PROKPROT_PATH_BIG(){return Shared.IGBVM ? PROKPROT_PATH_BIG_IGBVM : PROKPROT_PATH_BIG;}
	public static final String MITO_PATH(){return Shared.IGBVM ? MITO_PATH_IGBVM : Shared.AWS ? MITO_PATH_AWS : MITO_PATH;}
	public static final String FUNGI_PATH(){return Shared.IGBVM ? FUNGI_PATH_IGBVM : Shared.AWS ? FUNGI_PATH_AWS : FUNGI_PATH;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean useWhitelist() {return useWhitelist;}
	public static String blacklist() {return blacklist;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int CODING=A48;
	public static boolean deltaOut=true;
	public static boolean rcomp=true;

	public static final int defaultK=32;
	public static final int defaultK2=24;
	public static final int defaultKAmino=12;
	public static final int defaultK2Amino=9;
	public static int k=defaultK;
	public static int k2=defaultK2;
	public static int entropyK=3;
	public static boolean setK=false;
	public static boolean amino=false;
	public static boolean amino8=false;
	public static boolean translate=false;
	public static boolean sixframes=false;
	public static boolean processSSU=true;
	public static int min_SSU_len=1000;
	public static int HASH_VERSION=2;
	public static String pgmFile=null;
	public static GeneModel pgm=null;
	
	static boolean aminoOrTranslate(){return amino | translate;}
	
	private static int bitsPerCycle=8; //Optimal speed for K=31 is 9bpc (15% faster than 8). 9, 10, and 11 are similar.
	private static long cycleMask=~((-1L)<<bitsPerCycle);
	private static final long codeOrMask=1L;
	private static final long codeAndMask=~1L;
	private static int maxCycles=(64+bitsPerCycle-1)/bitsPerCycle;
	private static int codeIncrement=(int)(cycleMask+1);
	private static int codeMax;
	private static int codeMax2;
	
	private static long hashSeed=12345;
	private static long[][] codes;//=makeCodes(maxCycles, codeIncrement, hashSeed, false);
	private static long[] codes1D;//=makeCodes1D(codes);
	static boolean useToValue2=true;
	
	//These are set in postParse()
	private static int bitsPerBase=2;
	private static int basesPerCycle=bitsPerCycle/bitsPerBase;
	private static int hashCycles=64/bitsPerCycle; //Note: Needs to be variable based on k to make k2 codes compatible with k codes
	private static int hashCycles2;
	private static long k2mask;
	private static long k2submask;
	//Difference in length of k and k2, in bits
	private static int k2shift;
	private static long k2midmask;
	
	/** 
	 * Make the SketchHeap size this factor bigger,
	 * when minKeyOccuranceCount is used  
	 */
	public static double sketchHeapFactor=8;
//	static int minKeyOccuranceCount=1;
	
	public static int minSketchSize=3;
	public static int targetSketchSize=10000;
	public static boolean AUTOSIZE=true;
	public static float AUTOSIZE_FACTOR=1;
	public static boolean SET_AUTOSIZE_FACTOR=false;
	public static boolean SET_AUTOSIZE=false;
	public static boolean SET_TARGET_SIZE=false;
	public static boolean AUTOSIZE_LINEAR=false;
	public static double AUTOSIZE_LINEAR_DENSITY=0.001;
	public static float maxGenomeFraction=0.04f;//Was 0.015, but that is often too sparse
	public static float maxGenomeFractionSmall=0.10f;
	public static int smallSketchSize=150;
	public static boolean makeIndex=true;
	public static float prealloc=0;
	public static boolean allToAll=false;
	public static boolean compareSelf=false;
	public static boolean skipCompare=false;
	public static final int bitSetBits=2; //Needs to be 2 for unique counts.

	private static double keyFraction=0.16;
	private static double keyFraction2=keyFraction*1.2;
	public static long minHashValue=setMinHashValue();
	public static double keyFraction(){return keyFraction;}
	public static void setKeyFraction(double d){
		assert(d>0);
		keyFraction=Tools.mid(0.0001, d, 0.5);
		keyFraction2=Tools.mid(0.0001, keyFraction*1.2, 0.5);
	}
	public static long setMinHashValue(){
		double mult=1-2*keyFraction;
		minHashValue=(long)(mult*Long.MAX_VALUE);
		assert(minHashValue>=0 && minHashValue<Long.MAX_VALUE);
		return minHashValue;
	}
	
	public static int minFakeID=1900000000;
	static boolean hashNames=false;
	static boolean skipNonCanonical=true;
	static boolean useSizeEstimate=true;
	public static boolean allowMultithreadedFastq=false;
	static boolean forceDisableMultithreadedFastq=false;
	
	static boolean preferSSUMap=false;
	static boolean preferSSUMapForEuks=true;
	static boolean useSSUMapOnly=false;
	static boolean useSSUMapOnlyForEuks=false;
	static boolean ban18SForProks=true;
	static boolean ban16SForEuks=true;
	
	public static long sampleseed=-1L;
	
	public static TaxTree taxtree=null;
	private static String treefile=null;
	static String blacklist=null;
	static boolean useWhitelist=false;
	private static boolean postparsed=false;
	public static boolean KILL_OK=false;
	public static boolean EXACT_ANI=true;
	public static boolean useSSA=true;
	public static boolean useSSA3=false;
	
	//Needs to be last due to dependencies.
	public static DisplayParams defaultParams=new DisplayParams();
	
	static AlignmentThreadPool alignerPool=null;
	
	public static boolean verbose2=false;
	public static boolean LOADER2=true;
	
}
