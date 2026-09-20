package align2;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import bloom.BloomFilter;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.CoveragePileup;
import shared.Shared;
import parse.Parse;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.ReadStreamWriter;
import stream.SamLine;
import structures.ListNum;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Main BBMap alignment engine for mapping short reads to reference genomes.
 * Provides high-speed k-mer based alignment with configurable sensitivity modes.
 * Supports single and paired-end reads with comprehensive output format options.
 *
 * @author Brian Bushnell
 * @date Dec 22, 2012
 */
public final class BBMapS extends AbstractMapper implements Accumulator<BBMapS.ProcessThread> {
	private int hybridTipSearchCeiling; // Assigned during superclass construction; no field initializer.
	

	/**
	 * Program entry point for BBMap alignment.
	 * Initializes mapper, loads index, processes ambiguous mappings, and executes alignment.
	 * @param args Command-line arguments for alignment configuration
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		BBMapS mapper=new BBMapS(args);
		args=Tools.condenseStrict(args);
		if(!INDEX_LOADED){mapper.loadIndex();}
		if(Data.scaffoldPrefixes){mapper.processAmbig2();}
		mapper.testSpeed(args);
		ReadWrite.waitForWritingToFinish();
		t.stop();
		outstream.println("\nTotal time:     \t"+t);
		clearStatics();
		BBMapSplitterS.clearStatics();
	}
	
	/**
	 * Constructs BBMap instance with specified arguments.
	 * Inherits configuration parsing and validation from AbstractMapper.
	 * @param args Command-line arguments for mapper configuration
	 */
	public BBMapS(String[] args){
		super(args);
	}
	
	/**
	 * Sets BBMap-specific default values for alignment parameters.
	 * Configures compression, key density, alignment scoring, and output settings.
	 * Called during initialization to establish baseline configuration.
	 */
	@Override
	public void setDefaults(){
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_BGZIP=ReadWrite.USE_UNBGZIP=true;
		ReadWrite.PREFER_BGZIP=true;
		ReadWrite.ZIPLEVEL=2;
		MAKE_MATCH_STRING=true;
		keylen=13;
		
		MINIMUM_ALIGNMENT_SCORE_RATIO=0.56f;

		keyDensity=1.9f;//2.3f;
		maxKeyDensity=3f;//4f;
		minKeyDensity=1.5f;//1.8f;
		maxDesiredKeys=15;
		
		SLOW_ALIGN_PADDING=4;
		SLOW_RESCUE_PADDING=4+SLOW_ALIGN_PADDING;
		TIP_SEARCH_DIST=100;
		
		MSA_TYPE="MultiStateAligner11ts";
		MAX_SITESCORES_TO_PRINT=5;
		PRINT_SECONDARY_ALIGNMENTS=false;
		AbstractIndex.MIN_APPROX_HITS_TO_KEEP=1;

		// Keep normal multi-scaffold packing unchanged, but allow one wheat-scale
		// scaffold to occupy its own chromosome.  The limit stays below both the
		// Java byte[] ceiling and BBIndex's largest v4 site field (chrombits=1).
		FastaToChromArrays2.MAX_SINGLE_SCAFFOLD=Shared.MAX_ARRAY_LEN/2-200000;
	}
	
	/**
	 * Pre-processes arguments to apply speed/accuracy mode presets.
	 * Modifies key density, alignment strictness, and index parameters based on
	 * fast, slow, or vslow mode selection before main argument parsing.
	 *
	 * @param args Original command-line arguments
	 * @return Modified argument array with mode-specific parameters added
	 */
	@Override
	public String[] preparse(String[] args){
		boolean quantumOnlyRequested=false;
		boolean pseudoAlignRequested=false;
		boolean quantumHybridRequested=false;
		for(int i=0; i<args.length; i++){
			final String s=args[i];
			if(s==null){continue;}
			final int equals=s.indexOf('=');
			final String key=(equals<0 ? s : s.substring(0, equals));
			if(key.equalsIgnoreCase("quantumonebase")){
				final String value=(equals<0 ? null : s.substring(equals+1));
				System.setProperty("bbmap3.quantumTieredMutate",
						Boolean.toString(Parse.parseBoolean(value)));
				args[i]=null;
			}else if(key.equalsIgnoreCase("quantumonly")){
				final String value=(equals<0 ? null : s.substring(equals+1));
				quantumOnlyRequested=Parse.parseBoolean(value);
				System.setProperty("bbmap3.quantumOnly",
						Boolean.toString(quantumOnlyRequested));
				args[i]=null;
			}else if(key.equalsIgnoreCase("pseudoalign")){
				final String value=(equals<0 ? null : s.substring(equals+1));
				pseudoAlignRequested=Parse.parseBoolean(value);
				System.setProperty("bbmap3.pseudoAlign",
						Boolean.toString(pseudoAlignRequested));
				args[i]=null;
			}else if(key.equalsIgnoreCase("quantumhybrid")){
				final String value=(equals<0 ? null : s.substring(equals+1));
				quantumHybridRequested=Parse.parseBoolean(value);
				System.setProperty("bbmap3.quantumHybrid",
						Boolean.toString(quantumHybridRequested));
				args[i]=null;
			}
		}
		if(quantumHybridRequested && (quantumOnlyRequested || pseudoAlignRequested ||
				Boolean.getBoolean("bbmap3.quantumTieredMutate"))){
			throw new IllegalArgumentException("quantumhybrid=t is incompatible with quantumonly, pseudoalign, and quantumonebase");
		}
		BBIndex.EMIT_POLYCRYSTALLINE_MATCH=pseudoAlignRequested;
		if(quantumOnlyRequested || pseudoAlignRequested){
			final ArrayList<String> list=new ArrayList<String>();
			list.add("rescue=f");
			list.add("local=f");
			list.add("secondary=f");
			if(pseudoAlignRequested){list.add("match=f");}
			for(String s : args){if(s!=null){list.add(s);}}
			args=list.toArray(new String[list.size()]);
		}
		if(fast){
			ArrayList<String> list=new ArrayList<String>();
			list.add("tipsearch="+TIP_SEARCH_DIST/5);
			list.add("maxindel=80");
			list.add("minhits=2");
			list.add("bwr=0.18");
			list.add("bw=40");
			list.add("minratio=0.65");
			list.add("midpad=150");
			list.add("minscaf=50");
			list.add("quickmatch=t");
			list.add("rescuemismatches=15");
			list.add("rescuedist=800");
			list.add("maxsites=3");
			list.add("maxsites2=100");
//			list.add("k=13");
			
			//TODO:  Make these adjustable.
//			MIN_TRIM_SITES_TO_RETAIN_SINGLE
//			MIN_TRIM_SITES_TO_RETAIN_PAIRED
//			MAX_TRIM_SITES_TO_RETAIN
			//TODO:  Make trimLists adjustable via an offset or multiplier
			
			BBIndex.setFractionToExclude(BBIndex.FRACTION_GENOME_TO_EXCLUDE*1.25f);
			
			for(String s : args){if(s!=null){list.add(s);}}
			args=list.toArray(new String[list.size()]);
			
			keyDensity*=0.9f;
			maxKeyDensity*=0.9f;
			minKeyDensity*=0.9f;
		}else if(vslow){
			ArrayList<String> list=new ArrayList<String>();
			list.add("tipsearch="+(TIP_SEARCH_DIST*3)/2);
			list.add("minhits=1");
			list.add("minratio=0.22");
			list.add("usequality=f");
			list.add("rescuemismatches=50");
			list.add("rescuedist=2500");
			list.add("maxindel=100");
			
			BBIndex.setFractionToExclude(0);
			
			for(String s : args){if(s!=null){list.add(s);}}
			args=list.toArray(new String[list.size()]);
			
			SLOW_ALIGN_PADDING=SLOW_ALIGN_PADDING*2+8;
			SLOW_RESCUE_PADDING=SLOW_RESCUE_PADDING*2+2;

			AbstractIndex.SLOW=true;
			AbstractIndex.VSLOW=true;
			keyDensity*=2.5f;
			maxKeyDensity*=2.5f;
			minKeyDensity*=2.5f;
		}else if(slow){
			ArrayList<String> list=new ArrayList<String>();
			list.add("tipsearch="+(TIP_SEARCH_DIST*3)/2);
//			list.add("maxindel=80");
			list.add("minhits=1");
//			list.add("bwr=0.18");
//			list.add("bw=40");
			list.add("minratio=0.45");
//			list.add("midpad=150");
//			list.add("minscaf=50");
//			list.add("k=13");
			
			BBIndex.setFractionToExclude(BBIndex.FRACTION_GENOME_TO_EXCLUDE*0.4f);
			
			for(String s : args){if(s!=null){list.add(s);}}
			args=list.toArray(new String[list.size()]);
			
			AbstractIndex.SLOW=true;
			keyDensity*=1.2f;
			maxKeyDensity*=1.2f;
			minKeyDensity*=1.2f;
		}
		
		if(excludeFraction>=0){
			BBIndex.setFractionToExclude(excludeFraction);
		}
		return args;
	}
	
	/**
	 * Post-processes parsed arguments to finalize configuration.
	 * Applies bandwidth constraints, handles input file detection,
	 * configures ambiguous read handling, and validates parameter combinations.
	 * @param args Parsed command-line arguments
	 */
	@Override
	void postparse(String[] args){
		hybridTipSearchCeiling=TIP_SEARCH_DIST; // Preserve parsed user/preset ceiling before startup indel clamp.
		if((hybridPair || hybridMaxIndel) && !explicitIndelBoundSet){
			// Force the validated low-pass bound through the existing install/clamp logic below,
			// rather than assigning BBIndex.MAX_INDEL directly, so TIP_SEARCH_DIST and any later
			// bandwidth clamp derive from it exactly as they would for a real user maxindel= flag.
			maxIndel1=50; maxIndel2=100;
		}

		if(MSA.bandwidthRatio>0 && MSA.bandwidthRatio<.2){
			SLOW_ALIGN_PADDING=Tools.min(SLOW_ALIGN_PADDING, 3);
			SLOW_RESCUE_PADDING=Tools.min(SLOW_RESCUE_PADDING, 6);
		}
		
		if(maxIndel1>-1){
			TIP_SEARCH_DIST=Tools.min(TIP_SEARCH_DIST, maxIndel1);
			BBIndex.MAX_INDEL=maxIndel1;
		}
		if(maxIndel2>-1){
			BBIndex.MAX_INDEL2=maxIndel2;
		}
		
		if(minApproxHits>-1){
			BBIndex.MIN_APPROX_HITS_TO_KEEP=minApproxHits;
		}
		
		if(expectedSites>-1){
			BBMapThread.setExpectedSites(expectedSites);
			outstream.println("Set EXPECTED_SITES to "+expectedSites);
		}
		
		if(fractionGenomeToExclude>=0){
			BBIndex.setFractionToExclude(fractionGenomeToExclude);
		}
		
		{
			final String a=(args.length>0 ? args[0] : null);
			final String b=(args.length>1 ? args[1] : null);
			if(in1==null && a!=null && a.indexOf('=')<0 && (a.startsWith("stdin") || new File(a).exists())){in1=a;}
			if(in2==null && b!=null && b.indexOf('=')<0 && new File(b).exists()){in2=b;}
			if(ERROR_ON_NO_OUTPUT && !OUTPUT_READS && in1!=null){throw new RuntimeException("Error: no output file, and ERROR_ON_NO_OUTPUT="+ERROR_ON_NO_OUTPUT);}
		}

		assert(synthReadlen<BBMapThread.ALIGN_ROWS);
		
		if(MSA.bandwidth>0){
			int halfwidth=MSA.bandwidth/2;
			TIP_SEARCH_DIST=Tools.min(TIP_SEARCH_DIST, halfwidth/2);
			hybridTipSearchCeiling=Tools.min(hybridTipSearchCeiling, halfwidth/2);
			BBIndex.MAX_INDEL=Tools.min(BBIndex.MAX_INDEL, halfwidth/2);
			BBIndex.MAX_INDEL2=Tools.min(BBIndex.MAX_INDEL2, halfwidth);
			SLOW_ALIGN_PADDING=Tools.min(SLOW_ALIGN_PADDING, halfwidth/4);
			SLOW_RESCUE_PADDING=Tools.min(SLOW_RESCUE_PADDING, halfwidth/4);
		}
		
		if(PRINT_SECONDARY_ALIGNMENTS){
			REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
			BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
		}
		
		if(in1!=null){
			if(ambigMode==AMBIG_BEST){
				REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
				if(!PRINT_SECONDARY_ALIGNMENTS){BBIndex.QUIT_AFTER_TWO_PERFECTS=true;}
				outstream.println("Retaining first best site only for ambiguous mappings.");
			}else if(ambigMode==AMBIG_ALL){
				PRINT_SECONDARY_ALIGNMENTS=ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS=true;
				REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
				BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
				SamLine.MAKE_NH_TAG=true;
				ambiguousAll=true;
				outstream.println("Retaining all best sites for ambiguous mappings.");
			}else if(ambigMode==AMBIG_RANDOM){
				REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
				BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
				ambiguousRandom=true;
				outstream.println("Choosing a site randomly for ambiguous mappings.");
			}else if(ambigMode==AMBIG_TOSS){
				REMOVE_DUPLICATE_BEST_ALIGNMENTS=true;
				BBIndex.QUIT_AFTER_TWO_PERFECTS=true;
				outstream.println("Ambiguously mapped reads will be considered unmapped.");
			}else{
				throw new RuntimeException("Unknown ambiguous mapping mode: "+ambigMode);
			}
		}

		if(hybridPair){
			// Validate the FINAL effective state, after bandwidth clamps and ambigMode derivation
			// above -- not an early snapshot -- so a bandwidth clamp that reduces the forced/explicit
			// 50/100 below that, or ambig=all deriving PRINT_SECONDARY_ALIGNMENTS, are both caught.
			if(BBIndex.MAX_INDEL!=50 || BBIndex.MAX_INDEL2!=100){
				throw new RuntimeException("hybridpair requires effective maxindel=50/100 after all "
					+"derivations (explicit setting or bandwidth clamp may have changed this); got "
					+BBIndex.MAX_INDEL+"/"+BBIndex.MAX_INDEL2);
			}
			if(QUICK_MATCH_STRINGS){throw new RuntimeException("hybridpair is incompatible with quickmatch=t");}
			if(STRICT_MAX_INDEL){throw new RuntimeException("hybridpair is incompatible with strictmaxindel");}
			if(PRINT_SECONDARY_ALIGNMENTS){throw new RuntimeException("hybridpair is incompatible with secondary=t");}
			if(PERFECTMODE || SEMIPERFECTMODE){throw new RuntimeException("hybridpair is incompatible with perfectmode/semiperfectmode");}
		}
		if(hybridMaxIndel){
			if(hybridPair){throw new IllegalArgumentException("hybridmaxindel and hybridpair are separate modes");}
			if(QUICK_MATCH_STRINGS){throw new IllegalArgumentException("hybridmaxindel is incompatible with quickmatch=t");}
			if(STRICT_MAX_INDEL){throw new IllegalArgumentException("hybridmaxindel is incompatible with strictmaxindel");}
			if(PRINT_SECONDARY_ALIGNMENTS){throw new IllegalArgumentException("hybridmaxindel is incompatible with secondary=t");}
			if(PERFECTMODE || SEMIPERFECTMODE){throw new IllegalArgumentException("hybridmaxindel is incompatible with perfectmode/semiperfectmode");}
			if(Boolean.getBoolean("bbmap3.quantumOnly") || Boolean.getBoolean("bbmap3.quantumHybrid") ||
					Boolean.getBoolean("bbmap3.quantumTieredMutate") || Boolean.getBoolean("bbmap3.pseudoAlign")){
				throw new IllegalArgumentException("hybridmaxindel is incompatible with Quantum and pseudoalignment modes");
			}
			hybridMaxIndelConfig=new HybridMaxIndelConfig(BBIndex.MAX_INDEL,BBIndex.MAX_INDEL2,
				retryMaxIndel1,retryMaxIndel2,retryMinMapq);
		}
		if(Boolean.getBoolean("bbmap3.quantumOnly")){
			if(Boolean.getBoolean("bbmap3.quantumTieredMutate")){
				throw new RuntimeException("quantumonly=t is incompatible with quantumonebase=t");
			}
			if(hybridPair){throw new RuntimeException("quantumonly=t is incompatible with hybridpair=t");}
			if(RESCUE){throw new RuntimeException("quantumonly=t requires rescue=f");}
			if(LOCAL_ALIGN){throw new RuntimeException("quantumonly=t requires local=f");}
			if(PRINT_SECONDARY_ALIGNMENTS){
				throw new RuntimeException("quantumonly=t requires secondary=f and ambig!=all");
			}
		}
		if(Boolean.getBoolean("bbmap3.quantumHybrid")){
			if(hybridPair){throw new RuntimeException("quantumhybrid=t is incompatible with hybridpair=t");}
			if(PERFECTMODE || SEMIPERFECTMODE){
				throw new RuntimeException("quantumhybrid=t is incompatible with perfectmode/semiperfectmode");
			}
		}
		if(Boolean.getBoolean("bbmap3.pseudoAlign")){
			if(Boolean.getBoolean("bbmap3.quantumOnly") ||
					Boolean.getBoolean("bbmap3.quantumTieredMutate")){
				throw new RuntimeException("pseudoalign=t is incompatible with Quantum modes");
			}
			if(hybridPair){throw new RuntimeException("pseudoalign=t is incompatible with hybridpair=t");}
			if(RESCUE){throw new RuntimeException("pseudoalign=t requires rescue=f");}
			if(LOCAL_ALIGN){throw new RuntimeException("pseudoalign=t requires local=f");}
			if(PRINT_SECONDARY_ALIGNMENTS){throw new RuntimeException("pseudoalign=t requires secondary=f");}
			if(MAKE_MATCH_STRING){throw new RuntimeException("pseudoalign=t requires match=f");}
			if(STRICT_MAX_INDEL){throw new RuntimeException("pseudoalign=t is incompatible with strictmaxindel");}
			if(PERFECTMODE || SEMIPERFECTMODE){
				throw new RuntimeException("pseudoalign=t is incompatible with perfectmode/semiperfectmode");
			}
			if(MIN_IDFILTER>0 || AbstractMapThread.SUBFILTER>=0 ||
					AbstractMapThread.DELFILTER>=0 || AbstractMapThread.INSFILTER>=0 ||
					AbstractMapThread.INDELFILTER>=0 || AbstractMapThread.DELLENFILTER>=0 ||
					AbstractMapThread.INSLENFILTER>=0 || AbstractMapThread.EDITFILTER>=0 ||
					AbstractMapThread.NFILTER>=0){
				throw new RuntimeException("pseudoalign=t is incompatible with identity/edit filters");
			}
		}

	}
	
	/**
	 * Performs pre-alignment setup and validation.
	 * Configures minimum identity thresholds, output streams, blacklists,
	 * and validates required parameters like build number and reference.
	 */
	@Override
	public void setup(){
		
		assert(!useRandomReads || maxReads>0 || (in1!=null && in1.equals("sequential"))) : "Please specify number of reads to use.";
		
		if(minid!=-1){
			MINIMUM_ALIGNMENT_SCORE_RATIO=MSA.minIdToMinRatio(minid, MSA_TYPE);
			outstream.println("Set MINIMUM_ALIGNMENT_SCORE_RATIO to "+Tools.format("%.3f",MINIMUM_ALIGNMENT_SCORE_RATIO));
		}
		
		if(!setxs){SamLine.MAKE_XS_TAG=(SamLine.INTRON_LIMIT<1000000000);}
		if(setxs && !setintron){SamLine.INTRON_LIMIT=10;}
		
		if(outFile==null && outFile2==null && outFileM==null && outFileM2==null && outFileU==null && outFileU2==null
				&& outFileB==null && outFileB2==null && splitterOutputs==null && BBMapSplitterS.streamTable==null){
			outstream.println("No output file.");
			OUTPUT_READS=false;
		}else{
			OUTPUT_READS=true;
			if(bamscript!=null){
				BBMapSplitterS.makeBamScript(bamscript, splitterOutputs, outFile, outFile2, outFileM, outFileM2, outFileU, outFileU2, outFileB, outFileB2);
			}
		}
//		assert(false) : bamscript+", "+BBMapSplitterS.streamTable+", "+OUTPUT_READS;
		
		
		
		FastaReadInputStream.MIN_READ_LEN=Tools.max(keylen+2, FastaReadInputStream.MIN_READ_LEN);
		assert(FastaReadInputStream.settingsOK());
		
		if(build<0){throw new RuntimeException("Must specify a build number, e.g. build=1");}
		else{Data.GENOME_BUILD=build;}
		
		if(blacklist!=null && blacklist.size()>0){
			Timer t=new Timer();
			t.start();
			for(String s : blacklist){
				Blacklist.addToBlacklist(s);
			}
			t.stop();
			outstream.println("Created blacklist:\t"+t);
			t.start();
		}
		
		if(ziplevel!=-1){ReadWrite.ZIPLEVEL=ziplevel;}
		if(reference!=null){RefToIndex.makeIndex(reference, build, outstream, keylen);}
	}
	

	/**
	 * Configures handling of reads that map to multiple references.
	 * Sets parameters for splitting, first-reference assignment, random selection,
	 * or discarding based on AMBIGUOUS2_MODE setting.
	 */
	@Override
	void processAmbig2(){
		assert(Data.scaffoldPrefixes) : "Only process this block if there are multiple references.";
		if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_SPLIT){
			REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
			BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
			outstream.println("Reads that map to multiple references will be written to special output streams.");
		}else if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_FIRST){
			REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
			BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
			outstream.println("Reads that map to multiple references will be written to the first reference's stream only.");
		}else if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_TOSS){
			BBIndex.QUIT_AFTER_TWO_PERFECTS=true;
			outstream.println("Reads that map to multiple references will be considered unmapped.");
		}else if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_RANDOM){
			REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
			BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
			outstream.println("Reads that map to multiple references will be written to a random stream.");
		}else if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_ALL){
			REMOVE_DUPLICATE_BEST_ALIGNMENTS=false;
			BBIndex.QUIT_AFTER_TWO_PERFECTS=false;
			outstream.println("Reads that map to multiple references will be written to all relevant output streams.");
		}else{
			BBSplitter.AMBIGUOUS2_MODE=BBSplitter.AMBIGUOUS2_FIRST;
		}
	}
	
	/**
	 * Loads reference genome index and prepares for alignment.
	 * Initializes chromosome data structures, generates k-mer index,
	 * applies genome size optimizations, and optionally creates Bloom filter.
	 * Configures coverage analysis structures if requested.
	 */
	@Override
	void loadIndex(){
		Timer t=new Timer(outstream, true);
		
		if(build>-1){
			Data.setGenome(build);
			AbstractIndex.MINCHROM=1;
			AbstractIndex.MAXCHROM=Data.numChroms;
			if(minChrom<0){minChrom=1;}
			if(maxChrom<0 || maxChrom>Data.numChroms){maxChrom=Data.numChroms;}
			outstream.println("Set genome to "+Data.GENOME_BUILD);
			
			if(RefToIndex.AUTO_CHROMBITS){RefToIndex.chrombits=autoChromBits();}
			if(RefToIndex.chrombits!=-1){
				BBIndex.setChromBits(RefToIndex.chrombits);
				if(verbose_stats>0){outstream.println("Set CHROMBITS to "+RefToIndex.chrombits);}
			}
		}
		
		assert(minChrom>=AbstractIndex.MINCHROM && maxChrom<=AbstractIndex.MAXCHROM) :
			minChrom+", "+maxChrom+", "+AbstractIndex.MINCHROM+", "+AbstractIndex.MAXCHROM;
		AbstractIndex.MINCHROM=minChrom;
		AbstractIndex.MAXCHROM=maxChrom;
		
		if(targetGenomeSize>0){
			long bases=Data.numDefinedBases;
			long x=Tools.max(1, Math.round(0.25f+bases*1d/targetGenomeSize));
			BBMapThread.setExpectedSites((int)x);
			outstream.println("Set EXPECTED_SITES to "+x);
		}
		
		assert(!(PERFECTMODE && SEMIPERFECTMODE));
		if(PERFECTMODE){setPerfectMode();}
		if(SEMIPERFECTMODE){setSemiperfectMode();}
		
		//Optional section for discrete timing of chrom array loading
		if((SLOW_ALIGN && !Boolean.getBoolean("bbmap3.pseudoAlign")) ||
				AbstractIndex.USE_EXTENDED_SCORE || useRandomReads || MAKE_MATCH_STRING){
			outstream.println();
			if(INDEX_LOADED){
				//do nothing
			}else if(RefToIndex.chromlist==null){
				Data.loadChromosomes(minChrom, maxChrom);
			}else{
				assert(RefToIndex.chromlist.size()==maxChrom-minChrom+1) : RefToIndex.chromlist.size();
				for(ChromosomeArray cha : RefToIndex.chromlist){
					Data.chromosomePlusMatrix[cha.chromosome]=cha;
				}
			}
			t.stop();
			outstream.println("Loaded Reference:\t"+t);
			t.start();
		}
		RefToIndex.chromlist=null;
		
		t.start();
		BBIndex.loadIndex(minChrom, maxChrom, keylen, !RefToIndex.NODISK, RefToIndex.NODISK);
		
		{
			long len=Data.numDefinedBases;
			if(len<300000000){
				BBIndex.MAX_HITS_REDUCTION2+=1;
				BBIndex.MAXIMUM_MAX_HITS_REDUCTION+=1;
				if(len<30000000){
					BBIndex.setFractionToExclude(BBIndex.FRACTION_GENOME_TO_EXCLUDE*0.5f);
					BBIndex.MAXIMUM_MAX_HITS_REDUCTION+=1;
					BBIndex.HIT_REDUCTION_DIV=Tools.max(BBIndex.HIT_REDUCTION_DIV-1, 3);
				}else if(len<100000000){
					BBIndex.setFractionToExclude(BBIndex.FRACTION_GENOME_TO_EXCLUDE*0.6f);
				}else{
					BBIndex.setFractionToExclude(BBIndex.FRACTION_GENOME_TO_EXCLUDE*0.75f);
				}
			}
		}
		
		t.stop();
		outstream.println("Generated Index:\t"+t);
		t.start();
		
		if((!SLOW_ALIGN || Boolean.getBoolean("bbmap3.pseudoAlign")) &&
				!AbstractIndex.USE_EXTENDED_SCORE && !useRandomReads && !MAKE_MATCH_STRING){
			for(int chrom=minChrom; chrom<=maxChrom; chrom++){
				Data.unload(chrom, true);
			}
		}
		
		if(ReadWrite.countActiveThreads()>0){
			ReadWrite.waitForWritingToFinish();
			t.stop();
			outstream.println("Finished Writing:\t"+t);
			t.start();
		}
		
		if(coverageBinned!=null || coverageBase!=null || rangeCov!=null || coverageHist!=null || coverageStats!=null || coverageRPKM!=null || normcov!=null || normcovOverall!=null || calcCov){
			String[] cvargs=("covhist="+coverageHist+"\tcovstats="+coverageStats+"\tbasecov="+coverageBase+"\trangecov="+rangeCov+"\tbincov="+coverageBinned+"\tphyscov="+coveragePhysical+
					"\t32bit="+cov32bit+"\tnzo="+covNzo+"\ttwocolumn="+covTwocolumn+"\tsecondary="+PRINT_SECONDARY_ALIGNMENTS+"\tcovminscaf="+coverageMinScaf+
					"\tksb="+covKsb+"\tbinsize="+covBinSize+"\tk="+covK+"\tstartcov="+covStartOnly+"\tstopcov="+covStopOnly+"\tstrandedcov="+covStranded+"\trpkm="+coverageRPKM+
					"\tnormcov="+normcov+"\tnormcovo="+normcovOverall+(in1==null ? "" : "\tin1="+in1)+(in2==null ? "" : "\tin2="+in2)+
					(covSetbs ? ("\tbitset="+covBitset+"\tarrays="+covArrays) : "")).split("\t");
			pileup=new CoveragePileup(cvargs);
			pileup.createDataStructures();
			pileup.loadScaffoldsFromIndex(minChrom, maxChrom);
		}
		
		if(!forceanalyze && (in1==null || maxReads==0)){return;}
		
		BBIndex.analyzeIndex(minChrom, maxChrom, BBIndex.FRACTION_GENOME_TO_EXCLUDE, keylen);
		
		t.stop("Analyzed Index:   ");
		t.start();
		
		if(makeBloomFilter){
			String serialPath=RefToIndex.bloomLoc(build);
			File serialFile=new File(serialPath);
//			System.err.println(serialPath+", "+serialFile.exists()+", "+bloomSerial+", "+RefToIndex.NODISK);
			if(bloomSerial && !RefToIndex.NODISK && serialFile.exists()){
				bloomFilter=ReadWrite.read(BloomFilter.class, RefToIndex.bloomLoc(build), true);
				t.stop("Loaded Bloom Filter: ");
			}else{
				if(bloomSerial){System.out.println("Could not read "+serialPath+", generating filter from reference.");}
				bloomFilter=new BloomFilter(true, bloomFilterK, bloomFilterK, 1, bloomFilterHashes, bloomFilterMinHits, true);
				t.stop("Made Bloom Filter: ");
				if(bloomSerial && !RefToIndex.NODISK && !RefToIndex.FORCE_READ_ONLY){
//					 && serialFile.canWrite()
					try {
						ReadWrite.writeObjectInThread(bloomFilter, serialPath, true);
						outstream.println("Writing Bloom Filter.");
					} catch (Throwable e) {
						e.printStackTrace();
						outstream.println("Can't Write Bloom Filter.");
					}
				}
			}
			outstream.println(bloomFilter.filter.toShortString());
			t.start();
		}
//		assert(false) : makeBloomFilter;
//		assert(false) : RefToIndex.chrombits+", "+AbstractIndex.CHROMS_PER_BLOCK;
	}

	/**
	 * Selects the largest chromosome grouping that satisfies both packed-coordinate
	 * capacity and the per-block {@code int[]} allocation ceiling.  The first bound
	 * preserves BBIndex's signed-positive v4 coordinate layout; the second carries
	 * BBMap2's independent block-budget idea into the existing compatible layout.
	 */
	static int autoChromBits(){
		final int maxLength=Tools.max(Data.chromLengths);
		int chrombits=Tools.min(Integer.numberOfLeadingZeros(maxLength)-1, 16);
		while(chrombits>0 && maxBlockBases(chrombits)>Shared.MAX_ARRAY_LEN){chrombits--;}
		return chrombits;
	}

	/** Conservative upper bound for the sites array in any mask-aligned block. */
	static long maxBlockBases(int chrombits){
		final int low=(1<<chrombits)-1;
		long max=0;
		for(int i=1; i<=Data.numChroms;){
			final int a=Tools.max(1, i&~low);
			final int b=Tools.min(Data.numChroms, i|low);
			long sum=0;
			for(int c=a; c<=b; c++){sum+=Data.chromLengths[c];}
			max=Tools.max(max, sum);
			i=b+1;
		}
		return max;
	}
		
	/**
	 * Executes the alignment pipeline with Streamer/Writer factories.
	 * Each ProcessThread directly claims input lists and submits the same
	 * dense list IDs to the shared Writer routes, following the pattern in
	 * template/A_SampleStreamerMT.java.
	 */
	@Override
	public void testSpeed(String[] args){
		if(in1==null || maxReads==0){
			outstream.println("No reads to process; quitting.");
			return;
		}

		Timer t=new Timer();
		// Preserve BBMap's memory-based worker limit before allocating reader/mapper workers.
		if(Shared.USE_JNI){
			final int requestedThreads=Shared.threads();
			adjustThreadsforMemory(105);
			if(Shared.threads()<requestedThreads*0.9){
				outstream.println("Disabling JNI due to low system memory.");
				Shared.USE_JNI=false;
				Shared.setThreads(requestedThreads);
			}
		}
		if(!Shared.USE_JNI){adjustThreadsforMemory(65);}
		final int threads=Tools.max(1, Shared.threads());
		Read.VALIDATE_IN_CONSTRUCTOR=(threads<2);

		final FileFormat ffIn1=FileFormat.testInput(in1, FileFormat.FASTQ, 0, 0, true, true, false);
		final FileFormat ffIn2=FileFormat.testInput(in2, FileFormat.FASTQ, 0, 0, true, true, false);
		final Streamer streamer=StreamerFactory.getReadInputStream(maxReads, ffIn1.samOrBam(),
				ffIn1, ffIn2, qfin1, qfin2, threads);
		streamer.setSampleRate(samplerate, sampleseed);
		streamer.start();
		final boolean paired=streamer.paired();
		if(paired){BBIndex.QUIT_AFTER_TWO_PERFECTS=false;}
		if(hybridPair && !paired){
			// Real pairedness (interleaved=t included) is only known here, not in postparse();
			// close the just-started Streamer explicitly since the try/finally below that would
			// otherwise cover it does not begin until after writer/engine construction.
			ReadWrite.closeStream(streamer);
			throw new RuntimeException("hybridpair requires paired input; got single-ended reads.");
		}
		final int buff=(!ORDERED ? 12 : Tools.max(32, 2*threads));
		final Writer[] writers=openWriters(args, buff, paired);

		AbstractMapThread.CALC_STATISTICS=CALC_STATISTICS;
		final AbstractMapThread[] mtts=new AbstractMapThread[threads];
		final ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			final BBMapThread engine=new BBMapThread(paired, keylen,
					pileup, SLOW_ALIGN, CORRECT_THRESH, minChrom,
					maxChrom, keyDensity, maxKeyDensity, minKeyDensity, maxDesiredKeys, REMOVE_DUPLICATE_BEST_ALIGNMENTS,
					SAVE_AMBIGUOUS_XY, MINIMUM_ALIGNMENT_SCORE_RATIO, TRIM_LIST, MAKE_MATCH_STRING, QUICK_MATCH_STRINGS,
					null, null, null, null,
					SLOW_ALIGN_PADDING, SLOW_RESCUE_PADDING, OUTPUT_MAPPED_ONLY, DONT_OUTPUT_BLACKLISTED_READS, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS,
					REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, KILL_BAD_PAIRS, rcompMate,
					PERFECTMODE, SEMIPERFECTMODE, FORBID_SELF_MAPPING, TIP_SEARCH_DIST,
					ambiguousRandom, ambiguousAll, KFILTER, MIN_IDFILTER, qtrimLeft, qtrimRight, untrim, TRIM_QUALITY, minTrimLength,
					LOCAL_ALIGN, RESCUE, STRICT_MAX_INDEL, MSA_TYPE, bloomFilter, hybridTipSearchCeiling, hybridPair,
					hybridMaxIndelConfig);
			engine.idmodulo=idmodulo;
			if(verbose){
				engine.verbose=verbose;
				engine.index().verbose=verbose;
			}
			mtts[i]=engine;
			alpt.add(new ProcessThread(streamer, writers, engine, i));
		}

		boolean success=false;
		try{
			success=ThreadWaiter.startAndWait(alpt, this);
			if(success){
				for(Writer writer : writers){
					// ReadWrite's Writer closure checks both errorState and successful completion.
					if(writer!=null && (writer.poisonAndWait() || !writer.finishedSuccessfully())){
						success=false;
					}
				}
			}else{
				for(Writer writer : writers){
					if(writer!=null){writer.finishError();}
				}
			}
		}finally{
			// Some readers report malformed input through errorState rather than nextList().
			if(ReadWrite.closeStream(streamer)){
				errorStateS=true;
				success=false;
			}
			errorStateS|=closeSplitterStreams(!success);
		}

		t.stop();
		if(printStats){outstream.println("\n\n   ------------------   Results   ------------------   ");}
		if(Boolean.getBoolean("bbmap3.quantumTieredShadow") ||
				Boolean.getBoolean("bbmap3.quantumTieredMutate")){
			final TieredQuantumStats total=new TieredQuantumStats();
			for(AbstractMapThread mtt : mtts){
				final TieredQuantumStats worker=((BBMapThread)mtt).quantumTieredStats();
				assert(worker!=null) : "Tiered Quantum mode must create worker statistics";
				total.add(worker);
			}
			outstream.print(total.toTsv());
		}
		if(Boolean.getBoolean("bbmap3.quantumOnly")){
			final QuantumOnlyStats total=new QuantumOnlyStats();
			for(AbstractMapThread mtt : mtts){
				final QuantumOnlyStats worker=((BBMapThread)mtt).quantumOnlyStats();
				assert(worker!=null) : "Quantum-only mode must create worker statistics";
				total.add(worker);
			}
			outstream.print(total.toTsv());
		}
		if(Boolean.getBoolean("bbmap3.quantumHybrid")){
			final QuantumHybridStats total=new QuantumHybridStats();
			for(AbstractMapThread mtt : mtts){
				final QuantumHybridStats worker=((BBMapThread)mtt).quantumHybridStats();
				assert(worker!=null) : "Broad Quantum hybrid mode must create worker statistics";
				total.add(worker);
			}
			outstream.print(total.toTsv());
		}
		if(hybridMaxIndelConfig!=null){
			final HybridMaxIndelStats total=new HybridMaxIndelStats();
			for(AbstractMapThread mtt : mtts){
				final HybridMaxIndelStats worker=((BBMapThread)mtt).hybridMaxIndelStats();
				assert(worker!=null) : "Selective max-indel mode must create worker statistics";
				total.add(worker);
			}
			outstream.print(total.toTsv());
		}

		printOutput(mtts, t, keylen, paired, false, pileup, scafNzo, sortStats, statsOutputFile);
		if(!success || errorStateS){throw new RuntimeException("BBMapS terminated in an error state; the output may be corrupt.");}
	}

	private Writer[] openWriters(String[] args, int buff, boolean paired){
		final Writer[] writers=new Writer[4]; // A, M, U, B
		if(OUTPUT_READS){
			ReadStreamWriter.MINCHROM=minChrom;
			ReadStreamWriter.MAXCHROM=maxChrom;
			writers[0]=makeWriter(outFile, outFile2, qfout, qfout2, buff);
			writers[1]=makeWriter(outFileM, outFileM2, qfoutM, qfoutM2, buff);
			writers[2]=makeWriter(outFileU, outFileU2, qfoutU, qfoutU2, buff);
			writers[3]=Data.scaffoldPrefixes ? null : makeWriter(outFileB, outFileB2, qfoutB, qfoutB2, buff);
		}
		if(Data.scaffoldPrefixes){
			BBMapSplitterS.streamTable=BBMapSplitterS.makeOutputStreams(args, OUTPUT_READS, true, buff, paired, overwrite, append, false);
			if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_SPLIT){
				BBMapSplitterS.streamTableAmbiguous=BBMapSplitterS.makeOutputStreams(args, OUTPUT_READS, true, buff, paired, overwrite, append, true);
			}
		}else{
			BBSplitter.TRACK_SET_STATS=false;
		}
		if(BBSplitter.TRACK_SET_STATS){
			outstream.print("Creating ref-set statistics table: ");
			BBSplitter.makeSetCountTable();
			outstream.println("done.");
		}
		// Scaffold counting requires the shared table consumed by AbstractMapper reporting.
		if(BBSplitter.TRACK_SCAF_STATS && BBSplitter.scafCountTable==null){
			BBSplitter.makeScafCountTable();
		}
		return writers;
	}

	private Writer makeWriter(String file1, String file2, String qf1, String qf2, int buff){
		if(file1==null){return null;}
		final FileFormat ff1=FileFormat.testOutput(file1, DEFAULT_OUTPUT_FORMAT, 0, 0, true, overwrite, append, true);
		final FileFormat ff2=file2==null ? null : FileFormat.testOutput(file2, DEFAULT_OUTPUT_FORMAT, 0, 0, true, overwrite, append, true);
		AbstractMapThread.OUTPUT_SAM|=ff1.samOrBam();
		final Writer writer=WriterFactory.getStream(ff1, ff2, qf1, qf2, buff, null, false, Shared.threads());
		writer.start();
		return writer;
	}

	/** @return true on any splitter close error, including an aborted stream. */
	private boolean closeSplitterStreams(boolean error){
		boolean closeError=false;
		if(BBMapSplitterS.streamTable!=null){
			for(Writer ros : BBMapSplitterS.streamTable.values()){
				closeError|=closeSplitterStream(ros, error);
			}
		}
		if(BBMapSplitterS.streamTableAmbiguous!=null){
			for(Writer ros : BBMapSplitterS.streamTableAmbiguous.values()){
				closeError|=closeSplitterStream(ros, error);
			}
		}
		return closeError;
	}

	/** ReadWrite.closeStream folds errorState and unsuccessful writer completion. */
	private boolean closeSplitterStream(Writer ros, boolean error){
		if(ros==null){return false;}
		if(error){ros.finishError();return true;}
		return ros.poisonAndWait() || !ros.finishedSuccessfully();
	}

	/** Abandon invalid output and wake mapper workers blocked on missing ordered lists. */
	private void abortSplitterStreams(){
		if(BBMapSplitterS.streamTable!=null){
			for(Writer ros : BBMapSplitterS.streamTable.values()){
				if(ros!=null){ros.finishError();}
			}
		}
		if(BBMapSplitterS.streamTableAmbiguous!=null){
			for(Writer ros : BBMapSplitterS.streamTableAmbiguous.values()){
				if(ros!=null){ros.finishError();}
			}
		}
	}

	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			readsProcessedS+=pt.readsProcessedT;
			basesProcessedS+=pt.basesProcessedT;
			errorStateS|=!pt.success;
		}
	}

	@Override
	public final boolean success(){return !errorStateS;}

	@Override
	public final ReadWriteLock rwlock(){return rwlockS;}

	private final ReadWriteLock rwlockS=new ReentrantReadWriteLock();
	private long readsProcessedS=0;
	private long basesProcessedS=0;
	private volatile boolean errorStateS=false;

	final class ProcessThread extends Thread {
		ProcessThread(Streamer streamer_, Writer[] writers_, BBMapThread engine_, int tid_){
			streamer=streamer_;
			writers=writers_;
			engine=engine_;
			tid=tid_;
		}

		@Override
		public void run(){
			try{
				ListNum<Read> ln=streamer.nextList();
				// Sampling can empty a data batch; Streamer.nextList() uses null for EOF.
				// Process empty lists too so each ordered Writer receives every list ID.
				// TODO: A_SampleStreamerMT has the same empty-list EOF predicate; audit it separately.
				while(ln!=null){
					processList(ln);
					ln=streamer.nextList();
				}
				success=true;
			}catch(Throwable t){
				error=t;
				t.printStackTrace();
				// A failed worker may never submit its ordered list ID. Signal immediately:
				// waiting for startAndWait to return would strand siblings in output backpressure.
				// Writer.finishError is idempotent/nonblocking; CROS.abort wakes ordered add waits.
				for(Writer writer : writers){
					if(writer!=null){writer.finishError();}
				}
				abortSplitterStreams();
			}
		}

		private void processList(ListNum<Read> ln){
			// Skip original input IDs before shredding, statistics, or mapping. Keep ln.id
			// even for a wholly skipped batch so ordered output receives every list ID.
			if(AbstractMapThread.SKIP_INITIAL>0){
				final ArrayList<Read> input=ln.list;
				int kept=0;
				for(int i=0; i<input.size(); i++){
					final Read r=input.get(i);
					if(r.numericID>=AbstractMapThread.SKIP_INITIAL){input.set(kept++, r);}
				}
				input.subList(kept, input.size()).clear();
			}
			ArrayList<Read> readlist=engine.handleLongReads(ln.list);
			final LongList bloomBuffer=(engine.bloomFilter==null ? null : new LongList(150));

			for(Read r : readlist){
				final long startTime=AbstractMapThread.TIME_TAG ? System.nanoTime() : 0;
				readsProcessedT+=r.pairCount();
				basesProcessedT+=r.pairLength();
				engine.readsIn1++;
				engine.readsIn2+=r.mateCount();
				engine.basesIn1+=r.length();
				engine.basesIn2+=r.mateLength();
				final Read r2=r.mate;

				final boolean passesBloom=(engine.bloomFilter!=null && engine.bloomFilter.passes(r, r2, bloomBuffer, 1));
				if(passesBloom){
					engine.basesUsed1+=r.length();
					engine.basesUsed2+=r.mateLength();
					engine.readsPassedBloomFilter+=r.pairCount();
					engine.basesPassedBloomFilter+=r.pairLength();
					engine.readsUsed1++;
					engine.readsUsed2+=r.mateCount();
					continue;
				}

				if(r.synthetic()){
					engine.syntheticReads++;
					if(r.originalSite==null){r.makeOriginalSite();}
					r.clearSite();
					if(r2!=null){
						assert(r2.synthetic());
						if(r2.originalSite==null){r2.makeOriginalSite();}
						r2.clearSite();
					}
				}
				r.clearAnswers(true);
				assert(r.bases==null || r.length()<=engine.maxReadLength()) :
					"Read "+r.numericID+", length "+r.length()+" exceeds the limit of "+engine.maxReadLength();

				if(engine.readstats!=null){
					if(ReadStats.COLLECT_QUALITY_STATS){engine.readstats.addToQualityHistogram(r);}
					if(ReadStats.COLLECT_BASE_STATS){engine.readstats.addToBaseHistogram(r);}
					if(ReadStats.COLLECT_LENGTH_STATS){engine.readstats.addToLengthHistogram(r);}
					if(ReadStats.COLLECT_GC_STATS){engine.readstats.addToGCHistogram(r);}
				}
				if(engine.TRIM_LEFT || engine.TRIM_RIGHT){
					TrimRead.trim(r, engine.TRIM_LEFT, engine.TRIM_RIGHT, engine.TRIM_QUAL, engine.TRIM_ERROR_RATE, engine.TRIM_MIN_LENGTH);
					TrimRead.trim(r2, engine.TRIM_LEFT, engine.TRIM_RIGHT, engine.TRIM_QUAL, engine.TRIM_ERROR_RATE, engine.TRIM_MIN_LENGTH);
				}
				if(AbstractMapThread.RCOMP){r.reverseComplementFast();}

				if(r2==null){
					final byte[] basesM=AminoAcid.reverseComplementBases(r.bases);
					engine.basesUsed1+=(basesM==null ? 0 : basesM.length);
					engine.processRead(r, basesM);
					engine.capSiteList(r, engine.MAX_SITESCORES_TO_PRINT, engine.PRINT_SECONDARY_ALIGNMENTS);
					assert(Read.CHECKSITES(r, basesM));
				}else{
					if(engine.RCOMP_MATE!=AbstractMapThread.RCOMP){r2.reverseComplementFast();}
					final byte[] basesM1=AminoAcid.reverseComplementBases(r.bases);
					final byte[] basesM2=AminoAcid.reverseComplementBases(r2.bases);
					engine.basesUsed1+=(basesM1==null ? 0 : basesM1.length);
					engine.basesUsed2+=(basesM2==null ? 0 : basesM2.length);
					assert(r2.bases==null || r2.length()<=engine.maxReadLength()) :
						"Read "+r2.numericID+" exceeds the limit of "+engine.maxReadLength();
					engine.processReadPair(r, basesM1, basesM2);
					engine.capSiteList(r, engine.MAX_SITESCORES_TO_PRINT, engine.PRINT_SECONDARY_ALIGNMENTS);
					engine.capSiteList(r2, engine.MAX_SITESCORES_TO_PRINT, engine.PRINT_SECONDARY_ALIGNMENTS);
					assert(Read.CHECKSITES(r, basesM1));
					assert(Read.CHECKSITES(r2, basesM2));
				}

				if(engine.UNTRIM && (engine.TRIM_LEFT || engine.TRIM_RIGHT)){
					TrimRead.untrim(r);
					TrimRead.untrim(r2);
				}
				if(engine.readstats!=null){
					if(ReadStats.COLLECT_MATCH_STATS){engine.readstats.addToMatchHistogram(r);}
					if(ReadStats.COLLECT_INSERT_STATS && r.paired()){engine.readstats.addToInsertHistogram(r, (engine.SAME_STRAND_PAIRS || !engine.REQUIRE_CORRECT_STRANDS_PAIRS));}
					if(ReadStats.COLLECT_QUALITY_ACCURACY){engine.readstats.addToQualityAccuracy(r);}
					if(ReadStats.COLLECT_ERROR_STATS){engine.readstats.addToErrorHistogram(r);}
					if(ReadStats.COLLECT_INDEL_STATS){engine.readstats.addToIndelHistogram(r);}
					if(ReadStats.COLLECT_IDENTITY_STATS){engine.readstats.addToIdentityHistogram(r);}
				}
				if(AbstractMapThread.TIME_TAG){
					final Long elapsed=(System.nanoTime()-startTime+500)/1000;
					r.setObj(elapsed);
					if(r2!=null){r2.setObj(elapsed);}
					if(engine.readstats!=null && ReadStats.COLLECT_TIME_STATS){engine.readstats.addToTimeHistogram(r);}
				}
			}

			if(engine.RenameByInsert){
				final boolean ignoreStrand=(!engine.REQUIRE_CORRECT_STRANDS_PAIRS || engine.SAME_STRAND_PAIRS);
				for(Read r : readlist){
					if(r.mapped() && r.mateMapped() && r.paired()){
						final int insert=Read.insertSizeMapped(r, r.mate, ignoreStrand);
						final String s="insert="+insert;
						r.id=s+" 1:"+r.numericID;
						r.mate.id=s+" 2:"+r.numericID;
					}
				}
			}
			if(engine.pileup!=null){
				synchronized(engine.pileup){
					for(Read r : readlist){
						engine.pileup.processRead(r);
						if(r.mate!=null){engine.pileup.processRead(r.mate);}
					}
				}
			}
			emit(ln, readlist);
		}

		private void emit(ListNum<Read> ln, ArrayList<Read> readlist){
			final long id=ln.id;
			final boolean black=Blacklist.hasBlacklist();
			if(BBMapSplitterS.streamTable!=null || BBSplitter.TRACK_SET_STATS || BBSplitter.TRACK_SCAF_STATS){
				BBMapSplitterS.printReads(readlist, id, null, engine.CLEARZONE1());
			}
			final ArrayList<Read> mapped=new ArrayList<Read>(readlist.size());
			final ArrayList<Read> unmapped=new ArrayList<Read>(readlist.size());
			final ArrayList<Read> blacklisted=new ArrayList<Read>(readlist.size());
			for(Read r : readlist){
				if(r!=null){
					final Read r2=r.mate;
					final boolean isMapped=(r.mapped() || (r2!=null && r2.mapped()));
					if(isMapped){
						if(!black || !Blacklist.inBlacklist(r)){mapped.add(r);}
					}else{unmapped.add(r);}
					if(black && Blacklist.inBlacklist(r)){blacklisted.add(r);}
				}
			}
			if(writers[1]!=null){writers[1].addReads(new ListNum<Read>(mapped, id));}
			if(writers[3]!=null){writers[3].addReads(new ListNum<Read>(blacklisted, id));}
			if(writers[2]!=null){writers[2].addReads(new ListNum<Read>(unmapped, id));}
			if(writers[0]!=null){
				if(engine.OUTPUT_MAPPED_ONLY){AbstractMapThread.removeUnmapped(readlist);}
				if(engine.DONT_OUTPUT_BLACKLISTED_READS){AbstractMapThread.removeBlacklisted(readlist);}
				for(Read r : readlist){
					if(r!=null){
						if(AbstractMapThread.CLEAR_ATTACHMENT){r.nullifyObject();}
						assert(r.bases!=null);
						if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
					}
				}
				writers[0].addReads(new ListNum<Read>(readlist, id));
			}
		}

		long readsProcessedT=0;
		long basesProcessedT=0;
		boolean success=false;
		Throwable error=null;
		final Streamer streamer;
		final Writer[] writers;
		final BBMapThread engine;
		final int tid;
	}

	/** Only BBMapS wires the per-invocation hybrid controller through to BBMapThread. */
	@Override
	boolean supportsHybridPair(){return true;}
	@Override
	boolean supportsHybridMaxIndel(){return true;}

	/**
	 * Configures parameters for semi-perfect alignment mode.
	 * Reduces key density requirements and alignment score thresholds
	 * to allow alignments with small numbers of mismatches.
	 */
	@Override
	void setSemiperfectMode() {
		assert(SEMIPERFECTMODE);
		if(SEMIPERFECTMODE){
			TRIM_LIST=false;
			keyDensity/=2;
			maxKeyDensity/=2;
			minKeyDensity=1.1f;
			maxDesiredKeys/=2;
			MINIMUM_ALIGNMENT_SCORE_RATIO=0.45f;
			BBIndex.setSemiperfectMode();
		}
	}

	/**
	 * Configures parameters for perfect alignment mode.
	 * Sets maximum alignment score ratio to require exact matches
	 * and adjusts key density for perfect-match detection.
	 */
	@Override
	void setPerfectMode() {
		assert(PERFECTMODE);
		if(PERFECTMODE){
			TRIM_LIST=false;
			keyDensity/=2;
			maxKeyDensity/=2;
			minKeyDensity=1.1f;
			maxDesiredKeys/=2;
			MINIMUM_ALIGNMENT_SCORE_RATIO=1.0f;
			BBIndex.setPerfectMode();
		}
	}
	

	/**
	 * Prints current alignment configuration settings.
	 * Displays key density, index parameters, hit filtering settings,
	 * and other alignment options based on verbosity level.
	 * @param k K-mer length for alignment
	 */
	@Override
	void printSettings(int k){
		
		printSettings0(k, BBIndex.MAX_INDEL, MINIMUM_ALIGNMENT_SCORE_RATIO);
		
		if(verbose_stats>=2){
			outstream.println("Key Density:          \t"+keyDensity+" ("+minKeyDensity+" ~ "+maxKeyDensity+")");
			outstream.println("Max keys:             \t"+maxDesiredKeys);
			
			outstream.println("Block Subsections:     \t"+BBIndex.CHROMS_PER_BLOCK);
			outstream.println("Fraction To Remove:    \t"+Tools.format("%.4f", (BBIndex.REMOVE_FREQUENT_GENOME_FRACTION ? BBIndex.FRACTION_GENOME_TO_EXCLUDE : 0)));
			//		sysout.println("ADD_SCORE_Z:           \t"+Index4.ADD_SCORE_Z);
			outstream.println("Hits To Keep:          \t"+BBIndex.MIN_APPROX_HITS_TO_KEEP);
		}
		
		if(verbose_stats>=3){
			outstream.println("Remove Clumpy:         \t"+BBIndex.REMOVE_CLUMPY);
			if(BBIndex.REMOVE_CLUMPY){
				outstream.println("CLUMPY_MAX_DIST:       \t"+BBIndex.CLUMPY_MAX_DIST);
				outstream.println("CLUMPY_MIN_LENGTH:     \t"+BBIndex.CLUMPY_MIN_LENGTH_INDEX);
				outstream.println("CLUMPY_FRACTION:       \t"+BBIndex.CLUMPY_FRACTION);
			}
			outstream.println("Remove Long Lists:     \t"+BBIndex.TRIM_LONG_HIT_LISTS);
			if(BBIndex.TRIM_LONG_HIT_LISTS){
				outstream.println("HIT_FRACTION_TO_RETAIN:\t"+BBIndex.HIT_FRACTION_TO_RETAIN);
			}
			outstream.println("Trim By Greedy:        \t"+BBIndex.TRIM_BY_GREEDY);
			outstream.println("Trim By Total Sites:   \t"+BBIndex.TRIM_BY_TOTAL_SITE_COUNT);
			if(BBIndex.TRIM_BY_TOTAL_SITE_COUNT){
				outstream.println("MAX_AVG_SITES:         \t"+BBIndex.MAX_AVERAGE_LIST_TO_SEARCH);
				outstream.println("MAX_AVG_SITES_2:       \t"+BBIndex.MAX_AVERAGE_LIST_TO_SEARCH2);
				outstream.println("MAX_SHORTEST_SITE:     \t"+BBIndex.MAX_SHORTEST_LIST_TO_SEARCH);
			}
			outstream.println("Index Min Score:       \t"+BBIndex.MIN_SCORE_MULT);

			outstream.println("Dynamic Trim:          \t"+BBIndex.DYNAMICALLY_TRIM_LOW_SCORES);
			if(BBIndex.DYNAMICALLY_TRIM_LOW_SCORES){
				outstream.println("DYNAMIC_SCORE_THRESH:  \t"+BBIndex.DYNAMIC_SCORE_THRESH);
			}
		}
		
	}

}
