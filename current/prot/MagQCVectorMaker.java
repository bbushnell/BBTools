package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntLongHashMap;
import tax.TaxNode;
import tax.TaxTree;
import tracker.KmerTracker;

/**
 * Generates MAG-QC training vectors by synthesizing bins from the per-contig
 * precompute cache (Stage 1) and emitting each bin's feature vector with its
 * ACHIEVED completeness and contamination as continuous regression targets.
 *
 * <p>A synthetic bin = a random subset of one TARGET organism's contigs (to hit a
 * sampled completeness against the target's true genome size, plasmids included)
 * plus, optionally, contigs from one or more CONTAMINANT organisms (to hit a sampled
 * contamination). Labels are computed from the explicit selected target, not a
 * re-inferred winner: completeness = target_bp / genomeSize[target];
 * contamination = foreign_bp / (target_bp + foreign_bp). Foreign target bases are
 * c/(1-c)*cleanBases. The generator emits the ACHIEVED labels after whole-contig
 * selection.
 *
 * <p>Organisms are split globally BEFORE sampling; a held-out organism appears as
 * neither target nor contaminant in the training file (Barbara design review #3).
 * The feature vector uses the same cache the deployment path will produce, so the
 * aggregation here is the train==serve contract.
 *
 * <p>Vector layout: [family features] + [13 global stats] + [phylum one-hot], then
 * the two outputs completeness, contamination. The family-count encoding is selected
 * by enc= : ratio (default) N/(1+N); raw min(N,32)/32; log log2(1+N)/log2(65); two =
 * two columns per family (presence 0/1, then excess-copies min(N-1,16)/16); norm = 0 if
 * absent else N/avgCopyWhenPresent[family] (per-family baseline so present-at-typical is
 * ~1.0 and duplication is &gt;1.0, calibrated to each family's natural copy number). The
 * two/norm schemes keep the duplication signal that flags contamination uncompressed
 * (CheckM2 uses raw counts for exactly this reason).
 *
 * <p>Idiom rewrite (2026-08): I/O runs on {@link ByteFile}/{@link ByteStreamWriter}
 * with a reused {@link LineParser1} (byte-range parsing, zero per-field String
 * allocation except where a value must persist as a String or is on a load-once,
 * bounded-frequency path). Every tid-keyed lookup used inside the per-bin hot path
 * ({@code makeBin}) is a primitive {@link structures.IntHashMap}/{@link IntLongHashMap}
 * or a dense array (byTid's contig lists, since IntHashMap can't hold list values) -
 * no boxed {@code HashMap<Integer,...>} on the hot path. {@code writeSet}'s per-attempt
 * buffers ({@code fam}/{@code glob}/{@code labels}) and the per-bin {@link Agg}
 * accumulator are now hoisted out of the attempt loop and reused (cleared, not
 * reallocated). {@code fmt()}'s exact string representation (whole numbers print
 * without a decimal point, everything else at 6 fixed decimals) is replicated by
 * {@link #appendFmt(ByteBuilder, double)} for the hot output paths; the String-
 * returning {@link #fmt(double)} is KEPT for the aggregator's roundtrip-through-string
 * rounding (low frequency, and an actual String is what {@code Float.parseFloat} needs).
 * NOT threaded this pass (Brian's explicit constraint): {@code makeBin} draws from
 * one sequential {@link Random} stream per output set, and parallelizing bin synthesis
 * would reorder those draws and change the output. Algorithm, RNG call sequence, and
 * every output byte are UNCHANGED - see {@code MagQCVectorMakerTest}/
 * {@code MagQCAggVectorTest} for the pinned behavioral contract, verified against the
 * original by UMP45's differential (byte-identical across all fixture modes).
 */
public class MagQCVectorMaker {

	public static void main(String[] args){
		MagQCVectorMaker x=new MagQCVectorMaker(args);
		x.process();
	}

	public MagQCVectorMaker(String[] args){
		for(String arg : args){
			int eq=arg.indexOf('=');
			if(eq<0){continue;}
			String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("cache")){cacheFile=b;}
			else if(a.equals("sizemap")){sizemapFile=b;}
			else if(a.equals("familylist")){familyFile=b;}
			else if(a.equals("features")){featuresFile=b;}
			else if(a.equals("taxpgm")){taxpgmFile=b;}
			else if(a.equals("tree")){treeFile=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("outval")){outval=b;}
			else if(a.equals("n")){n=Long.parseLong(b);}
			else if(a.equals("valn")){valn=Long.parseLong(b);}
			else if(a.equals("valfrac")){valfrac=Double.parseDouble(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("minlen")){minlen=Integer.parseInt(b);}
			else if(a.equals("mixcomp")){mixComp=Double.parseDouble(b);}
			else if(a.equals("mixcont")){mixCont=Double.parseDouble(b);}
			else if(a.equals("cleanspike")){cleanSpike=Double.parseDouble(b);}
			else if(a.equals("multicontamprob")){multiContamProb=Double.parseDouble(b);}
			else if(a.equals("perfectfrac")){perfectFrac=Double.parseDouble(b);}
			else if(a.equals("nearperfectfrac")){nearPerfectFrac=Double.parseDouble(b);}
			else if(a.equals("benchtruth")){benchTruthFile=b; benchMode=true;}
			else if(a.equals("benchmanifest")){benchManifestFile=b; benchMode=true;}
			else if(a.equals("benchvec")){benchVecFile=b; benchMode=true;}
			else if(a.equals("benchbins")){benchBins=Long.parseLong(b);}
			else if(a.equals("samefamprob")){sameFamProb=Double.parseDouble(b);}
			else if(a.equals("enc")){enc=parseEnc(b);}
			else if(a.equals("subnet")){subnetName=b.toLowerCase();}
			else if(a.equals("subnetout")){subnetOut=b;}
			else if(a.equals("subnetvalout")){subnetValOut=b;}
			else if(a.equals("subsetfile")){subsetFile=b;}
			else if(a.equals("aggmanifest")){aggManifestFile=b;}
			else if(a.equals("aggout")){aggOut=b;}
			else if(a.equals("aggvalout")){aggValOut=b;}
			else if(a.equals("densehead")){denseHead=Integer.parseInt(b);}
			else if(a.equals("aggobs")){aggObsServe=parseAggObs(b);}
			else if(a.equals("poolmode")){poolMode=parsePoolMode(b);}
			else{System.err.println("Warning: unknown arg "+arg);}
		}
		if(cacheFile==null || sizemapFile==null || taxpgmFile==null || out==null){
			throw new RuntimeException("Required: cache= sizemap= taxpgm= out= [outval= familylist= tree= n= valn= ...]");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Per-contig          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * One frozen subnet participating in aggregator-vector emission: its subset
	 * definition (family ranks, or the ncRNA accessor for the special "ncrna" row),
	 * its loaded CellNet, and an exact-width input buffer (CellNet.applyInput asserts
	 * the width, which doubles as a wiring guard).
	 */
	static final class AggSubnet {
		String name;
		int numObs;
		int[] ranks;//null for the ncrna row
		ml.CellNet net;
		boolean ncrna;
		float[] buf;
	}

	/** One cached contig's sufficient statistics; family counts stored sparsely. */
	static final class Contig {
		int tid, length, gc, acgt, cds, mapped, coding, r16, r23, r5, rother, trna;
		long glenSum, glenSq;
		int[] famRank, famCount;
		int[] dimer;//16 dinucleotide counts (KmerTracker k=2 native index), cache field 18; null on pre-rebuild caches
		int[] antiRank, antiCount;//sparse tRNA anticodon code(0-63)->count, cache field 17; EMPTY if absent
		String name;//contig_id (cache field 0); populated only in benchmark mode, for the FASTA manifest
	}

	/*--------------------------------------------------------------*/
	/*----------------            Loaders           ----------------*/
	/*--------------------------------------------------------------*/

	private void loadAux(){
		// family list -> count of family columns
		if(familyFile!=null){
			final ByteFile bf=ByteFile.makeByteFile(familyFile, true);
			bf.nextLine();//header
			int c=0;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0){c++;}
			}
			bf.close();
			numFam=c;
		}
		// taxpgm: tid <tab> phylum <tab> pgm -> tid->phylum (staged), phylum vocabulary.
		// phylumIndex needs every distinct phylum name collected first (TreeSet, sorted
		// deterministic order), so tid->phylum is staged in parallel primitive/String
		// lists here and resolved to tid2phylumIdx (int->int) in a second pass below.
		final IntList stagedTid=new IntList();
		final ArrayList<String> stagedPhylum=new ArrayList<String>();
		{
			final ByteFile bf=ByteFile.makeByteFile(taxpgmFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			final TreeSet<String> phyla=new TreeSet<String>();
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				lp.set(line);
				if(lp.terms()>=2){
					final int tid=lp.parseInt(0);
					final String phy=lp.parseString(1);
					stagedTid.add(tid); stagedPhylum.add(phy);
					phyla.add(phy);
				}
			}
			bf.close();
			phylumList=new ArrayList<String>(phyla);
			phylumList.add("other");
			for(int i=0; i<phylumList.size(); i++){phylumIndex.put(phylumList.get(i), i);}
			numPhyla=phylumList.size();
		}
		for(int i=0; i<stagedTid.size(); i++){
			final Integer pi=phylumIndex.get(stagedPhylum.get(i));
			tid2phylumIdx.put(stagedTid.get(i), pi==null ? phylumIndex.get("other") : pi);
		}
		// sizemap: tid <tab> bp. Last row for a tid wins (matches HashMap.put's overwrite
		// semantics) - IntLongHashMap.put() does NOT overwrite, so remove-then-put.
		{
			final ByteFile bf=ByteFile.makeByteFile(sizemapFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				lp.set(line);
				if(lp.terms()>=2){
					final int tid=lp.parseInt(0);
					final long bp=lp.parseLong(1);
					genomeSize.remove(tid);
					genomeSize.put(tid, bp);
				}
			}
			bf.close();
		}
	}

	/** Loads the per-contig cache into per-tid contig lists (dense tid->index + parallel
	 *  ArrayList<Contig>[] - IntHashMap can't hold list values, so this is the one field
	 *  that keeps a per-tid object list rather than becoming a flat primitive map). */
	private void loadCache(){
		final ByteFile bf=ByteFile.makeByteFile(cacheFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final IntList rankBuf=new IntList(), countBuf=new IntList();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			final int terms=lp.terms();
			assert(terms>=17) : "Malformed cache row: "+terms+" fields (need >=17): "+new String(line);
			final Contig c=new Contig();
			c.tid=lp.parseInt(1);
			c.length=lp.parseInt(3);
			if(c.length<minlen){continue;}
			if(benchMode){c.name=lp.parseString(0);}//retain contig_id only when a benchmark manifest is needed
			c.gc=lp.parseInt(4);
			c.acgt=lp.parseInt(5);
			c.cds=lp.parseInt(6);
			c.mapped=lp.parseInt(7);
			c.glenSum=lp.parseLong(8);
			c.glenSq=lp.parseLong(9);
			c.coding=lp.parseInt(10);
			c.r16=lp.parseInt(11);
			c.r23=lp.parseInt(12);
			c.r5=lp.parseInt(13);
			c.rother=lp.parseInt(14);
			c.trna=lp.parseInt(15);
			final int flen=lp.length(16);
			if(flen==0){c.famRank=EMPTY; c.famCount=EMPTY;}
			else{
				rankBuf.clear(); countBuf.clear();
				parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
				c.famRank=rankBuf.toArray();
				c.famCount=countBuf.toArray();
			}
			// Rebuild fields (backward-compatible): field 17 = tRNA anticodon sparse code:count; field 18 = 16
			// dinucleotide counts (dense CSV). Absent on pre-rebuild 17-field caches -> anti EMPTY, dimer null.
			if(terms>=18){
				final int alen=lp.length(17);
				if(alen==0){c.antiRank=EMPTY; c.antiCount=EMPTY;}
				else{
					rankBuf.clear(); countBuf.clear();
					parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
					c.antiRank=rankBuf.toArray();
					c.antiCount=countBuf.toArray();
				}
			}else{c.antiRank=EMPTY; c.antiCount=EMPTY;}
			if(terms>=19){
				lp.length(18);
				c.dimer=parseDimers(lp.line(), lp.a(), lp.b());
			}
			int idx=tidToIdx.get(c.tid);
			if(idx<0){
				idx=contigLists.size();
				tidToIdx.put(c.tid, idx);
				contigLists.add(new ArrayList<Contig>());
				if(idx>=domainIdxArr.length){domainIdxArr=Arrays.copyOf(domainIdxArr, Math.max(idx+1, domainIdxArr.length*2));}
				domainIdxArr[idx]=domainIndex(lp.parseString(2));
			}
			contigLists.get(idx).add(c);
		}
		bf.close();
	}

	/** Parses "rank:count;rank:count;..." within [a,b) of line into the (cleared) output
	 *  lists, in order. Zero allocation beyond the two IntList's own backing-array growth. */
	private static void parseFamCounts(byte[] line, int a, int b, IntList outRank, IntList outCount){
		int start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==';'){
				if(i>start){
					int colon=-1;
					for(int j=start; j<i; j++){if(line[j]==':'){colon=j; break;}}
					assert(colon>start) : "Malformed famcounts field (missing ':' in a rank:count pair).";
					outRank.add(Parse.parseInt(line, start, colon));
					outCount.add(Parse.parseInt(line, colon+1, i));
				}
				start=i+1;
			}
		}
	}

	/** Parses exactly 16 comma-separated ints within [a,b) into a new int[16] (dense dinucleotide
	 *  counts, KmerTracker native index order). Zero-alloc beyond the returned array. */
	private static int[] parseDimers(byte[] line, int a, int b){
		final int[] out=new int[16];
		int idx=0, start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==','){
				if(i>start && idx<16){out[idx]=Parse.parseInt(line, start, i);}
				idx++;
				start=i+1;
			}
		}
		assert(idx==16) : "dimer field expected 16 comma-separated counts, got "+idx+": "+new String(line, a, b-a);
		return out;
	}

	/** Returns tid's contig list, or null if tid was never seen in the cache (mirrors
	 *  the original {@code byTid.get(tid)}'s null-on-absent contract exactly). */
	private ArrayList<Contig> getContigs(int tid){
		final int idx=tidToIdx.get(tid);
		return idx<0 ? null : contigLists.get(idx);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Process            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads a family-feature subset (one rank per line) for reduced-width vectors. */
	private int[] loadRanks(String file){
		final IntList l=new IntList();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			int a=0, b=line.length;
			while(a<b && line[a]<=' '){a++;}
			while(b>a && line[b-1]<=' '){b--;}
			if(b>a){l.add(Parse.parseInt(line, a, b));}
		}
		bf.close();
		final int[] a=l.toArray();
		System.err.println("feature subset: "+a.length+" family ranks kept");
		return a;
	}

	private static boolean parseAggObs(String s){
		s=s.toLowerCase();
		if(s.equals("serve") || s.equals("whole") || s.equals("wholebin")){return true;}
		if(s.equals("clean") || s.equals("target") || s.equals("targetonly")){return false;}
		throw new RuntimeException("Unknown aggobs="+s+" (serve|clean)");
	}

	private static int parsePoolMode(String s){
		s=s.toLowerCase();
		if(s.equals("trainval")){return POOL_TRAINVAL;}
		if(s.equals("valsplit")){return POOL_VALSPLIT;}
		if(s.equals("allbutc")){return POOL_ALLBUTC;}
		throw new RuntimeException("Unknown poolmode="+s+" (trainval|valsplit|allbutc)");
	}

	/**
	 * Loads the aggregator manifest: tab-separated rows
	 * base, numobs, numIn, subset_path, val_path, net_path (val_path unused here).
	 * base "ncrna" marks the special ncRNA row (subset_path "-", observed from the Agg
	 * ncRNA fields rather than fam[] ranks). Every net's input width is asserted against
	 * obs+phylum+SHARED_CONTEXT_WIDTH - the FROZEN 2026-08-24 vector-layout rebuild
	 * retired the old per-net representation flags (sndomain/snhhcaga/sngenelen/
	 * snbinscaled/sncodingaffine): every subnet now gets the identical standard context,
	 * so there is nothing left to override per row. #-lines are comments.
	 */
	private void loadAggManifest(String file){
		aggSubnets=new ArrayList<AggSubnet>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] lineB=bf.nextLine(); lineB!=null; lineB=bf.nextLine()){
			String line=new String(lineB).trim();
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			String[] p=line.split("\t");
			if(p.length<6){throw new RuntimeException("Manifest row needs 6 columns "
				+"(base numobs numIn subset_path val_path net_path): "+line);}
			AggSubnet s=new AggSubnet();
			s.name=p[0];
			s.numObs=Integer.parseInt(p[1]);
			s.ncrna="ncrna".equals(s.name);
			if(s.ncrna){
				if(s.numObs!=NCRNA_OBS){throw new RuntimeException("ncrna row numobs must be "+NCRNA_OBS);}
			}else{
				s.ranks=loadRanks(p[3]);
				if(s.ranks.length!=s.numObs){throw new RuntimeException("Subset "+s.name
					+": "+s.ranks.length+" ranks but numobs="+s.numObs);}
			}
			s.net=ml.CellNetParser.load(p[5]);
			if(s.net==null){throw new RuntimeException("Failed to load net "+p[5]);}
			final int expected=s.numObs+numPhyla+SHARED_CONTEXT_WIDTH;
			final int manifestIn=Integer.parseInt(p[2]);
			if(manifestIn>=0 && manifestIn!=expected){throw new RuntimeException("Subset "+s.name
				+": manifest numIn="+manifestIn+" but expected "+expected);}
			if(s.net.numInputs()!=expected){throw new RuntimeException("Subset "+s.name
				+": net takes "+s.net.numInputs()+" inputs but expected "+expected
				+" (this net predates the FROZEN vector-layout rebuild and must be retrained)");}
			s.buf=new float[expected];
			aggSubnets.add(s);
		}
		bf.close();
		if(aggSubnets.isEmpty()){throw new RuntimeException("Empty aggregator manifest: "+file);}
		System.err.println("aggregator manifest: "+aggSubnets.size()+" subnets loaded");
	}

	/**
	 * Selects the raw dense head: the K most-prevalent families (organism presence
	 * count over all usable orgs, ties broken by rank), a reference-DB constant like
	 * avgCopyWhenPresent. Their whole-bin counts feed the aggregator raw (enc=two
	 * presence+excess), giving it a direct low-missingness signal alongside the
	 * subnet summaries.
	 */
	private int[] computeDenseHead(IntList usable, int k){
		final int[] present=new int[numFam];
		final int[] orgCount=new int[numFam];
		for(int i=0; i<usable.size(); i++){
			Arrays.fill(orgCount, 0);
			for(Contig c : getContigs(usable.get(i))){
				for(int j=0; j<c.famRank.length; j++){orgCount[c.famRank[j]]+=c.famCount[j];}
			}
			for(int f=0; f<numFam; f++){if(orgCount[f]>0){present[f]++;}}
		}
		Integer[] order=new Integer[numFam];
		for(int i=0; i<numFam; i++){order[i]=i;}
		Arrays.sort(order, (a, b) -> (present[a]!=present[b] ? present[b]-present[a] : a-b));
		final int[] head=new int[Math.min(k, numFam)];
		for(int i=0; i<head.length; i++){head[i]=order[i];}
		return head;
	}

	void process(){
		loadAux();
		loadCache();

		// usable tids: present in cache + sizemap, with contigs
		final IntList usable=new IntList();
		final int[] allTids=tidToIdx.toArray();
		for(int tid : allTids){
			if(genomeSize.contains(tid) && !getContigs(tid).isEmpty()){usable.add(tid);}
		}
		usable.sort();
		System.err.println("usable orgs="+usable.size()+", numFam="+numFam+", numPhyla="+numPhyla);

		// avg~0.5 normalization scales (Brian 2026-08-24): corpus reference for the size/#genes/gene-length
		// context features. Computed now (layout-independent); WIRED INTO the emit methods in the one-pass
		// vector-layout rewrite. Currently compute-and-log only, so all existing output is byte-identical.
		computeNormScales(usable);

		// optional taxonomy for same-family contaminant bias
		if(treeFile!=null){
			tree=TaxTree.loadTaxTree(treeFile, System.err, false, false);
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				TaxNode fn=(tree==null ? null : tree.getNodeAtLevel(tid, TaxTree.FAMILY));
				tid2family.put(tid, fn==null ? -1 : fn.id);
			}
		}

		// global organism split BEFORE sampling. shuffleInPlace replicates
		// java.util.Collections.shuffle(List,Random)'s EXACT algorithm/call sequence -
		// IntList's own shuffle() draws from Shared.threadLocalRandom(), NOT this seeded
		// stream, and would silently desync the whole train/val split from the original.
		final Random split=new Random(seed);
		shuffleInPlace(usable, split);
		final int nVal=(int)Math.round(usable.size()*valfrac);
		IntList valTids=subList(usable, 0, nVal);
		IntList trainTids=subList(usable, nVal, usable.size());
		System.err.println("train orgs="+trainTids.size()+", val orgs="+valTids.size());

		// precompute recoverable bp per tid
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			long sum=0; for(Contig c : getContigs(tid)){sum+=c.length;}
			recoverable.put(tid, sum);
		}

		// precompute each organism's NATIVE ncRNA complement (the subnet denominator):
		// {r16,r23,r5,rother,trna} summed over ALL of the tid's contigs.
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			int nR16=0, nR23=0, nR5=0, nRother=0, nTrna=0;
			for(Contig c : getContigs(tid)){
				nR16+=c.r16; nR23+=c.r23; nR5+=c.r5; nRother+=c.rother; nTrna+=c.trna;
			}
			nativeNcR16.put(tid, nR16); nativeNcR23.put(tid, nR23); nativeNcR5.put(tid, nR5);
			nativeNcRother.put(tid, nRother); nativeNcTrna.put(tid, nTrna);
		}

		if(enc==ENC_NORM){
			avgCopy=computeAvgCopy(usable);
			System.err.println("enc=norm: computed avgCopyWhenPresent over "+usable.size()+" orgs");
		}
		if(featuresFile!=null){keptRanks=loadRanks(featuresFile);}
		precomputeNStrings();

		final int baseFam=(keptRanks!=null ? keptRanks.length : numFam);
		final int famCols=baseFam*(enc==ENC_TWO ? 2 : 1);
		numInputs=famCols+SHARED_CONTEXT_WIDTH+numPhyla;//formatRow's new context block, not the old raw NUM_GLOBALS dump
		final boolean subnetNcrna=("ncrna".equals(subnetName));
		subnetFamset=("famset".equals(subnetName));
		if(subnetName!=null && !subnetNcrna && !subnetFamset){
			throw new RuntimeException("Unknown subnet="+subnetName+" (ncrna|famset)");
		}
		final boolean subnet=subnetNcrna || subnetFamset;
		if(subnetFamset){
			if(subsetFile==null){throw new RuntimeException("subnet=famset requires subsetfile=<one rank per line>");}
			subsetRanks=loadRanks(subsetFile);
			subsetMask=new boolean[numFam];
			for(int r : subsetRanks){subsetMask[r]=true;}
			lastFamObs=new int[subsetRanks.length];
			// Native subset complement per organism: the subset families' counts summed over
			// ALL the tid's contigs (the famset subnet's denominator target).
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				int sum=0;
				for(Contig c : getContigs(tid)){
					for(int j=0; j<c.famRank.length; j++){if(subsetMask[c.famRank[j]]){sum+=c.famCount[j];}}
				}
				nativeFamTotal.put(tid, sum);
			}
		}
		// Subnet input width: obs block + phylum one-hot + shared context (FROZEN, standard on
		// every net - the old optional domain/hhcaga/genelen blocks are retired, see CTX_N above).
		final int obsCols=(subnetFamset ? subsetRanks.length : NCRNA_OBS);
		subnetInputs=obsCols+numPhyla+SHARED_CONTEXT_WIDTH;
		// Aggregator mode: load manifest, then dense head + buffers.
		if(aggManifestFile!=null){
			loadAggManifest(aggManifestFile);
			if(aggOut==null){throw new RuntimeException("aggmanifest= requires aggout=");}
			denseRanks=computeDenseHead(usable, denseHead);
			cleanFamBuf=new int[numFam];
			rowCtx=new double[CTX_N];
			// dense head(2x) + phylum + shared context + raw ncRNA two-channel(2x5) + mapped-fraction(1)
			numAggInputs=aggSubnets.size()*4+1+2*denseRanks.length+numPhyla
				+SHARED_CONTEXT_WIDTH+2*NCRNA_OBS+1;
			System.err.println("agg: "+aggSubnets.size()+" subnets, dense head "+denseRanks.length
				+", numAggInputs="+numAggInputs+", obs="+(aggObsServe ? "serve" : "clean"));
		}
		if(rowCtx==null){rowCtx=new double[CTX_N];}//non-aggregator runs still need it for formatRow/subnet rows
		// poolmode=valsplit: both output sets come from the ORIGINAL val orgs (never seen
		// by any subnet trained on the seed-matched train side): first half = aggregator-train
		// (B), second half = final-test (C). Stacking discipline for free (Barbara).
		// poolmode=allbutc: train pool = EVERY usable org except C (Brian 2026-08-11: "hold out
		// vectors, never organisms" - 49 orgs starved the aggregator; the same C stays out so the
		// novel-org readout remains comparable across modes). C is IDENTICAL to valsplit's C.
		if(poolMode==POOL_VALSPLIT || poolMode==POOL_ALLBUTC){
			final int half=valTids.size()/2;
			if(half<1){throw new RuntimeException("poolmode needs >=2 val orgs, have "+valTids.size());}
			final IntList b=subList(valTids, 0, half);
			final IntList c=subList(valTids, half, valTids.size());
			if(poolMode==POOL_VALSPLIT){
				trainTids=b;
			}else{
				trainTids.addAll(b);//all usable orgs except C
			}
			valTids=c;
			System.err.println("poolmode="+(poolMode==POOL_VALSPLIT ? "valsplit" : "allbutc")
				+": aggregator-train orgs="+trainTids.size()+", final-test orgs="+c.size());
		}
		if(benchMode){
			if(benchTruthFile==null || benchManifestFile==null){
				throw new RuntimeException("benchmark mode needs benchtruth= and benchmanifest=");
			}
			if(poolMode!=POOL_ALLBUTC && poolMode!=POOL_VALSPLIT){
				throw new RuntimeException("benchmark mode requires poolmode=allbutc (a defined held-out pool C)");
			}
			writeBench(valTids, benchBins, new Random(seed*3+7));
			System.err.println("done (benchmark).");
			return;
		}
		writeSet(out, (subnet ? subnetOut : null), aggOut, trainTids, n, new Random(seed*2+1), "train");
		if(outval!=null && valn>0 && !valTids.isEmpty()){
			writeSet(outval, (subnet ? subnetValOut : null), aggValOut, valTids, valn, new Random(seed*2+2), "val");
		}
		System.err.println("done.");
	}

	/** Replicates java.util.Collections.shuffle(List,Random)'s exact algorithm and RNG
	 *  call sequence (Fisher-Yates, i from size down to 2, swap(i-1, rnd.nextInt(i))) -
	 *  determinism-critical: the global train/val organism split depends on drawing the
	 *  SAME sequence of rnd.nextInt() calls the original produced. */
	private static void shuffleInPlace(IntList list, Random rnd){
		for(int i=list.size(); i>1; i--){
			final int j=rnd.nextInt(i);
			final int tmp=list.get(i-1);
			list.set(i-1, list.get(j));
			list.set(j, tmp);
		}
	}

	private static IntList subList(IntList src, int from, int to){
		final IntList out=new IntList(Math.max(1, to-from));
		for(int i=from; i<to; i++){out.add(src.get(i));}
		return out;
	}

	/** Builds family->tids index within a pool (for same-family contaminant selection). */
	private HashMap<Integer,IntList> familyIndex(IntList pool){
		final HashMap<Integer,IntList> m=new HashMap<Integer,IntList>();
		for(int i=0; i<pool.size(); i++){
			final int tid=pool.get(i);
			final int fam=tid2family.get(tid);//IntHashMap.get returns -1 when absent, matching getOrDefault(tid,-1)
			IntList l=m.get(fam);
			if(l==null){m.put(fam, l=new IntList());}
			l.add(tid);
		}
		return m;
	}

	private void writeSet(String file, String subnetFile, String aggFile, IntList pool, long count, Random rnd, String tag){
		final HashMap<Integer,IntList> famIdx=familyIndex(pool);
		final ByteStreamWriter bsw=new ByteStreamWriter(file, true, false, true);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(numInputs*4+64);
		bb.append("#dims\t").append(numInputs).append("\t2\t0").nl();
		bsw.print(bb); bb.clear();

		final ByteStreamWriter sbsw=(subnetFile==null ? null : new ByteStreamWriter(subnetFile, true, false, true));
		final ByteBuilder bb2=(sbsw==null ? null : new ByteBuilder(subnetInputs*4+64));
		if(sbsw!=null){
			sbsw.start();
			bb2.append("#dims\t").append(subnetInputs).append("\t1\t0").nl();
			sbsw.print(bb2); bb2.clear();
		}

		final ByteStreamWriter absw=(aggFile==null ? null : new ByteStreamWriter(aggFile, true, false, true));
		final ByteBuilder bb3=(absw==null ? null : new ByteBuilder(numAggInputs*4+64));
		if(absw!=null){
			absw.start();
			bb3.append("#dims\t").append(numAggInputs).append("\t2\t0").nl();
			absw.print(bb3); bb3.clear();
		}

		// Per-attempt buffers hoisted OUT of the loop and reused (S3): fam[]/glob[]/labels[]
		// are cleared, never reallocated. fam[] must be explicitly zeroed every attempt -
		// Agg.add() only INCREMENTS specific indices, it never zeroes the whole array, so a
		// stale value from a REJECTED prior attempt would otherwise bleed into the next one.
		// Agg itself (S4) is built ONCE and reset() between attempts, so its scratchArr/lens
		// buffers (already designed to reuse-if-large-enough) actually get to amortize.
		final double[] labels=new double[2];
		final int[] fam=new int[numFam];
		final double[] glob=new double[NUM_GLOBALS];
		final Agg agg=new Agg();
		agg.setFam(fam);

		long made=0, tries=0;
		while(made<count && tries<count*20+1000){
			tries++;
			Arrays.fill(fam, 0, numFam, 0);
			agg.reset();
			final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg);
			if(targetPhylumIdx<0){continue;}
			formatRow(bb, fam, glob, targetPhylumIdx, labels, agg);
			bsw.print(bb); bb.clear();
			if(sbsw!=null){
				if(subnetFamset){formatFamsetRow(bb2, glob, targetPhylumIdx, agg);}
				else{formatNcrnaRow(bb2, glob, targetPhylumIdx, agg);}
				sbsw.print(bb2); bb2.clear();
			}
			if(absw!=null){
				formatAggRow(bb3, fam, glob, targetPhylumIdx, labels, agg);
				absw.print(bb3); bb3.clear();
			}
			made++;
			if((made%50000)==0){System.err.println(tag+": "+made+"/"+count);}
		}
		bsw.poisonAndWait();
		if(sbsw!=null){sbsw.poisonAndWait(); System.err.println(tag+": wrote "+made+" ncRNA-subnet rows to "+subnetFile);}
		if(absw!=null){absw.poisonAndWait(); System.err.println(tag+": wrote "+made+" aggregator rows to "+aggFile);}
		System.err.println(tag+": wrote "+made+" rows to "+file+" (tries="+tries+")");
	}

	/** Benchmark generation: draws {@code count} synthetic bins from the held-out pool (C in allbutc
	 *  mode) using the SAME makeBin sampler as training, and emits (1) a truth table
	 *  (binID, tid, completeness, contamination, totalBp, nContigs, nForeign), (2) a contig-name
	 *  manifest (binID, contig, native|foreign) for splitting the shred FASTA into per-bin FASTAs,
	 *  and (3) optional aggregator vectors so our net scores the IDENTICAL bins CheckM does. The
	 *  truth labels are the makeBin achieved labels (completeness=cleanBp/gsize, contamination=
	 *  foreignBp/totalBp) - exact ground truth by construction. */
	private void writeBench(IntList pool, long count, Random rnd){
		assert(benchMode) : "writeBench called outside benchmark mode";
		assert(pool!=null && !pool.isEmpty()) : "benchmark pool empty (need held-out orgs; use poolmode=allbutc)";
		final HashMap<Integer,IntList> famIdx=familyIndex(pool);
		final ByteStreamWriter tsw=new ByteStreamWriter(benchTruthFile, true, false, true); tsw.start();
		final ByteStreamWriter msw=new ByteStreamWriter(benchManifestFile, true, false, true); msw.start();
		final ByteStreamWriter vsw=(benchVecFile==null ? null : new ByteStreamWriter(benchVecFile, true, false, true));
		final ByteBuilder tb=new ByteBuilder(256), mb=new ByteBuilder(256);
		final ByteBuilder vb=(vsw==null ? null : new ByteBuilder(numAggInputs*4+64));
		tb.append("#binID\ttid\tcompleteness\tcontamination\ttotalBp\tnContigs\tnForeign").nl(); tsw.print(tb); tb.clear();
		mb.append("#binID\tcontig\trole").nl(); msw.print(mb); mb.clear();
		if(vsw!=null){vsw.start(); vb.append("#dims\t").append(numAggInputs).append("\t2\t0").nl(); vsw.print(vb); vb.clear();}

		final double[] labels=new double[2];
		final int[] fam=new int[numFam];
		final double[] glob=new double[NUM_GLOBALS];
		final Agg agg=new Agg(); agg.setFam(fam);
		long made=0, tries=0;
		while(made<count && tries<count*20+1000){
			tries++;
			Arrays.fill(fam, 0, numFam, 0); agg.reset();
			benchNative.clear(); benchForeign.clear();
			final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg);
			if(targetPhylumIdx<0){continue;}
			final int nContigs=benchNative.size()+benchForeign.size();
			assert(nContigs>0) : "bin passed makeBin with 0 collected contigs (bin"+made+", tid "+lastTarget+")";
			assert(labels[0]>=0 && labels[0]<=1.0001) : "completeness out of range: "+labels[0];
			assert(labels[1]>=0 && labels[1]<=CONT_MAX+1e-9) : "contamination "+labels[1]+" > CONT_MAX "+CONT_MAX;
			final String binID="bin"+made;
			tb.append(binID).tab().append(lastTarget).tab();
			appendFmt(tb, labels[0]); tb.tab(); appendFmt(tb, labels[1]); tb.tab();
			tb.append(lastTotalBp).tab().append(nContigs).tab().append(benchForeign.size()).nl();
			tsw.print(tb); tb.clear();
			for(final Contig c : benchNative){mb.append(binID).tab().append(c.name).append("\tnative").nl();}
			for(final Contig c : benchForeign){mb.append(binID).tab().append(c.name).append("\tforeign").nl();}
			msw.print(mb); mb.clear();
			if(vsw!=null){formatAggRow(vb, fam, glob, targetPhylumIdx, labels, agg); vsw.print(vb); vb.clear();}
			made++;
		}
		tsw.poisonAndWait(); msw.poisonAndWait(); if(vsw!=null){vsw.poisonAndWait();}
		System.err.println("bench: wrote "+made+" bins (tries="+tries+") to "+benchTruthFile+" + "+benchManifestFile
			+(benchVecFile==null ? "" : " + "+benchVecFile));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Bin sampler         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Synthesizes one bin into the provided fam[]/glob[] accumulators and labels[].
	 * @return the target's phylum index, or -1 if the attempt failed (retry).
	 */
	private int makeBin(IntList pool, HashMap<Integer,IntList> famIdx,
			Random rnd, int[] fam, double[] glob, double[] labels, Agg agg){
		final int target=pool.get(rnd.nextInt(pool.size()));
		final long gsize=genomeSize.get(target);
		final long recov=recoverable.get(target);
		if(gsize<=0 || recov<=0){return -1;}

		// Bin-type mixture (isolate/high-quality spike, Brian 2026-08-24): a fraction of bins are drawn
		// PERFECT (comp=1, contam=0) or NEAR-PERFECT (comp 0.90-1.0, contam 0-0.05) so the net gets sharp at
		// the low-contamination/high-completeness end (where CheckM2 beats us; the isolate-QC regime). When
		// both fracs are 0 the roll is skipped, so the RNG sequence - and all existing output - is unchanged.
		final boolean spike=(perfectFrac>0 || nearPerfectFrac>0);
		final double binTypeRoll=(spike ? rnd.nextDouble() : 1.0);
		final boolean perfectBin=(binTypeRoll<perfectFrac);
		final boolean nearPerfectBin=(!perfectBin && binTypeRoll<perfectFrac+nearPerfectFrac);
		// sampled completeness (flat + sqrt-high mixture), or the spiked high-completeness box
		final double comp=perfectBin ? 1.0 : (nearPerfectBin ? 0.90+0.10*rnd.nextDouble() : sampleComp(rnd));
		final long targetBp=(long)(comp*gsize);
		final long cleanBp=selectContigs(getContigs(target), targetBp, rnd, agg, benchMode?benchNative:null);
		if(cleanBp<=0){return -1;}
		// Snapshot the TARGET's observed ncRNA (before contaminants are added) and its id, so a
		// subnet emitter can pair this bin's observed ncRNA with the target's native complement.
		// Read-only w.r.t. the global vector (agg is only inspected).
		lastTarget=target;
		lastNcObs[0]=agg.r16; lastNcObs[1]=agg.r23; lastNcObs[2]=agg.r5;
		lastNcObs[3]=agg.rother; lastNcObs[4]=agg.trna;
		// Same snapshot for a famset subnet: the subset families' observed counts, target-only
		// (fam[] holds only target contributions here; contaminants are added below).
		if(subsetRanks!=null){
			for(int i=0; i<subsetRanks.length; i++){lastFamObs[i]=fam[subsetRanks[i]];}
		}
		// Full target-only fam snapshot for aggregator aggobs=clean (fam[] gains
		// contaminant counts below; this preserves the pre-contaminant state).
		if(cleanFamBuf!=null){System.arraycopy(fam, 0, cleanFamBuf, 0, numFam);}

		// sampled contamination: 0 for perfect bins, U[0,0.05] for near-perfect, else clean-spike/square-low
		final double cont=perfectBin ? 0.0 : (nearPerfectBin ? 0.05*rnd.nextDouble()
			: (rnd.nextDouble()<cleanSpike ? 0.0 : sampleCont(rnd)));
		long foreignBp=0;
		if(cont>0){
			final long foreignTarget=(long)((cont/(1.0-cont))*cleanBp);
			if(foreignTarget>0){
				int nContam=1;
				if(rnd.nextDouble()<multiContamProb){nContam=2+rnd.nextInt(2);}//2 or 3
				final long per=Math.max(1, foreignTarget/nContam);
				for(int k=0; k<nContam; k++){
					final int ctid=pickContaminant(pool, famIdx, target, rnd);
					if(ctid<0){break;}
					foreignBp+=selectContigs(getContigs(ctid), per, rnd, agg, benchMode?benchForeign:null);
				}
			}
		}

		final long totalBp=cleanBp+foreignBp;
		if(totalBp<=0 || agg.contigs<=0){return -1;}
		lastTotalBp=totalBp;
		// Whole-bin (serve-faithful) ncRNA observed, contaminants included - what a
		// deployed bin actually shows. glob[] lacks rother, so capture all 5 here.
		lastNcServe[0]=agg.r16; lastNcServe[1]=agg.r23; lastNcServe[2]=agg.r5;
		lastNcServe[3]=agg.rother; lastNcServe[4]=agg.trna;

		// achieved labels from the explicit target
		labels[0]=Math.min(1.0, cleanBp/(double)gsize);      // completeness
		labels[1]=foreignBp/(double)totalBp;                 // contamination
		// Whole-contig overshoot can push foreign past clean, flipping the dominant organism
		// (Barbara #4). Such a bin's contamination relative to the chosen target is out of the
		// 0-50% spec and ambiguous vs a bin-grader, so reject and retry.
		if(labels[1]>CONT_MAX){return -1;}

		// global stats
		glob[0]=log2(totalBp);
		glob[1]=log2(agg.contigs);
		glob[2]=log2(agg.l50());
		glob[3]=agg.acgt>0 ? agg.gc/(double)agg.acgt : 0;
		glob[4]=totalBp>0 ? agg.coding/(double)totalBp : 0;
		glob[5]=agg.cds>0 ? agg.glenSum/(double)agg.cds : 0;
		glob[6]=agg.geneStd();
		glob[7]=agg.cds>0 ? agg.mapped/(double)agg.cds : 0;
		glob[8]=log2(agg.richness());
		glob[9]=agg.r16; glob[10]=agg.r23; glob[11]=agg.r5; glob[12]=agg.trna;

		final Integer pi=tid2phylumIdx0(target);
		return pi==null ? phylumIndex.get("other") : pi;
	}

	/** tid2phylumIdx.get() returns int (-1 sentinel); wraps as Integer only at this one
	 *  low-frequency (once-per-bin) call site to keep makeBin's return contract unchanged. */
	private Integer tid2phylumIdx0(int target){
		final int v=tid2phylumIdx.get(target);
		return v<0 ? null : Integer.valueOf(v);
	}

	/** Picks a contaminant tid != target, favoring the target's family. */
	private int pickContaminant(IntList pool, HashMap<Integer,IntList> famIdx,
			int target, Random rnd){
		if(!tid2family.isEmpty() && rnd.nextDouble()<sameFamProb){
			final int fam=tid2family.get(target);
			final IntList l=famIdx.get(fam);
			if(l!=null && l.size()>1){
				for(int t=0; t<8; t++){final int c=l.get(rnd.nextInt(l.size())); if(c!=target){return c;}}
			}
		}
		for(int t=0; t<8; t++){final int c=pool.get(rnd.nextInt(pool.size())); if(c!=target){return c;}}
		return -1;
	}

	/** Randomly adds contigs (by shuffled draw) until reaching targetBp; returns bp added.
	 *  When {@code collect!=null} (benchmark mode), the selected Contigs are appended to it so a
	 *  FASTA manifest can be emitted - the sampling itself is unchanged (collect is inspect-only). */
	private long selectContigs(ArrayList<Contig> contigs, long targetBp, Random rnd, Agg agg, ArrayList<Contig> collect){
		final int nc=contigs.size();
		final int[] order=agg.scratch(nc);
		for(int i=0; i<nc; i++){order[i]=i;}
		// partial Fisher-Yates: shuffle enough to draw without replacement
		long added=0;
		for(int i=0; i<nc && added<targetBp; i++){
			final int j=i+rnd.nextInt(nc-i);
			final int tmp=order[i]; order[i]=order[j]; order[j]=tmp;
			final Contig c=contigs.get(order[i]);
			agg.add(c);
			if(collect!=null){collect.add(c);}
			added+=c.length;
		}
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Aggregator          ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulates per-bin sufficient statistics over selected contigs. Reused (S4) across
	 *  attempts via reset() instead of being constructed fresh per makeBin() call - this is
	 *  what lets scratch()'s reuse-if-large-enough buffer actually amortize as designed. */
	static final class Agg {
		int[] fam;
		long gc, acgt, coding, glenSum, glenSq;
		int contigs, cds, mapped, r16, r23, r5, rother, trna;
		int[] lens=new int[64]; int nlens=0;
		final int[] dimer=new int[16];//summed dinucleotide counts -> bin-faithful HH/CAGA (null-dimer contigs skipped)
		int[] scratchArr;
		void setFam(int[] fam_){fam=fam_;}
		void reset(){
			gc=acgt=coding=glenSum=glenSq=0;
			contigs=cds=mapped=r16=r23=r5=rother=trna=0;
			nlens=0;
			Arrays.fill(dimer, 0);
		}
		int[] scratch(int need){
			if(scratchArr==null || scratchArr.length<need){scratchArr=new int[need];}
			return scratchArr;
		}
		void add(Contig c){
			contigs++; gc+=c.gc; acgt+=c.acgt; coding+=c.coding; cds+=c.cds; mapped+=c.mapped;
			glenSum+=c.glenSum; glenSq+=c.glenSq; r16+=c.r16; r23+=c.r23; r5+=c.r5; rother+=c.rother; trna+=c.trna;
			for(int i=0; i<c.famRank.length; i++){fam[c.famRank[i]]+=c.famCount[i];}
			if(c.dimer!=null){for(int i=0; i<16; i++){dimer[i]+=c.dimer[i];}}
			if(nlens>=lens.length){lens=Arrays.copyOf(lens, lens.length*2);}
			lens[nlens++]=c.length;
		}
		double geneStd(){
			if(cds<=0){return 0;}
			double mean=glenSum/(double)cds;
			double var=glenSq/(double)cds-mean*mean;
			return var>0 ? Math.sqrt(var) : 0;
		}
		/** Bin-faithful homopolymer-dimer fraction from the SUMMED dinucleotide counts (additive; never
		 *  average per-contig ratios). Zero when no dimer data is present (KmerTracker guards the denominator). */
		float hh(){return KmerTracker.HH(dimer);}
		/** Bin-faithful CAGA strand-skew metric from the summed dinucleotide counts. */
		float caga(){return KmerTracker.CAGA(dimer);}
		int richness(){int r=0; for(int v : fam){if(v>0){r++;}} return r;}
		long l50(){
			if(nlens<=0){return 1;}
			int[] a=Arrays.copyOf(lens, nlens);
			Arrays.sort(a);
			long total=0; for(int v : a){total+=v;}
			long half=total/2, cum=0;
			for(int i=a.length-1; i>=0; i--){cum+=a[i]; if(cum>=half){return a[i];}}
			return a[0];
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Distributions        ----------------*/
	/*--------------------------------------------------------------*/

	// completeness in [COMP_MIN,1]: flat + sqrt(flat) (extra density high)
	private double sampleComp(Random rnd){
		double u=rnd.nextDouble();
		double base=(rnd.nextDouble()<mixComp) ? Math.sqrt(u) : u;
		return COMP_MIN+(1.0-COMP_MIN)*base;
	}
	// contamination in [0,CONT_MAX]: flat + square(flat) (extra density low)
	private double sampleCont(Random rnd){
		double u=rnd.nextDouble();
		double base=(rnd.nextDouble()<mixCont) ? u*u : u;
		return CONT_MAX*base;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	private void precomputeNStrings(){
		nOverN1=new String[N_STR_MAX+1];
		//The agg dense head uses the two-channel encoding regardless of the global enc mode.
		if(enc==ENC_TWO || aggManifestFile!=null){excessArr=new String[N_STR_MAX+1];}
		double logCap=Math.log(1+LOG_CAP)/LOG2;
		for(int i=0; i<=N_STR_MAX; i++){
			nOverN1[i]=fmt(encodeCount(i, logCap));
			if(excessArr!=null){int ex=Math.min(Math.max(i-1, 0), EXC_CAP); excessArr[i]=fmt(ex/(double)EXC_CAP);}
		}
	}
	/** Single-column encoding of a family's summed count under the active enc mode. */
	private double encodeCount(int count, double logCap){
		if(enc==ENC_LOG){double v=(Math.log(1+count)/LOG2)/logCap; return v>1 ? 1 : v;}
		if(enc==ENC_RAW){return Math.min(count, RAW_CAP)/(double)RAW_CAP;}
		return count/(double)(1+count);//ENC_RATIO (ENC_TWO handles presence separately)
	}
	private void appendFamStr(ByteBuilder bb, int count){
		if(count<=N_STR_MAX){bb.append(nOverN1[count]);}
		else{appendFmt(bb, encodeCount(count, Math.log(1+LOG_CAP)/LOG2));}
	}
	/** Two-channel family: presence (0/1) then excess-copies min(count-1,cap)/cap. */
	private void appendTwo(ByteBuilder bb, int count){
		bb.append(count>0 ? '1' : '0'); bb.tab();
		if(count<=N_STR_MAX){bb.append(excessArr[count]);}
		else{final int ex=Math.min(count-1, EXC_CAP); appendFmt(bb, ex/(double)EXC_CAP);}
		bb.tab();
	}
	/** Fixed-notation float, no exponent (RegressionTrainer's fast parser rejects 'e').
	 *  String-returning form, KEPT for the aggregator's per-subnet context roundtrip (Float.parseFloat(fmt(v)))
	 *  which genuinely needs a String - low frequency (once per bin, not per family). */
	private static String fmt(double v){
		if(v==(long)v){return Long.toString((long)v);}
		return String.format("%.6f", v);
	}
	/** Zero-allocation equivalent of fmt(), appended directly - MUST replicate fmt()'s
	 *  whole-number shortcut exactly: ByteBuilder's own fast append(double,decimals) always
	 *  prints the requested decimals (String.format-style), but a naive append(v,6) would
	 *  differ from fmt() on whole numbers, which appendSlow's precise-but-slow path does NOT
	 *  special-case either - so the whole-number branch is replicated here explicitly. */
	private static void appendFmt(ByteBuilder bb, double v){
		if(v==(long)v){bb.append((long)v);}
		else{bb.appendSlow(v, 6);}
	}
	private static int parseEnc(String s){
		s=s.toLowerCase();
		if(s.equals("ratio")){return ENC_RATIO;}
		if(s.equals("log") || s.equals("log1p")){return ENC_LOG;}
		if(s.equals("raw") || s.equals("linear")){return ENC_RAW;}
		if(s.equals("two") || s.equals("twochannel")){return ENC_TWO;}
		if(s.equals("norm") || s.equals("avgnorm")){return ENC_NORM;}
		throw new RuntimeException("Unknown enc="+s+" (ratio|raw|log|two|norm)");
	}

	private static boolean parseBool(String s){
		return s==null || s.equals("t") || s.equals("true") || s.equals("1") || s.equals("yes");
	}

	/** Maps a domain string to the 8-way one-hot index [bact,arch,fungi,plant,animal,protist,virus,other]. */
	private static int domainIndex(String d){
		if(d==null){return DOMAIN_OTHER;}
		final String s=d.toLowerCase();
		if(s.startsWith("bacteri")){return 0;}
		if(s.startsWith("archae")){return 1;}
		if(s.startsWith("fung")){return 2;}
		if(s.contains("viridiplant") || s.startsWith("plant")){return 3;}
		if(s.startsWith("metazoa") || s.startsWith("animal")){return 4;}
		if(s.startsWith("protist")){return 5;}
		if(s.startsWith("vir")){return 6;}
		return DOMAIN_OTHER;//incl. bare "eukaryota" (subkingdom needs tax lineage)
	}


	/**
	 * Per-family expected copy number WHEN PRESENT, over the reference organisms:
	 * for each family, the mean of the organism-level count across the organisms where
	 * that family appears (count &gt; 0). A per-family baseline so enc=norm can encode a
	 * present-at-typical-copies family as ~1.0 and a duplicated family as &gt;1.0,
	 * calibrated to each family's own natural copy number (Brian 2026-08-06). This is a
	 * reference-DB constant (not a label), computed over ALL usable orgs for stability.
	 */
	private double[] computeAvgCopy(IntList usable){
		final long[] sum=new long[numFam];
		final int[] present=new int[numFam];
		final int[] orgCount=new int[numFam];
		for(int i=0; i<usable.size(); i++){
			Arrays.fill(orgCount, 0);
			for(Contig c : getContigs(usable.get(i))){
				for(int j=0; j<c.famRank.length; j++){orgCount[c.famRank[j]]+=c.famCount[j];}
			}
			for(int f=0; f<numFam; f++){if(orgCount[f]>0){sum[f]+=orgCount[f]; present[f]++;}}
		}
		final double[] avg=new double[numFam];
		for(int f=0; f<numFam; f++){avg[f]=present[f]>0 ? sum[f]/(double)present[f] : 1.0;}
		return avg;
	}
	/** enc=norm: 0 if absent, else count / avgCopyWhenPresent[rank]. */
	private void appendNormStr(ByteBuilder bb, int count, int rank){
		if(count==0){bb.append('0'); return;}
		final double a=avgCopy[rank];
		appendFmt(bb, a>0 ? count/a : count);
	}

	/**
	 * Corpus-scale reference for the avg~0.5 normalized context features (Brian 2026-08-24). Over the
	 * usable orgs, computes the mean of log2(genome bp), log2(1+total CDS), mean gene length, and
	 * gene-length std; each such feature is later emitted as {@code 0.5*value/mean} so the population
	 * centers near 0.5 (bounded, well-scaled net inputs instead of huge raw magnitudes). This is an
	 * ORG-level reference (stable and deterministic); partial bins center a little below 0.5, which is
	 * fine - the goal is bounded scaling, not exact 0.5. These 4 scales are PERSISTED with the trained
	 * nets so MagQCScore and the production tool reproduce the identical encoding (persistence added in
	 * the one-pass emit rewrite). size/#genes are scaled on log2 (wide dynamic range); gene-length
	 * mean/std on the raw value.
	 */
	private void computeNormScales(IntList usable){
		double sumLogBp=0, sumLogCds=0, sumGlen=0, sumGlenStd=0;
		int nBp=0, nCds=0, nGlen=0;
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			final long bp=genomeSize.get(tid);
			if(bp>0){sumLogBp+=log2(bp); nBp++;}
			long cds=0, glenSum=0, glenSq=0;
			for(Contig c : getContigs(tid)){cds+=c.cds; glenSum+=c.glenSum; glenSq+=c.glenSq;}
			sumLogCds+=log2(1+cds); nCds++;
			if(cds>0){
				final double mean=glenSum/(double)cds;
				final double var=glenSq/(double)cds-mean*mean;
				sumGlen+=mean; sumGlenStd+=(var>0 ? Math.sqrt(var) : 0); nGlen++;
			}
		}
		scaleLogBp=(nBp>0 && sumLogBp>0 ? sumLogBp/nBp : 1);
		scaleLogCds=(nCds>0 && sumLogCds>0 ? sumLogCds/nCds : 1);
		scaleGlen=(nGlen>0 && sumGlen>0 ? sumGlen/nGlen : 1);
		scaleGlenStd=(nGlen>0 && sumGlenStd>0 ? sumGlenStd/nGlen : 1);
		System.err.println("normScales (avg~0.5 reference): mean_log2bp="+scaleLogBp+" mean_log2cds="+scaleLogCds
			+" mean_glen="+scaleGlen+" mean_glenStd="+scaleGlenStd);
	}

	private void formatRow(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg){
		if(enc==ENC_TWO){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){appendTwo(bb, fam[keptRanks[i]]);}
			}else{
				for(int i=0; i<numFam; i++){appendTwo(bb, fam[i]);}
			}
		}else if(enc==ENC_NORM){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){final int r=keptRanks[i]; appendNormStr(bb, fam[r], r); bb.tab();}
			}else{
				for(int i=0; i<numFam; i++){appendNormStr(bb, fam[i], i); bb.tab();}
			}
		}else{
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){appendFamStr(bb, fam[keptRanks[i]]); bb.tab();}
			}else{
				for(int i=0; i<numFam; i++){appendFamStr(bb, fam[i]); bb.tab();}
			}
		}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		appendFmt(bb, labels[0]); bb.tab();
		appendFmt(bb, labels[1]); bb.nl();
	}

	/**
	 * Fills ctx[0..CTX_N) with the FROZEN shared-context scalar block for the current bin
	 * (magqc_rebuild_20260824.plan "FROZEN VECTOR LAYOUT", Brian-confirmed 2026-08-24):
	 * size/#genes/gene-length mean+std at the persisted avg~0.5 corpus scale (computeNormScales),
	 * GC/coding-fraction/log2(richness) unchanged, HH/CAGA BIN-FAITHFUL from agg's SUMMED
	 * dinucleotide counts (additive discipline - never averaged per-contig; KmerTracker.HH/CAGA
	 * guard the zero-denominator case). Shared by every emit method so every net - subnets,
	 * aggregator, and the Tier-A monolith - sees byte-identical context.
	 */
	private void computeContext(double[] ctx, double[] glob, Agg agg){
		ctx[CTX_SIZE]=0.5*glob[0]/scaleLogBp;
		ctx[CTX_GC]=glob[3];
		ctx[CTX_CODING]=glob[4];
		ctx[CTX_GENES]=0.5*log2(1+agg.cds)/scaleLogCds;
		ctx[CTX_L2RICH]=glob[8];
		ctx[CTX_GLEN]=0.5*glob[5]/scaleGlen;
		ctx[CTX_GLENSTD]=0.5*glob[6]/scaleGlenStd;
		ctx[CTX_HH]=agg.hh();
		ctx[CTX_CAGA]=agg.caga();
	}

	/** Appends the FROZEN shared-context block for tid: the CTX_N scalars (already computed
	 *  into ctx by {@link #computeContext}) followed by the domain one-hot(8) - shared-context
	 *  item 10, categorical so it isn't part of the scalar array. */
	private void appendContext(ByteBuilder bb, double[] ctx, int tid){
		for(int i=0; i<CTX_N; i++){appendFmt(bb, ctx[i]); bb.tab();}
		final int di=domainIdxOf(tid);
		for(int i=0; i<DOMAINS; i++){bb.append(i==di ? '1' : '0'); bb.tab();}
	}

	/** Two-channel encoding for an aggregator-direct raw count (ncRNA types): presence(0/1)
	 *  then encodeCount(count) under the active enc mode - distinct from the family two-channel
	 *  {@link #appendTwo} (presence+excess-copies), since these aren't duplication-flag counts. */
	private void appendCountTwo(ByteBuilder bb, int count){
		bb.append(count>0 ? '1' : '0'); bb.tab();
		appendFmt(bb, encodeCount(count, Math.log(1+LOG_CAP)/LOG2)); bb.tab();
	}

	/**
	 * Emits one ncRNA-subnet training row for the current bin: the target organism's
	 * observed ncRNA counts (r16,r23,r5,rother,trna) plus the FROZEN shared context (phylum
	 * one-hot + the standard context block) as inputs, and the target's NATIVE ncRNA
	 * complement (summed over its whole genome) as the single regression target. The
	 * subnet thus learns the EXPECTED denominator; completeness = observed/expected is
	 * derived downstream (Barbara's refinement of the subset-relative-label design).
	 */
	private void formatNcrnaRow(ByteBuilder bb, double[] glob, int phylumIdx, Agg agg){
		int obsTotal=0;
		for(int i=0; i<NCRNA_OBS; i++){bb.append(lastNcObs[i]); bb.tab(); obsTotal+=lastNcObs[i];}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		final int nativeTotal=nativeNcR16.get(lastTarget)+nativeNcR23.get(lastTarget)+nativeNcR5.get(lastTarget)
			+nativeNcRother.get(lastTarget)+nativeNcTrna.get(lastTarget);
		// The bin's target contigs are a subset of the genome, so observed<=native always.
		assert(obsTotal<=nativeTotal) : "ncRNA observed "+obsTotal+" > native "+nativeTotal+" (tid "+lastTarget+")";
		bb.append(nativeTotal).nl();
	}

	/**
	 * Emits one famset-subnet training row for the current bin: the subset families'
	 * observed counts (target organism only, snapshotted before contaminants) plus the
	 * SAME FROZEN shared context as the ncRNA subnet, and the target's NATIVE subset
	 * total (the subset families' counts over its whole genome) as the single
	 * regression target. The subset is defined by subsetfile= (one family rank per
	 * line), so per-phylum marker sets and co-occurrence modules train with identical
	 * machinery; evaluate with SubnetRatioScore numobs=(subset size).
	 */
	private void formatFamsetRow(ByteBuilder bb, double[] glob, int phylumIdx, Agg agg){
		int obsTotal=0;
		for(int i=0; i<lastFamObs.length; i++){bb.append(lastFamObs[i]); bb.tab(); obsTotal+=lastFamObs[i];}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		final int nativeTotal=nativeFamTotal.get(lastTarget);
		assert(nativeTotal>=0) : "nativeFamTotal missing for tid "+lastTarget;
		// The bin's target contigs are a subset of the genome, so observed<=native always.
		assert(obsTotal<=nativeTotal) : "famset observed "+obsTotal+" > native "+nativeTotal+" (tid "+lastTarget+")";
		bb.append(nativeTotal).nl();
	}

	/**
	 * Emits one aggregator training row for the current bin. For every manifest
	 * subnet: gathers its observed counts (whole-bin by default - serve-faithful,
	 * contaminants included; target-only under aggobs=clean), builds the subnet's
	 * input exactly as its training rows were built (same column order, same fmt()
	 * rounding so in-process values match what a file round-trip would deliver),
	 * feed-forwards the frozen net, and emits [ratio, log-obs, log-pred, zero-flag].
	 * Then the pooled ratio baseline, the raw dense head (enc=two presence+excess of
	 * the top-K prevalent families), phylum one-hot, the FROZEN shared context, the
	 * raw ncRNA two-channel(5x2) direct aggregator input, mapped-fraction, and the
	 * GLOBAL comp/contam targets.
	 */
	private void formatAggRow(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg){
		final int[] famArr=(aggObsServe ? fam : cleanFamBuf);
		final int[] ncArr=(aggObsServe ? lastNcServe : lastNcObs);
		// Per-bin context, computed ONCE (shared by every subnet AND the aggregator tail - every
		// net now sees the identical FROZEN context, the old per-net rep flags are retired).
		computeContext(rowCtx, glob, agg);
		final int domainIdx=domainIdxOf(lastTarget);

		long sumObs=0;
		double sumPred=0;
		for(AggSubnet s : aggSubnets){
			// Input layout mirrors formatFamsetRow/formatNcrnaRow exactly: obs, phylum one-hot,
			// shared context (CTX_N scalars + domain one-hot).
			int p=0;
			long obsTotal=0;
			if(s.ncrna){
				for(int i=0; i<NCRNA_OBS; i++){s.buf[p++]=ncArr[i]; obsTotal+=ncArr[i];}
			}else{
				for(int r : s.ranks){final int v=famArr[r]; s.buf[p++]=v; obsTotal+=v;}
			}
			for(int i=0; i<numPhyla; i++){s.buf[p++]=(i==phylumIdx ? 1 : 0);}
			// Round-trip through fmt()'s exact string representation, same as the training-file
			// write+reparse a subnet actually saw (appendContext below writes via appendFmt) -
			// keeps this in-process feed bit-identical to what the frozen net was trained on.
			for(int i=0; i<CTX_N; i++){s.buf[p++]=Float.parseFloat(fmt(rowCtx[i]));}
			for(int i=0; i<DOMAINS; i++){s.buf[p++]=(i==domainIdx ? 1 : 0);}
			assert(p==s.buf.length) : s.name+": filled "+p+" of "+s.buf.length;
			s.net.applyInput(s.buf);
			s.net.feedForward();
			final double pred=s.net.getOutput(0);
			final double ratio=Math.min(RATIO_CAP, obsTotal/Math.max(0.5, pred));
			appendFmt(bb, ratio); bb.tab();
			appendFmt(bb, log2(1+obsTotal)); bb.tab();
			appendFmt(bb, log2(1+Math.max(0, pred))); bb.tab();
			bb.append(obsTotal==0 ? '1' : '0'); bb.tab();
			sumObs+=obsTotal;
			sumPred+=Math.max(0, pred);
		}
		appendFmt(bb, Math.min(RATIO_CAP, sumObs/Math.max(1, sumPred))); bb.tab();
		// Dense head is ALWAYS whole-bin (the aggregator's direct deployment signal, not a
		// subnet input; aggobs= only A/Bs the subnet-obs question - settled with Eru 2026-08-11).
		for(int r : denseRanks){appendTwo(bb, fam[r]);}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		appendContext(bb, rowCtx, lastTarget);
		// Aggregator-only direct inputs (FROZEN VECTOR LAYOUT): raw ncRNA two-channel (5 types x
		// [presence, encodeCount]) then mapped-fraction. Whole-bin/target-only follows the same
		// aggobs= switch as the subnet inputs above (ncArr), for serve-faithfulness consistency.
		for(int i=0; i<NCRNA_OBS; i++){appendCountTwo(bb, ncArr[i]);}
		appendFmt(bb, glob[7]); bb.tab();
		appendFmt(bb, labels[0]); bb.tab();
		appendFmt(bb, labels[1]); bb.nl();
	}

	private int domainIdxOf(int tid){
		final int v=tidToIdx.get(tid);
		return v<0 || v>=domainIdxArr.length ? DOMAIN_OTHER : domainIdxArr[v];
	}

	private static double log2(double v){return v<=1 ? 0 : Math.log(v)/LOG2;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String cacheFile, sizemapFile, familyFile, taxpgmFile, treeFile, out, outval, featuresFile;
	private String subnetName, subnetOut, subnetValOut, subsetFile;
	private String aggManifestFile, aggOut, aggValOut;
	private int denseHead=100;
	private boolean aggObsServe=true;
	private int poolMode=POOL_TRAINVAL;
	private static final int POOL_TRAINVAL=0, POOL_VALSPLIT=1, POOL_ALLBUTC=2;
	private ArrayList<AggSubnet> aggSubnets;
	private int[] denseRanks;
	private int[] cleanFamBuf;
	// FROZEN shared-context scratch buffer (CTX_N scalars), reused across every format*Row call
	// for the current bin - computeContext() overwrites it fully each time, no partial-update risk.
	private double[] rowCtx;
	private int numAggInputs;
	private final int[] lastNcServe=new int[5];
	private boolean subnetFamset=false;
	private int[] subsetRanks;
	private boolean[] subsetMask;
	private int[] lastFamObs;
	private final IntHashMap nativeFamTotal=new IntHashMap();
	private int subnetInputs;
	private int lastTarget;
	private long lastTotalBp;
	private final int[] lastNcObs=new int[5];
	private int[] keptRanks;
	private long n=400000, valn=40000, seed=1;
	private double valfrac=0.10, mixComp=0.5, mixCont=0.5, cleanSpike=0.15, multiContamProb=0.15, sameFamProb=0.70;
	private int minlen=0;
	private int enc=ENC_RATIO;

	private int numFam=8000, numPhyla, numInputs;
	private TaxTree tree;

	// byTid replacement: dense tid->index (IntHashMap can't hold list values) + a
	// parallel growable list-of-lists indexed by that dense int.
	private final IntHashMap tidToIdx=new IntHashMap();
	private final ArrayList<ArrayList<Contig>> contigLists=new ArrayList<ArrayList<Contig>>();
	private int[] domainIdxArr=new int[64];//parallel to tidToIdx's dense index (was tid2domainIdx)

	private final IntLongHashMap genomeSize=new IntLongHashMap();
	// Benchmark mode (off unless benchtruth=/benchmanifest= given): emit held-out synthetic bins as a
	// truth table + contig-name manifest (+ optional aggregator vectors) so CheckM1/CheckM2 and our net
	// score the IDENTICAL bins. Additive - existing output paths and the differential gate are untouched.
	private boolean benchMode=false;
	private String benchManifestFile=null, benchTruthFile=null, benchVecFile=null;
	private long benchBins=500;
	private final ArrayList<Contig> benchNative=new ArrayList<Contig>();
	private final ArrayList<Contig> benchForeign=new ArrayList<Contig>();
	// Isolate/high-quality training spike (Brian 2026-08-24): fraction of bins drawn perfect / near-perfect.
	private double perfectFrac=0.0, nearPerfectFrac=0.0;
	private final IntLongHashMap recoverable=new IntLongHashMap();
	private final IntHashMap nativeNcR16=new IntHashMap(), nativeNcR23=new IntHashMap(), nativeNcR5=new IntHashMap(),
		nativeNcRother=new IntHashMap(), nativeNcTrna=new IntHashMap();
	private final IntHashMap tid2phylumIdx=new IntHashMap();
	private final IntHashMap tid2family=new IntHashMap();
	private final HashMap<String,Integer> phylumIndex=new HashMap<String,Integer>();
	private ArrayList<String> phylumList;
	private String[] nOverN1;
	private String[] excessArr;
	private double[] avgCopy;
	// avg~0.5 normalization scales (corpus reference over usable orgs; persisted with the nets). Wired into
	// the emit methods in the one-pass vector-layout rewrite; computed by computeNormScales().
	private double scaleLogBp=1, scaleLogCds=1, scaleGlen=1, scaleGlenStd=1;

	private static final int NUM_GLOBALS=13;
	private static final int NCRNA_OBS=5;
	private static final int DOMAINS=8, DOMAIN_OTHER=7;
	private static final int N_STR_MAX=4096;
	private static final int ENC_RATIO=0, ENC_LOG=1, ENC_RAW=2, ENC_TWO=3, ENC_NORM=4;
	private static final int RAW_CAP=32, LOG_CAP=64, EXC_CAP=16;
	private static final double COMP_MIN=0.10, CONT_MAX=0.50, LOG2=Math.log(2);
	/** Cap for the per-subnet and pooled obs/pred ratios (duplication saturates at 2x). */
	private static final double RATIO_CAP=2.0;
	/** Indices into the FROZEN shared-context scalar block (Brian 2026-08-24 vector-layout
	 *  rebuild, magqc_rebuild_20260824.plan "FROZEN VECTOR LAYOUT"): 9 scalars, standard on
	 *  EVERY net (subnets + aggregator + monolith) - the old per-net optional sndomain/
	 *  snhhcaga/sngenelen/snbinscaled/sncodingaffine representation flags are RETIRED, there
	 *  is now exactly one encoding per feature. Domain one-hot(8) is the 10th shared-context
	 *  item but is categorical, not scalar, so it is appended separately right after this
	 *  block (see appendContext) rather than living in this array. */
	private static final int CTX_SIZE=0, CTX_GC=1, CTX_CODING=2, CTX_GENES=3, CTX_L2RICH=4,
		CTX_GLEN=5, CTX_GLENSTD=6, CTX_HH=7, CTX_CAGA=8, CTX_N=9;
	/** Shared-context wire width: the 9 scalars above plus the domain one-hot(8). */
	private static final int SHARED_CONTEXT_WIDTH=CTX_N+DOMAINS;
	private static final int[] EMPTY=new int[0];
}
