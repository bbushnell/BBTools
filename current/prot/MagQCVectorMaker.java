package prot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeSet;

import tax.TaxNode;
import tax.TaxTree;

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
			else if(a.equals("samefamprob")){sameFamProb=Double.parseDouble(b);}
			else if(a.equals("enc")){enc=parseEnc(b);}
			else if(a.equals("subnet")){subnetName=b.toLowerCase();}
			else if(a.equals("subnetout")){subnetOut=b;}
			else if(a.equals("subnetvalout")){subnetValOut=b;}
			else if(a.equals("sncodingaffine")){snCodingAffine=parseBool(b);}
			else if(a.equals("snbinscaled")){snBinScaled=parseBool(b);}
			else if(a.equals("sndomain")){snDomain=parseBool(b);}
			else if(a.equals("snhhcaga")){snHHCAGA=parseBool(b);}
			else if(a.equals("sngenelen")){snGeneLen=parseBool(b);}
			else if(a.equals("kmerfile")){kmerFile=b;}
			else{System.err.println("Warning: unknown arg "+arg);}
		}
		if(cacheFile==null || sizemapFile==null || taxpgmFile==null || out==null){
			throw new RuntimeException("Required: cache= sizemap= taxpgm= out= [outval= familylist= tree= n= valn= ...]");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Per-contig          ----------------*/
	/*--------------------------------------------------------------*/

	/** One cached contig's sufficient statistics; family counts stored sparsely. */
	static final class Contig {
		int tid, length, gc, acgt, cds, mapped, coding, r16, r23, r5, rother, trna;
		long glenSum, glenSq;
		int[] famRank, famCount;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Loaders           ----------------*/
	/*--------------------------------------------------------------*/

	private void loadAux(){
		try{
			// family list -> count of family columns
			if(familyFile!=null){
				BufferedReader br=new BufferedReader(new FileReader(familyFile));
				String line=br.readLine();//header
				int c=0; while((line=br.readLine())!=null){if(line.length()>0){c++;}}
				br.close(); numFam=c;
			}
			// taxpgm: tid <tab> phylum <tab> pgm  -> tid->phylum, and phylum vocabulary
			BufferedReader br=new BufferedReader(new FileReader(taxpgmFile));
			TreeSet<String> phyla=new TreeSet<String>();
			String line;
			while((line=br.readLine())!=null){
				String[] p=line.split("\t");
				if(p.length>=2){tid2phylum.put(Integer.parseInt(p[0]), p[1]); phyla.add(p[1]);}
			}
			br.close();
			phylumList=new ArrayList<String>(phyla);
			phylumList.add("other");
			for(int i=0; i<phylumList.size(); i++){phylumIndex.put(phylumList.get(i), i);}
			numPhyla=phylumList.size();
			// sizemap: tid <tab> bp
			br=new BufferedReader(new FileReader(sizemapFile));
			while((line=br.readLine())!=null){
				String[] p=line.split("\t");
				if(p.length>=2){genomeSize.put(Integer.parseInt(p[0]), Long.parseLong(p[1]));}
			}
			br.close();
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Loads the per-contig cache into per-tid contig lists. */
	private void loadCache(){
		try{
			BufferedReader br=new BufferedReader(new FileReader(cacheFile));
			String line;
			while((line=br.readLine())!=null){
				if(line.length()==0 || line.charAt(0)=='#'){continue;}
				String[] f=line.split("\t", -1);
				if(f.length<17){continue;}
				Contig c=new Contig();
				c.tid=Integer.parseInt(f[1]);
				c.length=Integer.parseInt(f[3]);
				if(c.length<minlen){continue;}
				c.gc=Integer.parseInt(f[4]);
				c.acgt=Integer.parseInt(f[5]);
				c.cds=Integer.parseInt(f[6]);
				c.mapped=Integer.parseInt(f[7]);
				c.glenSum=Long.parseLong(f[8]);
				c.glenSq=Long.parseLong(f[9]);
				c.coding=Integer.parseInt(f[10]);
				c.r16=Integer.parseInt(f[11]);
				c.r23=Integer.parseInt(f[12]);
				c.r5=Integer.parseInt(f[13]);
				c.rother=Integer.parseInt(f[14]);
				c.trna=Integer.parseInt(f[15]);
				String fc=f[16];
				if(fc.length()>0){
					String[] pairs=fc.split(";");
					c.famRank=new int[pairs.length];
					c.famCount=new int[pairs.length];
					for(int i=0; i<pairs.length; i++){
						int colon=pairs[i].indexOf(':');
						c.famRank[i]=Integer.parseInt(pairs[i].substring(0, colon));
						c.famCount[i]=Integer.parseInt(pairs[i].substring(colon+1));
					}
				}else{c.famRank=EMPTY; c.famCount=EMPTY;}
				ArrayList<Contig> list=byTid.get(c.tid);
				if(list==null){byTid.put(c.tid, list=new ArrayList<Contig>());}
				list.add(c);
				if(!tid2domainIdx.containsKey(c.tid)){tid2domainIdx.put(c.tid, domainIndex(f[2]));}
			}
			br.close();
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Process            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads a family-feature subset (one rank per line) for reduced-width vectors. */
	private int[] loadRanks(String file){
		try{
			ArrayList<Integer> l=new ArrayList<Integer>();
			BufferedReader br=new BufferedReader(new FileReader(file));
			String s;
			while((s=br.readLine())!=null){s=s.trim(); if(s.length()>0){l.add(Integer.parseInt(s));}}
			br.close();
			int[] a=new int[l.size()];
			for(int i=0; i<a.length; i++){a[i]=l.get(i);}
			System.err.println("feature subset: "+a.length+" family ranks kept");
			return a;
		}catch(Exception e){throw new RuntimeException(e);}
	}

	void process(){
		loadAux();
		loadCache();

		// usable tids: present in cache + sizemap, with contigs
		ArrayList<Integer> usable=new ArrayList<Integer>();
		for(Integer tid : byTid.keySet()){
			if(genomeSize.containsKey(tid) && !byTid.get(tid).isEmpty()){usable.add(tid);}
		}
		java.util.Collections.sort(usable);
		System.err.println("usable orgs="+usable.size()+", numFam="+numFam+", numPhyla="+numPhyla);

		// optional taxonomy for same-family contaminant bias
		if(treeFile!=null){
			tree=TaxTree.loadTaxTree(treeFile, System.err, false, false);
			for(Integer tid : usable){
				TaxNode fn=(tree==null ? null : tree.getNodeAtLevel(tid, TaxTree.FAMILY));
				tid2family.put(tid, fn==null ? -1 : fn.id);
			}
		}

		// global organism split BEFORE sampling
		Random split=new Random(seed);
		java.util.Collections.shuffle(usable, split);
		int nVal=(int)Math.round(usable.size()*valfrac);
		ArrayList<Integer> valTids=new ArrayList<Integer>(usable.subList(0, nVal));
		ArrayList<Integer> trainTids=new ArrayList<Integer>(usable.subList(nVal, usable.size()));
		System.err.println("train orgs="+trainTids.size()+", val orgs="+valTids.size());

		// precompute recoverable bp per tid
		for(Integer tid : usable){
			long sum=0; for(Contig c : byTid.get(tid)){sum+=c.length;}
			recoverable.put(tid, sum);
		}

		// precompute each organism's NATIVE ncRNA complement (the subnet denominator):
		// {r16,r23,r5,rother,trna} summed over ALL of the tid's contigs.
		for(Integer tid : usable){
			int[] nc=new int[5];
			for(Contig c : byTid.get(tid)){
				nc[0]+=c.r16; nc[1]+=c.r23; nc[2]+=c.r5; nc[3]+=c.rother; nc[4]+=c.trna;
			}
			nativeNc.put(tid, nc);
		}

		if(enc==ENC_NORM){
			avgCopy=computeAvgCopy(usable);
			System.err.println("enc=norm: computed avgCopyWhenPresent over "+usable.size()+" orgs");
		}
		if(featuresFile!=null){keptRanks=loadRanks(featuresFile);}
		precomputeNStrings();

		final int baseFam=(keptRanks!=null ? keptRanks.length : numFam);
		final int famCols=baseFam*(enc==ENC_TWO ? 2 : 1);
		numInputs=famCols+NUM_GLOBALS+numPhyla;
		// ncRNA subnet input width: 5 obs + phylum one-hot + 5 context (+ optional blocks).
		subnetInputs=NCRNA_OBS+numPhyla+NCRNA_CONTEXT
			+(snDomain?DOMAINS:0)+(snHHCAGA?2:0)+(snGeneLen?2:0);
		final boolean subnet=("ncrna".equals(subnetName));
		if(subnetName!=null && !subnet){throw new RuntimeException("Unknown subnet="+subnetName+" (only ncrna)");}
		if(subnet && snHHCAGA){
			if(kmerFile==null){throw new RuntimeException("snhhcaga=t requires kmerfile=<tid HH CAGA>");}
			loadKmerFile(kmerFile);
		}
		writeSet(out, (subnet ? subnetOut : null), trainTids, n, new Random(seed*2+1), "train");
		if(outval!=null && valn>0 && !valTids.isEmpty()){
			writeSet(outval, (subnet ? subnetValOut : null), valTids, valn, new Random(seed*2+2), "val");
		}
		System.err.println("done.");
	}

	/** Builds family->tids index within a pool (for same-family contaminant selection). */
	private HashMap<Integer,ArrayList<Integer>> familyIndex(ArrayList<Integer> pool){
		HashMap<Integer,ArrayList<Integer>> m=new HashMap<Integer,ArrayList<Integer>>();
		for(Integer tid : pool){
			int fam=tid2family.getOrDefault(tid, -1);
			ArrayList<Integer> l=m.get(fam);
			if(l==null){m.put(fam, l=new ArrayList<Integer>());}
			l.add(tid);
		}
		return m;
	}

	private void writeSet(String file, String subnetFile, ArrayList<Integer> pool, long count, Random rnd, String tag){
		HashMap<Integer,ArrayList<Integer>> famIdx=familyIndex(pool);
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file), 1<<20);
			bw.write("#dims\t"+numInputs+"\t2\t0\n");
			BufferedWriter sbw=null;
			StringBuilder sb2=null;
			if(subnetFile!=null){
				sbw=new BufferedWriter(new FileWriter(subnetFile), 1<<20);
				sbw.write("#dims\t"+subnetInputs+"\t1\t0\n");
				sb2=new StringBuilder(subnetInputs*4);
			}
			long made=0, tries=0;
			StringBuilder sb=new StringBuilder(numInputs*4);
			while(made<count && tries<count*20+1000){
				tries++;
				double[] labels=new double[2];
				int[] fam=new int[numFam];
				double[] glob=new double[NUM_GLOBALS];
				int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels);
				if(targetPhylumIdx<0){continue;}
				formatRow(sb, fam, glob, targetPhylumIdx, labels);
				bw.write(sb.toString());
				sb.setLength(0);
				if(sbw!=null){
					formatNcrnaRow(sb2, glob, targetPhylumIdx);
					sbw.write(sb2.toString());
					sb2.setLength(0);
				}
				made++;
				if((made%50000)==0){System.err.println(tag+": "+made+"/"+count);}
			}
			bw.close();
			if(sbw!=null){sbw.close(); System.err.println(tag+": wrote "+made+" ncRNA-subnet rows to "+subnetFile);}
			System.err.println(tag+": wrote "+made+" rows to "+file+" (tries="+tries+")");
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Bin sampler         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Synthesizes one bin into the provided fam[]/glob[] accumulators and labels[].
	 * @return the target's phylum index, or -1 if the attempt failed (retry).
	 */
	private int makeBin(ArrayList<Integer> pool, HashMap<Integer,ArrayList<Integer>> famIdx,
			Random rnd, int[] fam, double[] glob, double[] labels){
		int target=pool.get(rnd.nextInt(pool.size()));
		long gsize=genomeSize.get(target);
		long recov=recoverable.get(target);
		if(gsize<=0 || recov<=0){return -1;}

		// sampled completeness (flat + sqrt-high mixture), capped by recoverable fraction
		double comp=sampleComp(rnd);
		long targetBp=(long)(comp*gsize);
		Agg agg=new Agg(fam);
		long cleanBp=selectContigs(byTid.get(target), targetBp, rnd, agg);
		if(cleanBp<=0){return -1;}
		// Snapshot the TARGET's observed ncRNA (before contaminants are added) and its id, so a
		// subnet emitter can pair this bin's observed ncRNA with the target's native complement.
		// Read-only w.r.t. the global vector (agg is only inspected).
		lastTarget=target;
		lastNcObs[0]=agg.r16; lastNcObs[1]=agg.r23; lastNcObs[2]=agg.r5;
		lastNcObs[3]=agg.rother; lastNcObs[4]=agg.trna;

		// sampled contamination (clean spike, else flat + square-low mixture)
		double cont=(rnd.nextDouble()<cleanSpike ? 0.0 : sampleCont(rnd));
		long foreignBp=0;
		if(cont>0){
			long foreignTarget=(long)((cont/(1.0-cont))*cleanBp);
			if(foreignTarget>0){
				int nContam=1;
				if(rnd.nextDouble()<multiContamProb){nContam=2+rnd.nextInt(2);}//2 or 3
				long per=Math.max(1, foreignTarget/nContam);
				for(int k=0; k<nContam; k++){
					int ctid=pickContaminant(pool, famIdx, target, rnd);
					if(ctid<0){break;}
					foreignBp+=selectContigs(byTid.get(ctid), per, rnd, agg);
				}
			}
		}

		long totalBp=cleanBp+foreignBp;
		if(totalBp<=0 || agg.contigs<=0){return -1;}
		lastTotalBp=totalBp;

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

		String phy=tid2phylum.getOrDefault(target, "other");
		Integer pi=phylumIndex.get(phy);
		return pi==null ? phylumIndex.get("other") : pi;
	}

	/** Picks a contaminant tid != target, favoring the target's family. */
	private int pickContaminant(ArrayList<Integer> pool, HashMap<Integer,ArrayList<Integer>> famIdx,
			int target, Random rnd){
		if(!tid2family.isEmpty() && rnd.nextDouble()<sameFamProb){
			int fam=tid2family.getOrDefault(target, -1);
			ArrayList<Integer> l=famIdx.get(fam);
			if(l!=null && l.size()>1){
				for(int t=0; t<8; t++){int c=l.get(rnd.nextInt(l.size())); if(c!=target){return c;}}
			}
		}
		for(int t=0; t<8; t++){int c=pool.get(rnd.nextInt(pool.size())); if(c!=target){return c;}}
		return -1;
	}

	/** Randomly adds contigs (by shuffled draw) until reaching targetBp; returns bp added. */
	private long selectContigs(ArrayList<Contig> contigs, long targetBp, Random rnd, Agg agg){
		int nc=contigs.size();
		int[] order=agg.scratch(nc);
		for(int i=0; i<nc; i++){order[i]=i;}
		// partial Fisher-Yates: shuffle enough to draw without replacement
		long added=0;
		for(int i=0; i<nc && added<targetBp; i++){
			int j=i+rnd.nextInt(nc-i);
			int tmp=order[i]; order[i]=order[j]; order[j]=tmp;
			Contig c=contigs.get(order[i]);
			agg.add(c);
			added+=c.length;
		}
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Aggregator          ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulates per-bin sufficient statistics over selected contigs. */
	final class Agg {
		final int[] fam;
		long gc, acgt, coding, glenSum, glenSq;
		int contigs, cds, mapped, r16, r23, r5, rother, trna;
		int[] lens=new int[64]; int nlens=0;
		int[] scratchArr;
		Agg(int[] fam){this.fam=fam;}
		int[] scratch(int need){
			if(scratchArr==null || scratchArr.length<need){scratchArr=new int[need];}
			return scratchArr;
		}
		void add(Contig c){
			contigs++; gc+=c.gc; acgt+=c.acgt; coding+=c.coding; cds+=c.cds; mapped+=c.mapped;
			glenSum+=c.glenSum; glenSq+=c.glenSq; r16+=c.r16; r23+=c.r23; r5+=c.r5; rother+=c.rother; trna+=c.trna;
			for(int i=0; i<c.famRank.length; i++){fam[c.famRank[i]]+=c.famCount[i];}
			if(nlens>=lens.length){lens=Arrays.copyOf(lens, lens.length*2);}
			lens[nlens++]=c.length;
		}
		double geneStd(){
			if(cds<=0){return 0;}
			double mean=glenSum/(double)cds;
			double var=glenSq/(double)cds-mean*mean;
			return var>0 ? Math.sqrt(var) : 0;
		}
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
		if(enc==ENC_TWO){excessArr=new String[N_STR_MAX+1];}
		double logCap=Math.log(1+LOG_CAP)/LOG2;
		for(int i=0; i<=N_STR_MAX; i++){
			nOverN1[i]=fmt(encodeCount(i, logCap));
			if(enc==ENC_TWO){int ex=Math.min(Math.max(i-1, 0), EXC_CAP); excessArr[i]=fmt(ex/(double)EXC_CAP);}
		}
	}
	/** Single-column encoding of a family's summed count under the active enc mode. */
	private double encodeCount(int count, double logCap){
		if(enc==ENC_LOG){double v=(Math.log(1+count)/LOG2)/logCap; return v>1 ? 1 : v;}
		if(enc==ENC_RAW){return Math.min(count, RAW_CAP)/(double)RAW_CAP;}
		return count/(double)(1+count);//ENC_RATIO (ENC_TWO handles presence separately)
	}
	private String famStr(int count){
		if(count<=N_STR_MAX){return nOverN1[count];}
		return fmt(encodeCount(count, Math.log(1+LOG_CAP)/LOG2));
	}
	/** Two-channel family: presence (0/1) then excess-copies min(count-1,cap)/cap. */
	private void appendTwo(StringBuilder sb, int count){
		sb.append(count>0 ? '1' : '0'); sb.append('\t');
		if(count<=N_STR_MAX){sb.append(excessArr[count]);}
		else{int ex=Math.min(count-1, EXC_CAP); sb.append(fmt(ex/(double)EXC_CAP));}
		sb.append('\t');
	}
	/** Fixed-notation float, no exponent (RegressionTrainer's fast parser rejects 'e'). */
	private static String fmt(double v){
		if(v==(long)v){return Long.toString((long)v);}
		return String.format("%.6f", v);
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

	/** Loads per-organism scaled HH/CAGA (tid&lt;TAB&gt;HH&lt;TAB&gt;CAGA) for the ncRNA subnet. */
	private void loadKmerFile(String file){
		try{
			tidKmer=new HashMap<Integer,float[]>();
			BufferedReader br=new BufferedReader(new FileReader(file));
			String line;
			while((line=br.readLine())!=null){
				if(line.length()==0 || line.charAt(0)=='#'){continue;}
				String[] p=line.split("\t");
				if(p.length>=3){tidKmer.put(Integer.parseInt(p[0]),
					new float[]{Float.parseFloat(p[1]), Float.parseFloat(p[2])});}
			}
			br.close();
			System.err.println("loaded HH/CAGA for "+tidKmer.size()+" orgs");
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/**
	 * Per-family expected copy number WHEN PRESENT, over the reference organisms:
	 * for each family, the mean of the organism-level count across the organisms where
	 * that family appears (count &gt; 0). A per-family baseline so enc=norm can encode a
	 * present-at-typical-copies family as ~1.0 and a duplicated family as &gt;1.0,
	 * calibrated to each family's own natural copy number (Brian 2026-08-06). This is a
	 * reference-DB constant (not a label), computed over ALL usable orgs for stability.
	 */
	private double[] computeAvgCopy(ArrayList<Integer> usable){
		long[] sum=new long[numFam];
		int[] present=new int[numFam];
		int[] orgCount=new int[numFam];
		for(Integer tid : usable){
			Arrays.fill(orgCount, 0);
			for(Contig c : byTid.get(tid)){
				for(int i=0; i<c.famRank.length; i++){orgCount[c.famRank[i]]+=c.famCount[i];}
			}
			for(int f=0; f<numFam; f++){if(orgCount[f]>0){sum[f]+=orgCount[f]; present[f]++;}}
		}
		double[] avg=new double[numFam];
		for(int f=0; f<numFam; f++){avg[f]=present[f]>0 ? sum[f]/(double)present[f] : 1.0;}
		return avg;
	}
	/** enc=norm: 0 if absent, else count / avgCopyWhenPresent[rank]. */
	private String normStr(int count, int rank){
		if(count==0){return "0";}
		double a=avgCopy[rank];
		return fmt(a>0 ? count/a : count);
	}

	private void formatRow(StringBuilder sb, int[] fam, double[] glob, int phylumIdx, double[] labels){
		if(enc==ENC_TWO){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){appendTwo(sb, fam[keptRanks[i]]);}
			}else{
				for(int i=0; i<numFam; i++){appendTwo(sb, fam[i]);}
			}
		}else if(enc==ENC_NORM){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){int r=keptRanks[i]; sb.append(normStr(fam[r], r)); sb.append('\t');}
			}else{
				for(int i=0; i<numFam; i++){sb.append(normStr(fam[i], i)); sb.append('\t');}
			}
		}else{
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){sb.append(famStr(fam[keptRanks[i]])); sb.append('\t');}
			}else{
				for(int i=0; i<numFam; i++){sb.append(famStr(fam[i])); sb.append('\t');}
			}
		}
		for(int i=0; i<NUM_GLOBALS; i++){sb.append(fmt(glob[i])); sb.append('\t');}
		for(int i=0; i<numPhyla; i++){sb.append(i==phylumIdx ? '1' : '0'); sb.append('\t');}
		sb.append(fmt(labels[0])); sb.append('\t');
		sb.append(fmt(labels[1])); sb.append('\n');
	}

	/**
	 * Emits one ncRNA-subnet training row for the current bin: the target organism's
	 * observed ncRNA counts (r16,r23,r5,rother,trna) plus a shared context block (phylum
	 * one-hot and size/composition globals) as inputs, and the target's NATIVE ncRNA
	 * complement (summed over its whole genome) as the single regression target. The
	 * subnet thus learns the EXPECTED denominator; completeness = observed/expected is
	 * derived downstream (Barbara's refinement of the subset-relative-label design).
	 */
	private void formatNcrnaRow(StringBuilder sb, double[] glob, int phylumIdx){
		int obsTotal=0;
		for(int i=0; i<NCRNA_OBS; i++){sb.append(lastNcObs[i]); sb.append('\t'); obsTotal+=lastNcObs[i];}
		for(int i=0; i<numPhyla; i++){sb.append(i==phylumIdx ? '1' : '0'); sb.append('\t');}
		if(snDomain){//domain one-hot (forward-infra; constant for a bacteria-only corpus)
			final int di=tid2domainIdx.getOrDefault(lastTarget, DOMAIN_OTHER);
			for(int i=0; i<DOMAINS; i++){sb.append(i==di ? '1' : '0'); sb.append('\t');}
		}
		// context block (5 columns; values may be transformed): bin size, GC, coding density,
		// log2(contigs), log2(richness). F2/F3 rescale for the weight-decay-regularized optimizer.
		final double binSize=(snBinScaled ? log2(1.0+lastTotalBp/2048.0)*0.0625 : glob[0]);
		final double coding=(snCodingAffine ? glob[4]*1.05-0.5 : glob[4]);
		sb.append(fmt(binSize)); sb.append('\t');
		sb.append(fmt(glob[3])); sb.append('\t');
		sb.append(fmt(coding)); sb.append('\t');
		sb.append(fmt(glob[1])); sb.append('\t');
		sb.append(fmt(glob[8])); sb.append('\t');
		if(snHHCAGA){//per-org HH, CAGA (min-max scaled over training orgs), from kmerfile
			final float[] hc=tidKmer.get(lastTarget);
			sb.append(fmt(hc==null ? 0 : hc[0])); sb.append('\t');
			sb.append(fmt(hc==null ? 0 : hc[1])); sb.append('\t');
		}
		if(snGeneLen){//mean gene length + gene-length stddev (cache-derived globals)
			sb.append(fmt(glob[5])); sb.append('\t');
			sb.append(fmt(glob[6])); sb.append('\t');
		}
		final int[] nc=nativeNc.get(lastTarget);
		int nativeTotal=0;
		for(int v : nc){nativeTotal+=v;}
		// The bin's target contigs are a subset of the genome, so observed<=native always.
		assert(obsTotal<=nativeTotal) : "ncRNA observed "+obsTotal+" > native "+nativeTotal+" (tid "+lastTarget+")";
		sb.append(nativeTotal); sb.append('\n');
	}

	private static double log2(double v){return v<=1 ? 0 : Math.log(v)/LOG2;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String cacheFile, sizemapFile, familyFile, taxpgmFile, treeFile, out, outval, featuresFile;
	private String subnetName, subnetOut, subnetValOut, kmerFile;
	private boolean snCodingAffine=false, snBinScaled=false, snDomain=false, snHHCAGA=false, snGeneLen=false;
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

	private final HashMap<Integer,ArrayList<Contig>> byTid=new HashMap<Integer,ArrayList<Contig>>();
	private final HashMap<Integer,Long> genomeSize=new HashMap<Integer,Long>();
	private final HashMap<Integer,Long> recoverable=new HashMap<Integer,Long>();
	private final HashMap<Integer,int[]> nativeNc=new HashMap<Integer,int[]>();
	private final HashMap<Integer,Integer> tid2domainIdx=new HashMap<Integer,Integer>();
	private HashMap<Integer,float[]> tidKmer;
	private final HashMap<Integer,String> tid2phylum=new HashMap<Integer,String>();
	private final HashMap<Integer,Integer> tid2family=new HashMap<Integer,Integer>();
	private final HashMap<String,Integer> phylumIndex=new HashMap<String,Integer>();
	private ArrayList<String> phylumList;
	private String[] nOverN1;
	private String[] excessArr;
	private double[] avgCopy;

	private static final int NUM_GLOBALS=13;
	private static final int NCRNA_OBS=5, NCRNA_CONTEXT=5;
	private static final int DOMAINS=8, DOMAIN_OTHER=7;
	private static final int N_STR_MAX=4096;
	private static final int ENC_RATIO=0, ENC_LOG=1, ENC_RAW=2, ENC_TWO=3, ENC_NORM=4;
	private static final int RAW_CAP=32, LOG_CAP=64, EXC_CAP=16;
	private static final double COMP_MIN=0.10, CONT_MAX=0.50, LOG2=Math.log(2);
	private static final int[] EMPTY=new int[0];
}
