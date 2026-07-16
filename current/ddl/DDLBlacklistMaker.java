package ddl;

import java.io.PrintStream;
import java.util.ArrayList;

import cardinality.DynamicDemiLog;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntHashSet;
import map.LongHashSet;
import map.LongObjectMap;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Builds a DDL kmer blacklist from pre-built DDL sketch files.
 * Reads sketches with kmer arrays (kmers=t), promotes each TaxID
 * to genus level, counts distinct genera per kmer, and outputs
 * kmers exceeding the threshold as FASTA with multi-level headers.
 *
 * Accepts multiple comma-separated input files; processes each
 * sequentially and discards sketch data after kmer extraction
 * to minimize memory usage.
 *
 * Condense mode (condense=t blacklist=X): condenses double-sized
 * sketches to half size while preferring non-blacklisted kmers,
 * then validates. Approximates full re-sketching without reading
 * raw genome data.
 *
 * @author Brian Bushnell, Noire
 * @date May 24, 2026
 */
public class DDLBlacklistMaker {

	public static void main(String[] args){
		Timer t=new Timer();
		DDLBlacklistMaker bm=new DDLBlacklistMaker(args);
		bm.process(t);
	}

	public DDLBlacklistMaker(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("mintaxcount") || a.equals("mintaxa") || a.equals("mincount")){
				minTaxCount=Integer.parseInt(b);
			}else if(a.equals("minentropy") || a.equals("entropy")){
				minEntropy=Float.parseFloat(b);
			}else if(a.equals("entropyk") || a.equals("ek")){
				entropyK=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("validate")){
				validate=Parse.parseBoolean(b);
			}else if(a.equals("samples") || a.equals("validatesamples")){
				validateSamples=Integer.parseInt(b);
			}else if(a.equals("exponent") || a.equals("ebits")){
				DynamicDemiLog.setExponent(Integer.parseInt(b));
			}else if(a.equals("blacklist") || a.equals("bl")){
				blacklistFile=b;
			}else if(a.equals("condense")){
				condense=Parse.parseBoolean(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}else if(a.equals("prefilterbits") || a.equals("pfbits")){
				prefilterBits=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			Parser.processQuality();
			overwrite=parser.overwrite;
			out=parser.out1;
			in=parser.in1;
		}

		assert(in!=null) : "No input file specified.";
		assert(out!=null || validate || condense) : "No output file or action specified.";
	}

	void process(Timer t){
		if(condense){
			processCondense(t);
		}else{
			processBlacklist(t);
		}
	}

	/*--------------------------------------------------------------*/
	/*-----------        Condense + Validate        ----------------*/
	/*--------------------------------------------------------------*/

	void processCondense(Timer t){
		final LongHashSet bl;
		if(blacklistFile!=null){
			DynamicDemiLog.loadBlacklist(blacklistFile, k);
			bl=getBlacklistSet();
		}else{
			bl=null;
		}

		final String[] inFiles=in.split(",");
		final ArrayList<DDLRecord> condensed=new ArrayList<DDLRecord>();
		long totalRecords=0;
		int blacklisted=0, kept=0;

		for(String inFile : inFiles){
			outstream.println("Loading "+inFile.trim());
			ArrayList<DDLRecord> records=DDLLoader.loadFile(inFile.trim(), k);
			outstream.println("  "+records.size()+" records");
			totalRecords+=records.size();

			for(DDLRecord rec : records){
				DynamicDemiLog ddl=rec.ddl;
				int oldBuckets=ddl.maxArray().length;
				if(oldBuckets<2){continue;}
				int newBuckets=oldBuckets/2;

				char[] oldMax=ddl.maxArray();
				long[] oldKmers=ddl.kmerArray();

				char[] newMax=new char[newBuckets];
				for(int i=0; i<newBuckets; i++){
					int j=i+newBuckets;
					char scoreA=oldMax[i], scoreB=oldMax[j];
					boolean blA=(bl!=null && oldKmers!=null && oldKmers[i]!=0 && bl.contains(oldKmers[i]));
					boolean blB=(bl!=null && oldKmers!=null && oldKmers[j]!=0 && bl.contains(oldKmers[j]));

					if(blA && !blB){
						newMax[i]=scoreB; blacklisted++;
					}else if(blB && !blA){
						newMax[i]=scoreA; blacklisted++;
					}else{
						newMax[i]=(scoreA>=scoreB ? scoreA : scoreB);
						if(blA && blB){kept++;}
					}
				}

				DynamicDemiLog cddl=DynamicDemiLog.fromArray(newMax, (int)rec.id, rec.name, k, -1);
				condensed.add(new DDLRecord(cddl, rec.id, rec.taxID, rec.name));
			}
			records.clear();
		}

		outstream.println("\nTotal records: "+totalRecords+", condensed: "+condensed.size());
		outstream.println("Blacklisted bucket winners replaced: "+blacklisted);
		outstream.println("Both candidates blacklisted (kept higher): "+kept);

		validate(condensed);

		t.stop();
		outstream.println("\nTime: \t"+t);
	}

	private static LongHashSet getBlacklistSet(){
		try{
			java.lang.reflect.Field f=DynamicDemiLog.class.getDeclaredField("blacklist");
			f.setAccessible(true);
			return (LongHashSet)f.get(null);
		}catch(Exception e){
			throw new RuntimeException("Cannot access blacklist field", e);
		}
	}

	/*--------------------------------------------------------------*/
	/*-----------        Blacklist Generation        ---------------*/
	/*--------------------------------------------------------------*/

	@SuppressWarnings("unchecked")
	void processBlacklist(Timer t){
		final TaxTree tree=TaxTree.sharedTree();
		final String[] inFiles=in.split(",");
		final int shards=Tools.mid(1, Shared.threads(), MAX_SHARDS);
		outstream.println("Sharding by bucket across "+shards+" thread(s).");
		final Timer pt=new Timer(outstream, true);//Phase timer, independent of the total-runtime timer t.

		//Optional pass 1: build a lossy counting table so the exact maps below only hold candidates.
		final byte[][] counts=(prefilter ? countKmers(inFiles, shards) : null);
		if(prefilter){pt.stopAndStart("Prefilter pass 1 total:");}

		//Disjoint per-shard maps.  A kmer's bucket = hash(kmer)&bucketMask is fixed across every sketch,
		//so kmerArray[i] always holds a kmer whose bucket is i; splitting buckets into contiguous ranges
		//therefore partitions kmer space exactly.  The per-shard maps never share a key, so they need no
		//locks and no merge, and no single map approaches LongObjectMap's array-backed capacity.
		final LongObjectMap<IntHashSet>[] maps=new LongObjectMap[shards];
		for(int s=0; s<shards; s++){maps[s]=new LongObjectMap<IntHashSet>(IntHashSet.class);}

		long totalRecords=0, totalWithKmers=0, prefiltered=0, loadNanos=0, mapNanos=0;
		int skippedTids=0;

		for(String inFile : inFiles){
			outstream.println("Loading "+inFile);
			long ts=System.nanoTime();
			ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(inFile.trim(), k);
			loadNanos+=System.nanoTime()-ts;
			outstream.println("  "+records.size()+" records");
			totalRecords+=records.size();

			final int buckets=bucketCount(records);
			final int[] genusTids=genusTids(records, tree);//Resolved once, not once per shard.

			int withKmers=0;
			for(int r=0; r<records.size(); r++){
				if(records.get(r).ddl.hasKmers()){
					withKmers++;
					if(genusTids[r]<0){skippedTids++;}
				}
			}
			totalWithKmers+=withKmers;
			if(buckets<1){records.clear(); continue;}

			ts=System.nanoTime();
			final ArrayList<MapThread> list=new ArrayList<MapThread>(shards);
			for(int s=0; s<shards; s++){
				list.add(new MapThread(records, genusTids, counts, maps[s],
					shardLo(s, buckets, shards), shardLo(s+1, buckets, shards)));
			}
			for(MapThread mt : list){mt.start();}
			for(MapThread mt : list){
				while(mt.getState()!=Thread.State.TERMINATED){
					try{mt.join();}catch(InterruptedException e){e.printStackTrace();}
				}
				prefiltered+=mt.prefiltered;
			}
			mapNanos+=System.nanoTime()-ts;
			records.clear();
			outstream.println("  Kmer map size: "+mapSize(maps));
		}
		if(counts!=null){
			outstream.println("Prefilter rejected "+prefiltered+" kmer instances before the exact map.");
		}
		outstream.println("  [pass 2 phases: load "+sec(loadNanos)+"s (MT load, reader-bound), accumulate "
			+sec(mapNanos)+"s ("+shards+" threads)]");
		pt.stopAndStart("Accumulation pass 2 total:");

		outstream.println("\nTotal records: "+totalRecords+" ("+totalWithKmers+" with kmers)");
		if(skippedTids>0){outstream.println("Skipped "+skippedTids+" records with out-of-range TaxIDs.");}
		outstream.println("Unique kmers across all files: "+mapSize(maps));

		if(validate && inFiles.length==1){
			outstream.println("\nReloading for validation...");
			ArrayList<DDLRecord> records=DDLLoader.loadFile(inFiles[0].trim(), k);
			validate(records);
		}else if(validate){
			outstream.println("Validate mode requires a single input file; skipping.");
		}

		if(out!=null && totalWithKmers>0){
			buildBlacklist(maps, t);
		}else if(out!=null){
			outstream.println("ERROR: No records have kmer arrays. Rebuild with kmers=t.");
		}
	}

	/*--------------------------------------------------------------*/

	/**
	 * Optional pass 1: count kmer occurrences in a lossy counting table, so that pass 2 only needs
	 * an exact map entry for kmers that could possibly survive.
	 *
	 * The table is counts[bucket][kmer & cellMask], saturating at 255.  The bucket index is free:
	 * a DDL sketch stores at most one kmer per bucket, so a kmer's array position IS its bucket, and
	 * a kmer's bucket is fixed (hash & bucketMask) across every sketch.
	 *
	 * Losslessness: a kmer in g distinct genera occupies g distinct taxids, hence at least g
	 * records, and its cell aggregates at least its own records, so
	 *     cellCount >= recordCount >= genusCount.
	 * No kmer that would pass minTaxCount can be filtered out here; collisions only INFLATE counts,
	 * admitting extra kmers to pass 2, which counts them exactly and discards them.  The blacklist is
	 * therefore identical to the one-pass build.
	 *
	 * This exists because the exact map is the memory bottleneck: LongObjectMap is array-backed and
	 * capped near 1.825 billion entries (SAFE_ARRAY_LEN * 0.85), which a 64k-bucket RefSeq build
	 * exceeds outright.
	 *
	 * @param inFiles Sketch files to scan
	 * @param shards Number of bucket ranges to count in parallel
	 * @return counts[bucket][cell], or null if no records had kmer arrays
	 */
	private byte[][] countKmers(String[] inFiles, final int shards){
		final int cells=1<<prefilterBits;
		byte[][] counts=null;
		long instances=0, loadNanos=0, countNanos=0;

		for(String inFile : inFiles){
			outstream.println("Prefilter pass 1: counting "+inFile.trim());
			long ts=System.nanoTime();
			ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(inFile.trim(), k);
			loadNanos+=System.nanoTime()-ts;
			final int buckets=bucketCount(records);
			if(buckets<1){records.clear(); continue;}
			if(counts==null){
				counts=new byte[buckets][cells];
				outstream.println("  Counting table: "+buckets+" buckets x "+cells+" cells = "
					+String.format("%.1f", buckets*(double)cells/(1024*1024*1024))+" GB");
			}
			assert(buckets==counts.length) : "Sketches disagree on bucket count: "+buckets+" vs "+counts.length;

			ts=System.nanoTime();
			final ArrayList<CountThread> list=new ArrayList<CountThread>(shards);
			for(int s=0; s<shards; s++){
				list.add(new CountThread(records, counts, shardLo(s, buckets, shards), shardLo(s+1, buckets, shards)));
			}
			for(CountThread ct : list){ct.start();}
			for(CountThread ct : list){
				while(ct.getState()!=Thread.State.TERMINATED){
					try{ct.join();}catch(InterruptedException e){e.printStackTrace();}
				}
				instances+=ct.instances;
			}
			countNanos+=System.nanoTime()-ts;
			records.clear();
		}
		outstream.println("Prefilter pass 1 done: "+instances+" kmer instances counted.");
		outstream.println("  [pass 1 phases: load "+sec(loadNanos)+"s (MT load, reader-bound), count "
			+sec(countNanos)+"s ("+shards+" threads)]");
		return counts;
	}

	/** Counts kmers for one contiguous range of buckets.  Threads own disjoint rows of the counting
	 * table (a kmer at array index i lands in counts[i]), so increments never race and no atomics
	 * are needed. */
	private class CountThread extends Thread {

		CountThread(ArrayList<DDLRecord> records_, byte[][] counts_, int lo_, int hi_){
			records=records_; counts=counts_; lo=lo_; hi=hi_;
		}

		@Override
		public void run(){
			final int cellMask=(1<<prefilterBits)-1;
			for(DDLRecord rec : records){
				final long[] kmers=rec.ddl.kmerArray();
				if(kmers==null){continue;}
				final int limit=Math.min(hi, kmers.length);
				for(int i=lo; i<limit; i++){
					final long kmer=kmers[i];
					if(kmer==0){continue;}
					final byte[] row=counts[i];
					final int cell=(int)(kmer&cellMask);
					if((row[cell]&0xFF)<255){row[cell]++;}//Saturate rather than wrap
					instances++;
				}
			}
		}

		private final ArrayList<DDLRecord> records;
		private final byte[][] counts;
		private final int lo, hi;
		long instances=0;
	}

	/** Accumulates kmer->genera for one contiguous bucket range into a map owned solely by this thread.
	 * Bucket ranges partition kmer space exactly (a kmer's bucket is fixed across every sketch), so the
	 * per-shard maps are disjoint and require no merge. */
	private class MapThread extends Thread {

		MapThread(ArrayList<DDLRecord> records_, int[] genusTids_, byte[][] counts_,
				LongObjectMap<IntHashSet> map_, int lo_, int hi_){
			records=records_; genusTids=genusTids_; counts=counts_; map=map_; lo=lo_; hi=hi_;
		}

		@Override
		public void run(){
			final int cellMask=(counts!=null ? (1<<prefilterBits)-1 : 0);
			for(int r=0; r<records.size(); r++){
				final int genusTid=genusTids[r];
				if(genusTid<0){continue;}//Unusable TaxID
				final long[] kmers=records.get(r).ddl.kmerArray();
				if(kmers==null){continue;}
				final int limit=Math.min(hi, kmers.length);
				for(int i=lo; i<limit; i++){
					final long kmer=kmers[i];
					if(kmer==0){continue;}
					if(counts!=null && !survivesPrefilter(counts[i][(int)(kmer&cellMask)], kmer)){
						prefiltered++;
						continue;
					}
					IntHashSet genera=map.get(kmer);
					if(genera==null){
						genera=new IntHashSet(4);
						map.put(kmer, genera);
					}
					genera.add(genusTid);
				}
			}
		}

		private final ArrayList<DDLRecord> records;
		private final int[] genusTids;
		private final byte[][] counts;
		private final LongObjectMap<IntHashSet> map;
		private final int lo, hi;
		long prefiltered=0;
	}

	/** First bucket owned by shard s, splitting [0,buckets) into contiguous ranges. */
	private static int shardLo(int s, int buckets, int shards){
		return (int)((long)buckets*s/shards);
	}

	/** Genus TaxID per record, or -1 where the TaxID is unusable.  Resolved once per file rather than
	 * once per shard, since every shard would otherwise repeat the same tree walks. */
	private static int[] genusTids(ArrayList<DDLRecord> records, TaxTree tree){
		final int[] out=new int[records.size()];
		for(int r=0; r<out.length; r++){
			final int tid=records.get(r).taxID;
			if(tid<1 || (tree!=null && tid>=tree.nodes.length)){out[r]=-1; continue;}
			int genusTid=tid;
			if(tree!=null){
				TaxNode gn=tree.getNodeAtLevel(tid, TaxTree.GENUS);
				if(gn!=null){genusTid=gn.id;}
			}
			out[r]=genusTid;
		}
		return out;
	}

	/** Total live entries across all per-shard maps. */
	private static long mapSize(LongObjectMap<IntHashSet>[] maps){
		long sum=0;
		for(LongObjectMap<IntHashSet> m : maps){sum+=m.size();}
		return sum;
	}

	/** Format nanoseconds as seconds with 2 decimals, for phase timing. */
	private static String sec(long nanos){return Tools.format("%.2f", nanos/1e9);}

	/** Bucket count of the first record carrying kmers; 0 if none do. */
	private static int bucketCount(ArrayList<DDLRecord> records){
		for(DDLRecord rec : records){
			final long[] kmers=rec.ddl.kmerArray();
			if(kmers!=null){return kmers.length;}
		}
		return 0;
	}

	/**
	 * Decides whether a kmer earns an exact map entry, given its (possibly inflated) cell count.
	 * A saturated cell always passes, so minTaxCount above 255 stays safe.  Low-entropy kmers bypass
	 * the count entirely, mirroring buildBlacklist's entropy criterion (entropy depends only on the
	 * kmer, so it needs no map).
	 */
	private boolean survivesPrefilter(final byte cellCount, final long kmer){
		final int count=cellCount&0xFF;
		if(count>=255 || count>=minTaxCount){return true;}
		return minEntropy>0 && kmerEntropy(kmer, k, entropyK)<minEntropy;
	}

	/** Gathers and rank-promotes the surviving kmers of ONE shard map into its own results list.
	 * The shard maps are disjoint and TaxTree lookups are read-only, so with a per-thread scratch set
	 * and per-thread results, the promotion parallelizes with nothing mutable shared. */
	private class BuildThread extends Thread {

		BuildThread(LongObjectMap<IntHashSet> map_, TaxTree tree_, int[] levels_){
			map=map_; tree=tree_; LEVELS=levels_;
		}

		@Override
		public void run(){
			final IntHashSet promoted=new IntHashSet(256);
			final long[] keys=map.keys();
			final long invalid=map.invalid();
			for(long kmer : keys){
				if(kmer==invalid){continue;}
				final IntHashSet genera=map.get(kmer);
				if(genera==null){continue;}
				final int genusCount=genera.size();
				final float ent=kmerEntropy(kmer, k, entropyK);
				final boolean passesTaxa=(genusCount>=minTaxCount);
				final boolean passesEntropy=(minEntropy>0 && ent<minEntropy);
				if(!passesTaxa && !passesEntropy){continue;}
				if(!passesTaxa){byEntropy++;}

				long[] entry=new long[2+LEVELS.length+1];
				entry[0]=kmer;
				entry[1]=genusCount;
				entry[2]=genusCount;
				entry[2+LEVELS.length]=(long)(ent*1000);

				if(tree!=null){
					int[] tids=genera.toArray();
					for(int lev=1; lev<LEVELS.length; lev++){
						promoted.clear();
						for(int tid : tids){
							TaxNode node=tree.getNodeAtLevel(tid, LEVELS[lev]);
							if(node!=null){promoted.add(node.id);}
						}
						entry[2+lev]=promoted.size();
					}
				}
				results.add(entry);
			}
		}

		private final LongObjectMap<IntHashSet> map;
		private final TaxTree tree;
		private final int[] LEVELS;
		final ArrayList<long[]> results=new ArrayList<long[]>();
		int byEntropy=0;
	}

	/*--------------------------------------------------------------*/

	void buildBlacklist(LongObjectMap<IntHashSet>[] maps, Timer t){
		final TaxTree tree=TaxTree.sharedTree();
		final int[] LEVELS={TaxTree.GENUS, TaxTree.FAMILY, TaxTree.ORDER,
			TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM};
		final String[] LCODES={"g", "f", "o", "c", "p", "k", "sk"};
		final Timer bt=new Timer(outstream, true);//Build-phase timer.

		//Parallel gather+promote: one thread per disjoint shard map.  TaxTree lookups are read-only and
		//each thread keeps its own scratch set + results list, so the threads share nothing mutable.
		final ArrayList<BuildThread> builders=new ArrayList<BuildThread>(maps.length);
		for(LongObjectMap<IntHashSet> m : maps){builders.add(new BuildThread(m, tree, LEVELS));}
		for(BuildThread b : builders){b.start();}
		final ArrayList<long[]> results=new ArrayList<long[]>();
		int byEntropy=0;
		for(BuildThread b : builders){
			while(b.getState()!=Thread.State.TERMINATED){
				try{b.join();}catch(InterruptedException e){e.printStackTrace();}
			}
			results.addAll(b.results);//Concatenated in map order -> deterministic input to the sort below.
			byEntropy+=b.byEntropy;
		}
		bt.stopAndStart("  build: gather+promote:");

		//Genus count descending; kmer ascending breaks ties so output is reproducible regardless of the
		//shard count (which only changes the order entries were discovered in, never the set).
		results.sort((a, b)->{int c=Long.compare(b[1], a[1]); return c!=0 ? c : Long.compare(a[0], b[0]);});
		bt.stopAndStart("  build: sort:");
		outstream.println("Kmers output: "+results.size()
			+" ("+(results.size()-byEntropy)+" by genus count >= "+minTaxCount
			+(byEntropy>0 ? ", "+byEntropy+" by low entropy < "+minEntropy : "")
			+")");

		final FileFormat ff=FileFormat.testOutput(out, FileFormat.FASTA, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(long[] entry : results){
			final long kmer=entry[0];
			StringBuilder sb=new StringBuilder(80);
			sb.append(">kmer_").append(Long.toHexString(kmer));
			sb.append(" raw=").append(entry[1]);
			if(tree!=null){
				for(int lev=0; lev<LEVELS.length; lev++){
					sb.append(' ').append(LCODES[lev]).append('=').append(entry[2+lev]);
				}
			}
			sb.append(String.format(" e=%.3f", entry[2+LEVELS.length]/1000.0));
			bsw.print(sb.append('\n').toString());
			bsw.print(DynamicDemiLog.unpackKmer(kmer, k)+"\n");
		}
		bsw.poisonAndWait();
		bt.stop("  build: write:");

		t.stop();
		outstream.println("Wrote "+results.size()+" kmers to "+out);
		outstream.println("Time: \t"+t);
	}

	/*--------------------------------------------------------------*/

	void validate(ArrayList<DDLRecord> records){
		final TaxTree tree=TaxTree.sharedTree();
		if(tree==null){
			outstream.println("WARNING: TaxTree not loaded; cannot validate.");
			return;
		}

		final int[] levels={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY,
			TaxTree.ORDER, TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM};
		final int LIFE_BUCKET=levels.length;
		final int NUM_BUCKETS=levels.length+1;
		final String[] levelNames={"species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "life"};

		final shared.Random rand=Shared.threadLocalRandom(42);
		final int n=records.size();
		final long[] pairCounts=new long[NUM_BUCKETS];
		final double[] sharedSums=new double[NUM_BUCKETS];
		final double[] wkidSums=new double[NUM_BUCKETS];
		final double[] aniSums=new double[NUM_BUCKETS];
		long skippedNoRank=0;

		outstream.println("\nValidating shared keys across taxonomic distances ("+validateSamples+" random pairs)...");
		for(int s=0; s<validateSamples; s++){
			int i=rand.nextInt(n), j=rand.nextInt(n);
			if(i==j){continue;}
			DDLRecord a=records.get(i), b=records.get(j);
			if(a.taxID<1 || b.taxID<1){continue;}
			if(a.taxID>=tree.nodes.length || b.taxID>=tree.nodes.length){continue;}

			TaxNode an=tree.getNode(a.taxID), bn=tree.getNode(b.taxID);
			if(an==null || bn==null){continue;}
			if(an.id==bn.id){continue;}
			TaxNode ca=tree.commonAncestor(an, bn);
			if(ca==null){continue;}

			int bucket=-1;
			TaxNode node=ca;
			while(node!=null && bucket<0){
				for(int lev=0; lev<levels.length; lev++){
					if(node.level==levels[lev]){bucket=lev; break;}
				}
				if(bucket>=0){break;}
				if(node.pid<0 || node.pid==node.id){break;}
				node=tree.getNode(node.pid);
			}
			if(bucket<0){bucket=LIFE_BUCKET;}

			int[] comp=a.ddl.compareToDetailed(b.ddl);
			int lower=comp[0], shared=comp[1], higher=comp[2];
			pairCounts[bucket]++;
			sharedSums[bucket]+=shared;

			int div=Math.min(shared+lower, shared+higher);
			if(div<6){div=6;}
			float wkid=(shared>0 ? Math.min(1f, (float)shared/div) : 0);
			float ani=(wkid>0 ? (float)Math.exp(Math.log(wkid)/k) : 0);
			wkidSums[bucket]+=wkid;
			aniSums[bucket]+=ani;
		}

		outstream.println("\nShared keys by common ancestor level:");
		outstream.println("Ancestor_level\tPairs\tAvg_shared_keys\tAvg_WKID\tAvg_ANI");
		for(int lev=0; lev<NUM_BUCKETS; lev++){
			if(pairCounts[lev]>0){
				outstream.println(levelNames[lev]+"\t"+pairCounts[lev]+"\t"+
					String.format("%.2f", sharedSums[lev]/pairCounts[lev])+"\t"+
					String.format("%.6f", wkidSums[lev]/pairCounts[lev])+"\t"+
					String.format("%.4f", aniSums[lev]/pairCounts[lev]));
			}
		}
		if(skippedNoRank>0){
			outstream.println("(skipped "+skippedNoRank+" pairs with no-rank common ancestor)");
		}
		long unrelatedPairs=0;
		double unrelatedShared=0;
		for(int lev=4; lev<NUM_BUCKETS; lev++){
			unrelatedPairs+=pairCounts[lev];
			unrelatedShared+=sharedSums[lev];
		}
		if(unrelatedPairs>0){
			outstream.println("(noise floor: "+String.format("%.2f", unrelatedShared/unrelatedPairs)+
				" avg shared keys at class+ distance)");
		}
	}

	/*--------------------------------------------------------------*/

	static float kmerEntropy(long kmer, int k, int ek){
		final int numSubkmers=k-ek+1;
		if(numSubkmers<1){return 1;}
		final int subkmerSpace=1<<(2*ek);
		final int subkmerMask=subkmerSpace-1;
		final int[] counts=new int[subkmerSpace];
		for(int i=0; i<numSubkmers; i++){
			int sub=(int)((kmer>>(2*i))&subkmerMask);
			counts[sub]++;
		}
		double entropy=0;
		final double logN=Math.log(numSubkmers);
		for(int c : counts){
			if(c>0){
				double p=(double)c/numSubkmers;
				entropy-=p*Math.log(p);
			}
		}
		return (float)(entropy/logN);
	}

	/*--------------------------------------------------------------*/

	private String in;
	private String out;
	private String blacklistFile;
	private int k=19;
	private int minTaxCount=20;
	private float minEntropy=0;
	private int entropyK=3;
	private boolean overwrite=false;
	private boolean verbose=false;
	private boolean validate=false;
	private boolean condense=false;
	/** Two-pass mode: count kmers in a lossy table first, so the exact map only holds candidates. */
	private boolean prefilter=false;
	/** Low kmer bits used to index the counting table; table is buckets x 2^prefilterBits bytes. */
	private int prefilterBits=18;
	/** Upper bound on bucket shards; past the thread count, more shards only add per-map overhead. */
	private static final int MAX_SHARDS=64;
	private int validateSamples=100000;

	private static final PrintStream outstream=System.err;
}
