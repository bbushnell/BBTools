package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import consensus.BaseGraph;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.ScrabbleAligner;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Builds per-anticodon consensus tRNA sequences using pivot-based
 * O(N) alignment with ScrabbleAligner.
 * @author Neptune, Brian Bushnell
 */
public class TrnaConsensusBuilder {

	public static void main(String[] args){
		Timer t=new Timer();
		TrnaConsensusBuilder x=new TrnaConsensusBuilder(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public TrnaConsensusBuilder(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("mingroupsize") || a.equals("mingroup")){
				minGroupSize=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("minidentity")){
				minIdentity=Float.parseFloat(b);
			}else if(a.equals("clusterid") || a.equals("clusteridentity")){
				clusterIdentity=Float.parseFloat(b);
			}else if(a.equals("cluster")){
				doClustering=Parse.parseBoolean(b);
			}else if(a.equals("minclustersize") || a.equals("mincluster")){
				minClusterSize=Integer.parseInt(b);
			}else if(a.equals("recruit")){
				doRecruit=Parse.parseBoolean(b);
			}else if(a.equals("recruitid") || a.equals("recruitidentity")){
				recruitIdentity=Float.parseFloat(b);
			}else if(a.equals("outmodel") || a.equals("model")){
				outModel=b;
			}else if(a.equals("endtrim") || a.equals("endtrimfrac")){
				endTrimFrac=Float.parseFloat(b);
			}else if(a.equals("lentilt") || a.equals("lengthtilt")){
				lenTilt=Float.parseFloat(b);
			}else if(a.equals("reassignrounds") || a.equals("rounds")){
				reassignRounds=Integer.parseInt(b);
			}else if(a.equals("census")){
				census=Parse.parseBoolean(b);
			}else if(a.equals("seedin") || a.equals("seed") || a.equals("seeds")){
				seedIn=b;
			}else if(a.equals("prunecount") || a.equals("prune")){
				pruneCount=Integer.parseInt(b);
			}else if(a.equals("prefix") || a.equals("consensusprefix")){
				consensusPrefix=b;
			}else if(a.equals("stochastic") || a.equals("stochasticconsensus") || a.equals("stochastictrace")){
				stochasticConsensus=Parse.parseBoolean(b);
			}else if(a.equals("stochasticseed")){
				stochasticSeed=Long.parseLong(b);
			}else if(a.equals("outassignments")){
				outAssignments=b;
			}else if(a.equals("family")){
				family=b;
			}else if(parser.parse(arg, a, b)){
				//handled
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in=parser.in1;
		out=parser.out1;

		if(in==null){throw new RuntimeException("Error - an input file is required.");}
		if(out==null){throw new RuntimeException("Error - an output file is required.");}
		//Citan, 2026-08-28: family namespace must be explicit whenever assignments are written,
		//so clusterKeys stay collision-free across a future joint multi-family manifest (kept
		//clusterKey embeds modelLabel, and modelLabel's format -- consensusPrefix+anticodon+
		//"_c"+fci (consensusPrefix defaults to "tRNA_consensus_") -- is NOT itself family-
		//namespaced, so two different families could otherwise produce byte-identical modelLabels).
		if(outAssignments!=null && (family==null || family.isEmpty())){
			throw new RuntimeException("Error - a non-empty family= is required when outassignments= is set.");
		}

		//Citan, 2026-08-28, round 5: resolved ONCE here, single-threaded, before any parallel
		//buildFromAlignments call can happen -- this is what makes the per-call local-RNG design
		//(see buildFromAlignments/deriveLocalSeed) schedule-independent even in "random seed"
		//mode: the ACTUAL seed driving every call in this run is fixed here, once, regardless of
		//how many threads later consume it or in what order. Zero cost when stochasticConsensus
		//is false (no Random ever constructed).
		resolvedRunSeed=(!stochasticConsensus ? 0L
			: (stochasticSeed>=0 ? stochasticSeed : new java.util.Random().nextLong()));

		overwrite=parser.overwrite;

		ffin=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		ArrayList<Read> reads=loadReads();
		outstream.println("Loaded "+reads.size()+" tRNA sequences.");
		if(outAssignments!=null){checkNoDuplicateRecordIds(reads);}

		LinkedHashMap<String, ArrayList<byte[]>> seedMap=null;
		if(seedIn!=null){
			seedMap=loadSeeds(seedIn);
		}

		LinkedHashMap<String, ArrayList<Read>> groups=groupByAnticodon(reads);
		outstream.println("Found "+groups.size()+" anticodon groups.");

		if(outAssignments!=null){assignmentMap=new java.util.IdentityHashMap<>();}

		ArrayList<Read> consensusList=new ArrayList<>();
		ArrayList<BaseGraph> modelList=(outModel!=null ? new ArrayList<>() : null);
		long num=0;
		int totalClusters=0;
		for(Entry<String, ArrayList<Read>> e : groups.entrySet()){
			String anticodon=e.getKey();
			ArrayList<Read> group=e.getValue();
			if(group.size()<minGroupSize){
				if(verbose){outstream.println("Skipping "+anticodon+" ("+group.size()+" sequences, below min="+minGroupSize+")");}
				//Citan, 2026-08-28: this whole group was NEVER sequence-clustered (skipped before
				//clusterSequences() is even called) -- there is no verified redundancy relationship
				//among its members. Sharing one small clusterKey here would INVENT redundancy and
				//could wrongly let a rare/unique sequence borrow another unrelated one's "K count"
				//downstream. Each member is recorded as its own singleton orphan instead.
				recordEachAsOrphan(group);
				continue;
			}

			ArrayList<byte[]> seeds=(seedMap!=null ? seedMap.get(anticodon) : null);

			if(doClustering && group.size()>1){
				final ArrayList<ArrayList<Read>> clusters=clusterSequences(group, seeds);
				recordDroppedOrphans(group, clusters);
				//Per-cluster consensus+HBM+census are independent: compute in
				//parallel, then assemble serially in cluster order so output
				//(labels, numericIDs, census lines) is deterministic.
				final byte[][] consensi=new byte[clusters.size()][];
				final BaseGraph[] graphs=new BaseGraph[clusters.size()];
				final String[] censusLines=new String[clusters.size()];
				final String[] labels=new String[clusters.size()];
				final String anticodonF=anticodon;
				ArrayList<Future<?>> futures=new ArrayList<>();
				for(int ci=0; ci<clusters.size(); ci++){
					final int fci=ci;
					final ArrayList<Read> cluster=clusters.get(ci);
					if(cluster.size()<minClusterSize){continue;}
					futures.add(pool().submit(new Runnable(){
						@Override
						public void run(){
							final String label=consensusPrefix+anticodonF+"_c"+fci+" n="+cluster.size();
							final byte[] consensus;
							if(outModel!=null){
								//Consensus and HBM built TOGETHER, from alignments to the SAME stabilized
								//reference (Brian's ruling, Aug 22 2026, via Noire) -- see
								//buildConsensusAndGraph's javadoc for why the old separate buildConsensus()+
								//buildBaseGraph() calls produced two DIFFERENT, coordinate-incompatible
								//reference sequences.
								final Object[] pair=buildConsensusAndGraph(cluster, label);
								if(pair==null || ((byte[])pair[0]).length<MIN_CONSENSUS_LEN){return;}
								consensus=(byte[])pair[0];
								graphs[fci]=(BaseGraph)pair[1];
							}else{
								final byte[] c=buildConsensus(cluster);
								if(c==null || c.length<MIN_CONSENSUS_LEN){return;}
								consensus=c;
							}
							consensi[fci]=consensus;
							labels[fci]=label;
							if(census){censusLines[fci]=censusString(cluster, consensus, label);}
						}
					}));
				}
				waitAll(futures);
				int kept=0;
				for(int ci=0; ci<clusters.size(); ci++){
					final ArrayList<Read> cluster=clusters.get(ci);
					if(consensi[ci]==null){
						//Either below minClusterSize (never attempted, `continue`d above) or a
						//size-eligible cluster whose consensus build failed (line ~188/193) --
						//both mean this cluster never became a real model.
						recordGroupAsSmall(cluster);
						continue;
					}
					if(censusLines[ci]!=null){outstream.println(censusLines[ci]);}
					Read r=new Read(consensi[ci], null, labels[ci], num);
					consensusList.add(r);
					if(modelList!=null){modelList.add(graphs[ci]);}
					recordClusterAsKept(cluster, labels[ci]);
					num++;
					kept++;
				}
				outstream.println(anticodon+": "+group.size()+" seqs -> "+clusters.size()+" clusters, "+kept+" kept");
				totalClusters+=clusters.size();
			}else{
				if(group.size()>=minClusterSize){
					final String label=consensusPrefix+anticodon+" n="+group.size();
					final byte[] consensus;
					BaseGraph bg=null;
					if(modelList!=null){
						//Consensus and HBM built TOGETHER -- see buildConsensusAndGraph's javadoc.
						final Object[] pair=buildConsensusAndGraph(group, label);
						consensus=(pair==null ? null : (byte[])pair[0]);
						if(pair!=null){bg=(BaseGraph)pair[1];}
					}else{
						consensus=buildConsensus(group);
					}
					if(consensus!=null && consensus.length>=MIN_CONSENSUS_LEN){
						Read r=new Read(consensus, null, label, num);
						consensusList.add(r);
						if(modelList!=null){modelList.add(bg);}
						recordClusterAsKept(group, label);
						num++;
					}else{
						recordGroupAsSmall(group);
					}
					outstream.println(anticodon+": "+group.size()+" seqs -> len "+(consensus==null ? 0 : consensus.length));
				}else{
					outstream.println(anticodon+": "+group.size()+" seqs (below minClusterSize="+minClusterSize+")");
					recordGroupAsSmall(group);
				}
				totalClusters++;
			}
		}

		if(pruneCount>0 && consensusList.size()>pruneCount){
			int[] sizes=new int[consensusList.size()];
			for(int i=0; i<sizes.length; i++){
				String id=consensusList.get(i).id;
				int nIdx=id.lastIndexOf("n=");
				sizes[i]=(nIdx>=0 ? Parse.parseInt(id, nIdx+2) : 0);
			}
			int[] sorted=sizes.clone();
			Arrays.sort(sorted);
			final int cutoff=sorted[Tools.min(pruneCount-1, sorted.length-1)];
			ArrayList<Read> kept=new ArrayList<>();
			ArrayList<BaseGraph> keptModels=(modelList!=null ? new ArrayList<>() : null);
			int removed=0;
			for(int i=0; i<consensusList.size(); i++){
				if(sizes[i]<=cutoff && removed<pruneCount){
					removed++;
					//A cluster removed here already went through recordClusterAsKept above --
					//demote its members back to small, since it did NOT actually ship. Matched by
					//clusterKey (family+"_kept_"+its label): assignmentMap has no reverse index
					//from label to members, so this is a scan, only ever run when pruneCount>0.
					demoteClusterKeyToSmall(family+"_kept_"+consensusList.get(i).id);
				}else{
					kept.add(consensusList.get(i));
					if(keptModels!=null){keptModels.add(modelList.get(i));}
				}
			}
			outstream.println("Pruned "+removed+" clusters (cutoff n<="+cutoff+"), "+kept.size()+" remaining.");
			consensusList=kept;
			modelList=keptModels;
		}

		writeAssignments(reads);

		if(ffout!=null){
			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, 4, null, false);
			ros.start();
			ros.add(consensusList, 0);
			errorState|=ReadWrite.closeStream(ros);
		}

		if(modelList!=null && outModel!=null){
			outstream.println("Writing "+modelList.size()+" BaseGraph models to "+outModel);
			if(outModel.endsWith(".txt") || outModel.endsWith(".hbm")){
				writeTextModels(modelList, outModel);
			}else{
				ReadWrite.writeObjectInThread(modelList, outModel, false);
			}
		}

		if(pool!=null){pool.shutdown();}

		t.stop();
		outstream.println();
		outstream.println("Input:  "+reads.size()+" sequences in "+groups.size()+" groups.");
		if(doClustering){outstream.println("Clusters: "+totalClusters+" total.");}
		outstream.println("Output: "+consensusList.size()+" consensus sequences.");
		if(modelList!=null){outstream.println("Models: "+modelList.size()+" BaseGraph HBMs.");}
		outstream.println("Time:   "+t);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private LinkedHashMap<String, ArrayList<byte[]>> loadSeeds(String fname){
		LinkedHashMap<String, ArrayList<byte[]>> map=new LinkedHashMap<>();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start();
		int count=0;
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					String ac=parseAnticodon(r.id);
					if(ac==null){ac="unknown";}
					ArrayList<byte[]> list=map.get(ac);
					if(list==null){list=new ArrayList<>(); map.put(ac, list);}
					list.add(r.bases);
					count++;
				}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		outstream.println("Loaded "+count+" seed consensus sequences in "+map.size()+" anticodon groups.");
		return map;
	}

	private ArrayList<Read> loadReads(){
		ArrayList<Read> reads=new ArrayList<>();
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ffin, null);
		cris.start();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					reads.add(r);
				}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		return reads;
	}

	static LinkedHashMap<String, ArrayList<Read>> groupByAnticodon(ArrayList<Read> reads){
		LinkedHashMap<String, ArrayList<Read>> groups=new LinkedHashMap<>();
		for(Read r : reads){
			String anticodon=parseAnticodon(r.id);
			if(anticodon==null){anticodon="unknown";}
			ArrayList<Read> list=groups.get(anticodon);
			if(list==null){
				list=new ArrayList<>();
				groups.put(anticodon, list);
			}
			list.add(r);
		}
		return groups;
	}

	static String parseAnticodon(String header){
		if(header==null){return null;}
		// Try Note=tRNA-Xxx(YYY) first (RefSeq format)
		int idx=header.indexOf("Note=tRNA-");
		if(idx>=0){
			int paren=header.indexOf('(', idx);
			if(paren>=0){
				int close=header.indexOf(')', paren);
				if(close>paren+1 && close-paren<=5){
					String ac=header.substring(paren+1, close);
					if(!ac.contains(":")){return ac;}
				}
			}
		}
		// Try product=tRNA-Xxx (tRNAscan-SE format) -> amino acid
		idx=header.indexOf("product=tRNA-");
		if(idx>=0){
			int start=idx+13;
			int end=start;
			while(end<header.length() && header.charAt(end)!=';'
					&& header.charAt(end)!=' ' && header.charAt(end)!='\t'){
				end++;
			}
			if(end>start){return header.substring(start, end);}
		}
		// Fallback: look for tRNA-Xxx(YYY) anywhere
		idx=header.indexOf("tRNA-");
		if(idx<0){idx=header.indexOf("trna-");}
		if(idx>=0){
			int paren=header.indexOf('(', idx);
			if(paren>=0){
				int close=header.indexOf(')', paren);
				if(close>paren+1 && close-paren<=5){
					String ac=header.substring(paren+1, close);
					if(!ac.contains(":")){return ac;}
				}
			}
			// No anticodon in parens — extract amino acid name
			int start=idx+5;
			int end=start;
			while(end<header.length() && Character.isLetterOrDigit(header.charAt(end))){
				end++;
			}
			if(end>start){return header.substring(start, end);}
		}
		return null;
	}

	ArrayList<ArrayList<Read>> clusterSequences(ArrayList<Read> group, ArrayList<byte[]> seeds){
		ArrayList<ArrayList<Read>> clusters=new ArrayList<>();
		ArrayList<byte[]> centroids=new ArrayList<>();

		group.sort((a, b)->b.length()-a.length());

		if(seeds!=null && !seeds.isEmpty()){
			// Seeded mode: assign reads to pre-existing centroids only
			centroids.addAll(seeds);
			for(int i=0; i<seeds.size(); i++){clusters.add(new ArrayList<>());}
			for(Read r : group){
				float bestId=0;
				int bestCluster=-1;
				for(int i=0; i<centroids.size(); i++){
					float id=ScrabbleAligner.alignStatic(r.bases, centroids.get(i), null);
					if(id>bestId){bestId=id; bestCluster=i;}
				}
				if(bestId>=recruitIdentity && bestCluster>=0){
					clusters.get(bestCluster).add(r);
				}
			}
		}else{
			// Greedy initial clustering. The per-read decision must stay sequential
			// because centroids accumulate as reads are assigned. The scan over the
			// current, frozen centroid list is independent and can be parallelized.
			final float[] matchOut=new float[2];
			for(Read r : group){
				byte[] seq=r.bases;
				bestCentroidMatch(seq, centroids, matchOut);
				final float bestId=matchOut[0];
				final int bestCluster=(int)matchOut[1];
				if(bestId>=clusterIdentity && bestCluster>=0){
					clusters.get(bestCluster).add(r);
				}else{
					ArrayList<Read> newCluster=new ArrayList<>();
					newCluster.add(r);
					clusters.add(newCluster);
					centroids.add(seq);
				}
			}
		}

		// Step 2: Build consensus for each cluster (optimal centroid)
		// Build consensus for ALL clusters, even small ones — they may
		// grow during reassignment. Size filter applied at final output.
		byte[][] consensusSeqs=new byte[clusters.size()][];
		rebuildConsensi(clusters, consensusSeqs);

		// Step 3: Iteratively reassign all sequences to the best-matching consensus.
		// Cluster choice gets a slight tilt toward the cluster whose median member
		// length is closest to the sequence's own; the pass/fail gate stays on raw
		// identity so the tilt only breaks near-ties, digging length-homogeneous
		// clusters (flush termini at birth) one sand grain at a time.
		HashMap<Read,Integer> assign=new HashMap<>();
		for(int ci=0; ci<clusters.size(); ci++){
			for(Read r : clusters.get(ci)){assign.put(r, ci);}
		}
		ArrayList<Read> orphans=new ArrayList<>();
		for(int round=0; round<reassignRounds; round++){
			if(round>0){//Round 0 uses the consensuses from Step 2
				rebuildConsensi(clusters, consensusSeqs);
			}
			final int[] medians=new int[clusters.size()];
			for(int i=0; i<clusters.size(); i++){medians[i]=medianLength(clusters.get(i));}

			//Score all reads against all consensuses in parallel (pure per-read
			//function of the frozen consensusSeqs/medians), then apply serially.
			final byte[][] consensiF=consensusSeqs;
			final int[] targets=new int[group.size()];
			final float[] bestIds=new float[group.size()];
			final ArrayList<Read> groupF=group;
			ArrayList<Future<?>> futures=new ArrayList<>();
			final int chunk=Tools.max(16, group.size()/(8*Tools.max(1, Shared.threads()))+1);
			for(int start=0; start<group.size(); start+=chunk){
				final int from=start, to=Tools.min(group.size(), start+chunk);
				futures.add(pool().submit(new Runnable(){
					@Override
					public void run(){
						for(int ri=from; ri<to; ri++){
							final Read r=groupF.get(ri);
							float bestScore=0, bestId=0;
							int bestTarget=-1;
							for(int ti=0; ti<consensiF.length; ti++){
								if(consensiF[ti]==null){continue;}
								float id=ScrabbleAligner.alignStatic(r.bases, consensiF[ti], null);
								float score=id;
								if(lenTilt>0 && medians[ti]>0){
									float lenSim=Tools.min(r.length(), medians[ti])/(float)Tools.max(r.length(), medians[ti]);
									score=id*(1-lenTilt*(1-lenSim));
								}
								if(score>bestScore){bestScore=score; bestId=id; bestTarget=ti;}
							}
							targets[ri]=bestTarget;
							bestIds[ri]=bestId;
						}
					}
				}));
			}
			waitAll(futures);

			ArrayList<ArrayList<Read>> newClusters=new ArrayList<>();
			for(int i=0; i<clusters.size(); i++){newClusters.add(new ArrayList<>());}
			orphans.clear();
			int moved=0;
			for(int ri=0; ri<group.size(); ri++){
				final Read r=group.get(ri);
				final Integer prev=assign.get(r);
				if(bestIds[ri]>=recruitIdentity && targets[ri]>=0){
					newClusters.get(targets[ri]).add(r);
					if(prev==null || prev.intValue()!=targets[ri]){moved++;}
					assign.put(r, targets[ri]);
				}else{
					orphans.add(r);
					if(prev!=null){moved++;}
					assign.remove(r);
				}
			}
			clusters=newClusters;
			if(verbose){outstream.println("  Round "+(round+1)+": "+moved+" moved, "+orphans.size()+" orphaned");}
			if(moved==0){break;}
		}

		// Step 4: Recruit orphans using k-mer index
		if(doRecruit && !orphans.isEmpty()){
			// Rebuild consensus after reassignment — include all clusters
			// so small ones can recruit orphans and potentially grow
			rebuildConsensi(clusters, consensusSeqs);

			int numKmers=1<<(2*RECRUIT_K);
			@SuppressWarnings("unchecked")
			ArrayList<Integer>[] lists=new ArrayList[numKmers];
			for(int i=0; i<numKmers; i++){lists[i]=new ArrayList<>();}
			for(int ci=0; ci<consensusSeqs.length; ci++){
				if(consensusSeqs[ci]==null){continue;}
				int kmer=0, len=0;
				for(int j=0; j<consensusSeqs[ci].length; j++){
					int x=AminoAcid.baseToNumber[consensusSeqs[ci][j]];
					if(x>=0){
						kmer=((kmer<<2)|x)&(numKmers-1);
						len++;
						if(len>=RECRUIT_K){lists[kmer].add(ci);}
					}else{len=0; kmer=0;}
				}
			}

			int recruited=0;
			for(Read r : orphans){
				int[] hits=new int[clusters.size()];
				int kmer=0, len=0;
				for(int j=0; j<r.length(); j++){
					int x=AminoAcid.baseToNumber[r.bases[j]];
					if(x>=0){
						kmer=((kmer<<2)|x)&(numKmers-1);
						len++;
						if(len>=RECRUIT_K){
							for(int ci : lists[kmer]){hits[ci]++;}
						}
					}else{len=0; kmer=0;}
				}

				float bestId=0;
				int bestTarget=-1;
				for(int ci=0; ci<hits.length; ci++){
					if(hits[ci]<RECRUIT_MIN_HITS || consensusSeqs[ci]==null){continue;}
					float id=ScrabbleAligner.alignStatic(r.bases, consensusSeqs[ci], null);
					if(id>bestId){bestId=id; bestTarget=ci;}
				}
				if(bestId>=recruitIdentity && bestTarget>=0){
					clusters.get(bestTarget).add(r);
					recruited++;
				}
				//TODO: Probable bug (found during outassignments review, G11/Citan, 2026-08-28) --
				//a Read that fails this recruit check is silently dropped: it is never added to
				//any cluster, never returned to the caller, and has zero accounting anywhere in
				//the pipeline (never reaches consensusList/modelList, never logged individually).
				//Not fixed here -- process()'s new outassignments= tracking now EXPOSES these as
				//orphan rows (one singleton per dropped Read) via a set-difference against the
				//original group, but deliberately does NOT change which Reads end up in the
				//returned `clusters` -- output membership (consensus.fa/model.hbm) is unchanged.
			}
			if(verbose){outstream.println("  Recruited "+recruited+" of "+orphans.size()+" orphans");}
		}

		return clusters;
	}

	/** Below this centroid count, scanning serially is cheaper than farming the
	 * current centroid list across the worker pool. */
	private static final int PARALLEL_CENTROID_THRESHOLD=200;

	/**
	 * Finds the exact same best centroid as the serial ascending-index scan.
	 * The caller appends to centroids only after this method returns, so the
	 * list is read-only while workers scan disjoint ranges. Chunk results are
	 * merged in ascending order with the same strict-'&gt;' comparison, preserving
	 * the serial first-index tie break.
	 */
	private void bestCentroidMatch(final byte[] seq, final ArrayList<byte[]> centroids,
			final float[] out){
		final int n=centroids.size();
		if(n<PARALLEL_CENTROID_THRESHOLD){
			float bestId=0;
			int bestCluster=-1;
			for(int i=0; i<n; i++){
				float id=ScrabbleAligner.alignStatic(seq, centroids.get(i), null);
				if(id>bestId){bestId=id; bestCluster=i;}
			}
			out[0]=bestId;
			out[1]=bestCluster;
			return;
		}

		final int nThreads=Tools.max(1, Shared.threads());
		final int chunk=Tools.max(1, (n+nThreads-1)/nThreads);
		final int nChunks=(n+chunk-1)/chunk;
		final float[] chunkBestId=new float[nChunks];
		final int[] chunkBestIdx=new int[nChunks];
		ArrayList<Future<?>> futures=new ArrayList<>(nChunks);
		for(int c=0; c<nChunks; c++){
			final int from=c*chunk, to=Tools.min(n, from+chunk);
			final int fc=c;
			futures.add(pool().submit(new Runnable(){
				@Override
				public void run(){
					float localBest=0;
					int localIdx=-1;
					for(int i=from; i<to; i++){
						float id=ScrabbleAligner.alignStatic(seq, centroids.get(i), null);
						if(id>localBest){localBest=id; localIdx=i;}
					}
					chunkBestId[fc]=localBest;
					chunkBestIdx[fc]=localIdx;
				}
			}));
		}
		waitAll(futures);

		float bestId=0;
		int bestCluster=-1;
		for(int c=0; c<nChunks; c++){
			if(chunkBestIdx[c]>=0 && chunkBestId[c]>bestId){
				bestId=chunkBestId[c];
				bestCluster=chunkBestIdx[c];
			}
		}
		out[0]=bestId;
		out[1]=bestCluster;
	}

	byte[] buildConsensus(ArrayList<Read> group){
		if(group.isEmpty()){return null;}
		if(group.size()==1){return group.get(0).bases.clone();}

		ArrayList<byte[]> seqs=new ArrayList<>(group.size());
		for(Read r : group){seqs.add(r.bases);}

		byte[] pivot=pickPivot(seqs);
		byte[] consensus=buildFromAlignments(pivot, seqs);
		if(consensus==null){return pivot.clone();}

		for(int pass=1; pass<passes; pass++){
			byte[] refined=buildFromAlignments(consensus, seqs);
			if(refined!=null){consensus=refined;}
		}

		return consensus;
	}

	/**
	 * Flushness census for one kept cluster: how many nt of member 5' sequence
	 * does the final model MISS?  Members longer than the model are aligned
	 * model-into-member; rStart of that alignment is the model's observable 5'
	 * deficit.  Members that fit inside the model cannot expose a deficit and
	 * count as 0.  Caveat: rStart inherits the fill tie-break's late-origin
	 * bias, so absolute deficits are overstated, but comparisons between
	 * libraries censused the same way are valid.
	 */
	private String censusString(ArrayList<Read> cluster, byte[] cons, String label){
		AlignmentStats stats=new AlignmentStats(true);
		int[] deficit=new int[8];//index=nt of missing model 5'; [7] means >=7
		int aligned=0, longer=0;
		for(Read r : cluster){
			stats.clear();
			stats.doTrace=true;
			if(r.length()<=cons.length){
				float id=ScrabbleAligner.alignAndTraceStatic(r.bases, cons, stats);
				if(stats.matchString==null || id<minIdentity){continue;}
				aligned++;
				//KNOWN LIMITATION: a member that fits inside the model cannot expose a
				//model deficit, so it counts as 0 even if the member itself is truncated.
				//This mildly understates deficits but is consistent across all libraries
				//censused with this method — do not change without re-censusing the series.
				deficit[0]++;
			}else{
				float id=ScrabbleAligner.alignAndTraceStatic(cons, r.bases, stats);
				if(stats.matchString==null || id<minIdentity){continue;}
				aligned++; longer++;
				deficit[Tools.min(Tools.max(stats.rStart, 0), 7)]++;
			}
		}
		int modal=0;
		for(int i=1; i<deficit.length; i++){
			if(deficit[i]>deficit[modal]){modal=i;}
		}
		return "CENSUS\t"+label+"\taligned="+aligned+"\tlonger="+longer
			+"\tmodalDeficit="+modal+"\tdeficitHist="+Arrays.toString(deficit);
	}

	/** Rebuilds each cluster's consensus in parallel; empty clusters get null. */
	private void rebuildConsensi(final ArrayList<ArrayList<Read>> clusters, final byte[][] consensusSeqs){
		ArrayList<Future<?>> futures=new ArrayList<>();
		for(int i=0; i<clusters.size(); i++){
			final int fi=i;
			futures.add(pool().submit(new Runnable(){
				@Override
				public void run(){
					consensusSeqs[fi]=(clusters.get(fi).isEmpty() ? null : buildConsensus(clusters.get(fi)));
				}
			}));
		}
		waitAll(futures);
	}

	private ExecutorService pool(){
		if(pool==null){pool=Executors.newFixedThreadPool(Tools.max(1, Shared.threads()));}
		return pool;
	}

	private static void waitAll(ArrayList<Future<?>> futures){
		for(Future<?> f : futures){
			try{f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
	}

	/** Median member length of a cluster; 0 if empty. */
	static int medianLength(ArrayList<Read> list){
		if(list.isEmpty()){return 0;}
		int[] lens=new int[list.size()];
		for(int i=0; i<lens.length; i++){lens[i]=list.get(i).length();}
		Arrays.sort(lens);
		return lens[lens.length/2];
	}

	static byte[] pickPivot(ArrayList<byte[]> seqs){
		//90th-percentile length rather than max: the longest members are exactly
		//the intron-bearing outliers, so max-length selection picks the worst
		//possible pivot in intron families.  Only the first pass is affected;
		//later passes align to the consensus, which is pivot-insensitive.
		int[] lens=new int[seqs.size()];
		for(int i=0; i<lens.length; i++){lens[i]=seqs.get(i).length;}
		int[] sorted=lens.clone();
		Arrays.sort(sorted);
		final int target=sorted[(int)(0.9f*(sorted.length-1))];
		byte[] best=seqs.get(0);
		int bestDif=Integer.MAX_VALUE;
		for(byte[] s : seqs){
			final int dif=Tools.absdif(s.length, target);
			if(dif<bestDif){bestDif=dif; best=s;}
		}
		return best;
	}

	byte[] buildFromAlignments(byte[] ref, ArrayList<byte[]> queries){
		final int refLen=ref.length;
		int[][] counts=new int[refLen][5];
		//Each member's base 0/last base IS its annotated tRNA boundary, so the modal
		//alignment start/stop in ref frame votes for the family's true termini.
		int[] startVotes=new int[refLen], endVotes=new int[refLen];

		AlignmentStats stats=new AlignmentStats(true);
		int aligned=0;
		for(byte[] query : queries){
			stats.clear();
			stats.doTrace=true;
			float identity=ScrabbleAligner.alignAndTraceStatic(query, ref, stats);

			if(stats.matchString==null || identity<minIdentity){continue;}
			aligned++;
			if(stats.rStart>=0 && stats.rStart<refLen){startVotes[stats.rStart]++;}
			if(stats.rStop>=0 && stats.rStop<refLen){endVotes[stats.rStop]++;}

			int refPos=stats.rStart;
			int queryPos=0;
			for(byte b : stats.matchString){
				switch(b){
					case 'm': case 'S': case 'N':
						if(refPos>=0 && refPos<refLen){
							int base=AminoAcid.baseToNumber[query[queryPos]];
							if(base>=0 && base<4){counts[refPos][base]++;}
						}
						refPos++; queryPos++;
						break;
					case 'D':
						if(refPos>=0 && refPos<refLen){counts[refPos][4]++;}
						refPos++;
						break;
					case 'I':
						queryPos++;
						break;
				}
			}
		}

		if(aligned<1){return null;}

		//Truncate the consensus at the modal member boundaries: a depth-fraction cliff
		//proved unreliable (ragged member alignment left ~half of all models 1-4nt
		//5'-truncated), but the members themselves carry the annotated termini, so the
		//most-voted start/stop columns are the family's true ends.
		int first=0, last=refLen-1;
		if(endTrimFrac>0){
			int bs=0, be=refLen-1;
			for(int i=0; i<refLen; i++){
				if(startVotes[i]>startVotes[bs]){bs=i;}
				if(endVotes[i]>endVotes[be]){be=i;}
			}
			if(bs<be){first=bs; last=be;}
		}

		//Citan, 2026-08-28, round 5: a LOCAL Random, not a shared field -- buildFromAlignments is
		//called from multiple threads concurrently (rebuildConsensi/the per-cluster consensus+HBM
		//loop each submit one Runnable per cluster to the pool), and a SHARED Random's internal
		//state gets consumed in THREAD-SCHEDULE order, not data order -- the same (ref,queries)
		//input could draw different random values on different runs depending purely on which
		//other clusters' Runnables happened to interleave first. Seeding this LOCAL instance from
		//a hash of (run seed, ref, ordered query bases) instead makes the draw a pure FUNCTION of
		//what this call is actually processing, so the output is invariant to thread count/
		//scheduling -- see deriveLocalSeed's javadoc for the exact mixing scheme and why length
		//prefixes (not just concatenated bytes) are required for unambiguous boundaries.
		final java.util.Random localRandom=(stochasticConsensus
			? new java.util.Random(deriveLocalSeed(resolvedRunSeed, ref, queries)) : null);

		ByteBuilder bb=new ByteBuilder(refLen);
		for(int i=first; i<=last; i++){
			int total=counts[i][0]+counts[i][1]+counts[i][2]+counts[i][3]+counts[i][4];
			if(total==0){bb.append(ref[i]); continue;}

			int gapCount=counts[i][4];
			int baseCount=total-gapCount;

			if(stochasticConsensus){
				//Proportional-weight trace (Brian, 2026-08-26, ported from Dori source
				//ca74a17062c1988d7219aa18d33dfb778361b7617f4b90206a91ae6e62eed3dc): gap-vs-base,
				//then which base, are both sampled by weight instead of taking the strict
				//majority, so a near-50/50 column becomes a genuine coin flip instead of always
				//resolving to whichever side the tie-break favors. This is what lets a cluster of
				//2-3 sequences produce a genuinely intermediate consensus instead of collapsing to
				//one member -- for a lopsided (near-unanimous) column it degenerates to the same
				//base the greedy path would pick, just probabilistically instead of by comparison.
				if(localRandom.nextInt(total)<gapCount){continue;}
				assert(baseCount>0) : "unreachable: baseCount=0 implies gapCount==total, "
					+"which the branch above always takes; i="+i+", gapCount="+gapCount+", total="+total;
				int r=localRandom.nextInt(baseCount);
				int base=0;
				for(; base<4; base++){
					r-=counts[i][base];
					if(r<0){break;}
				}
				assert(base<4) : "weighted base sample overran counts["+i+"]="+Arrays.toString(counts[i]);
				bb.append(AminoAcid.numberToBase[base]);
				continue;
			}

			if(gapCount>baseCount){continue;}

			int maxCount=0, maxBase=0;
			for(int b=0; b<4; b++){
				if(counts[i][b]>maxCount){
					maxCount=counts[i][b];
					maxBase=b;
				}
			}
			bb.append(AminoAcid.numberToBase[maxBase]);
		}

		return bb.toBytes();
	}

	/**
	 * Derives a per-call local RNG seed from the run seed and this call's DATA (ref + ordered
	 * query bases), not from thread/call-order identity -- this is what makes stochasticConsensus
	 * schedule-independent (Citan, 2026-08-28, round 5): the same (runSeed, ref, queries) always
	 * produces the same seed, regardless of which thread computes it or when.
	 *
	 * <p>FNV-1a-style byte mixing. Every length (ref.length, queries.size(), each query's own
	 * length) is folded in explicitly BEFORE its corresponding bytes -- "unambiguous boundaries":
	 * without an explicit length prefix, two different (ref,queries) splits that happen to
	 * concatenate to the same raw byte stream (e.g. queries=["AB","C"] vs queries=["A","BC"])
	 * could otherwise hash identically.
	 */
	static long deriveLocalSeed(long runSeed, byte[] ref, ArrayList<byte[]> queries){
		long h=0xcbf29ce484222325L;//FNV-1a 64-bit offset basis
		h=mixLong(h, runSeed);
		h=mixLong(h, ref.length);
		h=mixBytes(h, ref);
		h=mixLong(h, queries.size());
		for(byte[] q : queries){
			h=mixLong(h, q.length);
			h=mixBytes(h, q);
		}
		return h;
	}

	private static final long FNV_PRIME=0x100000001b3L;

	private static long mixByte(long h, byte b){
		h^=(b & 0xffL);
		h*=FNV_PRIME;
		return h;
	}

	private static long mixBytes(long h, byte[] arr){
		for(byte b : arr){h=mixByte(h, b);}
		return h;
	}

	private static long mixLong(long h, long v){
		for(int i=0; i<8; i++){h=mixByte(h, (byte)(v>>>(i*8)));}
		return h;
	}

	/**
	 * Builds a mapped Read from an alignment for BaseGraph consumption.
	 * Terminal insertions (unaligned query overhang, emitted by honest glocal
	 * traceback) are trimmed: BaseGraph.add cannot represent an alignment
	 * starting with I, and overhang should not contribute counts or scores.
	 * Terminal D is deliberately NOT stripped here (unlike invertToModelFrame):
	 * in this orientation the glocal traceback never ends in D, and D does not
	 * advance qpos, so no out-of-bounds read is possible.
	 * @return The aligned Read, or null if nothing aligned.
	 */
	static Read toAlignedRead(byte[] bases, AlignmentStats stats, String id, long numericID){
		byte[] match=stats.matchString;
		int lead=0, tail=0;
		while(lead<match.length && match[lead]=='I'){lead++;}
		while(tail<match.length-lead && match[match.length-1-tail]=='I'){tail++;}
		if(lead+tail>0){
			if(match.length-lead-tail<=0){return null;}
			match=Arrays.copyOfRange(match, lead, match.length-tail);
			bases=Arrays.copyOfRange(bases, lead, bases.length-tail);
		}
		Read aligned=new Read(bases, null, id, numericID);
		aligned.match=match;
		aligned.start=stats.rStart;
		aligned.stop=stats.rStop;
		aligned.setMapped(true);
		return aligned;
	}

	/**
	 * Builds a BaseGraph (HBM) from a cluster of tRNA sequences.
	 * Picks pivot as reference, aligns all others with traceback,
	 * and adds each alignment to the graph.
	 */
	static BaseGraph buildBaseGraph(ArrayList<Read> cluster, String name){
		if(cluster.isEmpty()){return null;}

		// Pick longest as pivot/reference
		Read pivotRead=cluster.get(0);
		for(Read r : cluster){
			if(r.length()>pivotRead.length()){pivotRead=r;}
		}
		byte[] pivotBases=pivotRead.bases;

		BaseGraph bg=new BaseGraph(name, pivotBases, null, 0, 0);

		AlignmentStats stats=new AlignmentStats(true);
		for(Read r : cluster){
			if(r.bases==pivotBases){continue;}
			stats.clear();
			stats.doTrace=true;
			float id=ScrabbleAligner.alignAndTraceStatic(r.bases, pivotBases, stats);
			if(stats.matchString==null || id<0.3f){continue;}

			Read aligned=toAlignedRead(r.bases, stats, r.id, r.numericID);
			if(aligned==null){continue;}
			bg.add(aligned);
		}
		consensus.BaseGraphHelper.initForScoring(bg);
		return bg;
	}

	/** Iteration cap for buildConsensusAndGraph's align-to-consensus/rebuild-HBM fixed-point
	 * loop -- a safety bound, not a target; converges in 1-2 iterations in practice for
	 * tRNA-scale clusters (bounded by the same identity/majority-vote dynamics buildConsensus's
	 * own `passes` refinement already exhibits). */
	private static final int MAX_GRAPH_ITERS=5;

	/**
	 * Builds the FINAL consensus AND its BaseGraph (HBM) together, from alignments to the SAME
	 * reference (Brian's ruling, Aug 22 2026, via Noire): the HBM must be built from alignments
	 * to the consensus, not the bootstrap pivot and not (the pre-fix bug) a separately-chosen
	 * "longest member" pivot -- otherwise model.original and the shipped consensus fasta entry
	 * are DIFFERENT sequences in DIFFERENT coordinate frames. Measured on the shipped library
	 * before this fix: only 20.4% of 4,474 models had matching length between
	 * trna_consensus.fa and trna_models.hbm, mean diff 5.98bp, max 202bp -- this caused a real
	 * downstream bug in the boundary-precision NN's fuzziness feature (see
	 * plans/hybrid_two_net_result_20260822.md, Neptune's bedroom).
	 *
	 * <p>Bootstraps from pickPivot() exactly like buildConsensus() (first pass only -- the
	 * pivot itself never appears in the output, per Brian's "pivots don't appear in the
	 * output"). Then iterates: align every cluster member to the CURRENT reference, rebuild the
	 * HBM from those alignments (buildGraphFromReference), and separately rebuild the consensus
	 * from the SAME alignments (buildFromAlignments); if the consensus text changed, the new
	 * consensus becomes next iteration's reference and the HBM is rebuilt again. Stops when the
	 * consensus stabilizes (byte-identical to the previous iteration) or MAX_GRAPH_ITERS is hit.
	 * @return {consensus (byte[]), graph (BaseGraph)}, or null if the cluster produced no usable
	 *   consensus (mirrors buildConsensus's own null case).
	 */
	Object[] buildConsensusAndGraph(ArrayList<Read> cluster, String name){
		if(cluster.isEmpty()){return null;}
		if(cluster.size()==1){
			final byte[] only=cluster.get(0).bases.clone();
			final BaseGraph bg0=new BaseGraph(name, only, null, 0, 0);//degenerate 1-member cluster: nothing to align
			consensus.BaseGraphHelper.initForScoring(bg0);
			return new Object[]{only, bg0};
		}

		final ArrayList<byte[]> seqs=new ArrayList<>(cluster.size());
		for(Read r : cluster){seqs.add(r.bases);}

		byte[] ref=pickPivot(seqs);//bootstrap only -- discarded once the loop below produces a consensus
		byte[] cons=buildFromAlignments(ref, seqs);
		if(cons==null){cons=ref.clone();}//nothing aligned even to the pivot: degenerate, mirrors buildConsensus's own fallback

		//Refine the CONSENSUS ALONE to a fixed point first (cheap -- no graph building per
		//iteration); build the HBM exactly ONCE at the end, from the final stabilized cons.
		//CAUGHT BY TESTING (Aug 22 2026, not by inspection): an earlier version of this loop
		//rebuilt the graph EVERY iteration from the CURRENT cons, then updated cons to the
		//newly-refined value for the NEXT iteration -- so whenever the loop exited by hitting
		//MAX_GRAPH_ITERS (rather than by convergence), the returned bg was built from the
		//SECOND-TO-LAST cons while the returned cons was the LAST (one iteration newer),
		//silently reintroducing exactly the coordinate-frame mismatch this method exists to
		//fix. Verified on a real 2,000-sequence smoke corpus: only 2/33 clusters had
		//byte-identical consensus/model sequences under the buggy version (the 2 that happened
		//to converge before hitting the cap); this version fixes that by building the graph
		//only after the consensus loop is fully done, so bg and cons can never desync.
		int iter=0;
		for(; iter<MAX_GRAPH_ITERS; iter++){
			final byte[] refined=buildFromAlignments(cons, seqs);
			if(refined==null || Arrays.equals(refined, cons)){break;}//stable, or nothing aligned this round
			cons=refined;
		}
		if(verbose && iter>=MAX_GRAPH_ITERS){
			outstream.println("WARNING: "+name+" did not stabilize within "+MAX_GRAPH_ITERS+" iterations -- "
				+"using the last computed consensus; the HBM below is still built to MATCH it exactly.");
		}
		final BaseGraph bg=buildGraphFromReference(cluster, cons, name);
		return new Object[]{cons, bg};
	}

	/** Builds a BaseGraph from alignments of every cluster member to the given EXTERNAL
	 * reference (the consensus, once buildConsensusAndGraph has one) -- the fixed-reference-vs-
	 * longest-member-as-its-own-pivot distinction is exactly the fix buildConsensusAndGraph
	 * exists for. Mirrors buildBaseGraph's own alignment/add loop exactly, parameterized by ref
	 * instead of internally picking the longest member as its own pivot. */
	private static BaseGraph buildGraphFromReference(ArrayList<Read> cluster, byte[] ref, String name){
		final BaseGraph bg=new BaseGraph(name, ref, null, 0, 0);
		final AlignmentStats stats=new AlignmentStats(true);
		for(Read r : cluster){
			stats.clear();
			stats.doTrace=true;
			final float id=ScrabbleAligner.alignAndTraceStatic(r.bases, ref, stats);
			if(stats.matchString==null || id<0.3f){continue;}
			final Read aligned=toAlignedRead(r.bases, stats, r.id, r.numericID);
			if(aligned==null){continue;}
			bg.add(aligned);
		}
		consensus.BaseGraphHelper.initForScoring(bg);
		return bg;
	}

	/**
	 * Scores a candidate sequence against a BaseGraph model.
	 * Aligns with traceback, then uses BaseGraph.score() for
	 * position-weighted likelihood.  The glocal aligner requires query<=ref,
	 * so a candidate longer than the model reference is aligned in the swapped
	 * orientation and the traceback inverted back to model coordinates.
	 */
	public static float scoreAgainstModel(byte[] candidate, BaseGraph model){
		if(model==null || candidate==null){return -999;}
		final byte[] cons=model.original;
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		final Read r;
		if(candidate.length<=cons.length){
			float id=ScrabbleAligner.alignAndTraceStatic(candidate, cons, stats);
			if(stats.matchString==null || id<0.2f){return -999;}
			r=toAlignedRead(candidate, stats, "candidate", 0);
		}else{
			float id=ScrabbleAligner.alignAndTraceStatic(cons, candidate, stats);
			if(stats.matchString==null || id<0.2f){return -999;}
			r=invertToModelFrame(candidate, stats);
		}
		if(r==null){return -999;}
		return model.score(r, false, true);
	}

	/**
	 * Converts a consensus-as-query alignment (candidate as the reference) into
	 * a candidate-as-query Read in model coordinates for BaseGraph scoring.
	 * I and D swap roles, the candidate span rStart..rStop supplies the bases,
	 * and terminal insertions the inversion exposes are trimmed (BaseGraph
	 * alignments cannot start with I).  The consensus is consumed fully by the
	 * glocal alignment, so model coverage always spans the whole model.
	 * @return The inverted Read, or null if the alignment is unusable.
	 */
	static Read invertToModelFrame(byte[] candidate, AlignmentStats stats){
		final byte[] match=stats.matchString;
		if(stats.rStart<0 || stats.rStop<stats.rStart || stats.rStop>=candidate.length){return null;}
		byte[] inv=new byte[match.length];
		for(int i=0; i<match.length; i++){
			final byte b=match[i];
			inv[i]=(b=='D' ? (byte)'I' : b=='I' ? (byte)'D' : b);
		}
		byte[] bases=Arrays.copyOfRange(candidate, stats.rStart, stats.rStop+1);
		//Strip terminal indels: terminal I is candidate overhang (consumes a base),
		//terminal D is model overhang the candidate never spans (advances the model
		//coordinate).  Neither should contribute counts or scores, and BaseGraph's
		//walk reads bases[qpos] on EVERY op, so a trailing D after the last base
		//would read past the end.  The first/last kept op is always m/S/N.
		int from=0, to=inv.length, bFrom=0, bTo=bases.length, startPos=0;
		while(from<to && (inv[from]=='I' || inv[from]=='D')){
			if(inv[from]=='I'){bFrom++;}else{startPos++;}
			from++;
		}
		while(to>from && (inv[to-1]=='I' || inv[to-1]=='D')){
			if(inv[to-1]=='I'){bTo--;}
			to--;
		}
		if(to-from<1 || bTo-bFrom<1){return null;}
		if(from>0 || to<inv.length){inv=Arrays.copyOfRange(inv, from, to);}
		if(bFrom>0 || bTo<bases.length){bases=Arrays.copyOfRange(bases, bFrom, bTo);}
		int modelSpan=0;
		for(byte b : inv){if(b!='I'){modelSpan++;}}
		if(modelSpan<1){return null;}
		Read aligned=new Read(bases, null, "candidate", 0);
		aligned.match=inv;
		aligned.start=startPos;
		aligned.stop=startPos+modelSpan-1;
		aligned.setMapped(true);
		return aligned;
	}

	/**
	 * Scores a candidate against all BaseGraph models and returns the best score.
	 */
	public static float bestModelScore(byte[] candidate, BaseGraph[] models){
		if(models==null || models.length==0 || candidate==null){return -999;}
		float best=-999;
		for(BaseGraph model : models){
			float score=scoreAgainstModel(candidate, model);
			if(score>best){best=score;}
		}
		return best;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@SuppressWarnings("unchecked")
	public static BaseGraph[] loadModels(String fname){
		if(fname==null){return null;}
		if(fname.endsWith(".txt") || fname.endsWith(".hbm")){
			return loadTextModels(fname);
		}
		ArrayList<BaseGraph> list=(ArrayList<BaseGraph>)ReadWrite.readObject(fname, false);
		return list.toArray(new BaseGraph[0]);
	}

	public static void writeTextModels(ArrayList<BaseGraph> models, String fname){
		ByteBuilder bb=new ByteBuilder();
		bb.append("#HBM\t1\n");
		for(BaseGraph bg : models){
			int len=consensus.BaseGraphHelper.length(bg);
			bb.append('>').append(bg.name).append('\t');
			bb.append("depth=").append(bg.ref[0].countSum+bg.del[0].countSum).nl();
			bb.append(bg.original).nl();
			for(int i=0; i<len; i++){
				int[] c=consensus.BaseGraphHelper.getCounts(bg, i);
				bb.append(c[0]).tab().append(c[1]).tab().append(c[2]).tab().append(c[3]).tab().append(c[4]).nl();
			}
		}
		ReadWrite.writeString(bb.toString(), fname);
	}

	public static BaseGraph[] loadTextModels(String fname){
		ArrayList<BaseGraph> list=new ArrayList<>();
		String[] lines=ReadWrite.readString(fname).split("\n");
		int i=0;
		while(i<lines.length){
			String line=lines[i];
			if(line.startsWith("#")){i++; continue;}
			if(!line.startsWith(">")){i++; continue;}
			String name=line.substring(1).split("\t")[0];
			i++;
			if(i>=lines.length){break;}
			byte[] bases=lines[i].getBytes();
			Tools.toUpperCase(bases);
			i++;
			BaseGraph bg=new BaseGraph(name, bases, null, list.size(), 0);
			for(int pos=0; pos<bases.length && i<lines.length; pos++, i++){
				String[] parts=lines[i].split("\t");
				if(parts.length<5){break;}
				int a=Integer.parseInt(parts[0]);
				int c=Integer.parseInt(parts[1]);
				int g=Integer.parseInt(parts[2]);
				int t=Integer.parseInt(parts[3]);
				int del=Integer.parseInt(parts[4]);
				consensus.BaseGraphHelper.setCounts(bg, pos, a, c, g, t, del);
			}
			consensus.BaseGraphHelper.initForScoring(bg);
			list.add(bg);
		}
		return list.toArray(new BaseGraph[0]);
	}

	public static String[] lastLoadedNames;

	public static byte[][] loadLibrary(String fname){
		if(fname==null){return null;}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start();
		ArrayList<byte[]> list=new ArrayList<>();
		ArrayList<String> names=new ArrayList<>();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					list.add(r.bases);
					names.add(r.id);
				}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		lastLoadedNames=names.toArray(new String[0]);
		return list.toArray(new byte[0][]);
	}

	/**
	 * Aligns a candidate tRNA sequence against all reference sequences
	 * and returns the best identity.
	 */
	public static float bestIdentity(byte[] candidate, byte[][] library){
		if(library==null || library.length==0 || candidate==null){return 0;}
		float best=0;
		for(byte[] ref : library){
			float id=ScrabbleAligner.alignStatic(candidate, ref, null);
			if(id>best){best=id;}
		}
		return best;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------   Assignment Tracking (a)    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Records (or, for the pruneCount post-process, OVERWRITES) one Read's clusterKey/
	 * modelLabel/status. No-op when outAssignments==null -- callers may call this
	 * unconditionally without their own null-check.
	 *
	 * <p>clusterKey scheme (Citan, 2026-08-28), every branch family-namespaced so two families'
	 * assignment files can be safely combined without collision:
	 * <ul>
	 * <li>kept: family+"_kept_"+modelLabel -- modelLabel itself is NOT family-namespaced
	 *     (TrnaConsensusBuilder.java's own label format is "tRNA_consensus_"+anticodon+"_c"+fci
	 *     regardless of family), so two different families could otherwise produce a
	 *     byte-identical modelLabel and collide.</li>
	 * <li>small: family+"_small_"+ the lexicographically smallest fullRecordId among that
	 *     cluster/group's ACTUAL members (content-derived, not an iteration index -- two
	 *     distinct clusters can never share a member, so this is collision-free by
	 *     construction).</li>
	 * <li>orphan: family+"_orphan_"+its own fullRecordId -- a true singleton, one per record,
	 *     never grouped with any other orphan (an orphan's "cluster" has exactly 1 distinct
	 *     source genome by definition, which is what correctly forces a genome carrying an
	 *     orphan gene to stay coverage-critical under any K>=1 redundancy gate downstream).</li>
	 * </ul>
	 * modelLabel column is the real shipped label for kept, or an explicit EMPTY string (never
	 * the text "null") for small/orphan.
	 */
	private void recordAssignment(Read r, String clusterKey, String modelLabel, String status){
		if(assignmentMap==null){return;}
		assignmentMap.put(r, new String[]{clusterKey, (modelLabel==null ? "" : modelLabel), status});
	}

	/** Records every member of `members` as SMALL, sharing one clusterKey derived from their own
	 * lexicographically-smallest fullRecordId (see recordAssignment). Correction (Citan,
	 * 2026-08-28): NOT used for the below-minGroupSize skip any more -- that case has no
	 * verified redundancy among its members at all (the whole group was never sequence-
	 * clustered) and now goes through recordEachAsOrphan instead. Actual call sites: (i) a
	 * cluster clusterSequences() ACTUALLY returned that fell below minClusterSize (genuine,
	 * algorithm-verified redundancy among its own members, just not enough of them); (ii) a
	 * size-eligible cluster (or, in no-clustering mode, the whole builderGroup treated as its
	 * own single build unit) whose consensus build failed despite meeting the size threshold;
	 * (iii) in no-clustering mode, a builderGroup below minClusterSize, treated as its own build
	 * unit since no-clustering mode never sub-clusters it in the first place. */
	private void recordGroupAsSmall(ArrayList<Read> members){
		if(assignmentMap==null || members.isEmpty()){return;}
		final String minId=minRecordId(members);
		final String clusterKey=family+"_small_"+minId;
		for(Read r : members){recordAssignment(r, clusterKey, null, "small");}
	}

	/** Records every member of one cluster as KEPT, sharing the cluster's real shipped label
	 * both as the modelLabel column and (family-namespaced) as its clusterKey. */
	private void recordClusterAsKept(ArrayList<Read> members, String modelLabel){
		if(assignmentMap==null){return;}
		final String clusterKey=family+"_kept_"+modelLabel;
		for(Read r : members){recordAssignment(r, clusterKey, modelLabel, "kept");}
	}

	/** Citan, 2026-08-28: fails loud before any assignment tracking begins if two different
	 * Read objects share the same fullRecordId -- the whole clusterKey scheme (orphan
	 * singletons, small-cluster min-id) assumes fullRecordId is a stable, collision-free
	 * per-record identity; two different records sharing one ID would silently make an
	 * orphan-singleton key non-unique (merging two unrelated records' redundancy under one
	 * key) or corrupt the deterministic min-id computation for a small cluster. */
	private void checkNoDuplicateRecordIds(ArrayList<Read> reads){
		final java.util.HashSet<String> seen=new java.util.HashSet<>();
		for(Read r : reads){
			if(!seen.add(r.id)){
				throw new RuntimeException("Error - duplicate fullRecordId found: '"+r.id+"' -- "
					+"outassignments= requires every record's id to be unique.");
			}
		}
	}

	/** Records every member of a group as its OWN singleton orphan -- used where no verified
	 * redundancy relationship exists among the members at all: specifically the below-
	 * minGroupSize skip, where the whole group was never sequence-clustered and never even
	 * became the no-clustering mode's own single build unit (Citan, 2026-08-28). Unlike
	 * recordGroupAsSmall, whose 3 documented call sites (see its own javadoc) all share one
	 * real, checked build unit -- either an actual clusterSequences() cluster or, in
	 * no-clustering mode, the builderGroup itself treated as its own attempted build. */
	private void recordEachAsOrphan(ArrayList<Read> members){
		if(assignmentMap==null){return;}
		for(Read r : members){recordAssignment(r, family+"_orphan_"+r.id, null, "orphan");}
	}

	/** Set-difference against everything clusterSequences() actually returned. Citan, 2026-08-28:
	 * this is a GENERAL catch-all, not specific to Step 4's k-mer recruit failure (the TODO:
	 * Probable bug comment there documents one concrete mechanism, but is not the only one this
	 * method covers) -- any Read absent from every returned cluster is caught here regardless of
	 * WHY, including: Step-4 k-mer-recruit failure (doRecruit=true), a Step-3 orphan when
	 * doRecruit=false (never gets a Step-4 attempt at all, dropped immediately), and seeded mode
	 * (clusterSequences' seeds!=null branch has no orphan bookkeeping whatsoever -- a Read that
	 * fails every seed's recruitIdentity check is simply never added anywhere). This is the only
	 * place that reconstructs which original group members those were, for any of these paths.
	 * Each gets recorded as its own singleton orphan (see recordAssignment) -- deliberately does
	 * NOT alter `clusters` or any consensus/model output. */
	private void recordDroppedOrphans(ArrayList<Read> group, ArrayList<ArrayList<Read>> clusters){
		if(assignmentMap==null){return;}
		final java.util.IdentityHashMap<Read,Boolean> placed=new java.util.IdentityHashMap<>();
		for(ArrayList<Read> cluster : clusters){
			for(Read r : cluster){placed.put(r, Boolean.TRUE);}
		}
		final ArrayList<Read> dropped=new ArrayList<>();
		for(Read r : group){
			if(!placed.containsKey(r)){dropped.add(r);}
		}
		recordEachAsOrphan(dropped);
	}

	/** Demotes every record currently assigned to the given (kept-style) clusterKey back to
	 * small -- used by process()'s pruneCount post-process, which can remove an already-built
	 * cluster from the final shipped output after recordClusterAsKept already ran for it. Citan,
	 * 2026-08-28: the demoted clusterKey MUST follow the declared small-key schema exactly
	 * (family+"_small_"+lexicographically-smallest member fullRecordId) -- derived from the
	 * matching assignmentMap entries' own Read.id fields, not a special "(pruned)" variant, so a
	 * downstream manifest reader never needs to special-case this path. */
	private void demoteClusterKeyToSmall(String prunedClusterKey){
		if(assignmentMap==null){return;}
		String minId=null;
		for(java.util.Map.Entry<Read,String[]> e : assignmentMap.entrySet()){
			if(e.getValue()[0].equals(prunedClusterKey)){
				final String id=e.getKey().id;
				if(minId==null || id.compareTo(minId)<0){minId=id;}
			}
		}
		assert(minId!=null) : "demoteClusterKeyToSmall called with a clusterKey matching no "
			+"assignmentMap entries: "+prunedClusterKey;
		final String demotedKey=family+"_small_"+minId;
		for(String[] a : assignmentMap.values()){
			if(a[0].equals(prunedClusterKey)){
				a[0]=demotedKey;
				a[1]="";
				a[2]="small";
			}
		}
	}

	private static String minRecordId(ArrayList<Read> members){
		String min=null;
		for(Read r : members){
			if(min==null || r.id.compareTo(min)<0){min=r.id;}
		}
		return min;
	}

	/** Writes the outAssignments TSV: fullRecordId, sourceTid, builderGroup, clusterKey,
	 * modelLabel, status. Iterates the ORIGINAL input `reads` list (not assignmentMap's own
	 * iteration order) so row order is deterministic and independent of IdentityHashMap
	 * internals. Asserts exactly one assignment per input record -- a gap here means a code
	 * path added a new way to lose or duplicate a Read that this method doesn't yet know about,
	 * and that must fail loud, not silently under- or over-count. Fails loud (not skip-and-
	 * count, unlike this file's tolerant flank/gc conventions elsewhere -- this is the record's
	 * whole downstream identity) on an unparseable tid or a literal tab in any written field
	 * (which would corrupt the TSV structure). */
	private void writeAssignments(ArrayList<Read> reads){
		if(outAssignments==null){return;}
		//Citan, 2026-08-28: an independent check on assignmentMap's own populated size BEFORE
		//iterating -- the loop's own written==reads.size() (kept below) only proves the loop ran
		//once per `reads` entry, which is tautological (the loop always does that unless it
		//throws first) and would NOT catch assignmentMap holding extra/stale entries beyond what
		//`reads` accounts for.
		assert(assignmentMap.size()==reads.size()) : "assignmentMap.size()="+assignmentMap.size()
			+" != reads.size()="+reads.size()+" -- some code path populated assignmentMap with a "
			+"different record count than was actually loaded.";
		final structures.ByteBuilder bb=new structures.ByteBuilder();
		bb.append("#fullRecordId").tab().append("sourceTid").tab().append("builderGroup").tab()
			.append("clusterKey").tab().append("modelLabel").tab().append("status").nl();
		int written=0;
		for(Read r : reads){
			final String[] a=assignmentMap.get(r);
			assert(a!=null) : "No assignment recorded for record '"+r.id+"' -- exactly one row per "
				+"input record is required; this means a code path can lose a Read without this "
				+"method knowing about it.";
			final String fullRecordId=r.id;
			final int tid=tax.TaxTree.parseTaxID(fullRecordId);
			if(tid<1){
				throw new RuntimeException("Error - could not parse a source tid from record '"
					+fullRecordId+"' -- outassignments= requires every record to carry a "
					+"parseable tid_NNNN (or tid|/ncbi|/ncbi_) token.");
			}
			final String builderGroup=parseAnticodon(fullRecordId);
			final String group=(builderGroup==null ? "unknown" : builderGroup);
			final String clusterKey=a[0], modelLabel=a[1], status=a[2];
			for(String field : new String[]{fullRecordId, group, clusterKey, modelLabel, status}){
				if(field.indexOf('\t')>=0){
					throw new RuntimeException("Error - field contains a literal tab, which would "
						+"corrupt the outassignments TSV: '"+field+"'");
				}
			}
			bb.append(fullRecordId).tab().append(tid).tab().append(group).tab()
				.append(clusterKey).tab().append(modelLabel).tab().append(status).nl();
			written++;
		}
		assert(written==reads.size()) : "written="+written+" != reads.size()="+reads.size();
		final fileIO.ByteStreamWriter bsw=new fileIO.ByteStreamWriter(outAssignments, overwrite, false, false);
		bsw.start();
		bsw.print(bb);
		errorState|=bsw.poisonAndWait();
		outstream.println("Wrote "+written+" assignment rows to "+outAssignments);
	}

	private String in;
	private String out;
	private final FileFormat ffin;
	private final FileFormat ffout;
	private boolean overwrite=true;
	private boolean verbose=false;
	private boolean errorState=false;
	private int passes=2;
	private int minGroupSize=1;
	private float minIdentity=0.3f;
	private float clusterIdentity=0.70f;
	private boolean doClustering=true;
	private boolean doRecruit=true;
	private int minClusterSize=3;
	private float recruitIdentity=0.70f;
	/** Truncate consensus at the modal member start/stop boundaries; 0 disables */
	private float endTrimFrac=0.3f;
	/** Slight reassignment incentive toward the cluster with the closest median
	 * member length; multiplies identity by (1-lenTilt*(1-lenSim)). 0 disables. */
	private float lenTilt=0.03f;
	/** Reassignment iterations, rebuilding consensus between rounds; early exit when no moves */
	private int reassignRounds=2;
	/** Print a per-cluster 5'-flushness census line for each kept model */
	private boolean census=false;
	/** Worker pool for per-cluster consensus/HBM/census tasks; sized by Shared.threads() */
	private ExecutorService pool;
	private String seedIn;
	private int pruneCount=0;
	private String outModel;
	/** Model-label prefix (default matches the pre-existing hardcoded literal exactly, so the
	 * default path is byte-invariant). */
	private String consensusPrefix="tRNA_consensus_";
	/** Proportional-weight (coin-flip) consensus trace instead of strict majority; see
	 * buildFromAlignments. Off by default -- byte-invariant to the deterministic majority-vote
	 * path when false. */
	private boolean stochasticConsensus=false;
	/** Explicit seed for stochasticConsensus; -1 (default) picks one random run seed once, at
	 * construction time (see resolvedRunSeed). */
	private long stochasticSeed=-1;
	/** The ACTUAL seed driving every buildFromAlignments call this run -- resolved exactly once,
	 * single-threaded, at construction time (Citan, 2026-08-28, round 5): either stochasticSeed
	 * itself (if &gt;=0) or one freshly chosen random value (if &lt;0). Every call derives its own
	 * local seed from this PLUS its own (ref,queries) content, so results are schedule-independent
	 * even in "random seed" mode -- the run-level randomness is chosen once, not per-call. */
	private long resolvedRunSeed;
	private PrintStream outstream=System.err;

	//Boundary-NN manifest instrumentation (Citan/Noire/Brian, 2026-08-28) -- opt-in, off by
	//default, byte-invariant to the existing consensus.fa/model.hbm/census outputs when null
	//(see TrnaConsensusBuilderAssignmentsTest's byte-invariance test). See recordAssignment's
	//javadoc for the clusterKey scheme.
	private String outAssignments=null;
	private String family=null;
	/** Read identity -> {clusterKey, modelLabel, status}, populated only when outAssignments!=
	 * null. IdentityHashMap deliberately: Read has no meaningful equals/hashCode override here
	 * and the SAME Read objects are threaded through group/cluster lists without being copied,
	 * so identity is the correct and only correct key. */
	private java.util.IdentityHashMap<Read, String[]> assignmentMap;

	private static final int MIN_CONSENSUS_LEN=50;
	private static final int RECRUIT_K=5;
	private static final int RECRUIT_MIN_HITS=3;
}
