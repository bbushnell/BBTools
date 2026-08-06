package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

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

		LinkedHashMap<String, ArrayList<Read>> groups=groupByAnticodon(reads);
		outstream.println("Found "+groups.size()+" anticodon groups.");

		ArrayList<Read> consensusList=new ArrayList<>();
		ArrayList<BaseGraph> modelList=(outModel!=null ? new ArrayList<>() : null);
		long num=0;
		int totalClusters=0;
		for(Entry<String, ArrayList<Read>> e : groups.entrySet()){
			String anticodon=e.getKey();
			ArrayList<Read> group=e.getValue();
			if(group.size()<minGroupSize){
				if(verbose){outstream.println("Skipping "+anticodon+" ("+group.size()+" sequences, below min="+minGroupSize+")");}
				continue;
			}

			if(doClustering && group.size()>1){
				ArrayList<ArrayList<Read>> clusters=clusterSequences(group);
				int kept=0;
				for(int ci=0; ci<clusters.size(); ci++){
					ArrayList<Read> cluster=clusters.get(ci);
					if(cluster.size()<minClusterSize){continue;}
					byte[] consensus=buildConsensus(cluster);
					if(consensus!=null && consensus.length>=MIN_CONSENSUS_LEN){
						String label="tRNA_consensus_"+anticodon+"_c"+ci+" n="+cluster.size();
						Read r=new Read(consensus, null, label, num);
						consensusList.add(r);
						if(modelList!=null){
							BaseGraph bg=buildBaseGraph(cluster, label);
							modelList.add(bg);
						}
						num++;
						kept++;
					}
				}
				outstream.println(anticodon+": "+group.size()+" seqs -> "+clusters.size()+" clusters, "+kept+" kept");
				totalClusters+=clusters.size();
			}else{
				byte[] consensus=buildConsensus(group);
				if(consensus!=null && consensus.length>=MIN_CONSENSUS_LEN){
					String label="tRNA_consensus_"+anticodon+" n="+group.size();
					Read r=new Read(consensus, null, label, num);
					consensusList.add(r);
					if(modelList!=null && group.size()>=minClusterSize){
						BaseGraph bg=buildBaseGraph(group, label);
						modelList.add(bg);
					}
					num++;
				}
				outstream.println(anticodon+": "+group.size()+" seqs -> len "+(consensus==null ? 0 : consensus.length));
				totalClusters++;
			}
		}

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
		if(idx<0){idx=header.indexOf("Note=tRNA-");}
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
		if(idx<0){idx=header.indexOf("product=tRNA-");}
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

	ArrayList<ArrayList<Read>> clusterSequences(ArrayList<Read> group){
		// Step 1: Greedy initial clustering
		ArrayList<ArrayList<Read>> clusters=new ArrayList<>();
		ArrayList<byte[]> centroids=new ArrayList<>();

		group.sort((a, b)->b.length()-a.length());

		for(Read r : group){
			byte[] seq=r.bases;
			float bestId=0;
			int bestCluster=-1;
			for(int i=0; i<centroids.size(); i++){
				float id=ScrabbleAligner.alignStatic(seq, centroids.get(i), null);
				if(id>bestId){
					bestId=id;
					bestCluster=i;
				}
			}
			if(bestId>=clusterIdentity && bestCluster>=0){
				clusters.get(bestCluster).add(r);
			}else{
				ArrayList<Read> newCluster=new ArrayList<>();
				newCluster.add(r);
				clusters.add(newCluster);
				centroids.add(seq);
			}
		}

		// Step 2: Build consensus for each cluster (optimal centroid)
		byte[][] consensusSeqs=new byte[clusters.size()][];
		for(int i=0; i<clusters.size(); i++){
			if(clusters.get(i).size()>=minClusterSize){
				consensusSeqs[i]=buildConsensus(clusters.get(i));
			}
		}

		// Step 3: Reassign all sequences to best-matching consensus
		ArrayList<ArrayList<Read>> newClusters=new ArrayList<>();
		for(int i=0; i<clusters.size(); i++){
			newClusters.add(new ArrayList<>());
		}
		ArrayList<Read> orphans=new ArrayList<>();
		int moved=0;

		for(int ci=0; ci<clusters.size(); ci++){
			for(Read r : clusters.get(ci)){
				float bestId=0;
				int bestTarget=-1;
				for(int ti=0; ti<consensusSeqs.length; ti++){
					if(consensusSeqs[ti]==null){continue;}
					float id=ScrabbleAligner.alignStatic(r.bases, consensusSeqs[ti], null);
					if(id>bestId){bestId=id; bestTarget=ti;}
				}
				if(bestId>=recruitIdentity && bestTarget>=0){
					newClusters.get(bestTarget).add(r);
					if(bestTarget!=ci){moved++;}
				}else{
					orphans.add(r);
				}
			}
		}
		clusters=newClusters;
		if(verbose){outstream.println("  Reassigned: "+moved+" moved, "+orphans.size()+" orphaned");}

		// Step 4: Recruit orphans using k-mer index
		if(doRecruit && !orphans.isEmpty()){
			// Rebuild consensus after reassignment
			for(int i=0; i<clusters.size(); i++){
				if(clusters.get(i).size()>=minClusterSize){
					consensusSeqs[i]=buildConsensus(clusters.get(i));
				}else{
					consensusSeqs[i]=null;
				}
			}

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
			}
			if(verbose){outstream.println("  Recruited "+recruited+" of "+orphans.size()+" orphans");}
		}

		return clusters;
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

	static byte[] pickPivot(ArrayList<byte[]> seqs){
		byte[] best=seqs.get(0);
		for(int i=1; i<seqs.size(); i++){
			if(seqs.get(i).length>best.length){
				best=seqs.get(i);
			}
		}
		return best;
	}

	byte[] buildFromAlignments(byte[] ref, ArrayList<byte[]> queries){
		final int refLen=ref.length;
		int[][] counts=new int[refLen][5];

		AlignmentStats stats=new AlignmentStats(true);
		int aligned=0;
		for(byte[] query : queries){
			stats.clear();
			stats.doTrace=true;
			float identity=ScrabbleAligner.alignAndTraceStatic(query, ref, stats);

			if(stats.matchString==null || identity<minIdentity){continue;}
			aligned++;

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

		ByteBuilder bb=new ByteBuilder(refLen);
		for(int i=0; i<refLen; i++){
			int total=counts[i][0]+counts[i][1]+counts[i][2]+counts[i][3]+counts[i][4];
			if(total==0){bb.append(ref[i]); continue;}

			int gapCount=counts[i][4];
			int baseCount=total-gapCount;
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

			Read aligned=new Read(r.bases, null, r.id, r.numericID);
			aligned.match=stats.matchString;
			aligned.start=stats.rStart;
			aligned.stop=stats.rStop;
			aligned.setMapped(true);
			bg.add(aligned);
		}
		consensus.BaseGraphHelper.initForScoring(bg);
		return bg;
	}

	/**
	 * Scores a candidate sequence against a BaseGraph model.
	 * Aligns with traceback, then uses BaseGraph.score() for
	 * position-weighted likelihood.
	 */
	public static float scoreAgainstModel(byte[] candidate, BaseGraph model){
		if(model==null || candidate==null){return -999;}
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		float id=ScrabbleAligner.alignAndTraceStatic(candidate, model.original, stats);
		if(stats.matchString==null || id<0.2f){return -999;}

		Read r=new Read(candidate, null, "candidate", 0);
		r.match=stats.matchString;
		r.start=stats.rStart;
		r.stop=stats.rStop;
		r.setMapped(true);
		try{
			return model.score(r, false, true);
		}catch(AssertionError e){
			return id;
		}
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
	private String outModel;
	private PrintStream outstream=System.err;

	private static final int MIN_CONSENSUS_LEN=50;
	private static final int RECRUIT_K=5;
	private static final int RECRUIT_MIN_HITS=3;
}
