package prok;

import java.util.ArrayList;

import consensus.BaseGraph;
import idaligner.AlignmentStats;
import idaligner.ScrabbleAligner;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Quick test: build HBMs for one anticodon, score E.coli tRNA candidates.
 * Usage: java prok.TestHBM trnas=all_trnas.fa anticodon=GUC genome=ecoli.fna.gz
 */
public class TestHBM {

	public static void main(String[] args) throws Exception {
		String trnaFile=null, anticodon=null, genomeFile=null;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("trnas") || a.equals("in")){trnaFile=b;}
			else if(a.equals("anticodon") || a.equals("ac")){anticodon=b;}
			else if(a.equals("genome")){genomeFile=b;}
		}
		if(trnaFile==null || anticodon==null || genomeFile==null){
			System.err.println("Usage: TestHBM trnas=<fa> anticodon=<XXX> genome=<fna.gz>");
			return;
		}

		Timer t=new Timer();

		// Load and filter tRNAs for this anticodon
		System.err.println("Loading tRNAs for anticodon "+anticodon+"...");
		ArrayList<Read> allReads=loadFasta(trnaFile);
		ArrayList<Read> filtered=new ArrayList<>();
		for(Read r : allReads){
			String ac=TrnaConsensusBuilder.parseAnticodon(r.id);
			if(anticodon.equals(ac)){
				Tools.toUpperCase(r.bases);
				filtered.add(r);
			}
		}
		System.err.println("Found "+filtered.size()+" sequences for "+anticodon);

		// Cluster
		System.err.println("Clustering at 70% identity...");
		TrnaConsensusBuilder tcb=new TrnaConsensusBuilder(
			new String[]{"in="+trnaFile, "out=/dev/null"});
		ArrayList<ArrayList<Read>> clusters=tcb.clusterSequences(filtered, null);
		System.err.println("Found "+clusters.size()+" clusters");

		// Build BaseGraph for each cluster with >=3 members
		ArrayList<BaseGraph> models=new ArrayList<>();
		for(int i=0; i<clusters.size(); i++){
			ArrayList<Read> cluster=clusters.get(i);
			if(cluster.size()<3){continue;}
			BaseGraph bg=TrnaConsensusBuilder.buildBaseGraph(cluster, anticodon+"_c"+i);
			if(bg!=null){
				models.add(bg);
				System.err.println("  Cluster "+i+": "+cluster.size()+" seqs, ref len="+bg.original.length);
			}
		}
		System.err.println("Built "+models.size()+" HBM models");
		BaseGraph[] modelArray=models.toArray(new BaseGraph[0]);

		// Also build flat consensus for comparison
		byte[][] consensusLib=new byte[models.size()][];
		for(int i=0; i<models.size(); i++){
			// Traverse the BaseGraph to get consensus
			Read consensus=models.get(i).traverse();
			consensusLib[i]=consensus.bases;
		}

		// Load genome and run TrnaCaller to get candidates
		// Actually, just score all E. coli tRNAs and some random CDS fragments
		System.err.println("\nLoading genome...");
		ArrayList<Read> genome=loadFasta(genomeFile);
		byte[] genomeBases=genome.get(0).bases;
		Tools.toUpperCase(genomeBases);

		// Extract E. coli tRNAs from the genome (known positions from earlier)
		// For now, just score the known GUC tRNAs and some random 76bp fragments
		System.err.println("Scoring "+anticodon+" tRNAs from training set (first 10)...");
		System.err.println("seq_len\tidentity\thbm_score\tlabel");
		int count=0;
		for(Read r : filtered){
			if(count>=10){break;}
			Tools.toUpperCase(r.bases);
			float id=TrnaConsensusBuilder.bestIdentity(r.bases, consensusLib);
			float hbm=TrnaConsensusBuilder.bestModelScore(r.bases, modelArray);
			System.out.println(r.length()+"\t"+String.format("%.3f", id)+"\t"+String.format("%.4f", hbm)+"\tTP_training");
			count++;
		}

		// Score random 76bp fragments from the genome as FP candidates
		System.err.println("Scoring random 76bp genome fragments (FP proxies)...");
		int step=genomeBases.length/20;
		for(int i=0; i<20 && i*step+76<genomeBases.length; i++){
			byte[] frag=new byte[76];
			System.arraycopy(genomeBases, i*step, frag, 0, 76);
			float id=TrnaConsensusBuilder.bestIdentity(frag, consensusLib);
			float hbm=TrnaConsensusBuilder.bestModelScore(frag, modelArray);
			System.out.println("76\t"+String.format("%.3f", id)+"\t"+String.format("%.4f", hbm)+"\tFP_random");
		}

		t.stop();
		System.err.println("\nTime: "+t);
	}

	static ArrayList<Read> loadFasta(String path){
		ArrayList<Read> reads=new ArrayList<>();
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null){reads.add(r);}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		return reads;
	}
}
