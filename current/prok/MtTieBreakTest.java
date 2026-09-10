package prok;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import shared.Shared;

/**
 * Adversarial tie-break fixture for TrnaConsensusBuilder.bestCentroidMatch's work-aware parallel gate.
 *
 * Forces GUARANTEED exact ties: byte-IDENTICAL centroids placed at indices that fall in DIFFERENT
 * parallel chunks, plus a query equal to them (identity to each is equal by construction). The serial
 * scan (t=1) and the parallel ascending-index strict-'>' merge (t=8) must BOTH return the LOWEST index
 * among the tied centroids. The centroid list is long enough that n*len >= PARALLEL_CENTROID_WORK, so
 * t=8 exercises the parallel branch while t=1 stays serial.
 *
 * Greedy clustering never produces duplicate centroids (it would merge them), so this exact-tie case is
 * only reachable by calling bestCentroidMatch directly — hence its package-private visibility.
 *
 * Run (local, ~instant): java -cp current --add-modules jdk.incubator.vector prok.MtTieBreakTest
 */
public class MtTieBreakTest{

	public static void main(String[] args) throws Exception{
		File in=File.createTempFile("mttie_in", ".fa");
		File out=File.createTempFile("mttie_out", ".fa");
		in.deleteOnExit(); out.deleteOnExit();
		try(PrintWriter pw=new PrintWriter(in)){pw.println(">a"); pw.println("ACGTACGTACGTACGTACGT");}

		final int N=40, L=3000;               // 40*3000 = 120000 >= PARALLEL_CENTROID_WORK(100000)
		final byte[] bases={'A','C','G','T'};
		final Random rng=new Random(12345);   // fixed seed -> fully deterministic fixture

		// tie1 duplicated at indices 0 and 39 (chunk 0 vs last chunk); tie2 at 5 and 33 (interior, cross-chunk).
		final byte[] tie1=randSeq(L, bases, rng);
		final byte[] tie2=randSeq(L, bases, rng);
		final ArrayList<byte[]> centroids=new ArrayList<>(N);
		for(int i=0; i<N; i++){
			if(i==0 || i==N-1){centroids.add(tie1.clone());}
			else if(i==5 || i==33){centroids.add(tie2.clone());}
			else{centroids.add(randSeq(L, bases, rng));}
		}

		int fails=0;
		fails+=runCase("tie1@{0,39}", tie1.clone(), 0, centroids, in, out);
		fails+=runCase("tie2@{5,33}", tie2.clone(), 5, centroids, in, out);

		System.out.println(fails==0 ? "TIEBREAK_ALL_PASS" : ("TIEBREAK_FAIL count="+fails));
		System.exit(fails==0 ? 0 : 1);
	}

	private static byte[] randSeq(int L, byte[] bases, Random rng){
		byte[] s=new byte[L];
		for(int j=0; j<L; j++){s[j]=bases[rng.nextInt(4)];}
		return s;
	}

	/** Runs the query at t=1 (serial) and t=8 (parallel); both must return expectIdx. */
	private static int runCase(String name, byte[] query, int expectIdx,
			ArrayList<byte[]> centroids, File in, File out){
		int fails=0;
		float id1=-1, id8=-1; int idx1=-1, idx8=-1;
		for(int t : new int[]{1, 8}){
			Shared.setThreads(t);
			TrnaConsensusBuilder b=new TrnaConsensusBuilder(
				new String[]{"in="+in.getAbsolutePath(), "out="+out.getAbsolutePath(), "ow=t"});
			float[] o=new float[2];
			b.bestCentroidMatch(query, centroids, o);
			if(t==1){id1=o[0]; idx1=(int)o[1];}else{id8=o[0]; idx8=(int)o[1];}
			boolean ok=((int)o[1]==expectIdx);
			System.out.println(name+" threads="+t+" -> bestId="+o[0]+" bestIdx="+(int)o[1]
				+"  (expect idx "+expectIdx+")  "+(ok?"PASS":"FAIL"));
			if(!ok){fails++;}
		}
		boolean same=(idx1==idx8 && id1==id8);
		System.out.println(name+" serial==parallel: idx "+idx1+"=="+idx8+" id "+id1+"=="+id8+"  "+(same?"IDENTICAL":"DIFFER"));
		if(!same){fails++;}
		return fails;
	}
}
