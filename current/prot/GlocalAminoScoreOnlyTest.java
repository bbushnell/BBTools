package prot;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

/**
 * Correctness (vs the full-matrix reference, random + real pairs) and timing (vs
 * AAAligner.align, the production aligner, same pairs, >=10 s floor each) for
 * {@link GlocalAminoScoreOnly}. Usage: java -ea prot.GlocalAminoScoreOnlyTest <queries.faa> <truth.tsv> <reps.fasta> [n]
 * @author UMP45
 */
public final class GlocalAminoScoreOnlyTest {

	public static void main(String[] args) throws Exception{
		//1. random pairs vs reference
		final Random rnd=new Random(1);
		final byte[] letters="ARNDCQEGHILKMFPSTWYVX".getBytes();
		final int[] out=new int[2];
		int checked=0;
		for(int t=0; t<2000; t++){
			final int ql=1+rnd.nextInt(400), rl=1+rnd.nextInt(400);
			final byte[] q=new byte[ql], r=new byte[rl];
			for(int i=0; i<ql; i++){q[i]=letters[rnd.nextInt(letters.length)];}
			for(int i=0; i<rl; i++){r[i]=letters[rnd.nextInt(letters.length)];}
			final byte[] qe=Blosum62.encode(q, "q"), re=Blosum62.encode(r, "r");
			GlocalAminoScoreOnly.scoreOnly(qe, re, out);
			final int[] ref=GlocalAminoScoreOnly.scoreOnlyReference(qe, re);
			if(out[0]!=ref[0] || out[1]!=ref[1]){throw new AssertionError("random pair "+t+": got "+out[0]+"/"+out[1]+" ref "+ref[0]+"/"+ref[1]);}
			checked++;
		}
		System.out.println("RANDOM_VS_REFERENCE PASS "+checked);

		if(args.length<3){System.out.println("no real pairs given; done"); return;}
		final int n=(args.length>3 ? Integer.parseInt(args[3]) : 1000);
		final HashMap<String,String> queries=loadFasta(args[0]), reps=loadFasta(args[2]);
		final List<String[]> truth=new ArrayList<String[]>();
		try(BufferedReader br=new BufferedReader(new FileReader(args[1]))){
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] x=line.split("\t"); truth.add(new String[]{x[0], x[1]});
			}
		}
		final int stride=Math.max(1, truth.size()/n);
		final ArrayList<byte[]> qRaw=new ArrayList<byte[]>(), rRaw=new ArrayList<byte[]>();
		for(int i=0; i<truth.size() && qRaw.size()<n; i+=stride){
			final String qs=queries.get(truth.get(i)[0]), rs=reps.get(truth.get(i)[1]);
			if(qs==null || rs==null){continue;}
			qRaw.add(qs.getBytes()); rRaw.add(rs.getBytes());
		}
		final int np=qRaw.size();
		final byte[][] qe=new byte[np][], re=new byte[np][];
		long cells=0;
		for(int i=0; i<np; i++){qe[i]=Blosum62.encode(qRaw.get(i), "q"); re[i]=Blosum62.encode(rRaw.get(i), "r"); cells+=(long)qe[i].length*re[i].length;}
		System.out.println("real pairs "+np+", mean cells/pair "+(cells/(double)np));

		//2. real pairs vs reference (the comparison with Eru's GlocalAminoBlosumFast.scoreOnly, branch eru/glocal-amino, was run
		//   separately: 990/1000 agree, the 10 differ only where a query overhang before the reference start is optimal -- his
		//   column-0 boundary is NEG, this class uses -i*GAP like idaligner.GlocalAligner; recorded 2026-09-03)
		for(int i=0; i<np; i++){
			GlocalAminoScoreOnly.scoreOnly(qe[i], re[i], out);
			final int[] ref=GlocalAminoScoreOnly.scoreOnlyReference(qe[i], re[i]);
			if(out[0]!=ref[0] || out[1]!=ref[1]){throw new AssertionError("real pair "+i+": got "+out[0]+"/"+out[1]+" ref "+ref[0]+"/"+ref[1]);}
		}
		System.out.println("REAL_VS_REFERENCE PASS "+np);

		//3. timing, >=10 s floor per aligner, same pairs, single thread
		final long FLOOR=10_000_000_000L;
		long sink=0;
		//warm
		for(int w=0; w<3; w++){for(int i=0; i<np; i++){GlocalAminoScoreOnly.scoreOnly(qe[i], re[i], out); sink+=out[0];}}
		double nsMine=time(FLOOR, cells, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoScoreOnly.scoreOnly(qe[i], re[i], out); s+=out[0];} return s;});
		//AAAligner.align takes Blosum62-ENCODED residues (asserts on raw letters)
		double nsAAA=time(FLOOR, cells, () -> {long s=0; for(int i=0; i<np; i++){final AAAlignment a=AAAligner.align(qe[i], re[i], false); s+=(a==null ? 0 : a.rawScore);} return s;});
		System.out.println("TIMING ns/cell: GlocalAminoScoreOnly(UMP45)="+fmt(nsMine)+"  AAAligner.align(production)="+fmt(nsAAA)
			+"  speedup vs production "+fmt(nsAAA/nsMine)+"x  (sink "+sink+")");
	}

	interface Work{long run();}

	private static double time(final long floor, final long cellsPerLoop, final Work w){
		long loops=0, sink=0; final long t0=System.nanoTime(); long el;
		do{sink+=w.run(); loops++;}while((el=System.nanoTime()-t0)<floor);
		if(sink==42){System.err.println("sink");}
		return el/(double)(loops*cellsPerLoop);
	}

	private static String fmt(double d){return String.format("%.3f", d);}

	private static HashMap<String,String> loadFasta(String path) throws Exception{
		final HashMap<String,String> m=new HashMap<String,String>();
		try(BufferedReader br=new BufferedReader(new FileReader(path))){
			String id=null; final StringBuilder sb=new StringBuilder();
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				if(line.startsWith(">")){
					if(id!=null){m.put(id, sb.toString());}
					id=line.substring(1).split("\\s")[0]; sb.setLength(0);
				}else{sb.append(line.trim());}
			}
			if(id!=null){m.put(id, sb.toString());}
		}
		return m;
	}
}
