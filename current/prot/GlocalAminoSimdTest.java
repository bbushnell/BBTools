package prot;

import java.util.Random;

import shared.Shared;

/**
 * Exactness + timing gate for {@link GlocalAminoSimd} against {@link GlocalAminoScoreOnly}
 * (the scalar contract) and its full-matrix reference. Every comparison is exact integer
 * equality on BOTH outputs (score and first-maximum reference end).
 *
 * <p>Cases: (1) the existing 2,000 random pairs (same generator as GlocalAminoScoreOnlyTest);
 * (2) lane/segment boundary lengths (q and r in {1,2,15,16,17,31,32,33,255,256,257});
 * (3) 300 homolog pairs: mutated copies with insertions/deletions, so lazy-F must run more
 * than one pass; (4) X-rich pairs; (5) long pairs inside the 16-bit bound and pairs beyond
 * it, asserting the fallback counter moves only for the latter; (6) the same suite with
 * {@code Shared.SIMD=false}, proving the scalar path is reached and identical; (7) timing,
 * >=10 s floor each, ns/cell SIMD vs scalar on 300x300 random pairs and on the homolog set.</p>
 *
 * <p>Run: {@code java -ea --add-modules jdk.incubator.vector -cp <build>:<BBTools/current> prot.GlocalAminoSimdTest [skiptiming]}
 * — without the module, {@code Shared.SIMD} is false and the suite proves the fallback only.</p>
 *
 * @author UMP45
 */
public final class GlocalAminoSimdTest {

	private static final byte[] LETTERS="ARNDCQEGHILKMFPSTWYVX".getBytes();

	public static void main(String[] args) throws Exception{
		final boolean skipTiming=(args.length>0 && "skiptiming".equals(args[0]));
		System.out.println("Shared.SIMD="+Shared.SIMD+" (SIMD kernel "+(Shared.SIMD ? "ACTIVE" : "INACTIVE: scalar fallback only")+")");
		final boolean simdOn=Shared.SIMD;
		runSuite("SIMD="+simdOn);
		if(simdOn){
			Shared.SIMD=false;
			try{runSuite("SIMD=false (forced)");}finally{Shared.SIMD=true;}
		}
		if(!skipTiming){timing(simdOn);}
		//Real pairs, same three files and stride sampling as GlocalAminoScoreOnlyTest:
		//  [skiptiming] <queries.faa> <truth.tsv> <reps.fasta> [n]   (skiptiming: exactness only, no 10 s floors)
		final String[] realArgs=(skipTiming ? java.util.Arrays.copyOfRange(args, 1, args.length) : args);
		if(realArgs.length>=3){realPairs(realArgs, simdOn, !skipTiming);}
		else{System.out.println("REAL_PAIRS 0 (no fixture files given)");}
		System.out.println("GlocalAminoSimdTest PASS");
	}

	/** SIMD-or-fallback vs scalar vs full-matrix reference on n stride-sampled real (query, assigned-rep) pairs; then timing. */
	private static void realPairs(final String[] args, final boolean simdOn, final boolean doTiming) throws Exception{
		final int n=(args.length>3 ? Integer.parseInt(args[3]) : 1000);
		final java.util.HashMap<String,String> queries=loadFasta(args[0]), reps=loadFasta(args[2]);
		final java.util.ArrayList<String[]> truth=new java.util.ArrayList<String[]>();
		try(java.io.BufferedReader br=new java.io.BufferedReader(new java.io.FileReader(args[1]))){
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] x=line.split("\t"); truth.add(new String[]{x[0], x[1]});
			}
		}
		final int stride=Math.max(1, truth.size()/n);
		final java.util.ArrayList<byte[]> qe=new java.util.ArrayList<byte[]>(), re=new java.util.ArrayList<byte[]>();
		for(int i=0; i<truth.size() && qe.size()<n; i+=stride){
			final String qs=queries.get(truth.get(i)[0]), rs=reps.get(truth.get(i)[1]);
			if(qs==null || rs==null){continue;}
			qe.add(Blosum62.encode(qs.getBytes(), "q")); re.add(Blosum62.encode(rs.getBytes(), "r"));
		}
		final int np=qe.size();
		check(np>0, "real pairs: none loaded from "+args[1]);
		final int[] a=new int[2], b=new int[2];
		long cells=0, simdBefore=GlocalAminoSimd.simdCalls, scalarBefore=GlocalAminoSimd.scalarCalls;
		final GlocalAminoSimd.Profile[] prof=new GlocalAminoSimd.Profile[np];
		for(int i=0; i<np; i++){
			prof[i]=GlocalAminoSimd.profile(qe.get(i));
			GlocalAminoSimd.scoreOnly(prof[i], re.get(i), a);
			GlocalAminoScoreOnly.scoreOnly(qe.get(i), re.get(i), b);
			check(a[0]==b[0] && a[1]==b[1], "real pair "+i+": simd "+a[0]+"/"+a[1]+" != scalar "+b[0]+"/"+b[1]+" (q="+qe.get(i).length+" r="+re.get(i).length+")");
			final int[] ref=GlocalAminoScoreOnly.scoreOnlyReference(qe.get(i), re.get(i));
			check(a[0]==ref[0] && a[1]==ref[1], "real pair "+i+": simd != reference");
			cells+=(long)qe.get(i).length*re.get(i).length;
		}
		System.out.println("REAL_PAIRS "+np+" EXACT_VS_SCALAR_AND_REFERENCE PASS (truth rows "+truth.size()+", stride "+stride+", mean cells/pair "
			+(cells/(double)np)+", simdCalls +"+(GlocalAminoSimd.simdCalls-simdBefore)+", scalarCalls +"+(GlocalAminoSimd.scalarCalls-scalarBefore)+")");
		if(!doTiming){return;}
		final int[] out=new int[2];
		final long FLOOR=10_000_000_000L;
		final long fc=cells;
		for(int w=0; w<3; w++){for(int i=0; i<np; i++){GlocalAminoSimd.scoreOnly(prof[i], re.get(i), out); GlocalAminoScoreOnly.scoreOnly(qe.get(i), re.get(i), out);}}
		final double nsSimd=time(FLOOR, fc, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoSimd.scoreOnly(prof[i], re.get(i), out); s+=out[0];} return s;});
		final double nsScalar=time(FLOOR, fc, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoScoreOnly.scoreOnly(qe.get(i), re.get(i), out); s+=out[0];} return s;});
		System.out.println("TIMING ns/cell real pairs("+np+"): GlocalAminoSimd="+fmt(nsSimd)+" GlocalAminoScoreOnly="+fmt(nsScalar)+" speedup="+fmt(nsScalar/nsSimd)+"x (Shared.SIMD="+simdOn+")");
	}

	private static java.util.HashMap<String,String> loadFasta(String path) throws Exception{
		final java.util.HashMap<String,String> m=new java.util.HashMap<String,String>();
		try(java.io.BufferedReader br=new java.io.BufferedReader(new java.io.FileReader(path))){
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

	private static void runSuite(final String label){
		final Random rnd=new Random(1);
		final int[] a=new int[2], b=new int[2];
		long cases=0;
		//(1) random pairs, the existing generator.
		for(int t=0; t<2000; t++){
			final byte[] q=random(rnd, 1+rnd.nextInt(400)), r=random(rnd, 1+rnd.nextInt(400));
			cases+=compare(q, r, a, b, "random "+t);
		}
		//(2) boundary lengths.
		final int[] lens={1,2,15,16,17,31,32,33,255,256,257};
		for(int ql : lens){for(int rl : lens){cases+=compare(random(rnd, ql), random(rnd, rl), a, b, "boundary q="+ql+" r="+rl);}}
		//(3) homologs with indels (high scores, multi-pass lazy-F).
		for(int t=0; t<300; t++){
			final byte[] q=random(rnd, 50+rnd.nextInt(400));
			cases+=compare(q, mutate(rnd, q), a, b, "homolog "+t);
			cases+=compare(mutate(rnd, q), q, a, b, "homolog-rev "+t);
		}
		//(4) X-rich.
		for(int t=0; t<100; t++){
			final byte[] q=random(rnd, 1+rnd.nextInt(60)), r=random(rnd, 1+rnd.nextInt(60));
			for(int i=0; i<q.length; i+=3){q[i]='X';}
			for(int i=1; i<r.length; i+=2){r[i]='X';}
			cases+=compare(q, r, a, b, "xrich "+t);
		}
		//(5) long pairs: inside the bound (SIMD when on) and beyond it (scalar fallback always).
		final long simdBefore=GlocalAminoSimd.simdCalls, scalarBefore=GlocalAminoSimd.scalarCalls;
		for(int t=0; t<20; t++){cases+=compare(random(rnd, 1200+rnd.nextInt(200)), random(rnd, 1200+rnd.nextInt(200)), a, b, "long-in "+t);}
		final long simdMid=GlocalAminoSimd.simdCalls, scalarMid=GlocalAminoSimd.scalarCalls;
		for(int t=0; t<5; t++){cases+=compare(random(rnd, 2400+rnd.nextInt(200)), random(rnd, 2400+rnd.nextInt(200)), a, b, "long-out "+t);}
		final long simdAfter=GlocalAminoSimd.simdCalls, scalarAfter=GlocalAminoSimd.scalarCalls;
		if(Shared.SIMD){
			check(simdMid-simdBefore==20 && scalarMid-scalarBefore==0, label+": in-bound long pairs must take the SIMD path ("+(simdMid-simdBefore)+"/"+(scalarMid-scalarBefore)+")");
			check(simdAfter-simdMid==0 && scalarAfter-scalarMid==5, label+": out-of-bound long pairs must fall back to scalar ("+(simdAfter-simdMid)+"/"+(scalarAfter-scalarMid)+")");
		}else{
			check(simdAfter-simdBefore==0, label+": SIMD path taken while Shared.SIMD=false");
		}
		//(7) ties: repetitive sequences make every last-row cell equal; refEnd must be the LEFTMOST maximum.
		for(int t=0; t<50; t++){
			final byte[] unit=random(rnd, 1+rnd.nextInt(4));
			final byte[] q=repeat(unit, 1+rnd.nextInt(6)), r=repeat(unit, 2+rnd.nextInt(40));
			cases+=compare(q, r, a, b, "tie "+t);
		}
		cases+=compare(new byte[]{'A'}, "AAAAAAAAAA".getBytes(), a, b, "tie A/A10");
		{
			GlocalAminoSimd.scoreOnly(Blosum62.encode(new byte[]{'A'}, "q"), Blosum62.encode("AAAAAAAAAA".getBytes(), "r"), a);
			check(a[1]==0, label+": tie A/A10 refEnd must be leftmost (0), got "+a[1]);
		}
		//(8) overflow bound edges: the largest admitted target takes the SIMD path and matches; one past it falls back;
		//    an all-W identical pair (max positive score 11/cell) and an all-gap-forcing pair (max negative) inside the bound match.
		{
			final byte[] q=Blosum62.encode(random(rnd, 500), "q");
			final GlocalAminoSimd.Profile p=GlocalAminoSimd.profile(q);
			final int edge=p.maxRLen;
			check(edge>500, label+": maxRLen for a 500-aa query is only "+edge);
			final byte[] rEdge=Blosum62.encode(random(rnd, edge), "r"), rOver=Blosum62.encode(random(rnd, edge+1), "r");
			final long s0=GlocalAminoSimd.simdCalls, c0=GlocalAminoSimd.scalarCalls;
			GlocalAminoSimd.scoreOnly(p, rEdge, a); GlocalAminoScoreOnly.scoreOnly(q, rEdge, b);
			check(a[0]==b[0] && a[1]==b[1], label+": bound-edge pair mismatch "+a[0]+"/"+a[1]+" vs "+b[0]+"/"+b[1]);
			GlocalAminoSimd.scoreOnly(p, rOver, a); GlocalAminoScoreOnly.scoreOnly(q, rOver, b);
			check(a[0]==b[0] && a[1]==b[1], label+": bound+1 pair mismatch");
			if(Shared.SIMD){
				check(GlocalAminoSimd.simdCalls-s0==1 && GlocalAminoSimd.scalarCalls-c0==1, label+": bound edge must be SIMD then scalar ("+(GlocalAminoSimd.simdCalls-s0)+"/"+(GlocalAminoSimd.scalarCalls-c0)+")");
			}
			cases+=2;
			final byte[] w=new byte[1500]; java.util.Arrays.fill(w, (byte)'W');
			cases+=compare(w, w, a, b, "max-positive W1500");//score 11*1500=16500, near the positive bound
			final byte[] c=new byte[1400]; java.util.Arrays.fill(c, (byte)'C');
			final byte[] g=new byte[1400]; java.util.Arrays.fill(g, (byte)'W');
			cases+=compare(c, g, a, b, "max-negative C1400/W1400");//C vs W is -2 everywhere; gaps compete
		}
		//(9) caller mutates its query array AFTER profile(): both paths must still score the ORIGINAL residues
		//    (the profile owns a snapshot; the scalar fallback must read the same snapshot, not the caller's array).
		{
			final byte[] orig=Blosum62.encode(random(rnd, 320), "q");
			final byte[] callerArray=orig.clone();
			final GlocalAminoSimd.Profile p=GlocalAminoSimd.profile(callerArray);
			for(int i=0; i<callerArray.length; i++){callerArray[i]=(byte)((callerArray[i]+7)%20);}//scramble every residue
			final byte[] rIn=Blosum62.encode(random(rnd, 300), "r");//inside the bound: SIMD path when Shared.SIMD
			final byte[] rOut=Blosum62.encode(random(rnd, p.maxRLen+1), "r");//beyond the bound: scalar fallback path
			GlocalAminoSimd.scoreOnly(p, rIn, a); GlocalAminoScoreOnly.scoreOnly(orig, rIn, b);
			check(a[0]==b[0] && a[1]==b[1], label+": after caller mutation, in-bound path scored the mutated query ("+a[0]+"/"+a[1]+" vs "+b[0]+"/"+b[1]+")");
			GlocalAminoSimd.scoreOnly(p, rOut, a); GlocalAminoScoreOnly.scoreOnly(orig, rOut, b);
			check(a[0]==b[0] && a[1]==b[1], label+": after caller mutation, fallback path scored the mutated query ("+a[0]+"/"+a[1]+" vs "+b[0]+"/"+b[1]+")");
			cases+=2;
		}
		//(10) shared profile, many threads, shuffled target order (design gate 2): every (target) result must equal the
		//     serial result exactly. Diagnostic counters are NOT asserted here (unsynchronized by design).
		{
			final byte[] q=Blosum62.encode(random(rnd, 280), "q");
			final GlocalAminoSimd.Profile shared=GlocalAminoSimd.profile(q);
			final int nt=64;
			final byte[][] targets=new byte[nt][];
			final int[][] serial=new int[nt][2];
			for(int t=0; t<nt; t++){targets[t]=Blosum62.encode(t%7==0 ? mutate(rnd, random(rnd, 280)) : random(rnd, 1+rnd.nextInt(500)), "r"); GlocalAminoSimd.scoreOnly(shared, targets[t], serial[t]);}
			final int threads=32;//design gate 2: t=1 (serial above) vs t=32
			final java.util.concurrent.ExecutorService pool=java.util.concurrent.Executors.newFixedThreadPool(threads);
			try{
				final java.util.List<java.util.concurrent.Future<int[][]>> fs=new java.util.ArrayList<java.util.concurrent.Future<int[][]>>();
				for(int w=0; w<threads*2; w++){
					final long seed=1000+w;
					fs.add(pool.submit(new java.util.concurrent.Callable<int[][]>(){public int[][] call(){
						final int[] order=new int[nt]; for(int i=0; i<nt; i++){order[i]=i;}
						final Random r=new Random(seed);
						for(int i=nt-1; i>0; i--){final int k=r.nextInt(i+1); final int tmp=order[i]; order[i]=order[k]; order[k]=tmp;}
						final int[][] got=new int[nt][2]; final int[] o=new int[2];
						for(int i=0; i<nt; i++){GlocalAminoSimd.scoreOnly(shared, targets[order[i]], o); got[order[i]][0]=o[0]; got[order[i]][1]=o[1];}
						return got;
					}}));
				}
				for(int w=0; w<fs.size(); w++){
					final int[][] got;
					try{got=fs.get(w).get(120, java.util.concurrent.TimeUnit.SECONDS);}
					catch(Exception e){throw new AssertionError(label+": threaded scoring task "+w+" failed or timed out", e);}
					for(int t=0; t<nt; t++){check(got[t][0]==serial[t][0] && got[t][1]==serial[t][1], label+": thread task "+w+" target "+t+": "+got[t][0]+"/"+got[t][1]+" != serial "+serial[t][0]+"/"+serial[t][1]);}
				}
			}finally{
				pool.shutdownNow();
			}
			cases+=nt;
		}
		//(6) profile reuse across many targets equals per-call profiles.
		final byte[] q=Blosum62.encode(random(rnd, 300), "q");
		final GlocalAminoSimd.Profile p=GlocalAminoSimd.profile(q);
		for(int t=0; t<200; t++){
			final byte[] r=Blosum62.encode(random(rnd, 1+rnd.nextInt(600)), "r");
			GlocalAminoSimd.scoreOnly(p, r, a);
			GlocalAminoScoreOnly.scoreOnly(q, r, b);
			check(a[0]==b[0] && a[1]==b[1], label+": profile reuse "+t+": "+a[0]+"/"+a[1]+" vs "+b[0]+"/"+b[1]);
			cases++;
		}
		System.out.println("EXACT_VS_SCALAR PASS ["+label+"] cases="+cases+" simdCalls="+GlocalAminoSimd.simdCalls+" scalarCalls="+GlocalAminoSimd.scalarCalls);
	}

	/** Compares SIMD-or-fallback vs scalar vs full-matrix reference; returns 1. */
	private static int compare(final byte[] qRaw, final byte[] rRaw, final int[] a, final int[] b, final String what){
		final byte[] q=Blosum62.encode(qRaw, "q"), r=Blosum62.encode(rRaw, "r");
		GlocalAminoSimd.scoreOnly(q, r, a);
		GlocalAminoScoreOnly.scoreOnly(q, r, b);
		check(a[0]==b[0] && a[1]==b[1], what+": simd "+a[0]+"/"+a[1]+" != scalar "+b[0]+"/"+b[1]+" (q="+q.length+" r="+r.length+")");
		if(q.length*r.length<=250000){//the O(n^2) reference is only cheap enough for the smaller pairs
			final int[] ref=GlocalAminoScoreOnly.scoreOnlyReference(q, r);
			check(a[0]==ref[0] && a[1]==ref[1], what+": simd "+a[0]+"/"+a[1]+" != reference "+ref[0]+"/"+ref[1]);
		}
		return 1;
	}

	private static byte[] repeat(final byte[] unit, final int times){
		final byte[] s=new byte[unit.length*times];
		for(int i=0; i<s.length; i++){s[i]=unit[i%unit.length];}
		return s;
	}

	private static byte[] random(final Random rnd, final int len){
		final byte[] s=new byte[len];
		for(int i=0; i<len; i++){s[i]=LETTERS[rnd.nextInt(LETTERS.length)];}
		return s;
	}

	/** ~10% substitutions, a few insertions and deletions, random flanks: a plausible homolog. */
	private static byte[] mutate(final Random rnd, final byte[] q){
		final java.io.ByteArrayOutputStream o=new java.io.ByteArrayOutputStream();
		final int lead=rnd.nextInt(20);
		for(int i=0; i<lead; i++){o.write(LETTERS[rnd.nextInt(LETTERS.length)]);}
		for(int i=0; i<q.length; i++){
			final double u=rnd.nextDouble();
			if(u<0.03){continue;}//deletion
			if(u<0.06){o.write(LETTERS[rnd.nextInt(LETTERS.length)]);}//insertion before
			o.write(rnd.nextDouble()<0.10 ? LETTERS[rnd.nextInt(LETTERS.length)] : q[i]);
		}
		final int tail=rnd.nextInt(20);
		for(int i=0; i<tail; i++){o.write(LETTERS[rnd.nextInt(LETTERS.length)]);}
		final byte[] m=o.toByteArray();
		return m.length==0 ? new byte[]{'A'} : m;
	}

	private static void timing(final boolean simdOn){
		final Random rnd=new Random(7);
		final int np=200;
		final byte[][] q=new byte[np][], r=new byte[np][], hq=new byte[np][], hr=new byte[np][];
		long cells=0, hcells=0;
		for(int i=0; i<np; i++){
			q[i]=Blosum62.encode(random(rnd, 300), "q"); r[i]=Blosum62.encode(random(rnd, 300), "r"); cells+=(long)q[i].length*r[i].length;
			final byte[] base=random(rnd, 100+rnd.nextInt(400));
			hq[i]=Blosum62.encode(base, "q"); hr[i]=Blosum62.encode(mutate(rnd, base), "r"); hcells+=(long)hq[i].length*hr[i].length;
		}
		final GlocalAminoSimd.Profile[] pq=new GlocalAminoSimd.Profile[np], phq=new GlocalAminoSimd.Profile[np];
		for(int i=0; i<np; i++){pq[i]=GlocalAminoSimd.profile(q[i]); phq[i]=GlocalAminoSimd.profile(hq[i]);}
		final int[] out=new int[2];
		final long FLOOR=10_000_000_000L;
		for(int w=0; w<5; w++){for(int i=0; i<np; i++){GlocalAminoSimd.scoreOnly(pq[i], r[i], out); GlocalAminoScoreOnly.scoreOnly(q[i], r[i], out);}}
		final double nsSimd=time(FLOOR, cells, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoSimd.scoreOnly(pq[i], r[i], out); s+=out[0];} return s;});
		final double nsScalar=time(FLOOR, cells, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoScoreOnly.scoreOnly(q[i], r[i], out); s+=out[0];} return s;});
		final double nsSimdH=time(FLOOR, hcells, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoSimd.scoreOnly(phq[i], hr[i], out); s+=out[0];} return s;});
		final double nsScalarH=time(FLOOR, hcells, () -> {long s=0; for(int i=0; i<np; i++){GlocalAminoScoreOnly.scoreOnly(hq[i], hr[i], out); s+=out[0];} return s;});
		System.out.println("TIMING ns/cell random300x300: GlocalAminoSimd="+fmt(nsSimd)+" GlocalAminoScoreOnly="+fmt(nsScalar)+" speedup="+fmt(nsScalar/nsSimd)+"x (Shared.SIMD="+simdOn+")");
		System.out.println("TIMING ns/cell homologs: GlocalAminoSimd="+fmt(nsSimdH)+" GlocalAminoScoreOnly="+fmt(nsScalarH)+" speedup="+fmt(nsScalarH/nsSimdH)+"x (Shared.SIMD="+simdOn+")");
		if(simdOn){check(nsScalar/nsSimd>=2.0, "SIMD kernel active but under 2x the scalar scorer: "+fmt(nsScalar/nsSimd)+"x (gate threshold is 5x in the receipt; 2x is the tripwire)");}
	}

	interface Work{long run();}
	private static double time(final long floor, final long cellsPerLoop, final Work w){
		long loops=0, sink=0; final long t0=System.nanoTime(); long el;
		do{sink+=w.run(); loops++;}while((el=System.nanoTime()-t0)<floor);
		if(sink==42){System.err.println("sink");}
		return el/(double)(loops*cellsPerLoop);
	}
	private static String fmt(double d){return String.format("%.3f", d);}
	private static void check(boolean ok, String msg){if(!ok){throw new AssertionError(msg);}}
}
