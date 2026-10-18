package test.prot;

import prot.Blosum62;
import prot.HbmPositionModel;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicReference;

/** Independent path enumeration and boundary fixtures for stop-facing profile tails. */
public final class HbmTailTest{
	public static void main(String[] args){
		try{test();}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}
	private static void test(){
		final Random random=new Random(20261017);
		for(int trial=0; trial<600; trial++){
			final int n=1+random.nextInt(4), m=1+random.nextInt(4), cap=1+random.nextInt(4);
			final int[][] scores=new int[22][n];
			for(int a=0; a<22; a++){for(int j=0; j<n; j++){scores[a][j]=(random.nextInt(17)-6)*64;}}
			final byte[] query=new byte[m];
			for(int i=0; i<m; i++){query[i]=(byte)random.nextInt(20);}
			int best=0;
			for(int q=Math.max(0, m-cap); q<m; q++){
				for(int r=Math.max(0, n-cap); r<n; r++){
					best=Math.max(best, enumerate(query, scores, q, r, 0));
				}
			}
			final HbmPositionModel model=HbmPositionModel.fromScores(scores);
			final HbmPositionModel.TailResult hit=model.alignTail(query, cap, true);
			final HbmPositionModel.TailResult scoreOnly=model.alignTail(query, cap, false);
			equal(legacyTail(scores, query, cap, true), hit, "enumerated recorded trial="+trial);
			equal(legacyTail(scores, query, cap, false), scoreOnly, "enumerated unrecorded trial="+trial);
			check((hit==null)==(best==0), "Positive suffix existence differs from exhaustive paths");
			check((hit==null)==(scoreOnly==null), "Trace recording changed suffix existence");
			if(hit!=null){
				check(hit.score==best && scoreOnly.score==best, "Tail score differs from exhaustive optimum");
				check(hit.queryStart==scoreOnly.queryStart && hit.start==scoreOnly.start, "Trace recording changed clipping");
				check(hit.queryStart>=Math.max(0, m-cap) && hit.start>=Math.max(0, n-cap), "Alignment escaped capped window");
				check(hit.queryEnd==m-1 && hit.end==n-1, "Tail alignment clipped a stop-facing end");
				int q=hit.queryStart, r=hit.start, sum=0;
				for(byte op : hit.path){
					if(op=='m'){sum+=scores[query[q++]][r++];}
					else if(op=='I'){q++; sum-=HbmPositionModel.GAP;}
					else if(op=='D'){r++; sum-=HbmPositionModel.GAP;}
					else{throw new AssertionError("Unknown tail path operation: "+op);}
				}
				check(q==m && r==n && sum==best, "Recorded path must consume both ends and reproduce its score");
			}
		}
		final int[][] scores=new int[22][7];
		for(int[] row : scores){Arrays.fill(row, -6*64);}
		for(int j=4; j<7; j++){scores[0][j]=5*64;}
		final HbmPositionModel model=HbmPositionModel.fromScores(scores);
		final byte[] query={1,1,1,1,1,0,0,0};
		final HbmPositionModel.TailResult hit=model.alignTail(query, 500, true);
		check(hit.score==15*64 && hit.queryStart==5 && hit.start==4, "Both unrelated left prefixes must clip freely");
		check(model.alignTail(query, 2, true).score==10*64, "Cap must bound both query and model windows");
		final byte[] wrongTail={0,0,0,1,1,1,1,1,1};
		check(model.alignTail(wrongTail, 500, true)==null, "An internal matching segment must not receive free right clipping");
		boolean rejected=false;
		try{model.alignTail(new byte[]{20}, 500, true);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Encoded stop cannot enter amino-acid alignment");
		final int tiny=tieTables();
		literalTies(); rectangularAndReset(); concurrent(); errors();
		System.err.println("HBM_TAIL_TEST_PASS exhaustive=600 free_left=true right_consumed=true capped=true path_scores=true"
			+" exact_tiny_trials="+tiny+" exact_random_trials=300 ties=true reset=true detached=true threads=4 errors=true");
	}
	/** Exhausts all two-residue tables of width1..2 over {-G,0,G,2G}, queries1..3 and caps1..3. */
	private static int tieTables(){
		final int g=HbmPositionModel.GAP; int trials=0;
		for(int n=1; n<=2; n++){
			for(int code=0; code<(1<<(4*n)); code++){
				final int[][] table=new int[22][n]; int x=code;
				for(int a=0; a<2; a++){for(int j=0; j<n; j++){table[a][j]=((x&3)-1)*g; x>>>=2;}}
				final HbmPositionModel model=HbmPositionModel.fromScores(table);
				for(int m=1; m<=3; m++){
					for(int bits=0; bits<(1<<m); bits++){
						final byte[] query=new byte[m];
						for(int i=0; i<m; i++){query[i]=(byte)((bits>>>i)&1);}
						for(int cap=1; cap<=3; cap++){compare(model, table, query, cap); trials++;}
					}
				}
			}
		}
		check(trials==11424, "Tiny-table enumeration must cover every declared table/query/cap combination");
		return trials;
	}
	private static void literalTies(){
		final int g=HbmPositionModel.GAP;
		int[][] table=new int[22][1]; table[0][0]=2*g; table[1][0]=g;
		literal(table, new byte[]{0, 1}, g, 1, 0, "m", "Diagonal must win an up tie");
		table=new int[22][2]; table[0][0]=2*g; table[0][1]=g;
		literal(table, new byte[]{0}, g, 0, 1, "m", "Diagonal must win a left tie");
		table=new int[22][2]; table[0][0]=3*g; table[1][1]=-3*g;
		literal(table, new byte[]{0, 1}, g, 0, 0, "mDI", "Up must win a left tie when diagonal loses");
		table[1][1]=-2*g;
		literal(table, new byte[]{0, 1}, g, 0, 0, "mm", "Diagonal must win a three-way positive tie");
		table=new int[22][2]; table[1][1]=g;
		literal(table, new byte[]{0, 1}, g, 1, 1, "m", "A zero-scoring diagonal must not extend traceback past restart");
		table=new int[22][2];
		check(HbmPositionModel.fromScores(table).alignTail(new byte[]{0, 1}, 2, true)==null, "All-zero endpoint must remain no hit");
		for(int[] row : table){Arrays.fill(row, -g);}
		check(HbmPositionModel.fromScores(table).alignTail(new byte[]{0, 1}, 2, false)==null, "All-negative endpoint must remain no hit");
	}
	private static void literal(int[][] table, byte[] query, int score, int qs, int rs, String path, String why){
		final HbmPositionModel model=HbmPositionModel.fromScores(table);
		final HbmPositionModel.TailResult hit=model.alignTail(query, 500, true);
		check(hit!=null && hit.score==score && hit.queryStart==qs && hit.start==rs
			&& Arrays.equals(hit.path, path.getBytes(java.nio.charset.StandardCharsets.US_ASCII)), why);
		compare(model, table, query, 500);
	}
	private static void rectangularAndReset(){
		final Random random=new Random(2026101809);
		for(int trial=0; trial<300; trial++){
			final int[][] table=randomTable(random, 1+random.nextInt(29));
			compare(HbmPositionModel.fromScores(table), table, query(random, 1+random.nextInt(23)), 1+random.nextInt(31));
		}
		//Vary both stride and height on the same thread, with independent caps/offsets and X residues.
		for(int[] shape : new int[][]{{1, 31}, {33, 1}, {3, 17}, {19, 2}, {499, 501}, {501, 499}, {500, 500}, {513, 701}, {2, 3}}){
			final int[][] table=randomTable(random, shape[1]); final byte[] query=query(random, shape[0]);
			final HbmPositionModel model=HbmPositionModel.fromScores(table);
			for(int cap : new int[]{1, 2, 30, 499, 500, 501}){compare(model, table, query, cap);}
		}
		final int[][] table=new int[22][9]; Arrays.fill(table[0], HbmPositionModel.GAP);
		final HbmPositionModel model=HbmPositionModel.fromScores(table);
		final byte[] query=new byte[9];
		final HbmPositionModel.TailResult saved=model.alignTail(query, 9, true);
		final byte[] path=saved.path.clone();
		final HbmPositionModel.Result ordinary=model.align(query, true);
		compare(model, table, new byte[]{1}, 1); compare(model, table, new byte[15], 7);
		compare(model, table, query, Integer.MAX_VALUE);
		final HbmPositionModel.Result after=model.align(query, true);
		check(Arrays.equals(saved.path, path), "Reusing tail scratch must not mutate a detached earlier path");
		equal(legacyTail(table, query, 9, true), saved, "Saved hit after scratch reuse");
		check(ordinary.score==after.score && ordinary.start==after.start && ordinary.end==after.end
			&& Arrays.equals(ordinary.path, after.path), "Tail scratch must not change ordinary alignment");
	}
	private static void concurrent(){
		final int[][] table=randomTable(new Random(17), 23);
		final HbmPositionModel model=HbmPositionModel.fromScores(table);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread[] threads=new Thread[4];
		for(int t=0; t<threads.length; t++){
			final int seed=t;
			threads[t]=new Thread(()->{
				try{
					final Random random=new Random(seed);
					for(int i=0; i<32; i++){compare(model, table, query(random, 1+random.nextInt(37)), 1+random.nextInt(29));}
				}catch(Throwable problem){failure.compareAndSet(null, problem);}
			});
			threads[t].start();
		}
		for(Thread thread : threads){
			try{thread.join();}catch(InterruptedException interrupted){Thread.currentThread().interrupt(); throw new AssertionError(interrupted);}
		}
		if(failure.get()!=null){throw new AssertionError("Concurrent models must use independent thread-local tail scratch", failure.get());}
	}
	private static void errors(){
		final int[][] table=new int[22][2]; final HbmPositionModel model=HbmPositionModel.fromScores(table);
		for(byte[] query : new byte[][]{null, new byte[0], {20}, {-1}, {22}}){
			for(boolean record : new boolean[]{false, true}){
				invalid(()->model.alignTail(query, 2, record)); invalid(()->legacyTail(table, query, 2, record));
			}
		}
		for(int cap : new int[]{0, -1}){invalid(()->model.alignTail(new byte[1], cap, true));}
		//Only the capped suffix was validated by b16; an invalid clipped prefix stays ignored.
		compare(model, table, new byte[]{20, 0}, 1);
		invalid(()->model.alignTail(new byte[170000], Integer.MAX_VALUE, false));
		final HbmPositionModel wide=HbmPositionModel.fromScores(new int[22][46340]);
		invalid(()->wide.alignTail(new byte[46340], Integer.MAX_VALUE, true));
	}
	private static void invalid(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid tail dimensions/residues must retain the existing failure");
	}
	private static int[][] randomTable(Random random, int length){
		final int[][] table=new int[22][length];
		for(int[] row : table){for(int j=0; j<length; j++){row[j]=random.nextInt(12801)-6400;}}
		return table;
	}
	private static byte[] query(Random random, int length){
		final byte[] query=new byte[length];
		for(int i=0; i<length; i++){final int a=random.nextInt(21); query[i]=(byte)(a==20 ? Blosum62.X_CODE : a);}
		return query;
	}
	private static void compare(HbmPositionModel model, int[][] table, byte[] query, int cap){
		for(boolean record : new boolean[]{false, true}){
			equal(legacyTail(table, query, cap, record), model.alignTail(query, cap, record),
				"m="+query.length+" n="+table[0].length+" cap="+cap+" record="+record);
		}
	}
	private static void equal(HbmPositionModel.TailResult expected, HbmPositionModel.TailResult actual, String why){
		check((expected==null)==(actual==null), "Tail null parity: "+why);
		if(expected==null){return;}
		check(expected.score==actual.score && expected.queryStart==actual.queryStart && expected.queryEnd==actual.queryEnd
			&& expected.start==actual.start && expected.end==actual.end && Arrays.equals(expected.path, actual.path),
			"Tail exact score/endpoints/path parity: "+why);
	}

	/** Frozen b16ca995 full-matrix tail routine. Test-only independent score/trace oracle.
	 * Only model fields and its private scratch are supplied here by the fixture. */
	static HbmPositionModel.TailResult legacyTail(int[][] scores, byte[] query, int maxLength, boolean record){
		if(query==null || query.length<1 || maxLength<1){
			throw new IllegalArgumentException("Tail alignment requires a nonempty query and positive window limit");
		}
		final int length=scores[0].length;
		final int m=Math.min(query.length, maxLength), n=Math.min(length, maxLength);
		final int qOffset=query.length-m, rOffset=length-n, stride=n+1;
		if((long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Oversized tail alignment: "+m+" x "+n);
		}
		final LegacyScratch scratch=LEGACY_SCRATCH.get(); scratch.ensure((m+1)*stride);
		final int[] h=scratch.matrix;
		Arrays.fill(h, 0, stride, 0);
		for(int i=1; i<=m; i++){
			final byte a=query[qOffset+i-1];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded tail residue: "+a);}
			final int[] row=scores[a];
			final int cur=i*stride, prev=cur-stride;
			int diag=h[prev], left=h[cur]=0;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=Math.max(0, diag+row[rOffset+j-1]);
				value=Math.max(value, up-GAP);
				value=Math.max(value, left-GAP);
				h[cur+j]=value; left=value; diag=up;
			}
		}
		final int score=h[m*stride+n];
		if(score<=0){return null;}
		int i=m, j=n;
		final byte[] ops=record ? new byte[m+n] : null;
		int at=record ? ops.length : 0;
		while(i>0 && j>0 && h[i*stride+j]>0){
			if(h[i*stride+j]==h[(i-1)*stride+j-1]+scores[query[qOffset+i-1]][rOffset+j-1]){
				if(record){ops[--at]='m';} i--; j--;
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-GAP){
				if(record){ops[--at]='I';} i--;
			}else{
				if(h[i*stride+j]!=h[i*stride+j-1]-GAP){throw new AssertionError("Tail traceback has no valid predecessor");}
				if(record){ops[--at]='D';} j--;
			}
		}
		assert(i<m && j<n) : "Positive tail evidence must consume residues from both sequences";
		return new HbmPositionModel.TailResult(score, qOffset+i, query.length-1, rOffset+j, length-1,
			record ? Arrays.copyOfRange(ops, at, ops.length) : null);
	}
	private static final class LegacyScratch{
		int[] matrix=new int[0];
		void ensure(int n){if(matrix.length<n){matrix=new int[n];}}
	}
	private static final ThreadLocal<LegacyScratch> LEGACY_SCRATCH=ThreadLocal.withInitial(LegacyScratch::new);
	private static final int GAP=HbmPositionModel.GAP;
	/** Enumerates every path from a chosen start through both terminal coordinates. */
	private static int enumerate(byte[] query, int[][] scores, int q, int r, int score){
		if(q==query.length && r==scores[0].length){return score;}
		int best=Integer.MIN_VALUE;
		if(q<query.length){best=Math.max(best, enumerate(query, scores, q+1, r, score-HbmPositionModel.GAP));}
		if(r<scores[0].length){best=Math.max(best, enumerate(query, scores, q, r+1, score-HbmPositionModel.GAP));}
		if(q<query.length && r<scores[0].length){
			best=Math.max(best, enumerate(query, scores, q+1, r+1, score+scores[query[q]][r]));
		}
		return best;
	}
	private static void check(boolean value, String reason){if(!value){throw new AssertionError(reason);}}
}
