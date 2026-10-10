package prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

/** Exact gate fixtures and differential selection against recording every candidate. @author Keqing */
public final class HbmCompetitiveGateTest {
	public static void main(String[] args){
		final HbmCompetitiveGate.Metrics metric=new HbmCompetitiveGate.Metrics();
		HbmCompetitiveGate.metrics(enc("ACDE"), enc("ADFE"), result(0, 3, "mImDm"), 1, 2, metric);
		check(metric.identity==60f && metric.paired==1 && metric.coverageQ==0.25f && metric.coverageT==0.5f && metric.coverage==0.25f,
			"Gap columns inflated identity or paired core coverage");
		HbmCompetitiveGate.metrics(enc("ACDX"), enc("ACDX"), result(0, 3, "mmmm"), 0, 3, metric);
		check(metric.identity==75f && metric.coverage==1f, "X identities must not count");
		metric.identity=HbmCompetitiveGate.MIN_IDENTITY; metric.coverage=HbmCompetitiveGate.MIN_COVERAGE;
		check(HbmCompetitiveGate.passes(metric), "Inclusive float32 threshold was rejected");
		metric.identity=Math.nextDown(metric.identity);
		check(!HbmCompetitiveGate.passes(metric), "Identity below float32 threshold passed");
		metric.identity=100; metric.coverage=Math.nextDown(HbmCompetitiveGate.MIN_COVERAGE);
		check(!HbmCompetitiveGate.passes(metric), "Coverage below float32 threshold passed");
		check(HbmCompetitiveGate.lengthCanPass(4, 5) && HbmCompetitiveGate.lengthCanPass(5, 4) && !HbmCompetitiveGate.lengthCanPass(3, 5), "Coverage upper bound differs");
		reject(()->HbmCompetitiveGate.metrics(enc("AA"), enc("AA"), result(0, 0, "m"), 0, 1, metric), "consume");
		reject(()->HbmCompetitiveGate.metrics(enc("AA"), enc("AA"), result(0, 1, "mZ"), 0, 1, metric), "Unknown");

		final double[] background=new double[20]; Arrays.fill(background, 0.05);
		final byte[][] refs={enc("CCCCC"), enc("AAAAC"), enc("AAAAC")};
		final AAGraph misleading=new AAGraph(refs[0], 0);
		for(int i=0; i<100; i++){misleading.addTrace(enc("AAAAA"), 0, ops("mmmmm"));}
		final HbmPositionModel[] profiles={profile(misleading, background), profile(new AAGraph(refs[1], 0), background), profile(new AAGraph(refs[2], 0), background)};
		final String[] ids={"misleading", "z", "a"}; final int[] first={0, 0, 0}, last={4, 4, 4};
		final HbmCompetitiveGate gate=new HbmCompetitiveGate(ids, refs, profiles, first, last);
		final HbmCompetitiveGate.Scratch work=new HbmCompetitiveGate.Scratch(3, 3);
		check(profiles[0].align(enc("AAAAA"), false).score>profiles[1].align(enc("AAAAA"), false).score, "Fixture must have an ineligible highest score");
		check(gate.select(enc("AAAAA"), new int[]{1, 0, 2}, 3, work)==2 && work.recorded==2 && work.metrics.identity==80f,
			"Selection failed to fall through the ineligible winner or break ties by family ID");
		check(gate.select(enc("AAAAA"), new int[]{2, 1, 0}, 3, work)==2, "Shortlist ordering changed the winner");
		check(gate.select(enc("WWWWWWWWWW"), new int[]{0, 1, 2}, 3, work)==-1 && work.scored==0 && work.lengthRejected==3 && Float.isNaN(work.metrics.identity),
			"Impossible lengths were aligned or stale winner evidence survived rejection");
		reject(()->gate.select(enc("AAAAA"), new int[]{0, 0}, 2, work), "Duplicate");
		reject(()->gate.select(enc("AAAAA"), new int[]{3}, 1, work), "invalid");

		final Random random=new Random(20261010);
		int compared=0, accepted=0;
		for(int trial=0; trial<600; trial++){
			final int n=6; final byte[][] r=new byte[n][]; final HbmPositionModel[] p=new HbmPositionModel[n];
			final String[] names={"f5", "f4", "f3", "f2", "f1", "f0"}; final int[] starts=new int[n], ends=new int[n], order=new int[n];
			for(int i=0; i<n; i++){
				r[i]=sequence(random, 3+random.nextInt(13));
				final AAGraph g=new AAGraph(r[i], 0);
				for(int j=0; j<3; j++){
					final byte[] member=r[i].clone(); member[random.nextInt(member.length)]=(byte)random.nextInt(20);
					g.addTrace(member, 0, filled(member.length));
				}
				p[i]=profile(g, background); starts[i]=random.nextInt(2); ends[i]=r[i].length-1-random.nextInt(2); order[i]=i;
			}
			final byte[] query=trial%2==0 ? r[trial%n].clone() : sequence(random, 3+random.nextInt(16));
			if(trial%3==0){query[random.nextInt(query.length)]=Blosum62.X_CODE;}
			final HbmCompetitiveGate current=new HbmCompetitiveGate(names, r, p, starts, ends);
			final HbmCompetitiveGate.Scratch scratch=new HbmCompetitiveGate.Scratch(n, n);
			final int expected=exhaustive(query, r, p, names, starts, ends);
			final int actual=current.select(query, order, n, scratch);
			check(expected==actual, "Deferred-trace winner differs from exhaustive production-metric oracle at trial "+trial);
			if(actual>=0){accepted++;}
			for(int i=n-1; i>0; i--){final int j=random.nextInt(i+1), tmp=order[i]; order[i]=order[j]; order[j]=tmp;}
			check(current.select(query, order, n, scratch)==expected, "Permuting candidates changed selection"); compared++;
		}
		check(accepted>0 && accepted<compared, "Differential fixture must exercise both accepted and rejected queries");
		System.err.println("HBM_COMPETITIVE_GATE_TEST_PASS trials="+compared+" accepted="+accepted+" paired_core=true float32=true fallback=true ties=true permutation=true");
	}

	/** Uses the existing construction coverage routine and AAAlignment identity, recording ALL candidates. */
	private static int exhaustive(byte[] q, byte[][] refs, HbmPositionModel[] models, String[] ids, int[] first, int[] last){
		int winner=-1, score=Integer.MIN_VALUE;
		for(int i=0; i<models.length; i++){
			final HbmPositionModel.Result r=models[i].align(q, true);
			int qi=0, ri=r.start, identities=0;
			for(byte op : r.path){
				if(op=='m'){if(q[qi]==refs[i][ri] && q[qi]<20){identities++;} qi++; ri++;}
				else if(op=='I'){qi++;}else if(op=='D'){ri++;}else{throw new AssertionError("Oracle received an unknown operation");}
			}
			final AAAlignment a=new AAAlignment(r.score, 0, q.length-1, r.start, r.end, identities, 0, 0, r.path.length, r.path);
			final ProteinSearcher.D55Metrics m=new ProteinSearcher.D55Metrics();
			ProteinSearcher.constructionCoreMetrics(a, q.length, first[i], last[i], last[i]-first[i]+1, m);
			if((float)a.pident()<40.612846f || (float)m.mutualOverlap<0.8f){continue;}
			if(winner<0 || r.score>score || (r.score==score && ids[i].compareTo(ids[winner])<0)){winner=i; score=r.score;}
		}
		return winner;
	}
	private static HbmPositionModel profile(AAGraph graph, double[] bg){return HbmPositionModel.deriveOne(graph, "logodds", .01, true, bg, -4);}
	private static byte[] sequence(Random r, int n){final byte[] out=new byte[n]; for(int i=0; i<n; i++){out[i]=(byte)r.nextInt(20);} return out;}
	private static byte[] filled(int n){final byte[] out=new byte[n]; Arrays.fill(out, (byte)'m'); return out;}
	private static byte[] enc(String s){return Blosum62.encode(s.getBytes(StandardCharsets.US_ASCII), s);}
	private static byte[] ops(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static HbmPositionModel.Result result(int start, int end, String path){return new HbmPositionModel.Result(0, start, end, ops(path));}
	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
	private static void reject(Runnable action, String text){
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(text), "Unexpected rejection: "+e); return;}
		throw new AssertionError("Malformed competitive input was accepted: "+text);
	}
}
