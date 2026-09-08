package prok;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Method;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import json.JsonObject;
import shared.Shared;
import stream.Read;

/** Focused checks for opt-in ordered rRNA consensus fallback. */
public class CallGenesRrnaFallbackConfigTest {

	public static void main(String[] args) throws Exception {
		testDefaultOff();
		testParserRoundTripAndIsolation();
		testHandMutatedUniversalFailFallbackPass();
		testHandMutatedUniversalPassStops();
		testAllFailPenaltyOnce();
		testR18SFallbackUsesFileOrder();
		testAttemptCountHistogramGrowth();
		testStatsOutputFormatting();
		System.out.println("PASS CallGenesRrnaFallbackConfigTest");
	}

	private static void testDefaultOff(){
		if(ProkObject.rrnaFallback){throw new AssertionError("rrnafallback is unexpectedly enabled by default");}
	}

	private static void testParserRoundTripAndIsolation(){
		final ConfigSnapshot before=new ConfigSnapshot();
		try{
			check(ProkObject.parse("rrnafallback=t", "rrnafallback", "t"), "rrnafallback=t was not parsed");
			check(ProkObject.rrnaFallback, "rrnafallback=t did not enable fallback");
			before.assertUnchangedExceptFallback();
			check(ProkObject.parse("rrnafallback=f", "rrnafallback", "f"), "rrnafallback=f was not parsed");
			check(!ProkObject.rrnaFallback, "rrnafallback=f did not disable fallback");
			check(!ProkObject.parse("notarealflag=t", "notarealflag", "t"), "Unknown flag was accepted");
			before.assertUnchangedExceptFallback();
		}finally{before.restore();}
	}

	private static void testHandMutatedUniversalFailFallbackPass(){
		for(boolean useIDA : new boolean[]{true, false}){
			runFixtureCase("rescue", useIDA, false, QUERY, RESCUE_ENTRY0, RESCUE_ENTRY1, false, 1, false);
			runFixtureCase("rescue", useIDA, true, QUERY, RESCUE_ENTRY0, RESCUE_ENTRY1, true, 2, true);
		}
	}

	private static void testHandMutatedUniversalPassStops(){
		for(boolean useIDA : new boolean[]{true, false}){
			runFixtureCase("pass", useIDA, false, QUERY, PASS_ENTRY0, PASS_ENTRY1_SENTINEL, true, 1, true);
			runFixtureCase("pass", useIDA, true, QUERY, PASS_ENTRY0, PASS_ENTRY1_SENTINEL, true, 1, true);
		}
	}

	/** Both entries reject (real production measurements, both aligners, see
	 * rrnafallback_hand_mutated_fixture_manifest_20260902.md's Allfail section) --
	 * exercises the once-only all-fail rrnaRejectedCandidates/orfScore-penalty path that
	 * neither the rescue nor pass fixture reaches (Citan's assignment, round 4/5,
	 * 2026-09-02). Counter assertions happen inside runFixtureCase itself. */
	private static void testAllFailPenaltyOnce(){
		for(boolean useIDA : new boolean[]{true, false}){
			runFixtureCase("allfail", useIDA, true, QUERY, ALLFAIL_ENTRY0, ALLFAIL_ENTRY1, false, 2, false);
		}
	}

	/** Forces a real histogram-row GROWTH, not just initial allocation (Citan's
	 * correction, 2026-09-02): runs the PASS fixture (1 attempt) first on a SHARED
	 * GeneCaller instance, then the RESCUE fixture (2 attempts) on the SAME instance/type
	 * -- the second call must grow rrnaAttemptCountHist[r16S] from size 2 (indices 0,1)
	 * to size 3 (indices 0,1,2), preserving the existing bucket-1 count while adding
	 * bucket 2, not silently losing it. Uses only already-built fixtures, no new ones. */
	private static void testAttemptCountHistogramGrowth(){
		final boolean oldIDA=GeneCaller.useIDAligner, oldFallback=ProkObject.rrnaFallback;
		final Read[] oldR16S=ProkObject.r16SSequence;
		try{
			GeneCaller.useIDAligner=true;
			final GeneCaller gc=new GeneCaller(1, 1, 1, 0f, 0f, 0f, 0f, 0f, new GeneModel(false));
			final StatsContainer sc=new StatsContainer(ProkObject.r16S);

			//First: PASS fixture, 1 attempt -- row starts null, grows to size 2 ([0]=0,[1]=1).
			ProkObject.rrnaFallback=false;
			ProkObject.r16SSequence=new Read[]{new Read(PASS_ENTRY0, null, "entry0", 0), new Read(PASS_ENTRY1_SENTINEL, null, "entry1", 1)};
			Orf orf=new Orf("synthetic", 50, 249, Shared.PLUS, 0, QUERY, false, ProkObject.r16S);
			check(gc.refineByAlignment(orf, QUERY, Shared.PLUS, sc), "pass fixture unexpectedly rejected");
			check(gc.rrnaAttemptCountHist[ProkObject.r16S].length==2, "expected row length 2 after 1-attempt candidate, got "
				+gc.rrnaAttemptCountHist[ProkObject.r16S].length);
			check(gc.rrnaAttemptCountHist[ProkObject.r16S][1]==1, "expected bucket[1]==1 after pass fixture");

			//Second: RESCUE fixture, 2 attempts -- row must GROW from size 2 to size 3,
			//preserving bucket[1]==1 (Citan's explicit growth-test requirement) while adding
			//a new bucket[2]==1.
			ProkObject.rrnaFallback=true;
			ProkObject.r16SSequence=new Read[]{new Read(RESCUE_ENTRY0, null, "entry0", 0), new Read(RESCUE_ENTRY1, null, "entry1", 1)};
			orf=new Orf("synthetic", 50, 249, Shared.PLUS, 0, QUERY, false, ProkObject.r16S);
			check(gc.refineByAlignment(orf, QUERY, Shared.PLUS, sc), "rescue fixture unexpectedly rejected");
			check(gc.rrnaAttemptCountHist[ProkObject.r16S].length==3, "expected row to GROW to length 3 after a 2-attempt candidate, got "
				+gc.rrnaAttemptCountHist[ProkObject.r16S].length);
			check(gc.rrnaAttemptCountHist[ProkObject.r16S][1]==1, "bucket[1] did not survive the resize -- growth lost existing data");
			check(gc.rrnaAttemptCountHist[ProkObject.r16S][2]==1, "bucket[2] was not correctly added after growth");
		}finally{
			GeneCaller.useIDAligner=oldIDA;
			ProkObject.rrnaFallback=oldFallback;
			ProkObject.r16SSequence=oldR16S;
		}
	}

	/** The opt-in ordered policy applies to 18S too; the prior first-only special case is retired. */
	private static void testR18SFallbackUsesFileOrder(){
		final boolean oldIDA=GeneCaller.useIDAligner, oldFallback=ProkObject.rrnaFallback;
		final Read[] oldR18S=ProkObject.r18SSequence;
		final float oldMin18S=ProkObject.min18SIdentity;
		try{
			GeneCaller.useIDAligner=true;
			ProkObject.rrnaFallback=true;
			ProkObject.min18SIdentity=0.62f;
			ProkObject.r18SSequence=new Read[]{new Read(RESCUE_ENTRY0, null, "entry0", 0), new Read(RESCUE_ENTRY1, null, "entry1", 1)};
			final ArrayList<Attempt> attempts=new ArrayList<Attempt>();
			final GeneCaller gc=new GeneCaller(1, 1, 1, 0f, 0f, 0f, 0f, 0f, new GeneModel(false));
			gc.setAttemptSink(new RefinementAttemptSink(){
				@Override
				public void onAttempt(int index, int total, int start, int stop, float score, boolean accepted, int postStart, int postStop, int strand, String scaffold, String consensusLabel, int consensusLength, float identity, String reason, String candidateId){
					attempts.add(new Attempt(index, total, start, stop, score, accepted));
				}
			});
			final Orf orf=new Orf("synthetic", 50, 249, Shared.PLUS, 0, QUERY, false, ProkObject.r18S);
			check(gc.refineByAlignment(orf, QUERY, Shared.PLUS, new StatsContainer(ProkObject.r18S)), "18S fallback did not accept entry1");
			check(attempts.size()==2 && attempts.get(0).index==0 && !attempts.get(0).accepted
					&& attempts.get(1).index==1 && attempts.get(1).accepted, "18S did not attempt ordered fallback");
		}finally{
			GeneCaller.useIDAligner=oldIDA;
			ProkObject.rrnaFallback=oldFallback;
			ProkObject.r18SSequence=oldR18S;
			ProkObject.min18SIdentity=oldMin18S;
		}
	}

	private static void runFixtureCase(String label, boolean useIDA, boolean fallback, byte[] query,
			byte[] entry0, byte[] entry1, boolean expectedResult, int expectedAttempts, boolean expectedLastAccepted){
		final boolean oldIDA=GeneCaller.useIDAligner, oldFallback=ProkObject.rrnaFallback;
		final Read[] oldR16S=ProkObject.r16SSequence;
		final float oldMin16S=ProkObject.min16SIdentity;
		try{
			check(oldMin16S==0.62f, "Fixture identities require min16SIdentity=0.62, got "+oldMin16S);
			GeneCaller.useIDAligner=useIDA;
			ProkObject.rrnaFallback=fallback;
			ProkObject.r16SSequence=new Read[]{new Read(entry0, null, "entry0", 0), new Read(entry1, null, "entry1", 1)};
			final ArrayList<Attempt> attempts=new ArrayList<Attempt>();
			final GeneCaller gc=new GeneCaller(1, 1, 1, 0f, 0f, 0f, 0f, 0f, new GeneModel(false));
			gc.setAttemptSink(new RefinementAttemptSink(){
				@Override
				public void onAttempt(int index, int total, int start, int stop, float score, boolean accepted, int postStart, int postStop, int strand, String scaffold, String consensusLabel, int consensusLength, float identity, String reason, String candidateId){
					attempts.add(new Attempt(index, total, start, stop, score, accepted));
				}
			});
			final StatsContainer sc=new StatsContainer(ProkObject.r16S);
			final Orf orf=new Orf("synthetic", 50, 249, Shared.PLUS, 0, query, false, ProkObject.r16S);
			final int start0=orf.start, stop0=orf.stop;
			final float score0=orf.orfScore;
			final boolean result=gc.refineByAlignment(orf, query, Shared.PLUS, sc);
			check(result==expectedResult, label+" useIDA="+useIDA+" fallback="+fallback+": result="+result);
			check(attempts.size()==expectedAttempts, label+": expected "+expectedAttempts+" attempts, got "+attempts.size());
			for(int i=0; i<attempts.size(); i++){
				Attempt a=attempts.get(i);
				check(a.index==i && a.total==2, label+": unexpected attempt "+i+" state");
			}
			check(attempts.get(attempts.size()-1).accepted==expectedLastAccepted, label+": unexpected final acceptance");
			if(attempts.size()>1){
				Attempt a=attempts.get(1);
				check(a.start==start0 && a.stop==stop0 && a.score==score0, label+": fallback saw mutated candidate state");
			}

			//Counter assertions (Citan-authorized, 2026-09-02): a fresh gc per call means
			//these are absolute post-call values, not deltas. Derived directly from the
			//already-known expectedAttempts/expectedResult/expectedLastAccepted, so every
			//existing runFixtureCase caller gets counter coverage for free.
			final int type=ProkObject.r16S;
			check(gc.rrnaCandidatesEnteringRefinement[type]==1, label+": expected 1 candidate entering refinement");
			check(gc.rrnaUniversalAttempts[type]==1, label+": expected exactly 1 universal attempt");
			check(gc.rrnaFallbackAttempts[type]==(expectedAttempts==2 ? 1 : 0), label+": unexpected fallback-attempt count");
			final boolean rescued=(expectedResult && expectedAttempts==2);
			check(gc.rrnaFallbackOnlyRescues[type]==(rescued ? 1 : 0), label+": unexpected fallback-only-rescue count");
			check(gc.rrnaRejectedCandidates[type]==(expectedResult ? 0 : 1), label+": unexpected rejected-candidate count");
			check(gc.rrnaAttemptCountHist[type].length>expectedAttempts, label+": histogram row too short for "+expectedAttempts+" attempts");
			check(gc.rrnaAttemptCountHist[type][expectedAttempts]==1, label+": expected histogram bucket["+expectedAttempts+"]==1");
			if(!expectedResult){
				//Once-only all-fail penalty (GeneCaller.java:1083): applied to whatever
				//orfScore was at loop exit, which the per-attempt restore-on-failure logic
				//guarantees is still score0 (no attempt in these fixtures ever accepts, so
				//nothing after score0 ever legitimately changes it).
				final float expectedPenalized=Math.min(-999f, score0-9999f);
				check(orf.orfScore==expectedPenalized, label+": all-fail orfScore penalty not applied exactly once "
					+"(expected "+expectedPenalized+", got "+orf.orfScore+")");
			}
		}finally{
			GeneCaller.useIDAligner=oldIDA;
			ProkObject.rrnaFallback=oldFallback;
			ProkObject.r16SSequence=oldR16S;
			ProkObject.min16SIdentity=oldMin16S;
		}
	}

	/** Focused output assertion for the two private CallGenes stats helpers (Citan's
	 * review correction, round 2, 2026-09-02: the histogram was collected/merged but
	 * never emitted -- this both proves the fix and gives the two formatting methods
	 * their first automated coverage at all, via reflection since they're private).
	 * Covers all 5 scalar fields plus the histogram string, in BOTH text and JSON. */
	private static void testStatsOutputFormatting() throws Exception {
		final long[] generated=new long[8], entered=new long[8], universal=new long[8], fallback=new long[8],
				rescues=new long[8], rejected=new long[8];
		//generated was added to the stats signature after this test's first version
		//(G11 port maintenance, 2026-09-05): keep it distinct from entered to prove
		//the two fields aren't cross-wired in either output path.
		generated[ProkObject.r16S]=4;
		entered[ProkObject.r16S]=3; universal[ProkObject.r16S]=3; fallback[ProkObject.r16S]=1;
		rescues[ProkObject.r16S]=1; rejected[ProkObject.r16S]=1;
		final long[][] hist=new long[8][];
		hist[ProkObject.r16S]=new long[]{0, 2, 1}; //bucket1=2, bucket2=1, bucket0 deliberately zero (must be omitted from the sparse string)

		//Direct coverage of the shared formatter first (package-visible, no reflection needed).
		check(CallGenes.formatAttemptHist(hist[ProkObject.r16S]).equals("1:2,2:1"),
			"formatAttemptHist produced wrong sparse string: "+CallGenes.formatAttemptHist(hist[ProkObject.r16S]));
		check(CallGenes.formatAttemptHist(null).equals(""), "formatAttemptHist(null) should be empty string");

		final Method printMethod=CallGenes.class.getDeclaredMethod("printRrnaFallbackStats",
				ByteStreamWriter.class, String.class, int.class,
				long[].class, long[].class, long[].class, long[].class, long[].class, long[].class, long[][].class);
		printMethod.setAccessible(true);
		final File tmp=new File(System.getProperty("java.io.tmpdir"),
				"CallGenesRrnaFallbackConfigTest_stats_"+System.nanoTime()+".txt");
		try{
			final ByteStreamWriter bsw=new ByteStreamWriter(tmp.getAbsolutePath(), true, false, false);
			bsw.start();
			printMethod.invoke(null, bsw, "16S", ProkObject.r16S, generated, entered, universal, fallback, rescues, rejected, hist);
			bsw.poisonAndWait();
			final String text=readFile(tmp);
			check(text.contains("entered=3"), "text output missing entered=3: "+text);
			check(text.contains("generated=4"), "text output missing generated=4: "+text);
			check(text.contains("universal=3"), "text output missing universal=3: "+text);
			check(text.contains("fallback=1"), "text output missing fallback=1: "+text);
			check(text.contains("rescues=1"), "text output missing rescues=1: "+text);
			check(text.contains("rejected=1"), "text output missing rejected=1: "+text);
			check(text.contains("attemptHist=1:2,2:1"), "text output missing correct histogram string: "+text);
		}finally{
			tmp.delete();
		}

		final Method jsonMethod=CallGenes.class.getDeclaredMethod("addRrnaFallbackStatsJson",
				JsonObject.class, String.class, int.class,
				long[].class, long[].class, long[].class, long[].class, long[].class, long[].class, long[][].class);
		jsonMethod.setAccessible(true);
		final JsonObject jo=new JsonObject();
		jsonMethod.invoke(null, jo, "16S", ProkObject.r16S, generated, entered, universal, fallback, rescues, rejected, hist);
		final String jsonText=jo.toString();
		check(jsonText.contains("\"16S Fallback Entered\": 3"), "JSON missing entered field: "+jsonText);
		check(jsonText.contains("\"16S Fallback Generated\": 4"), "JSON missing generated field: "+jsonText);
		check(jsonText.contains("\"16S Fallback Universal\": 3"), "JSON missing universal field: "+jsonText);
		check(jsonText.contains("\"16S Fallback Attempts\": 1"), "JSON missing fallback-attempts field: "+jsonText);
		check(jsonText.contains("\"16S Fallback Rescues\": 1"), "JSON missing rescues field: "+jsonText);
		check(jsonText.contains("\"16S Fallback Rejected\": 1"), "JSON missing rejected field: "+jsonText);
		check(jsonText.contains("\"16S Fallback AttemptHist\": \"1:2,2:1\""), "JSON missing correct histogram string: "+jsonText);
	}

	private static String readFile(File f) throws Exception {
		final StringBuilder sb=new StringBuilder();
		final BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		while((line=br.readLine())!=null){sb.append(line).append("\n");}
		br.close();
		return sb.toString();
	}

	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}

	private static final class Attempt {
		Attempt(int index_, int total_, int start_, int stop_, float score_, boolean accepted_){index=index_; total=total_; start=start_; stop=stop_; score=score_; accepted=accepted_;}
		final int index, total, start, stop;
		final float score;
		final boolean accepted;
	}

	private static final class ConfigSnapshot {
		ConfigSnapshot(){
			fallback=ProkObject.rrnaFallback;
			call16=ProkObject.call16S; call18=ProkObject.call18S; call23=ProkObject.call23S; call5=ProkObject.call5S;
			load16=ProkObject.load16SSequence; load18=ProkObject.load18SSequence; load23=ProkObject.load23SSequence; load5=ProkObject.load5SSequence;
			id16=ProkObject.min16SIdentity; id18=ProkObject.min18SIdentity; id23=ProkObject.min23SIdentity; id5=ProkObject.min5SIdentity;
			ssuStart=ProkObject.ssuStartSlop; ssuStop=ProkObject.ssuStopSlop; lsuStart=ProkObject.lsuStartSlop; lsuStop=ProkObject.lsuStopSlop; r5Start=ProkObject.r5SStartSlop; r5Stop=ProkObject.r5SStopSlop;
		}
		void assertUnchangedExceptFallback(){
			check(call16==ProkObject.call16S && call18==ProkObject.call18S && call23==ProkObject.call23S && call5==ProkObject.call5S, "rrnafallback changed call flags");
			check(load16==ProkObject.load16SSequence && load18==ProkObject.load18SSequence && load23==ProkObject.load23SSequence && load5==ProkObject.load5SSequence, "rrnafallback changed consensus-load flags");
			check(id16==ProkObject.min16SIdentity && id18==ProkObject.min18SIdentity && id23==ProkObject.min23SIdentity && id5==ProkObject.min5SIdentity, "rrnafallback changed identity thresholds");
			check(ssuStart==ProkObject.ssuStartSlop && ssuStop==ProkObject.ssuStopSlop && lsuStart==ProkObject.lsuStartSlop && lsuStop==ProkObject.lsuStopSlop && r5Start==ProkObject.r5SStartSlop && r5Stop==ProkObject.r5SStopSlop, "rrnafallback changed slop settings");
		}
		void restore(){
			ProkObject.rrnaFallback=fallback;
			ProkObject.call16S=call16; ProkObject.call18S=call18; ProkObject.call23S=call23; ProkObject.call5S=call5;
			ProkObject.load16SSequence=load16; ProkObject.load18SSequence=load18; ProkObject.load23SSequence=load23; ProkObject.load5SSequence=load5;
			ProkObject.min16SIdentity=id16; ProkObject.min18SIdentity=id18; ProkObject.min23SIdentity=id23; ProkObject.min5SIdentity=id5;
			ProkObject.ssuStartSlop=ssuStart; ProkObject.ssuStopSlop=ssuStop; ProkObject.lsuStartSlop=lsuStart; ProkObject.lsuStopSlop=lsuStop; ProkObject.r5SStartSlop=r5Start; ProkObject.r5SStopSlop=r5Stop;
		}
		final boolean fallback, call16, call18, call23, call5, load16, load18, load23, load5;
		final float id16, id18, id23, id5;
		final int ssuStart, ssuStop, lsuStart, lsuStop, r5Start, r5Stop;
	}

	private static final byte[] QUERY="GGGTTTCCTTTAGCACCCTGTGCAACTGGAAAGGGACTTAAAACCGGGTGGAGCTAATCCCAGCGGGATCGTTTCGAAGTTCTCTGTCACCCCGTAATTGAGATGATCCGGATCCGCTCTCCCTCCCATTTCAGTCGAAAGTGCCTAACAAAGGAAGATCGAGGGATAGTAGATTCCTGACACTTTAAAGGTGGGGTCCCGATGGTATTCGCCGAGTGCATATTACCGATTAGTACTCACTCATCATGTTTGTGCGCGGAGTGCCTTGTCGTGTGATGCAGTAGGAACGGTTGTATTTTG".getBytes();
	private static final byte[] RESCUE_ENTRY0="GGACTGAGATCAGCAGGATCGCCTGGCCTAGCTCTGCGGGTCCTGAATTTCGGCTATCCTGACCCAATTTTCCTGCCATGTAACGTGAAACGGGCACACAGGGTAAGCGCGCGGGCCAACCGACTGCCGACAACTTACTTGCAGCTTCCATCGAGTGTCGGCACAGTACATAAACTTGATTACGTCAAACGCACGGCGTT".getBytes();
	private static final byte[] RESCUE_ENTRY1="GAGCCAGTCCCAGCGGGATAGTTTCGAAGTTGTCTGTCACCCCGTAATTGAGATGATCCGTATCCGCTCTCCCTCCCATTTCAGTCGAAAGTGCCTAACAAAGGAAGATCGAGGGATAGTAGATTCCTGACACTTTAAAGGTCGGTTCCCGATGGTATTCGCCGAGGGCATATTACCGATAAGTACTCACTCATCATGTT".getBytes();
	private static final byte[] PASS_ENTRY0="GAGCTAATCCCAGCGGGATCGTTTCGAAGTTCTCTGTCACCCCGTAATTGAAATGATCCGGATCCGCTCTCCCTCCCATTTCAGTCGAAAGTGCCTAACAAAGGAAGATCGAGGGATAGTAGATTCCTGACACTTTAAAGGTGGGGTCCCGATGGTATTCGCCGAGTGCATATTACCGATTAGTACTCACTCATCATGTT".getBytes();
	private static final byte[] PASS_ENTRY1_SENTINEL="CTCGGCCTGGTGCACCATATCAAAACGGTTAGTCCGCAAATTTTGGGCCATCTTTGGCCGGGCTCGAAACGGGCACAGGGGTAGAGGACCCGGGTCCACAAGGAGCAGCCACCCATCTAAATAAGCAAAACGATGTTAAGAACTAATTTATTTAAACCGTAGTTTTATGCATGAGTACTCTTGGGAAGCTCACCGTGTCG".getBytes();
	private static final byte[] ALLFAIL_ENTRY0="CGGCACGACATAAAGGTTAGGCTCCGGAGTTCTCGGGCTCCACCTCCCCGAATCATATGATTCCCTCTAACGTAGGTAAAAAAGTCGCTCATGGTTAATTAAGGACAACCGATGCAGTTTGGGCAATAAACTCTTCCATGGTCCCGTACGGAAGCCTGACGCGCTGTTAAACTATTATGCTATCTCTAGCCGATGCTGGC".getBytes();
	private static final byte[] ALLFAIL_ENTRY1="CGGCACGACATAGCCTGTAGGCTCCGGAGTTCTCGGGCTCCACCTCCCTGTGTTACTATAGTTCGGCTCTAACCCTTCTCTTATACAGAAGTGGCCATCCATGGATTATCGAGATCTAGGCGGTTTCAGAGAATCATTTGACCCGGTCGGGGGGATCGTCACGCTAGCCTTAAGAAGTATATTATGTGAGTTATTATGCT".getBytes();
}
