package test.icecream;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import icecream.ReformatPacBio;
import icecream.ZMW;
import stream.Read;

/** Literal contributor/strand checks against an unstarted Reformat instance.
 * @author Brian Bushnell, Sucrose
 */
public final class PacBioConsensusContributorTest{

	/** Requires an existing tiny input only for constructor validation; starts no workers.
	 * @param args Exactly in=path
	 */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("in="), "Use in=existing-tiny-input");
		require(ReformatPacBio.class.desiredAssertionStatus(), "Production assertions must be enabled");
		final ReformatPacBio program=new ReformatPacBio(new String[]{args[0], "threads=1", "jni=f", "simd=f", "trim=f"});
		System.out.println("case\tcheck\texpected\tobserved\tresult");
		final int failures=observeCase(program, false)+observeCase(program, true);
		require(failures==0, failures+" consensus contributor observations differ from literal expectations");
	}

	/** Keeps original indices1/3 as the same-strand eligible reads in the discarded case. */
	private static int observeCase(ReformatPacBio program, boolean discard){
		final String label=discard ? "discarded_contributors" : "all_kept_alternating";
		final ZMW well=new ZMW(6);
		final Read[] original=new Read[6];
		final byte[][] initialBases=new byte[6][];
		final boolean[] initialDiscarded=new boolean[6];
		for(int i=0; i<6; i++){
			final byte[] bases=discard && i==5 ? null :
					(i%2==0 ? REVERSE : TEMPLATE).getBytes(StandardCharsets.US_ASCII);
			final Read read=new Read(bases, null, "m/42/"+(i*32)+"_"+((i+1)*32), i);
			initialDiscarded[i]=discard && i!=1 && i!=3;
			read.setDiscarded(initialDiscarded[i]);
			original[i]=read;
			initialBases[i]=bases==null ? null : bases.clone();
			well.add(read);
		}
		Read consensus=null;
		String completion="completed";
		try{consensus=program.makeConsensus(well);}
		catch(RuntimeException failure){
			completion=failure.getClass().getSimpleName();
			failure.printStackTrace(System.err);
		}
		int failures=0;
		failures+=observe(label, "completion", "completed", completion);
		failures+=observe(label, "consensus_bases", TEMPLATE,
				consensus==null || consensus.bases==null ? null : new String(consensus.bases, StandardCharsets.US_ASCII));
		failures+=observe(label, "average_identity", Float.valueOf(1.0f), consensus==null ? null : consensus.obj);
		boolean listPreserved=well.size()==6, basesPreserved=true, flagsPreserved=true;
		final StringBuilder strands=new StringBuilder();
		for(int i=0; i<6; i++){
			listPreserved&=i<well.size() && well.get(i)==original[i];
			basesPreserved&=Arrays.equals(initialBases[i], original[i].bases);
			flagsPreserved&=original[i].discarded()==initialDiscarded[i];
			strands.append(original[i].strand());
		}
		failures+=observe(label, "list_identity_and_order", Boolean.TRUE, listPreserved);
		failures+=observe(label, "original_index_strands", "010101", strands.toString());
		failures+=observe(label, "restored_input_bases", Boolean.TRUE, basesPreserved);
		failures+=observe(label, "discard_flags", Boolean.TRUE, flagsPreserved);
		return failures;
	}

	/** Prints every mismatch before the final failure so baseline controls remain visible. */
	private static int observe(String label, String check, Object expected, Object actual){
		final boolean pass=expected==null ? actual==null : expected.equals(actual);
		System.out.println(label+"\t"+check+"\t"+expected+"\t"+actual+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Keeps fixture setup/oracle checks active independently of Java assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Non-palindromic literal template and its independently written reverse complement. */
	private static final String TEMPLATE="AAAACCCCAAAACCCCAAAACCCCAAAACCCC";
	private static final String REVERSE="GGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTT";

	/** Static-only fixture. */
	private PacBioConsensusContributorTest(){}
}
