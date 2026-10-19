package test.consensus;

import java.util.Arrays;

import consensus.BaseGraph;
import consensus.ConsensusObject;
import stream.Read;

/** Exact graph-weight regressions for internal and terminal deletions. @author Fischl */
public final class BaseGraphDeletionQualityTest {

	/** Exercises every deletion position with unequal flanking qualities. */
	public static void main(final String[] args){
		final boolean oldMapq=ConsensusObject.useMapq, oldInvert=ConsensusObject.invertIdentity;
		try{
			ConsensusObject.useMapq=false;
			ConsensusObject.invertIdentity=false;
			if(args.length==1 && args[0].equals("legacy")){
				boolean failed=false;
				try{fixture(8, 7, 1, true);}
				catch(ArrayIndexOutOfBoundsException expected){failed=true;}
				check(failed, "Legacy add must reproduce the quality overrun before the last query base.");
				System.out.println("DELETION_QUALITY_LEGACY_REPRODUCED");
				return;
			}
			check(args.length==0, "Use no arguments for regression tests, or legacy for the old-code reproducer.");
			for(int length : new int[]{1, 2, 8}){
				for(int pos=0; pos<length; pos++){
					for(int deleted=1; deleted<=3; deleted++){
						fixture(length, pos, deleted, true);
						fixture(length, pos, deleted, false);
					}
				}
			}
			fixture(8, 4, 0, true);
			fixture(8, 4, 0, false);
			boolean rejected=false;
			try{fixture(8, 8, 1, true);}
			catch(AssertionError expected){
				if(!"D".equals(expected.getMessage())){throw expected;}
				rejected=true;
			}
			check(rejected, "Read.identityFlat's existing terminal-D rejection must remain explicit, not silently bypassed.");
			System.out.println("DELETION_QUALITY_TEST_PASS cases="+cases+" checks="+checks);
		}finally{
			ConsensusObject.useMapq=oldMapq;
			ConsensusObject.invertIdentity=oldInvert;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs exact reference-oriented match symbols without a mapper or quality pseudocount. */
	private static void fixture(final int length, final int pos, final int deleted, final boolean withQualities){
		assert(length>0 && pos>=0 && pos<=length && deleted>=0) : "The deletion lies between valid query bases or at one end.";
		final byte[] bases=new byte[length], qualities=new byte[length];
		final byte[] reference=new byte[length+deleted], match=new byte[reference.length];
		for(int i=0; i<length; i++){
			bases[i]=(byte)"ACGT".charAt(i%4);
			qualities[i]=(byte)(3+5*i);
		}
		System.arraycopy(bases, 0, reference, 0, pos);
		Arrays.fill(reference, pos, pos+deleted, (byte)'A');
		System.arraycopy(bases, pos, reference, pos+deleted, length-pos);
		Arrays.fill(match, (byte)'m');
		Arrays.fill(match, pos, pos+deleted, (byte)'D');
		final Read read=new Read(bases.clone(), withQualities ? qualities.clone() : null, "deletion", 0);
		read.match=match;
		read.start=0;
		read.stop=reference.length-1;
		read.setMapped(true);
		read.setShortMatch(false);
		final BaseGraph graph=new BaseGraph("reference", reference, null, 0, 0);
		graph.add(read);
		final int q=withQualities ? (qualities[Math.max(0, pos-1)]+qualities[Math.min(pos, length-1)])/2 : BaseGraph.fakeQuality;
		for(int i=0; i<reference.length; i++){
			final boolean deletion=i>=pos && i<pos+deleted;
			check(graph.del[i].countSum==(deletion ? 1 : 0), "Deletion observations must cover exactly the reference gap.");
			check(graph.del[i].weightSum==(deletion ? q+1 : 0), "Deletion weight must use its actual flanks, or its sole terminal flank.");
			if(!deletion){
				final int query=i<pos ? i : i-deleted;
				final int expected=withQualities ? qualities[query]+2 : BaseGraph.fakeQuality+2;
				check(graph.ref[i].countSum==2 && graph.ref[i].weightSum==expected,
						"An ordinary match retains its reference pseudocount and own query quality.");
			}
		}
		check(graph.delTotal==deleted && graph.msTotal==length && graph.baseTotal==length,
				"A deletion consumes reference positions but never consumes query bases.");
		check(Arrays.equals(read.bases, bases) && (!withQualities || Arrays.equals(read.quality, qualities)),
				"Graph accumulation must leave caller-owned bases and qualities unchanged.");
		for(int i=0; i<reference.length; i++){
			graph.ref[i].acgtProb=new float[4];
			graph.ref[i].acgtProb["ACGT".indexOf(reference[i])]=1;
		}
		check(Float.isFinite(graph.score(read, false, false)) && Float.isFinite(graph.scoreOld(read, false, false)),
				"Both graph scorers must accept the supported deletion geometry used by accumulation.");
		cases++;
	}

	/** Keeps test expectations active even if an external caller disables assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static int cases, checks;
}
