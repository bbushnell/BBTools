package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import shared.Shared;
import stream.Read;

/** Compares the selected packaged resources with explicit copies through production APIs.
 * @author Brian, Raiden */
public final class Euk18sFamilyResourceTest {
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("baseline=")){
			throw new IllegalArgumentException("Supply baseline=/selected/explicit/resources");
		}
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;
		final File baseline=new File(args[0].substring(9));
		final Euk18sRuntimeConfig defaults=config();
		final NcrnaFamily actual=load(defaults);
		final Euk18sRuntimeConfig explicit=config();
		explicit.parse("euk18sconsensus", new File(baseline, "euk18S_consensus_sequence.fa").toString());
		explicit.parse("euk18sseeds", new File(baseline, "euk18S_17mers.fa.gz").toString());
		explicit.parse("euk18smodelcutoffs", new File(baseline, "euk18S_model_cutoffs.tsv").toString());
		explicit.parse("euk18svoteoffsets", new File(baseline, "euk18S_seed_offsets.tsv").toString());
		explicit.parse("euk18svoteseedsha80", "47bd7e395c8398204083");
		explicit.parse("euk18slivepositions", new File(baseline, "euk18S_seed_positions.tsv").toString());
		final NcrnaFamily old=load(explicit);
		assert(Arrays.deepEquals(actual.library, old.library) && Arrays.equals(actual.modelNames, old.modelNames))
			: "The moved modern library must retain the exact ordered sequences and names";
		assert(actual.library.length==10 && actual.models==null && actual.rrnaEndpointFeatures==null)
			: "The modern package has ten models and no default HBM or endpoint network";
		final long[] keys=actual.kmerSet.toArray(), oldKeys=old.kmerSet.toArray();Arrays.sort(keys);Arrays.sort(oldKeys);
		assert(Arrays.equals(keys, oldKeys) && keys.length==395) : "Shared scanning must retain all395 selected caller-visible forward17mers";
		assert(actual.kLong==17 && actual.indexK==9 && actual.indexTopN==10 && actual.fixedMinHits==96
			&& actual.adaptive && actual.adaptFloor==96 && actual.adaptTopFrac==.8f && actual.adaptQFrac==0)
			: "Migration must preserve the measured 18S shortlist policy";
		assert(actual.minLen==916 && actual.maxLen==4000 && actual.windowPad==3109 && actual.seedMinHits==2
			&& actual.seedDistinct && actual.collapseFrac==.9f && actual.scoreA==0 && actual.scoreB==1)
			: "Migration must preserve candidate geometry, seed support and path scoring";
		assert(actual.pacBioConsensusAlignment && actual.pacBioRolling && !actual.modelEndClipping
			&& !actual.modelClipRescue && Float.isNaN(actual.modelClipRescueId) && actual.joinedSideSupport==1
			&& actual.joinedPositions.modelCount()==10 && actual.joinedPositions.seedCount()==keys.length
			&& actual.voteWindows && actual.voteBeforePadding && !actual.voteEnds && actual.voteSlack==60 && actual.voteEndsMaxSd==5f)
			: "Resource loading preserves rolling/joined/voting policies; interior rescue defaults off before physical-edge registration";
		NcrnaFamilyRegistry.enablePhysicalEdges(actual);
		assert(!actual.modelClipRescue && actual.modelClipRescueId==.66f && actual.physicalEdgeMinLen==actual.minLen)
			: "Default physical-edge registration retains the tested18S rescue cutoff without granting interior rescue";
		for(int i=0; i<actual.library.length; i++){
			assert(actual.modelThresholds.pass(i)==old.modelThresholds.pass(i)
				&& actual.modelThresholds.borderline(i)==old.modelThresholds.borderline(i))
				: "Every model must retain its old pass and borderline thresholds";
		}
		for(long key:keys){
			final SeedOffsetTable.KmerInfo a=actual.voteTable.get(key), b=old.voteTable.get(key);
			assert(a!=null && b!=null && a.leftOffset==b.leftOffset && a.rightOffset==b.rightOffset
				&& a.leftSD==b.leftSD && a.rightSD==b.rightSD && a.count==b.count)
				: "Moved voting offsets must preserve every trained key and uncertainty";
		}
		final Read[] legacy=ProkObject.loadConsensusSequenceType("18S", false, false);
		final byte[][] oldLegacy=TrnaConsensusBuilder.loadLibrary(new File(baseline, "18S_consensus_sequence.fa").toString());
		assert(legacy.length==oldLegacy.length && legacy.length!=actual.library.length)
			: "Legacy SSU consumers must receive their distinct library, never the modern ten models";
		for(int i=0; i<legacy.length; i++){
			assert(Arrays.equals(legacy[i].bases, oldLegacy[i])) : "Every ordered legacy18S sequence must survive the move";
		}
		final Euk18sRuntimeConfig missing=config();
		missing.parse("euk18svoteoffsets", new File(baseline, "definitely_missing_18s_offsets.tsv").toString());
		final ArrayList<NcrnaFamily> registered=new ArrayList<NcrnaFamily>();
		try{missing.apply(registered);throw new AssertionError("Missing required votes must fail before registration");}
		catch(IllegalArgumentException expected){
			assert(registered.isEmpty() && expected.getMessage().contains("vote_offsets"))
				: "Missing vote resources must fail at the expected binding without partial registration";
		}
		System.out.println("PASS modern18S resources, unchanged policies, distinct legacy library and atomic missing-offset failure");
	}
	private static Euk18sRuntimeConfig config(){
		final Euk18sRuntimeConfig c=new Euk18sRuntimeConfig();
		assert(c.parse("euk18s", "t")) : "The actual public18S flag must select the tested factory";return c;
	}
	private static NcrnaFamily load(Euk18sRuntimeConfig c){
		final ArrayList<NcrnaFamily> list=new ArrayList<NcrnaFamily>();c.apply(list);
		assert(list.size()==1) : "A complete18S load must publish exactly one family";return list.get(0);
	}
}
