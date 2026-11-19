package test.prok;

import java.util.Arrays;

import prok.NcrnaDdlProfile;
import prok.NcrnaDdlQueryKeys;

/** Finite checks of plurality, nonempty-key ties, bucket matching and ownership.
 * Expected votes are literal small examples, independent of sequence hashing.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlProfileTest{

	public static void main(final String[] args){
		testPluralityAndTies();
		testWideProfile();
		testOwnership();
		testConsensusWeights();
		testInvalidInputs();
		System.err.println("NCRNA_DDL_PROFILE_TEST_PASS");
	}

	/** Two of five votes can win; equal support prefers key5 over key7. */
	private static void testPluralityAndTies(){
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(32);
		final char[] member=new char[32];
		final int[] first={7, 5, 5, 7, 11};
		final int[] second={4, 4, 9, 0, 0};
		for(int i=0; i<first.length; i++){
			member[0]=(char)first[i];
			member[1]=(char)second[i];
			votes.add(member);
		}
		final NcrnaDdlProfile profile=votes.finish();
		assert(profile.members==5 && profile.keys()[0]==5)
			: "A two-of-five plurality must remain eligible, with smaller-key tie resolution";
		assert(profile.keys()[1]==4)
			: "Empty member buckets cannot beat observed key4";
		assert(profile.keys()[2]==0)
			: "An entirely empty bucket must not invent a voted key or support";
		final NcrnaDdlProfile highest=votes.finish(null, 0);
		assert(highest.keys()[0]==7 && highest.members==5)
			: "The explicit highest-key baseline must preserve member support while leaving legacy ties unchanged";
		Arrays.fill(member, (char)0);
		member[0]=7;
		member[1]=9;
		final int nonwinning=profile.equalBuckets(member);
		assert(nonwinning==0) : "Only the plurality key matches, not another observed member key";
		member[0]=5;
		member[1]=4;
		final int winning=profile.equalBuckets(member);
		assert(winning==2) : "Each matching nonempty profile bucket contributes one match";
		Arrays.fill(member, (char)0);
		final int empty=profile.equalBuckets(member);
		assert(empty==0) : "Empty query buckets never match empty profile buckets";
	}

	/** The widest experiment geometry must retain the full unsigned encoded key. */
	private static void testWideProfile(){
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(2048);
		final char[] member=new char[2048];
		member[2047]=(char)65535;
		votes.add(member);
		final NcrnaDdlProfile profile=votes.finish();
		final int matches=profile.equalBuckets(member);
		assert(profile.buckets()==2048 && profile.keys()[2047]==65535 && profile.members==1 && matches==1)
			: "Bucket2047 and encoded key65535 must survive packed primitive vote storage";
	}

	/** Neither a reused source sketch nor later votes may mutate a finished profile. */
	private static void testOwnership(){
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(32);
		final char[] member=new char[32];
		member[0]=9;
		votes.add(member);
		final NcrnaDdlProfile saved=votes.finish();
		member[0]=4;
		votes.add(member);
		votes.add(member);
		final NcrnaDdlProfile changed=votes.finish();
		assert(saved.members==1 && saved.keys()[0]==9)
			: "Finished profiles must own their keys independently of accumulator updates";
		assert(changed.members==3 && changed.keys()[0]==4)
			: "Reusing a member buffer must still contribute exactly one vote per addition";
	}

	/** Real colliding consensus words distinguish raw-pool weighting from first
	 * reducing consensus to its maximum. Votes below are literal counts, not
	 * recomputed by a second selector implementation.
	 */
	private static void testConsensusWeights(){
		final byte[] bases=new byte[512];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		long state=17;
		for(int i=0; i<bases.length; i++){
			state^=state<<13;
			state^=state>>>7;
			state^=state<<17;
			bases[i]=alphabet[(int)state&3];
		}
		final NcrnaDdlQueryKeys consensus=new NcrnaDdlQueryKeys(32, 10);
		consensus.build(bases, null, false, true);
		final int[] raw=consensus.packedKeys();
		int bucket=-1, low=0, high=0;
		for(int i=0; i<raw.length && bucket<0; i++){
			for(int j=i+1; j<raw.length; j++){
				if((raw[i]>>>16)==(raw[j]>>>16)){
					bucket=raw[i]>>>16;
					low=Math.min(raw[i]&65535, raw[j]&65535);
					high=Math.max(raw[i]&65535, raw[j]&65535);
					break;
				}
			}
		}
		assert(bucket>=0 && low>0 && high>low) : "The raw-pool fixture needs two distinct consensus keys in one bin";
		int outside=1;
		while(consensus.contains((bucket<<16)|outside)){outside++;}
		assert(outside<=65535) : "The finite consensus cannot contain every encoded key in the witness bin";
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(32);
		final char[] member=new char[32];
		member[bucket]=(char)low;
		for(int i=0; i<3; i++){votes.add(member);}
		member[bucket]=(char)high;
		votes.add(member);
		member[bucket]=(char)outside;
		for(int i=0; i<6; i++){votes.add(member);}
		final NcrnaDdlProfile baseline=votes.finish(consensus, 0);
		final NcrnaDdlProfile one=votes.finish(consensus, 1);
		final NcrnaDdlProfile four=votes.finish(consensus, 4);
		final NcrnaDdlProfile sixteen=votes.finish(consensus, 16);
		final NcrnaDdlProfile always=votes.finish(consensus, NcrnaDdlProfile.CONSENSUS_WINS);
		assert(baseline.keys()[bucket]==outside && one.keys()[bucket]==outside)
			: "Six real votes must beat three consensus votes with weight0 or1";
		assert(four.keys()[bucket]==low && sixteen.keys()[bucket]==low && always.keys()[bucket]==low)
			: "Boost every raw consensus key: lower key with three member votes beats the consensus maximum with one";
		assert(always.members==10 && baseline.members==10)
			: "Consensus weights must never change the number of real construction members";
		int emptyMemberBin=-1;
		for(final int packed:raw){
			if((packed>>>16)!=bucket){emptyMemberBin=packed>>>16; break;}
		}
		assert(emptyMemberBin>=0 && always.keys()[emptyMemberBin]!=0
			&& baseline.keys()[emptyMemberBin]==0)
			: "A consensus-only key may fill an empty member bin without creating additional voters";
		final NcrnaDdlProfile again=votes.finish(consensus, 0);
		assert(Arrays.equals(baseline.keys(), again.keys())) : "Selecting weighted variants cannot mutate original vote counts";
		consensus.build(new byte[0], null, false, true);
		final NcrnaDdlProfile empty=votes.finish(consensus, NcrnaDdlProfile.CONSENSUS_WINS);
		assert(Arrays.equals(empty.keys(), baseline.keys()) && always.keys()[bucket]==low)
			: "Empty consensus leaves member selection alone; later buffer reuse cannot alter an earlier result";
	}

	/** A missing member set or mismatched width cannot become a plausible profile. */
	private static void testInvalidInputs(){
		boolean rejected=false;
		try{new NcrnaDdlProfile.Votes(33);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "Non-power-of-two bucket counts cannot share native DDL geometry";
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(32);
		rejected=false;
		try{votes.finish();}catch(IllegalStateException expected){rejected=true;}
		assert(rejected) : "Every model profile must retain at least one actual assigned member";
		rejected=false;
		try{votes.add(new char[64]);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "A member sketch from another width must fail before votes are counted";
		votes.add(new char[32]);
		rejected=false;
		try{votes.finish(null, 1);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "A positive consensus weight requires the actual consensus keys";
		rejected=false;
		try{votes.finish(new NcrnaDdlQueryKeys(64, 10), 1);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "Consensus and member bucket geometries cannot be mixed";
		rejected=false;
		try{votes.finish(null, -2);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "Only the documented always-win sentinel can stand for a negative weight";
	}
}
