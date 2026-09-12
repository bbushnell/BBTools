package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import ukmer.Kmer;

/** Deterministic original-scan skip and whole-read transaction regressions.
 * Test-only string-key counts are independent of the production table backend.
 * @author Fischl */
public final class LocalEditReadGuardTest {
	public static void main(final String[] args){
		final boolean oldMask=Kmer.MASK_CORE,oldPacked=Kmer.PACKED;
		try{
			Kmer.MASK_CORE=false;Kmer.PACKED=true;
			for(int k:new int[]{31,62}){
				for(boolean reverse:new boolean[]{false,true}){
					checkCase(k,4*k,1,reverse,1); // Initial independent burden exceeds cap.
					checkCase(k,k+5,1,reverse,2); // Overlapping contexts undercount; discovery rolls back.
					checkCase(k,4*k,2,reverse,0); // Exactly at cap, exhausted, retain both edits.
				}
				indelRollback(k,false);indelRollback(k,true);
				preflightCost(k);depthGeometry(k);edgeAndUndefinedEstimate(k);
				mixedContext(k,false);mixedContext(k,true);ordinaryAndPair(k);
				unsupported(k);exceptionRestoration(k);pairedRejected(k);
			}
			System.out.println("LOCAL_EDIT_READ_GUARD_TEST_OK checks="+checks);
		}finally{Kmer.MASK_CORE=oldMask;Kmer.PACKED=oldPacked;}
	}
	private static void edgeAndUndefinedEstimate(final int k){
		final byte[] truth=sequence(12*k,k);final Counts counts=new Counts(k,truth);
		for(int mode=0;mode<3;mode++){
			final byte[] noisy=truth.clone();
			if(mode==0){noisy[0]=other(noisy[0]);noisy[noisy.length-1]=other(noisy[noisy.length-1]);}
			else{noisy[3*k]='N';noisy[7*k]='N';}
			if(mode==2){noisy[5*k]=other(noisy[5*k]);}
			final LocalEditCorrector corrector=new LocalEditCorrector(k,new HomopolymerIndelProposal.CountLookup(){
				@Override public int count(final Kmer key){return counts.count(key);}
			});
			final Read read=read(noisy);final byte[] original=read.bases;
			check(corrector.correctOne(read,false,1)==(mode==2 ? 1 : 0),"Edge or N-spanning troughs must not block an unrelated supported repair.");
			check(corrector.initialEstimatedEdits==(mode==2 ? 1 : 0),"Initial estimate excludes both read ends and complete N-spanning contexts.");
			check(corrector.callStatus!=LocalEditCorrector.CallStatus.INITIAL_LIMIT,"Excluded edge/N regions cannot force an initial skip.");
			if(mode<2){check(read.bases==original,"Edge/N-only profiles leave original arrays installed.");}
		}
	}
	private static void preflightCost(final int k){
		final byte[] truth=sequence(12*k,k),noisy=truth.clone();
		noisy[3*k]=other(noisy[3*k]);noisy[7*k]=other(noisy[7*k]);
		final Counts counts=new Counts(k,truth);
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return counts.count(key);}
		};
		final LocalEditCorrector guarded=new LocalEditCorrector(k,lookup),ordinary=new LocalEditCorrector(k,lookup);
		final Read rejected=read(noisy);final byte[] original=rejected.bases,quality=rejected.quality;
		check(guarded.correctOne(rejected,false,1)==0 && guarded.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT,"Separated original troughs must reject before editing.");
		check(guarded.initialEstimatedEdits==2 && guarded.profileQueries==noisy.length-k+1,"Preflight must reuse one depth fill and estimate two separated errors.");
		check(guarded.probeQueries==0 && guarded.verificationQueries==0 && guarded.pairQueries==0,"Depth-only rejection must not probe or verify any mutant, including pair witnesses.");
		check(rejected.bases==original && rejected.quality==quality,"Preflight cannot replace or mutate read arrays.");
		final Read a=read(noisy),b=read(noisy);
		check(guarded.correctOne(a,false,2)==1 && ordinary.correctOne(b)==1,"Nonbinding preflight must enter the original first-edit path.");
		check(Arrays.equals(a.bases,b.bases) && Arrays.equals(a.quality,b.quality),"Nonbinding preflight preserves edit order and quality output.");
		check(guarded.profileQueries==ordinary.profileQueries && guarded.probeQueries==ordinary.probeQueries && guarded.verificationQueries==ordinary.verificationQueries,"Nonbinding preflight must add zero count-table accesses, unlike repair-validation preflight.");
	}
	private static void depthGeometry(final int k){
		// Artificial profiles deliberately have no supported mutant keys. Their
		// K-wide dips still count as suspected errors, not verified repairs.
		for(int width:new int[]{k-2,k-1,k,k+1}){
			for(boolean weakFlanks:new boolean[]{false,true}){
				final byte[] bases=sequence(12*k,8123+k);final Counts counts=new Counts(k,bases);
				final Kmer key=new Kmer(k);
				for(int j=0;j<bases.length;j++){
					key.addRight(bases[j]);if(key.len()<k){continue;}
					final int start=j-k+1;
					final boolean low=(start>=2*k && start<2*k+width) || (start>=7*k && start<7*k+width);
					counts.map.put(Arrays.toString(key.key()),low ? 2 : weakFlanks ? 8 : 12);
				}
				final LocalEditCorrector corrector=new LocalEditCorrector(k,new HomopolymerIndelProposal.CountLookup(){
					@Override public int count(final Kmer candidate){return counts.count(candidate);}
				});
				final Read read=read(bases);final byte[] original=read.bases;
				check(corrector.correctOne(read,false,1)==0 && read.bases==original,"Unsupported artificial profile may not create an accepted repair.");
				final boolean reject=!weakFlanks && (width==k-1 || width==k);
				check((corrector.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT)==reject,"Only near-K troughs with strong contrast enter the initial burden estimate.");
				check(corrector.initialEstimatedEdits==(reject ? 2 : 0),"Narrow/merged/weak-flank troughs are not counted as independent single errors.");
				if(reject){check(corrector.probeQueries==0 && corrector.verificationQueries==0,"Suspected-error rejection needs no repair evidence lookups.");}
			}
		}
	}
	private static void checkCase(final int k,final int separation,final int cap,final boolean reverse,final int outcome){
		final byte[] truth=sequence(12*k,k),noisy=truth.clone();
		final int first=3*k,second=first+separation;
		noisy[first]=other(noisy[first]);noisy[second]=other(noisy[second]);
		final Counts counts=new Counts(k,truth);final Read read=read(reverse ? AminoAcid.reverseComplementBases(noisy) : noisy);
		final byte[] original=read.bases,originalQ=read.quality,before=original.clone(),beforeQ=originalQ.clone();
		final LocalEditEngine engine=new LocalEditEngine(k,counts,1);
		final int edits=engine.correct(read,cap);
		if(outcome==0){
			check(edits==2 && engine.substitutions==2 && engine.changedReads==1,"Exactly-max must commit only after a no-edit discovery pass.");
			check(Arrays.equals(read.bases,reverse ? AminoAcid.reverseComplementBases(truth) : truth),"Accepted two-edit output must equal independent truth.");
			check(engine.cappedReads==0 && engine.initialSkippedReads==0 && engine.rolledBackReads==0,"Exactly-max exhaustion is not over-limit.");
		}else{
			check(edits==0 && read.bases==original && read.quality==originalQ,"Rejected reads must retain exact original array identities.");
			check(engine.substitutions==0 && engine.insertions==0 && engine.deletions==0 && engine.changedReads==0,"Tentative/rejected edits must not enter committed statistics.");
			check(engine.initialSkippedReads==(outcome==1 ? 1 : 0) && engine.rolledBackReads==(outcome==2 ? 1 : 0),"Initial skip and rollback must be distinct outcomes.");
			check(engine.attemptedEdits==(outcome==1 ? 0 : 2),"Initial skip must precede any splice; rollback discovers the max+1 edit.");
		}
		check(Arrays.equals(original,before) && Arrays.equals(originalQ,beforeQ),"Original source buffers may never be mutated, including on success.");
		check(read.id.equals("guard") && read.numericID==7 && read.mate==null,"Correction transaction may not alter record identity or pair metadata.");
	}
	private static void unsupported(final int k){
		final Counts counts=new Counts(k,sequence(12*k,k));
		for(byte[] bases:new byte[][]{sequence(12*k,123456+k),new byte[k-1],new byte[6*k]}){
			if(bases.length!=12*k){Arrays.fill(bases,(byte)'N');}
			final Read read=read(bases);final byte[] original=read.bases,quality=read.quality;
			final LocalEditEngine engine=new LocalEditEngine(k,counts,1);
			check(engine.correct(read,1)==0 && engine.initialSkippedReads==0,"Low/undefined windows without supported repair proposals must not inflate initial burden.");
			check(read.bases==original && read.quality==quality,"Unsupported/short/N reads must remain untouched.");
		}
	}
	private static void indelRollback(final int k,final boolean extraBase){
		final byte[] truth=sequence(12*k,k);
		String noisy=new String(truth,StandardCharsets.US_ASCII);
		for(int p:new int[]{4*k+5,3*k}){
			noisy=extraBase ? noisy.substring(0,p)+noisy.charAt(p)+noisy.substring(p) : noisy.substring(0,p)+noisy.substring(p+1);
		}
		final Counts counts=new Counts(k,truth);final LocalEditEngine engine=new LocalEditEngine(k,counts,1);
		final Read rejected=read(noisy.getBytes(StandardCharsets.US_ASCII));final byte[] old=rejected.bases,oldQ=rejected.quality;
		check(engine.correct(rejected,1)==0 && engine.rolledBackReads==1,"Two nearby indels must exceed cap1 after initial underestimation.");
		check(rejected.bases==old && rejected.quality==oldQ,"Indel rollback restores original lengths, bases and qualities by identity.");
		check(engine.insertions==0 && engine.deletions==0 && engine.changedReads==0,"Rolled-back indels must not count as accepted repairs.");
		final Read accepted=read(noisy.getBytes(StandardCharsets.US_ASCII));
		check(engine.correct(accepted,2)==2 && Arrays.equals(accepted.bases,truth),"Same reusable engine must accept and exhaust the two-indel read at cap2.");
		check(engine.insertions==(extraBase ? 0 : 2) && engine.deletions==(extraBase ? 2 : 0) && engine.changedReads==1,"Only the second, committed transaction contributes indel/read counters.");
	}
	private static void exceptionRestoration(final int k){
		final byte[] truth=sequence(12*k,k),noisy=truth.clone();noisy[3*k]=other(noisy[3*k]);noisy[7*k]=other(noisy[7*k]);
		final Read read=read(noisy);final byte[] original=read.bases,quality=read.quality;
		final Counts counts=new Counts(k,truth);
		final LocalEditEngine engine=new LocalEditEngine(k,new LocalEditEngine.CountLookup(){
			@Override public int count(final Kmer key){if(read.bases!=original){throw new IntendedFailure();}return counts.count(key);}
		},1);
		boolean threw=false;try{engine.correct(read,2);}catch(IntendedFailure expected){threw=true;}
		check(threw && read.bases==original && read.quality==quality,"An exception after a tentative edit must restore original arrays and propagate.");
		check(engine.substitutions==0 && engine.changedReads==0,"Failed transactions cannot publish committed edits.");
	}
	private static void mixedContext(final int k,final boolean undefined){
		final byte[] truth=sequence(12*k,k),noisy=truth.clone();
		noisy[3*k]=other(noisy[3*k]);
		for(int i=7*k;i<9*k;i++){noisy[i]=undefined ? (byte)'N' : other(noisy[i]);}
		final byte[] expected=noisy.clone();expected[3*k]=truth[3*k];
		final LocalEditEngine engine=new LocalEditEngine(k,new Counts(k,truth),1);
		final Read read=read(noisy);
		check(engine.correct(read,1)==1 && Arrays.equals(read.bases,expected),"An unrelated unsupported or N-separated component must not block a valid single repair.");
		check(engine.initialSkippedReads==0 && engine.rolledBackReads==0 && engine.substitutions==1,"Unsupported/N components contribute zero to binding initial burden.");
	}
	private static void ordinaryAndPair(final int k){
		final byte[] truth=sequence(12*k,k),noisy=truth.clone();
		for(int p:new int[]{3*k,7*k,7*k+7}){noisy[p]=other(noisy[p]);}
		final Counts counts=new Counts(k,truth);
		for(int cap:new int[]{2,3}){
			final LocalEditEngine engine=new LocalEditEngine(k,counts,1);final Read read=read(noisy);
			final byte[] original=read.bases,quality=read.quality;final int edits=engine.correct(read,cap,true);
			if(cap==2){check(edits==0 && engine.rolledBackReads==1 && read.bases==original && read.quality==quality,"Ordinary plus pair-witness correction must not bypass the max+1 rollback.");}
			else{check(edits==3 && engine.substitutions==3 && engine.cappedReads==0 && Arrays.equals(read.bases,truth),"Ordinary plus pair witness may commit exactly max only after exhaustion.");}
		}
	}
	private static void pairedRejected(final int k){
		final byte[] truth=sequence(12*k,k);final Read read=read(truth),mate=read(truth);read.mate=mate;mate.mate=read;
		final LocalEditEngine engine=new LocalEditEngine(k,new Counts(k,truth),1);
		boolean threw=false;try{engine.correct(read,1);}catch(IllegalArgumentException expected){threw=true;}
		check(threw && read.mate==mate && mate.mate==read,"Unpaired-only precondition must precede any transaction, preserving links.");
	}
	private static byte[] sequence(final int length,final long seed){
		final Random random=new Random(seed);final byte[] out=new byte[length],alphabet={'A','C','G','T'};
		for(int i=0;i<length;i++){out[i]=alphabet[random.nextInt(4)];}return out;
	}
	private static byte other(final byte b){return b=='A' ? (byte)'C' : (byte)'A';}
	private static Read read(final byte[] bases){
		final byte[] quality=new byte[bases.length];for(int i=0;i<quality.length;i++){quality[i]=(byte)(i%41);}
		return new Read(bases.clone(),quality,"guard",7,false);
	}
	private static final class Counts implements LocalEditEngine.CountLookup {
		Counts(final int k,final byte[] truth){
			final Kmer key=new Kmer(k);for(byte b:truth){key.addRight(b);if(key.len()>=k){map.put(Arrays.toString(key.key()),12);}}
		}
		@Override public int count(final Kmer key){final Integer depth=map.get(Arrays.toString(key.key()));return depth==null ? 0 : depth;}
		final HashMap<String,Integer> map=new HashMap<String,Integer>();
	}
	private static final class IntendedFailure extends RuntimeException {private static final long serialVersionUID=1L;}
	private static void check(final boolean pass,final String why){checks++;if(!pass){throw new AssertionError(why);}}
	private static long checks;
}
