package assemble;

import java.util.Arrays;
import java.util.Random;
import java.io.File;
import fileIO.ByteStreamWriter;
import dna.AminoAcid;
import stream.Read;
import structures.IntList;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Controlled ideal-support truth tests for the integrated experimental caller.
 * Not a held-out library benchmark: support is explicitly forty clean copies.
 * @author Fischl */
public final class LocalEditCorrectorTest {
	public static void main(final String[] args){
		if(args.length==2 && args[0].equals("pairfixture")){pairFixture(new File(args[1]));return;}
		if(args.length==1){fixture(new File(args[0]));return;}
		if(args.length!=0){throw new IllegalArgumentException("Usage: [new_CLI_fixture_directory]");}
		final boolean packed=Kmer.PACKED,core=Kmer.MASK_CORE,quality=Read.CHANGE_QUALITY;
		try{
			Kmer.PACKED=true;Kmer.MASK_CORE=false;Read.CHANGE_QUALITY=false;
			minimalityOracle();
			for(int k:new int[]{31,62,63}){
				for(int seed=0;seed<4;seed++){
					for(int run:new int[]{1,4,10}){for(int kind=0;kind<3;kind++){repair(k,seed,run,kind);}}
				}
				controls(k);multiple(k);diagnostics(k);pairWitnesses(k);pairExclusions(k);
			}
		}finally{Kmer.PACKED=packed;Kmer.MASK_CORE=core;Read.CHANGE_QUALITY=quality;}
		System.out.println("LOCAL_EDIT_CORRECTOR_TEST_OK checks="+checks+" repairs="+repairs+"; isolated substitution/extra/missing-base truth including HP1/4/10, F/RC byte+quality, clean/3bp/ambiguous/selfRC controls, repeated-call coordinates. Ideal support only, no native CLI/performance claim.");
	}
	private static void repair(final int k,final int seed,final int run,final int kind){
		final String truth=sequence(seed,run);final int position=kind==0 ? 120+run/2 : 120;
		final String query=kind==0 ? truth.substring(0,position)+"T"+truth.substring(position+1) :
			kind==1 ? truth.substring(0,position)+(run==1 ? "T" : "A")+truth.substring(position) :
			truth.substring(0,position)+truth.substring(position+1);
		final Counts counts=new Counts(k);counts.add(truth,40);counts.add(query,1);
		final byte[] source=bytes(query),q=quality(source.length),expectedQ=repairQuality(q,position,kind);
		Read forward=null;
		for(boolean rc:new boolean[]{false,true}){
			final byte[] input=rc ? AminoAcid.reverseComplementBases(source) : source.clone(),quals=rc ? reverse(q) : q.clone();
			final byte[] snapshot=input.clone(),snapshotQ=quals.clone();
			final Read read=new Read(input,quals,"controlled",0,false);final LocalEditCorrector corrector=new LocalEditCorrector(k,counts);
			check(corrector.correctOne(read)==1,"Expected isolated repair K="+k+" seed="+seed+" run="+run+" kind="+kind+" rc="+rc);
			check(corrector.callStatus==LocalEditCorrector.CallStatus.APPLIED && corrector.acceptedTroughs>0 &&
				corrector.supportedCandidates>corrector.verificationRejectedCandidates,"Applied diagnostic must include a verified candidate from this call.");
			check(Arrays.equals(read.bases,rc ? AminoAcid.reverseComplementBases(bytes(truth)) : bytes(truth)),"Corrected sequence differs from explicit pre-error truth.");
			check(Arrays.equals(read.quality,rc ? reverse(expectedQ) : expectedQ),"Corrected qualities differ from independent original-position oracle.");
			check(Arrays.equals(input,snapshot) && Arrays.equals(quals,snapshotQ),"Integrated caller changed original input arrays.");
			if(!rc){forward=read;}else{
				check(Arrays.equals(forward.bases,AminoAcid.reverseComplementBases(read.bases)) && Arrays.equals(forward.quality,reverse(read.quality)),"Full selected correction must commute with RC.");
			}
			check(corrector.correctOne(read)==0,"Repaired supported read must not receive another correction.");
			check(corrector.callStatus==LocalEditCorrector.CallStatus.EXHAUSTED && corrector.lowRegions==0 &&
				corrector.supportedCandidates==0 && corrector.verificationRejectedCandidates==0,"Terminal zero after repair must reset prior successful-call diagnostics.");repairs++;
		}
	}
	private static void controls(final int k){
		final String truth=sequence(19,5);
		final Counts clean=new Counts(k);clean.add(truth,40);unchanged(k,truth,clean);
		for(int delta:new int[]{-3,3}){
			final Counts variants=new Counts(k);variants.add(sequence(19,5+delta),40);variants.add(truth,1);unchanged(k,truth,variants);
		}
		final Counts ambiguous=new Counts(k);ambiguous.add(truth,1);
		ambiguous.add(truth.substring(0,122)+"C"+truth.substring(123),40);
		ambiguous.add(truth.substring(0,122)+"G"+truth.substring(123),40);unchanged(k,truth,ambiguous);
		final String half=sequence(20,4);final String palindrome=half+new String(AminoAcid.reverseComplementBases(bytes(half)),java.nio.charset.StandardCharsets.US_ASCII);
		final Counts pal=new Counts(k);pal.add(palindrome,1);unchanged(k,palindrome,pal);
	}
	private static void multiple(final int k){
		final String truth=sequence(29,4);
		final int left=70,right=170;
		assert(left>=k && truth.length()-right>k && right-left>k) :
			"This two-error fixture specifically requires two disjoint low regions with complete supporting flanks at every tested K.";
		final char x=truth.charAt(left)=='A' ? 'C' : 'A',y=truth.charAt(right)=='A' ? 'C' : 'A';
		final String query=truth.substring(0,left)+x+truth.substring(left+1,right)+y+truth.substring(right+1);
		final Counts counts=new Counts(k);counts.add(truth,40);counts.add(query,1);
		final LocalEditCorrector corrector=new LocalEditCorrector(k,counts);
		final Read read=new Read(bytes(query),quality(query.length()),"two-errors",0,false);
		check(corrector.correctOne(read)==1 && corrector.correctOne(read)==1,"Separated errors need two calls with fresh coordinates/depths.");
		check(Arrays.equals(read.bases,bytes(truth)) && corrector.correctOne(read)==0,"Two successive corrections must restore exact truth and stop.");
	}
	/** Artificial count streams isolate control-flow accounting, not correction accuracy. */
	private static void diagnostics(final int k){
		final String truth=sequence(39,4),query=truth.substring(0,121)+"T"+truth.substring(122);
		final Counts counts=new Counts(k);counts.add(truth,40);counts.add(query,1);
		final LocalEditCorrector c=new LocalEditCorrector(k,counts);
		final String half=sequence(40,4);
		final String[] early={"ACG",truth.substring(0,90)+"X"+truth.substring(91),half+new String(AminoAcid.reverseComplementBases(bytes(half)),java.nio.charset.StandardCharsets.US_ASCII)};
		final LocalEditCorrector.CallStatus[] status={LocalEditCorrector.CallStatus.SHORT_READ,LocalEditCorrector.CallStatus.UNSUPPORTED_BASE,LocalEditCorrector.CallStatus.SELF_RC};
		for(int i=0;i<early.length;i++){
			check(c.correctOne(new Read(bytes(query),null,"seed",0,false))==1 && c.supportedCandidates>0,"Seed real successful-call state before checking early-return reset.");
			final Read read=new Read(bytes(early[i]),null,"early",0,false);final byte[] original=read.bases;
			check(c.correctOne(read)==0 && c.callStatus==status[i] && read.bases==original,"Early diagnostic must report its own reason without editing.");
			check(c.lowRegions==0 && c.acceptedTroughs==0 && c.skippedEdge==0 && c.skippedWide==0 && c.skippedUndefined==0 &&
				c.ambiguousLoci==0 && c.noSupportedCandidateLoci==0 && c.allCandidatesRejectedLoci==0 &&
				c.supportedCandidates==0 && c.verificationRejectedCandidates==0 && c.profileQueries==0 && c.probeQueries==0 && c.verificationQueries==0,
				"Early returns must clear all per-call diagnostics, not leak the previous successful read.");
		}
		final Read undefined=new Read(bytes(truth.substring(0,121)+"N"+truth.substring(122)),null,"undefined",0,false);
		final byte[] undefinedOriginal=undefined.bases;
		check(c.correctOne(undefined)==0 && undefined.bases==undefinedOriginal && c.callStatus==LocalEditCorrector.CallStatus.EXHAUSTED &&
			c.lowRegions==1 && c.skippedUndefined==1 && c.acceptedTroughs==0 && c.probeQueries==0,
			"An N-spanning low region must propagate the locator undefined skip, not a candidate failure.");
		final Counts ambiguity=new Counts(k);ambiguity.add(truth,1);
		ambiguity.add(truth.substring(0,121)+"C"+truth.substring(122),40);
		ambiguity.add(truth.substring(0,121)+"G"+truth.substring(122),40);
		final LocalEditCorrector ambiguous=new LocalEditCorrector(k,ambiguity);
		final Read unresolved=new Read(bytes(truth),null,"ambiguous",0,false);final byte[] unresolvedOriginal=unresolved.bases;
		check(ambiguous.correctOne(unresolved)==0 && unresolved.bases==unresolvedOriginal && ambiguous.callStatus==LocalEditCorrector.CallStatus.EXHAUSTED &&
			ambiguous.acceptedTroughs==1 && ambiguous.ambiguousLoci==1 && ambiguous.noSupportedCandidateLoci==0 && ambiguous.allCandidatesRejectedLoci==0,
			"Multiple supported alternatives must be counted as ambiguity, not absence of evidence or verification failure.");
		for(int kind=0;kind<4;kind++){
			final int mode=kind,windows=truth.length()-k+1;
			final HomopolymerIndelProposal.CountLookup artificial=new HomopolymerIndelProposal.CountLookup(){
				int calls;
				public int count(final Kmer key){
					assert(key!=null) : "Diagnostic stream substitutes count values, not null canonical keys.";
					final int j=calls++;
					if(j<windows){final int start=mode==2 ? 0 : 70,end=start+k+(mode==3 ? 1 : 0);return j>=start && j<end ? 1 : 40;}
					return mode==1 && j==windows ? 40 : 0;
				}
			};
			final LocalEditCorrector d=new LocalEditCorrector(k,artificial);
			final Read read=new Read(bytes(truth),null,"diagnostic",0,false);final byte[] original=read.bases;
			check(d.correctOne(read)==0 && read.bases==original && d.callStatus==LocalEditCorrector.CallStatus.EXHAUSTED && d.lowRegions==1,
				"Synthetic count stream must produce one exhausted low region without editing.");
			if(mode<2){
				check(d.acceptedTroughs==1 && d.noSupportedCandidateLoci==(mode==0 ? 1 : 0) && d.allCandidatesRejectedLoci==(mode==1 ? 1 : 0) &&
					d.supportedCandidates==mode && d.verificationRejectedCandidates==mode && d.ambiguousLoci==0,
					"No supported probe and one supported-but-rejected candidate are different terminal reasons.");
			}else{check(d.acceptedTroughs==0 && d.skippedEdge==(mode==2 ? 1 : 0) && d.skippedWide==(mode==3 ? 1 : 0) && d.probeQueries==0,
				"Edge/wide geometry skips must not be attributed to unsuccessful candidate probes.");}
		}
	}
	private static void unchanged(final int k,final String text,final Counts counts){
		for(boolean rc:new boolean[]{false,true}){
			final byte[] seq=rc ? AminoAcid.reverseComplementBases(bytes(text)) : bytes(text),q=quality(seq.length);
			final Read read=new Read(seq,q,"control",0,false);final LocalEditCorrector corrector=new LocalEditCorrector(k,counts);
			check(corrector.correctOne(read)==0 && read.bases==seq && read.quality==q,"Clean/variant/ambiguous control must remain unchanged, including array identities.");
			check(corrector.correctOne(read,true)==0 && read.bases==seq && read.quality==q,"Opt-in pair fallback must preserve clean/3bp/ambiguous/selfRC controls as well.");
		}
	}
	/** Controlled mixed pairs: successful witnesses must restore the complete truth. */
	private static void pairWitnesses(final int k){
		for(int spacing:new int[]{5,15,25}){
			for(int first=0;first<3;first++){for(int last=0;last<3;last++){
				final String truth=sequence(400+first*3+last,4);final int p=100,q=p+spacing;
				assertPair(k,truth,inject(inject(truth,q,last),p,first),"spacing="+spacing+" first="+first+" last="+last);
			}}
		}
		for(int run:new int[]{4,10}){for(int first:new int[]{1,2}){for(int last=0;last<3;last++){
			final String truth=sequence(600+run*10+first*3+last,run),right=inject(truth,120+run+6,last);
			final String query=first==1 ? right.substring(0,120)+"A"+right.substring(120) : right.substring(0,120)+right.substring(121);
			assertPair(k,truth,query,"pairedHP run="+run+" first="+first+" last="+last);
		}}}
		final String branchSource=sequence(731,4),branchTruth="A"+branchSource.substring(1,branchSource.length()-1)+"C";
		final String branchRight=inject(branchTruth,140,0),branchQuery=branchRight.substring(0,120)+"A"+branchRight.substring(120);
		final Counts branchCounts=new Counts(k);branchCounts.add(branchTruth,40);branchCounts.add(branchQuery,1);
		final LocalEditPairLookahead branch=new LocalEditPairLookahead(k,branchCounts);
		check(branch.propose(bytes(branchQuery),profile(k,bytes(branchQuery),branchCounts)) && branch.firstOperation==LocalSingleBaseEdit.Operation.DELETION &&
			branch.firstPosition==124 && Arrays.equals(branch.witness,bytes(branchTruth)),
			"Explicit canonical paired-HP fixture must use predecessor deletion at124, before first bad-window endpoint125.");
		final String truth=sequence(500,4);final Counts counts=new Counts(k);counts.add(truth,40);
		for(int delta:new int[]{-3,3}){
			final String variant=delta<0 ? truth.substring(0,100)+truth.substring(103) : truth.substring(0,100)+"CGT"+truth.substring(100);
			counts.add(variant,1);final LocalEditPairLookahead pair=new LocalEditPairLookahead(k,counts);
			check(!pair.propose(bytes(variant),profile(k,bytes(variant),counts)),"A three-base variant must not be accepted as a complete two-edit witness.");
		}
		for(int kind=0;kind<3;kind++){
			final String triple=inject(inject(inject(truth,110,kind),105,kind),100,kind);
			final Counts three=new Counts(k);three.add(truth,40);three.add(triple,1);
			check(!new LocalEditPairLookahead(k,three).propose(bytes(triple),profile(k,bytes(triple),three)),"Two-edit witness must not hide a third unresolved error.");
		}
	}
	private static void assertPair(final int k,final String truth,final String query,final String label){
		final Counts counts=new Counts(k);counts.add(truth,40);counts.add(query,1);
		Read forward=null;
		for(boolean rc:new boolean[]{false,true}){
			final byte[] input=rc ? AminoAcid.reverseComplementBases(bytes(query)) : bytes(query);
			final byte[] opposite=AminoAcid.reverseComplementBases(input);
			final boolean reverse=new String(input,java.nio.charset.StandardCharsets.US_ASCII).compareTo(new String(opposite,java.nio.charset.StandardCharsets.US_ASCII))>0;
			final byte[] canonical=reverse ? opposite : input,original=canonical.clone();
			final IntList depths=profile(k,canonical,counts);
			final int[] depthSnapshot=Arrays.copyOf(depths.array,depths.size);
			final LocalEditPairLookahead pair=new LocalEditPairLookahead(k,counts);
			final boolean found=pair.propose(canonical,depths);
			check(found,"Expected bounded mixed-pair witness K="+k+" "+label+" rc="+rc+
				" regions="+pair.regions+" candidates="+pair.candidates+" second="+pair.localSecondEdits+" verified="+pair.verifiedPairs+" queries="+pair.queries);
			final byte[] expected=rc ? AminoAcid.reverseComplementBases(bytes(truth)) : bytes(truth);
			check(Arrays.equals(reverse ? AminoAcid.reverseComplementBases(pair.witness) : pair.witness,expected),"Pair witness must restore full independent truth, not merely the probe word.");
			check(Arrays.equals(canonical,original) && pair.firstOperation!=null && pair.secondOperation!=null && pair.verifiedPairs>0,"Evidence-only lookahead must preserve original input and retain two verified operations.");
			check(Arrays.equals(depthSnapshot,Arrays.copyOf(depths.array,depths.size)),"Pair lookahead must not rewrite the caller original-depth profile.");
			check(pair.witnessContextEnd==pair.contextEnd+pair.witness.length-canonical.length && pair.contextStart>=0 && pair.witnessContextEnd<=pair.witness.length,
				"Original and two-edit witness context endpoints must differ by the total sequence-length delta.");
			check(!pair.propose(bytes(truth),profile(k,bytes(truth),counts)) && pair.witness==null && pair.firstOperation==null,"Clean subsequent call must reset previous pair witness.");
			final byte[] initialQ=rc ? reverse(quality(input.length)) : quality(input.length);
			final Read actual=new Read(input,initialQ,"pair-optin",0,false);
			final LocalEditCorrector integrated=new LocalEditCorrector(k,counts);
			check(distance(actual.bases,expected)==2,"Controlled pair must contain two global-distance errors, not merely two injected operations.");
			for(int step=0;step<2;step++){
				final byte[] before=actual.bases,beforeQ=actual.quality,snapshot=before.clone(),snapshotQ=beforeQ.clone();
				check(integrated.correctOne(actual,true)==1,"Opt-in must commit exactly one pair/ordinary edit per call: "+label);
				final int kind=integrated.lastOperation==LocalSingleBaseEdit.Operation.SUBSTITUTION ? 0 : integrated.lastOperation==LocalSingleBaseEdit.Operation.DELETION ? 1 : 2;
				check(Arrays.equals(actual.quality,repairQuality(beforeQ,integrated.lastPosition,kind)),"Opt-in quality changes must match an independent round-coordinate splice oracle.");
				check(Arrays.equals(before,snapshot) && Arrays.equals(beforeQ,snapshotQ),"Opt-in must not mutate prior-round arrays.");
				check(distance(actual.bases,expected)==1-step,"Each committed first/second edit must independently reduce truth distance by one.");
				if(integrated.usedPairLookahead){check(integrated.callStatus==LocalEditCorrector.CallStatus.APPLIED_PAIR && integrated.pairQueries>0,"A pair-first application must expose its source and lookup cost.");}
			}
			check(Arrays.equals(actual.bases,expected) && integrated.correctOne(actual,true)==0 && !integrated.usedPairLookahead && integrated.pairQueries==0,
				"Pair recovery must stop at full truth and reset pair diagnostics on the terminal call.");
			if(!rc){forward=actual;}else{check(Arrays.equals(forward.quality,reverse(actual.quality)),"Full pair correction qualities must commute with reverse complement.");}
			final Read implicit=new Read(input,initialQ,"default",0,false),explicit=new Read(input,initialQ,"explicit-off",0,false);
			final LocalEditCorrector off1=new LocalEditCorrector(k,counts),off2=new LocalEditCorrector(k,counts);
			check(off1.correctOne(implicit)==off2.correctOne(explicit,false) && Arrays.equals(implicit.bases,explicit.bases) && Arrays.equals(implicit.quality,explicit.quality) && off1.pairQueries==0 && off2.pairQueries==0,
				"Unspecified and explicit-off pair policies must agree and never perform lookahead queries.");
		}
	}
	/** Exhaustive short binary strings cover repeats, end skips and post-skip mismatches. */
	private static void minimalityOracle(){
		for(int na=0;na<=5;na++){for(int va=0;va<(1<<na);va++){
			final byte[] a=new byte[na];for(int i=0;i<na;i++){a[i]=(byte)(((va>>>i)&1)==0 ? 'A' : 'C');}
			for(int nb=0;nb<=5;nb++){for(int vb=0;vb<(1<<nb);vb++){
				final byte[] b=new byte[nb];for(int i=0;i<nb;i++){b[i]=(byte)(((vb>>>i)&1)==0 ? 'A' : 'C');}
				check(LocalEditPairLookahead.withinOneEdit(a,b)==(distance(a,b)<=1),"Minimality predicate must equal independent dynamic programming for lengths "+na+","+nb+" patterns "+va+","+vb);
			}}
		}}
	}
	/** Ordinary ambiguous locus must not block a disjoint pair or survive a call reset. */
	private static void pairExclusions(final int k){
		final String left=sequence(910,4),right=sequence(400,4),truth=left+right;
		final String query=left+inject(inject(right,115,0),100,0);
		final Counts counts=new Counts(k);counts.add(query,1);
		counts.add(truth.substring(0,121)+"C"+truth.substring(122),40);
		counts.add(truth.substring(0,121)+"G"+truth.substring(122),40);
		for(boolean rc:new boolean[]{false,true}){
			final byte[] raw=rc ? AminoAcid.reverseComplementBases(bytes(query)) : bytes(query);
			final byte[] expected=rc ? AminoAcid.reverseComplementBases(bytes(truth)) : bytes(truth),q=quality(raw.length);
			final Read read=new Read(raw,q,"ambiguous-and-pair",0,false);
			final LocalEditCorrector c=new LocalEditCorrector(k,counts);
			check(c.correctOne(read)==0 && c.ambiguousLoci==1,"Fixture must expose an ordinary ambiguous locus with an unresolved distant pair.");
			for(int step=0;step<2;step++){
				final byte[] old=read.bases,oldQ=read.quality,oldCopy=old.clone(),qCopy=oldQ.clone();
				check(c.correctOne(read,true)==1,"Per-locus exclusion must allow the disjoint recoverable pair K="+k+" rc="+rc);
				check(distance(read.bases,expected)==1-step,"Only the unrelated pair may improve; the ambiguous original allele must remain.");
				check(Arrays.equals(old,oldCopy) && Arrays.equals(oldQ,qCopy),"Exclusion handling must preserve original sequence and quality arrays.");
				check(Arrays.equals(read.quality,repairQuality(oldQ,c.lastPosition,0)),"Distant substitution repair must preserve independent quality splice semantics.");
			}
			final byte[] done=read.bases,doneQ=read.quality;
			check(Arrays.equals(done,expected) && c.correctOne(read,true)==0 && c.ambiguousLoci==1 && read.bases==done && read.quality==doneQ && c.pairQueries==0,
				"Ambiguous locus alone must abstain without probing it through pair lookahead.");
			final Read shortRead=new Read(bytes("AC"),null,"short-reset",0,false);
			check(c.correctOne(shortRead,true)==0 && c.ambiguousLoci==0 && c.pairQueries==0,"Early return must reset per-call exclusion diagnostics.");
			final Read again=new Read(raw,q,"pair-after-reset",0,false);
			check(c.correctOne(again,true)==1 && distance(again.bases,expected)==1,"Reused caller must recover the independent pair again after early reset.");
		}
		final Counts simple=new Counts(k);simple.add(right,40);
		final byte[] pair=bytes(inject(inject(right,115,0),100,0));simple.add(new String(pair,java.nio.charset.StandardCharsets.US_ASCII),1);
		final IntList depths=profile(k,pair,simple),excluded=new IntList();
		final LocalEditPairLookahead helper=new LocalEditPairLookahead(k,simple);
		check(helper.propose(pair,depths),"Exclusion reset fixture must initially have a pair witness.");
		excluded.add(helper.contextStart+1);
		check(!helper.propose(pair,depths,excluded) && helper.witness==null && helper.queries==0,"Explicit excluded region must be skipped without stale witness or lookups.");
		check(helper.propose(pair,depths),"A later call without exclusions must not inherit the prior call exclusions.");
		final String single=inject(right,100,1);final Counts oneCounts=new Counts(k);oneCounts.add(right,40);oneCounts.add(single,1);
		for(boolean rc:new boolean[]{false,true}){
			final byte[] one=rc ? AminoAcid.reverseComplementBases(bytes(single)) : bytes(single);
			final LocalEditPairLookahead detour=new LocalEditPairLookahead(k,oneCounts);
			check(!detour.propose(one,profile(k,one,oneCounts)) && detour.localSecondEdits>0 && detour.verifiedPairs==0,
				"A one-insertion query must reject locally completed SUB+DEL detours at the actual witness gate K="+k+" rc="+rc);
		}
	}
	/** Small independent unit-cost edit-distance oracle for controlled short fixtures. */
	private static int distance(final byte[] a,final byte[] b){
		assert(a!=null && b!=null) : "Truth-distance test oracle requires complete called sequences.";
		int[] previous=new int[b.length+1],current=new int[b.length+1];
		for(int j=0;j<=b.length;j++){previous[j]=j;}
		for(int i=1;i<=a.length;i++){
			current[0]=i;
			for(int j=1;j<=b.length;j++){current[j]=Math.min(previous[j]+1,Math.min(current[j-1]+1,previous[j-1]+(a[i-1]==b[j-1] ? 0 : 1)));}
			final int[] swap=previous;previous=current;current=swap;
		}
		return previous[b.length];
	}
	/** 0 substitution, 1 extra base, 2 missing base; applied right-to-left in tests. */
	private static String inject(final String input,final int p,final int kind){
		assert(p>0 && p<input.length()) : "Controlled error injection needs an internal original position.";
		char alternate='A';for(char b:new char[]{'A','C','G','T'}){if(b!=input.charAt(p) && b!=input.charAt(p-1)){alternate=b;break;}}
		if(kind==0){return input.substring(0,p)+alternate+input.substring(p+1);}
		if(kind==1){return input.substring(0,p)+alternate+input.substring(p);}
		return input.substring(0,p)+input.substring(p+1);
	}
	private static IntList profile(final int k,final byte[] bases,final Counts counts){
		final Kmer key=new Kmer(k);final IntList profile=new IntList();
		for(int i=0;i<bases.length;i++){key.addRight(bases[i]);if(i>=k-1){profile.add(counts.count(key));}}
		check(profile.size==bases.length-k+1,"Pair fixture requires one depth per original kmer start.");return profile;
	}
	private static String sequence(final int seed,final int run){
		assert(run>0) : "Controlled truth run must remain nonempty.";
		final Random random=new Random(2026091010L+seed);final StringBuilder b=new StringBuilder();
		for(int i=0;i<120;i++){b.append("ACGT".charAt(random.nextInt(4)));}b.setCharAt(119,'C');
		for(int i=0;i<run;i++){b.append('A');}b.append('G');
		for(int i=1;i<120;i++){b.append("ACGT".charAt(random.nextInt(4)));}return b.toString();
	}
	private static byte[] repairQuality(final byte[] q,final int p,final int kind){
		assert(p>=0 && p<q.length) : "Expected repair must address an original input position/boundary.";
		final int delta=kind==0 ? 0 : kind==1 ? -1 : 1;final byte[] out=new byte[q.length+delta];
		for(int i=0;i<out.length;i++){
			if(kind==0){out[i]=i==p ? 0 : q[i];}
			else if(kind==1){out[i]=q[i+(i>=p ? 1 : 0)];}
			else{out[i]=i==p ? 0 : q[i-(i>p ? 1 : 0)];}
		}return out;
	}
	private static byte[] quality(final int n){final byte[] q=new byte[n],values={0,1,60,93};for(int i=0;i<n;i++){q[i]=values[i%4];}return q;}
	private static byte[] bytes(final String s){return s.getBytes(java.nio.charset.StandardCharsets.US_ASCII);}
	private static byte[] reverse(final byte[] b){final byte[] out=b.clone();for(int i=0;i<b.length;i++){out[i]=b[b.length-1-i];}return out;}
	/** Files for actual CLI tests; all expected outputs are built without the caller. */
	private static void pairFixture(final File dir){
		if(dir.exists() || !dir.mkdirs()){throw new IllegalArgumentException("Pair fixture directory must be new: "+dir);}
		final ByteStreamWriter input=writer(dir,"queries.fq"),counts=writer(dir,"counts.fa"),expected=writer(dir,"expected.fa"),supportsOnly=writer(dir,"support-only.fa");
		int records=0;
		for(int kind=0;kind<13;kind++){
			final int run=kind==10 ? 10 : 4;
			final String truth=sequence(kind<9 ? 400+kind : 900+kind,run);String query=truth,wanted=truth;
			if(kind<9){query=inject(inject(truth,115,kind%3),100,kind/3);}
			else if(kind==9 || kind==10){
				final String right=inject(truth,120+run+6,0);
				query=kind==9 ? right.substring(0,120)+"A"+right.substring(120) : right.substring(0,120)+right.substring(121);
			}else if(kind==12){query=truth.substring(0,100)+"CGT"+truth.substring(100);wanted=query;}
			for(int copy=0;copy<40;copy++){
				fasta(counts,"support_"+kind+"_"+copy,bytes(truth));
				fasta(supportsOnly,"support_"+kind+"_"+copy,bytes(truth));
			}
			final byte[] q=quality(query.length());
			for(boolean rc:new boolean[]{false,true}){
				final String name="pair"+kind+"_"+(rc ? "RC" : "F");
				final byte[] bases=rc ? AminoAcid.reverseComplementBases(bytes(query)) : bytes(query);
				fastq(input,name,bases,rc ? reverse(q) : q);fasta(counts,name,bases);
				fasta(expected,name,rc ? AminoAcid.reverseComplementBases(bytes(wanted)) : bytes(wanted));records++;
			}
		}
		for(ByteStreamWriter w:new ByteStreamWriter[]{input,counts,expected,supportsOnly}){if(w.poisonAndWait()){throw new IllegalStateException("Pair CLI fixture write failed.");}}
		check(records==26,"Pair CLI panel has22 two-error queries and4 clean/3bp preservation controls.");
	}
	/** Files for actual CLI tests; all expected outputs are built without the caller. */
	private static void fixture(final File dir){
		if(dir.exists() || !dir.mkdirs()){throw new IllegalArgumentException("Fixture directory must be new: "+dir);}
		final ByteStreamWriter input=writer(dir,"queries.fq"),counts=writer(dir,"counts.fa"),expected=writer(dir,"expected.fq"),cap1=writer(dir,"expected-cap1.fq"),supportsOnly=writer(dir,"support-only.fa");
		int records=0,supports=0;
		for(int kind=0;kind<9;kind++){
			final int run=kind<3 ? 1 : 4;final String truth=sequence(13000+kind,run);
			String query=truth;
			if(kind==0){query=truth.substring(0,120)+"T"+truth.substring(121);}
			else if(kind==1 || kind==3){query=truth.substring(0,120)+(kind==1 ? "T" : "A")+truth.substring(120);}
			else if(kind==2 || kind==4){query=truth.substring(0,120)+truth.substring(121);}
			else if(kind==8){
				final char a=truth.charAt(70)=='A' ? 'C' : 'A',b=truth.charAt(170)=='A' ? 'C' : 'A';
				query=truth.substring(0,70)+a+truth.substring(71,170)+b+truth.substring(171);
			}
			final byte[] raw=bytes(query),q=quality(raw.length);
			byte[] outQ=q.clone(),capQ=q.clone();String capped=truth;
			if(kind<=4){outQ=repairQuality(q,120,kind==0 ? 0 : kind==1 || kind==3 ? 1 : 2);capQ=outQ.clone();}
			if(kind==8){
				outQ[70]=outQ[170]=0;
				final String reverseQuery=new String(AminoAcid.reverseComplementBases(raw),java.nio.charset.StandardCharsets.US_ASCII);
				final int selected=query.compareTo(reverseQuery)<0 ? 70 : 170;
				capped=query.substring(0,selected)+truth.charAt(selected)+query.substring(selected+1);capQ[selected]=0;
			}
			for(int copy=0;copy<40;copy++){
				if(kind==7){
					fasta(counts,"majorC_"+copy,bytes(truth.substring(0,122)+"C"+truth.substring(123)));
					fasta(counts,"majorG_"+copy,bytes(truth.substring(0,122)+"G"+truth.substring(123)));supports+=2;
					fasta(supportsOnly,"majorC_"+copy,bytes(truth.substring(0,122)+"C"+truth.substring(123)));
					fasta(supportsOnly,"majorG_"+copy,bytes(truth.substring(0,122)+"G"+truth.substring(123)));
				}else{
					fasta(counts,"major_"+kind+"_"+copy,bytes(kind==6 ? sequence(13000+kind,7) : truth));
					fasta(supportsOnly,"major_"+kind+"_"+copy,bytes(kind==6 ? sequence(13000+kind,7) : truth));supports++;
				}
			}
			for(boolean rc:new boolean[]{false,true}){
				final String name="case"+kind+"_"+(rc ? "RC" : "F");
				final byte[] oriented=rc ? AminoAcid.reverseComplementBases(raw) : raw;
				fastq(input,name,oriented,rc ? reverse(q) : q);
				fastq(expected,name,rc ? AminoAcid.reverseComplementBases(bytes(truth)) : bytes(truth),rc ? reverse(outQ) : outQ);
				fastq(cap1,name,rc ? AminoAcid.reverseComplementBases(bytes(capped)) : bytes(capped),rc ? reverse(capQ) : capQ);
				fasta(counts,name,oriented);supports++;records++;
			}
		}
		for(ByteStreamWriter w:new ByteStreamWriter[]{input,counts,expected,cap1,supportsOnly}){if(w.poisonAndWait()){throw new IllegalStateException("CLI fixture write failed.");}}
		check(records==18 && supports==418,"Pinned CLI panel must contain eighteen queries and418 count records.");
		System.out.println("LOCAL_EDIT_CLI_FIXTURE_OK queries=18 counts=418 changed=12 edits=14; maxedits1 gives12edits; no caller-derived output oracle.");
	}
	private static ByteStreamWriter writer(final File dir,final String name){final ByteStreamWriter w=new ByteStreamWriter(new File(dir,name).getPath(),false,false,false);w.start();return w;}
	private static void fasta(final ByteStreamWriter w,final String name,final byte[] bases){w.print('>').println(name).println(bases);}
	private static void fastq(final ByteStreamWriter w,final String name,final byte[] bases,final byte[] q){
		assert(bases.length==q.length) : "Independent CLI FASTQ oracle must retain one quality per base.";
		final byte[] encoded=new byte[q.length];for(int i=0;i<q.length;i++){encoded[i]=(byte)(q[i]+33);}
		w.print('@').println(name).println(bases).println("+").println(encoded);
	}
	private static final class Counts implements HomopolymerIndelProposal.CountLookup{
		Counts(final int k_){k=k_;final Kmer key=new Kmer(k);table=new HashArrayU1D(new int[]{2003},key.k,k);}
		void add(final String s,final int copies){
			assert(copies>0) : "Controlled count evidence uses explicit positive copy counts.";
			final Kmer key=new Kmer(k);for(int c=0;c<copies;c++){key.clearFast();for(int i=0;i<s.length();i++){key.addRight(s.charAt(i));if(key.len()>=k){table.increment(key);}}}
		}
		@Override public int count(final Kmer key){return table.getValue(key);}
		final int k;final HashArrayU1D table;
	}
	private static void check(final boolean ok,final String message){checks++;if(!ok){throw new AssertionError(message);}}
	private static long checks,repairs;
}
