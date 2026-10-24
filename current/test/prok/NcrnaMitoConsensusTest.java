package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import prok.NcrnaMitoConsensus;
import prok.NcrnaMitoConsensus.Data;
import prok.NcrnaMitoConsensus.Result;
import stream.Read;

/** Independent sequence expectations for copy weighting, indels and final trim.
 * @author Brian, Raiden */
public final class NcrnaMitoConsensusTest {
	public static void main(String[] args){
		final byte[] original=dna(240),changed=original.clone();changed[110]=different(changed[110]);
		final ArrayList<Read> votes=new ArrayList<Read>();add(votes,original,"a",3);add(votes,changed,"b",4);
		final Data data=new Data(votes);check(data.unique.size()==2,"Fixture must distinguish unique alignment count from seven copy votes");
		final Result serial=NcrnaMitoConsensus.refine(data,original,"weighted",20,0.5f,1,null);
		final Result parallel=NcrnaMitoConsensus.refine(data,original,"weighted",20,0.5f,3,null);
		check(Arrays.equals(serial.consensus,changed),"Four duplicate mutant copies must beat three reference copies");
		check(Arrays.equals(serial.consensus,parallel.consensus) && serial.graph.readTotal==7 && parallel.graph.readTotal==7,"Thread count must preserve sequence and every copy vote");
		check(Arrays.equals(serial.graph.original,serial.consensus),"HBM must use final FASTA coordinates");
		final byte[] inserted=new byte[original.length+3];System.arraycopy(original,0,inserted,0,120);System.arraycopy(new byte[]{'A','C','G'},0,inserted,120,3);System.arraycopy(original,120,inserted,123,120);
		final ArrayList<Read> indels=new ArrayList<Read>();add(indels,original,"a",3);add(indels,inserted,"b",5);
		final Result insertion=NcrnaMitoConsensus.refine(new Data(indels),original,"insertion",20,0.5f,2,null);
		check(Arrays.equals(insertion.consensus,inserted),"Graph refinement must retain the majority insertion omitted by base-count-only refinement");
		final ArrayList<Read> deletions=new ArrayList<Read>();add(deletions,original,"a",5);add(deletions,inserted,"b",3);
		final Result deletion=NcrnaMitoConsensus.refine(new Data(deletions),inserted,"deletion",20,0.5f,2,null);
		check(Arrays.equals(deletion.consensus,original),"Graph refinement must remove a majority-voted deletion and preserve long-member voting in final model frame");
		final byte[] extended=("ACGTCGATGACT"+new String(original,StandardCharsets.US_ASCII)+"TCAGCATGATCG").getBytes(StandardCharsets.US_ASCII);
		final ArrayList<Read> bounds=new ArrayList<Read>();add(bounds,original,"a",6);add(bounds,extended,"b",2);
		add(bounds,Arrays.copyOfRange(original,30,210),"partial",2);
		final Result endpoints=NcrnaMitoConsensus.refine(new Data(bounds),original,"endpoints",20,0.5f,2,null);
		check(endpoints.broad.length>endpoints.consensus.length && Arrays.equals(endpoints.consensus,original),"Low-depth real extensions must survive broad iteration and disappear only at final half-depth trim");
		check(endpoints.graph.readTotal==10,"Long and partial members must still vote after final trimming");
		check(endpoints.leftCovered==8 && endpoints.rightCovered==8,"Endpoint coverage counts actual original copies and excludes the two internal fragments, never reference pseudocounts");
		final Result eachRound=NcrnaMitoConsensus.refine(new Data(bounds),original,"each_round",20,0.5f,0.5f,2,null);
		check(Arrays.equals(eachRound.broad,original) && Arrays.equals(eachRound.consensus,original) && eachRound.graph.readTotal==10,"Explicit per-round half-depth trimming must suppress the minority extensions during iteration while retaining all ten votes");
		final ArrayList<Read> distant=new ArrayList<Read>();final byte[] longGene=dna(8000);add(distant,longGene,"long",2);
		final Result longResult=NcrnaMitoConsensus.refine(new Data(distant),longGene,"long_padding",5,0.5f,2,null);
		check(Arrays.equals(longGene,longResult.consensus) && longResult.graph.readTotal==2,"An8000bp seed after8000Ns must remain reachable; origin-anchored bands miss it");
		final ArrayList<Read> seeds=new ArrayList<Read>();seeds.add(new Read(original,null,"small n=3",0));seeds.add(new Read(changed,null,"largest n=4",1));
		check(NcrnaMitoConsensus.largestSeed(seeds)==seeds.get(1),"Seed choice follows largest cluster copy count");
		final Read trace=new Read("ACAN".getBytes(StandardCharsets.US_ASCII),null,"trace",0);trace.match="NIDmN".getBytes(StandardCharsets.US_ASCII);trace.start=0;trace.stop=3;
		final long counts=NcrnaMitoConsensus.referenceBaseCounts(trace,"NNAC".getBytes(StandardCharsets.US_ASCII));
		check((counts>>>32)==2 && (counts&0xffffffffL)==1,"Anchoring counts diagonal votes only: two non-N bases, one N base; insertion and deletion do not count as diagonal votes");
		int bad=0;
		try{NcrnaMitoConsensus.refine(data,original,"cap",1,0.5f,1,null);}catch(IllegalStateException e){bad++;}
		final ArrayList<Read> dup=new ArrayList<Read>(votes);dup.add(votes.get(0));try{new Data(dup);}catch(IllegalArgumentException e){bad++;}
		final ArrayList<Read> wrong=new ArrayList<Read>();add(wrong,original,"a",1);add(wrong,changed,"a",1);try{new Data(wrong);}catch(IllegalArgumentException e){bad++;}
		check(bad==3,"Nonconvergence, duplicate loci and conflicting unique IDs must fail loudly");
		final byte[] longOriginal=dna(2000),longChanged=longOriginal.clone();longChanged[1000]=different(longChanged[1000]);
		final ArrayList<Read> nearVotes=new ArrayList<Read>();add(nearVotes,longOriginal,"a",3);add(nearVotes,longChanged,"b",4);
		final Data nearData=new Data(nearVotes);
		final Result capped=NcrnaMitoConsensus.refine(nearData,longOriginal,"tolerance",1,0.5f,0.5f,2,null,true);
		check(!capped.exactStable && capped.toleranceStable && capped.exactTerminal.equals("MAX_ROUNDS") && Arrays.equals(capped.consensus,longChanged),"At the exact cap, one changed base in2000 may emit the last round under explicit tolerance without claiming exact stability");
		check(Arrays.equals(capped.consensus,capped.graph.original) && capped.graph.readTotal==7,"Tolerance FASTA/HBM coordinates and all seven original votes must agree");
		final Result continued=NcrnaMitoConsensus.refine(nearData,longOriginal,"continued",20,0.5f,0.5f,2,null,true);
		check(continued.exactStable && continued.toleranceStable && continued.rounds==2 && continued.exactTerminal.equals("EXACT"),"Tolerance on round1 must not stop the unchanged exact experiment before round2");
		check(NcrnaMitoConsensus.withinTolerance(longOriginal,longChanged,0.999f),"The declared0.999 threshold is inclusive");
		check(!NcrnaMitoConsensus.withinTolerance(original,changed,NcrnaMitoConsensus.consecutiveIdentity(original,changed)),"One substitution in240 bases is below0.999 despite unchanged length");
		check(!NcrnaMitoConsensus.withinTolerance(longOriginal,Arrays.copyOfRange(longOriginal,3,2000),1),"Three removed terminal bases fail the length criterion even when free-end identity is1");
		check(NcrnaMitoConsensus.withinTolerance(longOriginal,Arrays.copyOfRange(longOriginal,2,2000),1),"The declared two-base length tolerance is inclusive");
		System.out.println("MITO_CONSENSUS_TEST_PASS copy_weighting=true parallel_exact=true insertion=true deletion=true final_only_trim=true all_votes=true largest_seed=true long_padding=true endpoint_coverage=true roundtrim=true anchoring=true negatives=3");
		System.out.println("MITO_CONSENSUS_TOLERANCE_TEST_PASS terminal_only=true exact_preserved=true last_round=true all_votes=true identity_threshold=true length_threshold=true");
	}
	private static void add(ArrayList<Read> list,byte[] bases,String unique,int count){
		assert(count>0):"Fixture groups need real copy votes";
		for(int i=0;i<count;i++){final int n=list.size();list.add(new Read(bases,null,"locus"+n+" unique="+unique+" tid|1 partial=false pseudo=false",n));}
	}
	private static byte[] dna(int length){final Random r=new Random(8721);final byte[] b=new byte[length],alphabet={'A','C','G','T'};for(int i=0;i<length;i++){b[i]=alphabet[r.nextInt(4)];}return b;}
	private static byte different(byte b){return b=='A'?(byte)'C':(byte)'A';}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
