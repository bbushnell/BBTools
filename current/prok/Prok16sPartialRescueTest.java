package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import map.LongHashSet;

/** Actual scavenger fixtures for physical-edge length exceptions and endpoint isolation.
 * @author Brian, Raiden */
public final class Prok16sPartialRescueTest {
	public static void main(String[] args){
		final byte[] model=random(1800);
		for(int mode=1; mode<=3; mode++){
			final int from=mode==2 ? 0 : 400, to=mode==1 ? 1800 : 1650;
			final byte[] part=Arrays.copyOfRange(model, from, to);
			final byte[] bases=embed(part, mode==2 ? 80 : 0, mode==1 ? 80 : 0);
			for(int strand=0; strand<2; strand++){
				final NcrnaScavenger off=caller(model, 0), on=caller(model, 1200);
				final Recorder recorder=observe(on, model);
				check(off.scavenge("edge", bases, strand, new ArrayList<int[]>()).isEmpty(), "Strict primary must reject this clipped fixture before rescue");
				final ArrayList<Orf> calls=on.scavenge("edge", bases, strand, new ArrayList<int[]>());
				check(calls.size()==1, "One genuine clipped locus must reach the caller through seed/window/length gates");
				final Orf call=calls.get(0);
				check(call.ncrnaPartial==mode && call.accepted16sIdentity==1f && recorder.count==0,
					"Clipped rescue must retain verified identity and bypass complete-gene endpoint features/nets");
				check(call.start==(mode==2 ? 80 : 0) && call.stop==bases.length-1-(mode==1 ? 80 : 0),
					"Partial rescue must preserve raw clipped alignment ends");
				if(strand==1){call.flip();}
				final int mask=strand==0 ? mode : ((mode&1)<<1)|((mode&2)>>>1);
				check(call.toGff().contains(";partial="+(mask==1 ? "10" : mask==2 ? "01" : "11")),
					"Existing formatter must report genomic left/right clipping on either strand");
			}
		}
		final byte[] shortPart=Arrays.copyOfRange(model, 350, 1650);
		final NcrnaScavenger shortOn=caller(model, 1200);
		check(shortOn.scavenge("short_contig", shortPart, 0, new ArrayList<int[]>()).size()==1,
			"A whole contig below ordinary 1400 must reach the explicit partial route");
		check(caller(model, 1200).scavenge("below_floor", Arrays.copyOfRange(model, 400, 1599), 0, new ArrayList<int[]>()).isEmpty(),
			"The explicit partial floor remains binding at physical edges");
		check(caller(shortPart, 1200).scavenge("full_short_model", shortPart, 0, new ArrayList<int[]>()).isEmpty(),
			"Touching both edges cannot bypass 1400 unless query/model clipping actually occurred");
		final byte[] full=embed(model, 80, 80);
		final NcrnaScavenger fullOff=caller(model, 0), fullOn=caller(model, 1200);
		final Recorder completeObserver=observe(fullOn, model);
		final ArrayList<Orf> a=fullOff.scavenge("complete", full, 0, new ArrayList<int[]>()),
			b=fullOn.scavenge("complete", full, 0, new ArrayList<int[]>());
		check(a.size()==1 && calls(a).equals(calls(b))
			&& fullOff.alignmentCount()==fullOn.alignmentCount() && completeObserver.count==1,
			"Accepted complete primary must preserve calls, endpoint dispatch, and alignment work");
		for(boolean gap:new boolean[]{false, true}){
			final byte[] internal=embed(shortPart, 4000, 4000);
			if(gap){Arrays.fill(internal, 0, 4000, (byte)'N');Arrays.fill(internal, 5300, internal.length, (byte)'N');}
			final NcrnaScavenger x=caller(model, 0), y=caller(model, 1200);
			check(x.scavenge("internal", internal, 0, new ArrayList<int[]>()).isEmpty()
				&& y.scavenge("internal", internal, 0, new ArrayList<int[]>()).isEmpty()
				&& x.alignmentCount()==y.alignmentCount(), "Internal fragments and N gaps must add neither rescues nor extra alignments");
		}
		final NcrnaScavenger claimed=caller(model, 1200);
		final ArrayList<int[]> claims=new ArrayList<int[]>();claims.add(new int[]{0, 9});
		check(claimed.scavenge("claimed_left", Arrays.copyOfRange(model,400,1800), 0, claims).isEmpty(),
			"An artificial left boundary created by claim subtraction must not permit prefix clipping");
		config();
		System.out.println("PASS Prok16sPartialRescueTest actual_scavenger physical_edges strands floor actual_query_clipping primary_parity NN_bypass internal_N_and_claim_guards");
	}
	static void config(){
		final boolean enabled=CallGenes.RRNA17_ENABLED;final String profile=CallGenes.RRNA17_PROFILE;
		try{
			CallGenes.setRrna17Profile("16s");
			final NcrnaFamily f=family(random(1800), 0);
			final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();families.add(f);
			final Prok16sWindowConfig config=new Prok16sWindowConfig();
			check(config.parse("rrna16spartialminlen", "1200"), "The public spelling must be consumed by the run-local parser");
			config.parse("rrna16swindow", "containment");config.apply(families);
			check(scavenger(f).prok16sPartialMinLen==1200, "The configured floor must reach the actual worker");
			for(String value:new String[]{"-1", "1401", "NaN"}){
				boolean rejected=false;try{new Prok16sWindowConfig().parse("rrna16spartialminlen", value);}catch(IllegalArgumentException e){rejected=true;}
				check(rejected, "Invalid floors must fail before work begins");
			}
		}finally{CallGenes.RRNA17_ENABLED=enabled;CallGenes.RRNA17_PROFILE=profile;}
	}
	static NcrnaScavenger caller(byte[] model, int floor){return scavenger(family(model, floor));}
	static NcrnaScavenger scavenger(NcrnaFamily f){
		final NcrnaScavenger result=new NcrnaScavenger(f.library,f.models,f.modelNames,f.kmerSet,
			f.kLong,f.minLen,f.windowPad,f.indexK,f.indexTopN,f.adaptive,f.adaptFloor,
			f.adaptTopFrac,f.adaptQFrac,f.fixedMinHits,f.scoreA,f.scoreB,f.idPass,f.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(result,f);return result;
	}
	static String calls(ArrayList<Orf> calls){
		final StringBuilder text=new StringBuilder();for(Orf call:calls){text.append(call.toGff()).append('\n');}return text.toString();
	}
	static final class Recorder implements RrnaEndpointCallerFeatures.Sink{
		@Override public void capture(String contig,int strand,int model,int start,int stop,float[] five,boolean fu,float[] three,boolean tu){count++;}
		int count;
	}
	static NcrnaFamily family(byte[] model, int floor){
		assert(model.length>=17) : "Seeded fixtures need at least one complete discovery word";
		final LongHashSet seeds=new LongHashSet();
		for(int i=0; i+17<=model.length; i++){
			long key=0;for(int j=0; j<17; j++){key=(key<<2)|AminoAcid.baseToNumber[model[i+j]];}seeds.add(key);
		}
		final NcrnaFamily f=new NcrnaFamily("16S",new byte[][]{model},null,new String[]{"bacteria"},seeds,17,1400,3109,
			9,1,false,1f,0f,0f,1,0f,1f,1f,1f);
		f.outputType=ProkObject.r16S;f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;f.scavengePass2=false;
		f.preserveContainingWindows=false;f.prok16sPartialMinLen=floor;return f;
	}
	static Recorder observe(NcrnaScavenger caller, byte[] model){
		final RrnaPositionalKmerTable.Table five=new RrnaPositionalKmerTable.Table("bacteria","5prime",6,0,.5),
			three=new RrnaPositionalKmerTable.Table("bacteria","3prime",6,-5,.5);
		final byte[] training=embed(model,40,40);
		five.add(training,40);three.add(training,40+model.length-1);five.finish();three.finish();
		final Recorder recorder=new Recorder();
		caller.setRrnaEndpointFeatures(new RrnaEndpointCallerFeatures.Resources("16S",new String[]{"bacteria"},new byte[][]{model},
			new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{three},null,null,null,true),recorder);
		return recorder;
	}
	static byte[] random(int length){return random(length,419);}
	static byte[] random(int length, int seed){
		assert(length>0) : "Fixture DNA must be nonempty";
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};final Random random=new Random(seed);
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}return bases;
	}
	static byte[] embed(byte[] fragment, int left, int right){
		assert(left>=0 && right>=0) : "Flanks must preserve physical-end geometry";
		final byte[] bases=random(left+fragment.length+right,9281);System.arraycopy(fragment,0,bases,left,fragment.length);return bases;
	}
	static void check(boolean value, String why){if(!value){throw new AssertionError(why);}}
}
