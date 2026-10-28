package test.prok;

import java.util.ArrayList;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import parse.LineParser1;
import prok.CallGenes;
import prok.ConservedRnaSeedIndex;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.RrnaResourceIO;
import prok.TrnaKmerIndex;
import structures.ByteBuilder;

/** Exercises default-off loading and family-scoped controls on a real sealed bundle.
 * @author Brian Bushnell, Raiden
 */
public final class CallGenesCobalaminConfigTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=4){throw new IllegalArgumentException("Supply group directory, three comma-separated model counts, expected seed keys and a fresh fixture directory");}
		final String[] counts=args[1].split(",");check(counts.length==3, "All three B12 model counts must be explicit");
		final int[] models={Integer.parseInt(counts[0]), Integer.parseInt(counts[1]), Integer.parseInt(counts[2])};
		final int seeds=Integer.parseInt(args[2]);
		check(!CallGenes.COBALAMIN_ENABLED, "Measured B12 CPU overhead exceeds the default-on budget; selection must remain opt-in");
		final ArrayList<NcrnaFamily> previous=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		final boolean umbrella=CallGenes.NCRNA_FAMILIES_ENABLED;
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.parseCobalaminFlag("cobalaminroot", args[0]+"/absent");
			CallGenes.loadCobalaminResources();
			check(GeneCaller.ncrnaFamilies.isEmpty(), "Disabled cobalamin must neither open missing resources nor register a slot");
			CallGenes.parseCobalaminFlag("cobalamin", "t");
			boolean failed=false;
			try{CallGenes.loadCobalaminResources();}catch(IllegalArgumentException expected){failed=true;}
			check(failed && GeneCaller.ncrnaFamilies.isEmpty(), "Missing selected resources must fail before publishing a family");
			CallGenes.parseCobalaminFlag("cobalaminroot", args[0]);
			CallGenes.loadCobalaminResources();
			check(GeneCaller.ncrnaFamilies.size()==3, "One selected B12 group must register all three shared-index slots");
			final NcrnaFamily[] baseline=GeneCaller.ncrnaFamilies.toArray(new NcrnaFamily[3]);
			final int[] floors={109, 109, 66};
			final float[] identities={.64f, .66f, .66f};
			for(int i=0; i<baseline.length; i++){
				final NcrnaFamily f=baseline[i];
				check(f.name.equals("cobalamin_"+MEMBERS[i].toLowerCase(java.util.Locale.ROOT)) && f.kLong==17
					&& f.library.length==models[i] && f.models.length==models[i] && f.kmerSet.size()==seeds,
					"Group loading must conserve each subtype's model roster and shared seed denominator");
				check(f.indexK==9 && f.indexTopN==100 && f.minLen==floors[i] && f.hbmPass>1
					&& f.idPass==identities[i] && f.idBorderline==f.idPass,
					"Bare opt-in loading must reproduce the measured scalar identities, K9/top100 and independent length floors without HBM rescue");
				check(f.seedMinHits==2 && f.seedDistinct, "The selected seed resource requires two different keys before alignment");
				check(f.boundary5NetTemplate==null && f.boundary3NetTemplate==null,
					"Detection calibration precedes separately measured boundary networks");
			}
			checkSharedSeeds(baseline);
			final String saved=baseline[2].modelNames[0];
			try{
				baseline[2].modelNames[0]="cobalamin_RF01482_wrong_member";
				failed=false;try{CallGenes.requireCobalaminGroup(baseline);}catch(IllegalArgumentException expected){failed=true;}
				check(failed, "A wrong-directory model prefix must fail before subtype attribution can be corrupted");
			}finally{baseline[2].modelNames[0]=saved;}
			failed=false;try{CallGenes.loadCobalaminResources();}catch(IllegalStateException expected){failed=true;}
			check(failed && GeneCaller.ncrnaFamilies.size()==3, "Duplicate group loading must preserve the existing three callers");
			GeneCaller.ncrnaFamilies.clear();
			final Path fixture=Paths.get(args[3]);missingThirdConsensus(Paths.get(args[0]), fixture);
			CallGenes.parseCobalaminFlag("cobalaminroot", fixture.toString());
			failed=false;try{CallGenes.loadCobalaminResources();}catch(IllegalArgumentException expected){failed=true;}
			check(failed && GeneCaller.ncrnaFamilies.isEmpty(), "A missing third member must not publish the first two callers");
			CallGenes.parseCobalaminFlag("cobalaminroot", args[0]);
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILY_FILTER=CallGenes.parseNcrnaFamily("cobalamin");
			CallGenes.NCRNA_INDEX_K_OVERRIDE=9;CallGenes.NCRNA_INDEX_TOP_N_OVERRIDE=3;CallGenes.NCRNA_INDEX_MIN_HITS_OVERRIDE=3;
			CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=2;
			CallGenes.NCRNA_ID_PASS_OVERRIDE=.70f;CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=.60f;CallGenes.NCRNA_HBM_PASS_OVERRIDE=.75f;
			CallGenes.parseCobalaminFlag("cobalaminminlen", "120");
			CallGenes.validateNcrnaSweepOverrides();CallGenes.loadCobalaminResources();
			for(NcrnaFamily f:GeneCaller.ncrnaFamilies){
				check(f.indexK==9 && f.indexTopN==3 && f.fixedMinHits==3 && f.minLen==120
					&& f.idPass==.70f && f.idBorderline==.60f && f.hbmPass==.75f && f.seedMinHits==2 && f.seedDistinct,
					"The cobalamin group selector must apply sweep controls to all three members");
			}
			final NcrnaFamily changed=GeneCaller.ncrnaFamilies.get(0);
			final TrnaKmerIndex index=new TrnaKmerIndex(changed.library, changed.indexK, changed.adaptive,
				changed.adaptFloor, changed.adaptTopFrac, changed.adaptQFrac, changed.fixedMinHits);
			final int[] shortlist=index.shortlist(changed.library[0], changed.indexTopN);
			check(shortlist.length>0 && shortlist.length<=3 && index.lastSharedCount(shortlist[0])>=3,
				"The loaded index settings must produce a bounded real-model shortlist");
			checkSeedGate(changed, false);checkSeedGate(changed, true);
			GeneCaller.ncrnaFamilies.clear();CallGenes.COBALAMIN_MIN_LEN=-1;
			CallGenes.NCRNA_FAMILY_FILTER=CallGenes.parseNcrnaFamily("cobalamin_rf01482");
			CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=1;
			CallGenes.validateNcrnaSweepOverrides();CallGenes.loadCobalaminResources();
			for(int i=0; i<3; i++){
				final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(i);
				check(f.indexK==9 && f.indexTopN==(i==1?3:100) && f.minLen==floors[i]
					&& f.idPass==(i==1?.70f:identities[i]) && f.seedMinHits==(i==1?1:2) && f.seedDistinct==(i!=1),
					"An explicit subtype sweep must leave both sibling policies unchanged");
			}
			CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=0;failed=false;
			try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "An invalid seed threshold must fail before loading or calling");CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=2;
			CallGenes.NCRNA_MAX_LEN_OVERRIDE=300;failed=false;
			try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "B12 seed-gate support must not silently enable unrelated euk5S-only controls");CallGenes.NCRNA_MAX_LEN_OVERRIDE=-1;
			CallGenes.COBALAMIN_ENABLED=false;
			failed=false;try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "A cobalamin sweep cannot silently run with its family disabled");
			failed=false;try{CallGenes.parseCobalaminFlag("cobalaminminlen", "0");}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Nonpositive minimum length must fail before resource access");
		}finally{
			GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.addAll(previous);
			CallGenes.COBALAMIN_ENABLED=false;CallGenes.COBALAMIN_RESOURCE_ROOT=null;CallGenes.COBALAMIN_MIN_LEN=-1;
			CallGenes.NCRNA_FAMILIES_ENABLED=umbrella;CallGenes.NCRNA_FAMILY_FILTER=null;
			CallGenes.NCRNA_INDEX_K_OVERRIDE=CallGenes.NCRNA_INDEX_TOP_N_OVERRIDE=CallGenes.NCRNA_INDEX_MIN_HITS_OVERRIDE=-1;
			CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=-1;CallGenes.NCRNA_MAX_LEN_OVERRIDE=-1;
			CallGenes.NCRNA_ID_PASS_OVERRIDE=CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=CallGenes.NCRNA_HBM_PASS_OVERRIDE=Float.NaN;
		}
		System.out.println("CALLGENES_COBALAMIN_CONFIG_TEST_PASS");
	}
	/** Repeated occurrences of one real seed must fail; two different real keys pass the gate. */
	private static void checkSeedGate(NcrnaFamily f, boolean distinct){
		final long[] keys=f.kmerSet.toArray();check(keys.length>=2, "The two-key fixture requires two real library seeds");
		final byte[] bases=new byte[300];java.util.Arrays.fill(bases, (byte)'N');
		putSeed(bases, 100, keys[0]);putSeed(bases, 150, keys[distinct?1:0]);
		final prok.NcrnaScavenger s=new prok.NcrnaScavenger(f.library, f.models, f.modelNames, f.kmerSet, f.kLong, f.minLen, f.windowPad,
			f.indexK, f.indexTopN, f.adaptive, f.adaptFloor, f.adaptTopFrac, f.adaptQFrac, f.fixedMinHits, f.scoreA, f.scoreB, f.idPass, f.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(s, f);final int[] passed={0};
		s.setWorkloadSink(new prok.NcrnaWorkloadInstrumentSink(){
			@Override public void seedHits(String name, int strand, int[] positions){}
			@Override public void scheduledWindow(String name, int strand, int pass, int start, int stop){}
			@Override public void postSeedWindow(String name, int strand, int pass, int start, int stop, int count){
				check(count==2, "Both inputs contain exactly two occurrences; the effective distinct gate decides which proceeds");passed[0]++;
			}
		});
		s.scavenge("B12_seed_gate", bases, 0, new ArrayList<int[]>());
		check(s.kmerHitCount()==2 && (distinct?passed[0]>0:passed[0]==0 && s.alignmentCount()==0),
			"Loaded B12 controls must reject a repeated key before alignment and admit two independent keys");
	}
	private static void putSeed(byte[] bases, int start, long key){
		final byte[] alphabet={'A','C','G','T'};for(int i=0; i<17; i++){bases[start+i]=alphabet[(int)((key>>>(2*(16-i)))&3)];}
	}
	private static void checkSharedSeeds(final NcrnaFamily[] families){
		final long key=families[0].kmerSet.toArray()[0];final byte[] bases=new byte[17], alphabet={'A','C','G','T'};
		for(int i=0; i<bases.length; i++){bases[i]=alphabet[(int)((key>>>(2*(16-i)))&3)];}
		final ConservedRnaSeedIndex index=new ConservedRnaSeedIndex(new LongHashSet[]{families[0].kmerSet, families[1].kmerSet, families[2].kmerSet});
		final ConservedRnaSeedIndex.ScanResult scan=index.scan(bases);
		for(int i=0; i<3; i++){check(scan.hits(i).length==1 && scan.hits(i)[0]==8,
			"One forward17-mer scan must fan a shared B12 seed out to all three caller slots");}
	}
	/** Copy only the declared small resource payloads; break the final lookup without touching originals. */
	private static void missingThirdConsensus(Path source, Path fixture)throws Exception{
		Files.createDirectory(fixture);Files.copy(source.resolve("seeds.fa.gz"), fixture.resolve("seeds.fa.gz"));
		for(String member:MEMBERS){
			Files.createDirectory(fixture.resolve(member));
			for(String file:new String[]{"consensus.fa", "models.hbmc.gz"}){
				Files.copy(source.resolve(member).resolve(file), fixture.resolve(member).resolve(file));
			}
			final Path cutoffs=source.resolve(member).resolve("cutoffs.tsv");
			if(Files.exists(cutoffs)){Files.copy(cutoffs, fixture.resolve(member).resolve("cutoffs.tsv"));}
		}
		final ByteFile in=RrnaResourceIO.open(source.resolve("manifest.tsv").toString());
		final ByteStreamWriter out=RrnaResourceIO.writer(fixture.resolve("manifest.tsv").toString());
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();int changed=0;
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);
			if(p.terms()==8 && p.termEquals("RF01689", 0) && p.termEquals("consensus", 1)){
				b.clear();for(int i=0; i<8; i++){if(i>0){b.tab();}b.append(i==4?"absent.fa":p.parseString(i));}out.print(b.nl());changed++;
			}else{out.println(row);}
		}}finally{check(!in.close(), "Fixture manifest read failed");RrnaResourceIO.finish(out);}
		check(changed==1, "The negative fixture must alter exactly the third member's consensus lookup");
	}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
