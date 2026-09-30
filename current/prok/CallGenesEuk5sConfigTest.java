package prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import stream.Read;
import structures.ByteBuilder;

/** Standalone registration, rollback, and seed-cutoff fixtures; run in a fresh JVM.
 * @author Raiden
 */
public class CallGenesEuk5sConfigTest {
	public static void main(String[] args) throws Exception{
		check(!CallGenes.EUK5S_ENABLED, "euk5S must be off until release calibration is accepted");
		check(!CallGenes.NCRNA_FAMILIES_ENABLED && !CallGenes.RRNA17_ENABLED, "Fixture requires fresh process defaults");
		final boolean old16=ProkObject.call16S, old23=ProkObject.call23S;
		final boolean old5=ProkObject.call5S, old18=ProkObject.call18S;
		check(CallGenes.parseEuk5sFlag("EuK5S", "t"), "Case-insensitive standalone flag not parsed");
		CallGenes.validateRrna17GateCombo();
		CallGenes.validateNcrnaSweepOverrides();
		CallGenes.loadEuk5sResources();
		check(GeneCaller.ncrnaFamilies.size()==1, "Standalone euk5S must register exactly one family");
		final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(0);
		check("euk5S".equals(family.name) && family.seedMinHits==2 && family.kLong==17
			&& family.outputType==ProkObject.r5S && family.library.length>0, "Incomplete standalone euk5S registration");
		check(!CallGenes.NCRNA_FAMILIES_ENABLED && !CallGenes.RRNA17_ENABLED && CallGenes.RRNA17_PROFILE==null
			&& ProkObject.call16S==old16 && ProkObject.call23S==old23
			&& ProkObject.call5S==old5 && ProkObject.call18S==old18, "Standalone flag changed another calling gate");
		//Remove tRNA's optional shared-index slot to isolate the family registration.
		ProkObject.calltRNA=false;
		GeneCaller.initializeConservedRnaSeedIndex();
		check(GeneCaller.conservedRnaSeedSlotCount()==1, "Standalone euk5S did not enter the shared index");
		reject(()->CallGenes.loadEuk5sResources(), "Duplicate euk5S");
		check(GeneCaller.ncrnaFamilies.size()==1 && GeneCaller.ncrnaFamilies.get(0)==family,
			"Duplicate rejection must preserve the first registration");
		for(String profile : new String[]{"euk", "euk5s"}){
			CallGenes.RRNA17_ENABLED=true; CallGenes.RRNA17_PROFILE=profile;
			reject(()->CallGenes.validateRrna17GateCombo(), "Duplicate euk5S");
		}
		CallGenes.RRNA17_ENABLED=false; CallGenes.RRNA17_PROFILE=null;
		CallGenes.NCRNA_FAMILY_FILTER="euk5S";
		CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=1;
		reject(()->CallGenes.validateNcrnaSweepOverrides(), "ncrnaseedminhits>=2");
		CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=2;
		CallGenes.validateNcrnaSweepOverrides();
		CallGenes.NCRNA_SEED_MIN_HITS_OVERRIDE=-1; CallGenes.NCRNA_FAMILY_FILTER=null;
		testSeedCutoff(family);
		testDpCoexistence();
		GeneCaller.ncrnaFamilies.clear();
		testResourceFailures();
		check(CallGenes.parseEuk5sFlag("euk5s", "f") && !CallGenes.EUK5S_ENABLED,
			"Explicit false must disable the independent gate");
		check(!CallGenes.parseEuk5sFlag("rrna17", "t"), "Standalone parser consumed another flag");
		System.out.println("PASS CallGenesEuk5sConfigTest");
	}

	private static void testSeedCutoff(NcrnaFamily family){
		assert(family!=null) : "Seed-cutoff fixture needs the real registered family controls";
		final byte[] bases=family.library[0].clone();
		final LongHashSet seeds=new LongHashSet(2);
		seeds.add(seed(bases, 0));
		final NcrnaScavenger s=new NcrnaScavenger(family.library, family.models, family.modelNames,
			seeds, family.kLong, family.minLen, family.windowPad, family.indexK, family.indexTopN,
			family.adaptive, family.adaptFloor, family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits,
			family.scoreA, family.scoreB, family.idPass, family.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(s, family);
		final int[] oneHit=new ConservedRnaSeedIndex(new LongHashSet[]{seeds}).scan(bases).hits(0);
		check(oneHit.length==1, "Fixture's first seed must occur exactly once");
		check(s.scavenge("fixture", bases, 0, new ArrayList<int[]>(), oneHit).isEmpty()
			&& s.alignmentCount()==0, "One-hit window must be rejected before alignment");
		seeds.add(seed(bases, bases.length-17));
		final int[] twoHits=new ConservedRnaSeedIndex(new LongHashSet[]{seeds}).scan(bases).hits(0);
		check(twoHits.length==2, "Fixture must contain exactly two seed occurrences");
		check(!s.scavenge("fixture", bases, 0, new ArrayList<int[]>(), twoHits).isEmpty()
			&& s.alignmentCount()>0, "Two-hit exact consensus window must reach alignment and be called");
	}
	private static long seed(byte[] bases, int start){
		assert(start>=0 && start+17<=bases.length) : "Seed fixture needs a complete 17-mer";
		long key=0;
		for(int i=start; i<start+17; i++){
			final int b=AminoAcid.baseToNumber[bases[i]];
			check(b>=0, "Seed fixture must contain unambiguous bases"); key=(key<<2)|b;
		}
		return key;
	}

	/** Inject only the candidate-generation stage; use the real callGenes merge and DP.
	 * Both labels must win when better scored, and separate calls must coexist. */
	private static void testDpCoexistence(){
		final byte[] bases=new byte[1000]; java.util.Arrays.fill(bases, (byte)'A');
		for(int strand=0; strand<2; strand++){
			for(boolean eukWins : new boolean[]{false, true}){
				final Orf legacy=candidate(bases, strand, 100, null, eukWins ? 100f : 200f);
				final Orf euk=candidate(bases, strand, 100, "euk5S", eukWins ? 200f : 100f);
				final CandidateCaller caller=new CandidateCaller(legacy, euk);
				final ArrayList<Orf> selected=caller.callGenes(new Read(bases, null, "fixture", 0), false);
				check(selected.size()==1 && selected.get(0)==(eukWins ? euk : legacy),
					"Shared DP must select by score, not family priority; strand="+strand+" eukWins="+eukWins);
				check(legacy.pathScore()>-999999 && euk.pathScore()>-999999,
					"Both candidates must be scored by DP, including the rejected overlap");
			}
			final Orf legacy=candidate(bases, strand, 100, null, 200f);
			final Orf euk=candidate(bases, strand, 500, "euk5S", 200f);
			final ArrayList<Orf> selected=new CandidateCaller(legacy, euk).callGenes(new Read(bases, null, "fixture", 0), false);
			check(selected.size()==2 && selected.contains(legacy) && selected.contains(euk),
				"Nonoverlapping legacy/euk5S candidates must both survive the DP; strand="+strand);
		}
	}
	private static Orf candidate(byte[] bases, int strand, int start, String family, float score){
		assert(start>=0 && start+120<=bases.length) : "DP fixture candidates must lie inside the contig";
		final Orf orf=new Orf("fixture", start, start+119, strand, 0, bases, false, ProkObject.r5S);
		orf.ncrnaFamily=family; orf.orfScore=score;
		return orf;
	}
	private static final class CandidateCaller extends GeneCaller{
		CandidateCaller(Orf legacy_, Orf euk_){
			super(1001, 0, 0, 0f, 0f, 0f, 0f, 0f, new GeneModel(false));
			legacy=legacy_; euk=euk_;
		}
		@Override ArrayList<Orf>[] makeRnas(String name, byte[] bases){
			assert(legacy.strand==euk.strand) : "DP fixture compares both families on the same strand";
			@SuppressWarnings("unchecked")
			final ArrayList<Orf>[] lists=new ArrayList[]{new ArrayList<Orf>(), new ArrayList<Orf>()};
			lists[legacy.strand].add(legacy); lists[euk.strand].add(euk);
			return lists;
		}
		final Orf legacy, euk;
	}

	private static void testResourceFailures() throws Exception{
		final Path dir=Files.createTempDirectory("euk5s-config-fixture-");
		final String missing=dir.resolve("missing.fa").toString();
		final String valid=write(dir, "valid.fa", ">model\nACGTTGCAAGTCGATCGTACGATGCTTAGGCATCGAT\n");
		final String shortSeed=write(dir, "short.fa", ">short\nACGTTGCA\n");
		final String empty=write(dir, "empty.fa", "");
		final String duplicate=write(dir, "duplicate.fa", ">same\nACGTTGCAAGTCGATCGTACGATGC\n>same\nTCGTTGCAAGTCGATCGTACGATGC\n");
		try{
			reject(()->CallGenes.addEuk5sFamily(0, missing, valid), "Missing euk5S");
			reject(()->CallGenes.addEuk5sFamily(0, valid, missing), "Missing euk5S");
			reject(()->CallGenes.addEuk5sFamily(0, valid, shortSeed), "Incompatible euk5S seed");
			reject(()->CallGenes.addEuk5sFamily(0, shortSeed, valid), "Incompatible euk5S consensus");
			reject(()->CallGenes.addEuk5sFamily(0, empty, valid), "Incompatible euk5S consensus");
			reject(()->CallGenes.addEuk5sFamily(0, duplicate, valid), "Incompatible euk5S consensus");
			check(GeneCaller.ncrnaFamilies.isEmpty(), "Rejected resources left a partial family registration");
			CallGenes.EUK5S_CONSENSUS_OVERRIDE=valid;
			CallGenes.loadEuk5sResources();
			check(GeneCaller.ncrnaFamilies.size()==1 && GeneCaller.ncrnaFamilies.get(0).library.length==1
				&& GeneCaller.ncrnaFamilies.get(0).modelNames[0].equals("model"), "Explicit consensus must replace only the euk5S library");
			CallGenes.EUK5S_ENABLED=false;
			reject(()->CallGenes.validateRrna17GateCombo(), "euk5sconsensus requires");
			CallGenes.EUK5S_ENABLED=true;
		}finally{
			CallGenes.EUK5S_CONSENSUS_OVERRIDE=null; GeneCaller.ncrnaFamilies.clear();
			for(String file : new String[]{valid, shortSeed, empty, duplicate}){Files.deleteIfExists(java.nio.file.Paths.get(file));}
			Files.delete(dir);
		}
	}

	private static String write(Path dir, String name, String text){
		assert(dir!=null && name!=null) : "Fixture output must stay in its generated temporary directory";
		final String path=dir.resolve(name).toString();
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(new ByteBuilder(text));
		if(writer.poisonAndWait()){throw new IllegalStateException("Could not write fixture "+path);}
		return path;
	}
	private static void reject(Runnable action, String message){
		try{action.run();}
		catch(IllegalArgumentException e){check(e.getMessage().contains(message), "Unexpected rejection: "+e); return;}
		throw new AssertionError("Expected rejection containing: "+message);
	}
	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
