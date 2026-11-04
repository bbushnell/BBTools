package test.dna;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import dna.Exon;
import dna.Gene;
import dna.GeneSet;
import fileIO.ReadWrite;
import shared.Shared;
import structures.Range;

/** Fresh-JVM regression for active-build gene caches surviving Data.unloadAll.
 * Generates tiny Build1 and Build2 resources under work=, warms Build1, verifies
 * unloadAll dropped every public active-build cache that can be inspected before
 * calling a selected first Build2 accessor, then verifies all active Build2 values.
 * case=control proves the generated Build2 resources are loadable without reset.
 * The fixed gene_names_37.1.txt dictionary is checked separately.
 * This standalone test intentionally leaves Data globals dirty; each case runs
 * in a fresh JVM and relies on process exit to confine residue.
 * @author Jean
 */
public final class DataUnloadGeneCacheTest{

	/** Runs one fresh-JVM case with generated resources retained under work=. */
	public static void main(String[] args){
		String mode=null;
		File work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=new File(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		require(validMode(mode) && work!=null && work.isAbsolute(),
				"Require case=control|genes|sets|names|transcripts|ids|coding-range and absolute work=");
		require(checkPublicRangeSlots("fresh-start", false) & checkActiveIdSlot("fresh-start", false),
				"Fixture cases must start before the public active-build gene caches are populated");
		install(work);
		Data.setPath(work.getAbsolutePath()+File.separator);
		Data.GENE_MAP=GENE_MAP;
		final boolean ok=("control".equals(mode) ? control() : reset(mode));
		require(ok, "case="+mode+" must expose only Build2 after setup");
		System.out.println("PASS "+mode);
	}

	/** Verifies fresh Build2 resources independently of any unloadAll reset. */
	private static boolean control(){
		Data.GENOME_BUILD=2;
		return checkActiveBuild("fresh-build2", 2) & checkFixedDictionary();
	}

	/** Warms Build1, resets, checks exposed cleared slots, then uses the selected first accessor. */
	private static boolean reset(String mode){
		Data.GENOME_BUILD=1;
		require(checkActiveBuild("build1-warm", 1), "Build1 fixture must load before unloadAll");
		require(checkPublicRangeSlots("build1-warm", true) && checkActiveIdSlot("build1-warm", true),
				"Build1 warmup must populate the public slots used as unload sentinels");
		require(checkFixedDictionary(), "Fixed gene_names_37.1.txt dictionary must be readable");

		Data.unloadAll();
		Data.GENE_MAP=GENE_MAP;
		Data.GENOME_BUILD=2;

		final boolean unloaded=checkPublicRangeSlots("after-unload", false) &
				checkActiveIdSlot("after-unload", false);
		final boolean first=checkFirstAccessor("first-"+mode, mode, 2);
		final boolean remaining=checkActiveBuild("remaining-"+mode, 2);
		final boolean fixed=checkFixedDictionary();
		return unloaded & first & remaining & fixed;
	}

	/** Creates the fixed dictionary and one serialized one-exon gene per build/chromosome. */
	private static void install(File work){
		require(work.mkdir(), "Could not create fresh work directory "+work);
		final File geneRoot=new File(new File(work, "ref"), "genes");
		require(geneRoot.mkdirs(), "Could not create gene root "+geneRoot);
		ReadWrite.writeString("#source\tid\tname\nfixed\t3700\tFixedName\n",
				new File(geneRoot, "gene_names_37.1.txt").getPath(), false);
		for(int build=1; build<=2; build++){
			final File dir=new File(new File(geneRoot, "Build"+build), GENE_MAP);
			require(dir.mkdirs(), "Could not create "+dir);
			for(int chrom=1; chrom<=25; chrom++){
				final Gene[] genes={gene(build, chrom)};
				ReadWrite.writeAsync(genes, new File(dir, "chr"+chrom+".ga").getPath(), false);
			}
		}
	}

	/** Builds one valid, fully coding, one-exon transcript with build-distinct labels and IDs. */
	private static Gene gene(final int build, final int chrom){
		final int start=start(build, chrom), stop=stop(build, chrom);
		final Exon[] exons={new Exon(start, stop, (byte)chrom, Shared.PLUS, false, true)};
		return new Gene((byte)chrom, Shared.PLUS, start, stop, start, stop, id(build, chrom),
				symbol(build, chrom), transcript(build, chrom), null, null, null, exons,
				false, false, true, "seqGene", "fixture", "fixture");
	}

	/** Checks representative accessors across the active Build1/Build2 annotation caches. */
	private static boolean checkActiveBuild(String label, int build){
		boolean ok=true;
		ok&=checkGenes(label+"-genes", Data.getGenes(CHROM), build);
		ok&=checkGeneSets(label+"-sets", Data.getGeneSets(CHROM), build);
		ok&=checkGeneNameTable(label+"-names", Data.geneNameTable(CHROM), build);
		ok&=checkTranscriptNameTable(label+"-transcripts", Data.transcriptNameTable(CHROM), build);
		ok&=checkGeneIDTable(label+"-ids", Data.getGeneIDTable(), build);
		ok&=checkCodingRange(label+"-coding-range", Data.geneCodeAndExonRangeMatrix(CHROM), build);
		return ok;
	}

	/** Calls one family before any other post-reset getter can repopulate or mask it. */
	private static boolean checkFirstAccessor(String label, String mode, int build){
		if("genes".equals(mode)){return checkGenes(label, Data.getGenes(CHROM), build);}
		if("sets".equals(mode)){return checkGeneSets(label, Data.getGeneSets(CHROM), build);}
		if("names".equals(mode)){return checkGeneNameTable(label, Data.geneNameTable(CHROM), build);}
		if("transcripts".equals(mode)){return checkTranscriptNameTable(label, Data.transcriptNameTable(CHROM), build);}
		if("ids".equals(mode)){return checkGeneIDTable(label, Data.getGeneIDTable(), build);}
		if("coding-range".equals(mode)){
			return checkCodingRange(label, Data.geneCodeAndExonRangeMatrix(CHROM), build);
		}
		throw new AssertionError(mode);
	}

	/** Requires the single fixture gene from the expected build. */
	private static boolean checkGenes(String label, Gene[] genes, int build){
		final int expected=id(build, CHROM), stale=id(other(build), CHROM);
		final boolean ok=genes!=null && genes.length==1 && genes[0].id==expected &&
				symbol(build, CHROM).equals(genes[0].symbol) &&
				transcript(build, CHROM).equals(genes[0].mrnaAcc);
		System.out.println(label+" expected="+expected+" stale="+stale+
				" actual="+geneLabels(genes)+" pass="+ok);
		return ok;
	}

	/** Requires the single fixture gene set from the expected build. */
	private static boolean checkGeneSets(String label, GeneSet[] sets, int build){
		final int expected=id(build, CHROM), stale=id(other(build), CHROM);
		final boolean ok=sets!=null && sets.length==1 && sets[0].id==expected &&
				symbol(build, CHROM).equals(sets[0].name);
		System.out.println(label+" expected="+expected+" stale="+stale+
				" actual="+geneSetLabels(sets)+" pass="+ok);
		return ok;
	}

	/** Requires the symbol map to expose only the expected build on the checked chromosome. */
	private static boolean checkGeneNameTable(String label, HashMap<String, GeneSet> table, int build){
		final String expected=symbol(build, CHROM), stale=symbol(other(build), CHROM);
		final GeneSet found=table.get(expected);
		final boolean ok=found!=null && found.id==id(build, CHROM) && !table.containsKey(stale);
		System.out.println(label+" expected="+expected+" stale="+stale+
				" keys="+table.keySet()+" pass="+ok);
		return ok;
	}

	/** Requires the transcript map to expose only the expected build on the checked chromosome. */
	private static boolean checkTranscriptNameTable(String label, HashMap<String, Gene> table, int build){
		final String expected=transcript(build, CHROM), stale=transcript(other(build), CHROM);
		final Gene found=table.get(expected);
		final boolean ok=found!=null && found.id==id(build, CHROM) && !table.containsKey(stale);
		System.out.println(label+" expected="+expected+" stale="+stale+
				" keys="+table.keySet()+" pass="+ok);
		return ok;
	}

	/** Requires the active ID scan to have reached all 25 chromosomes from the expected build. */
	private static boolean checkGeneIDTable(String label, HashMap<Integer, ArrayList<GeneSet>> table, int build){
		final ArrayList<Integer> missing=new ArrayList<Integer>();
		final ArrayList<Integer> stale=new ArrayList<Integer>();
		for(int chrom=1; chrom<=25; chrom++){
			final int expectedId=id(build, chrom);
			final int staleId=id(other(build), chrom);
			if(!containsId(table, expectedId)){missing.add(Integer.valueOf(expectedId));}
			if(containsId(table, staleId)){stale.add(Integer.valueOf(staleId));}
		}
		final int size=(table==null ? -1 : table.size());
		final boolean ok=missing.isEmpty() && stale.isEmpty() && size==25;
		System.out.println(label+" expectedPresent="+(25-missing.size())+"/25 stalePresent="+
				stale.size()+"/25 size="+size+" missing="+missing+" stale="+stale+" pass="+ok);
		return ok;
	}

	/** Requires the coding-exon range to carry the expected build's coordinates. */
	private static boolean checkCodingRange(String label, Range[] ranges, int build){
		final int expectedStart=start(build, CHROM), expectedStop=stop(build, CHROM);
		final int staleStart=start(other(build), CHROM);
		final boolean ok=ranges!=null && ranges.length==1 && ranges[0].a==expectedStart &&
				ranges[0].b==expectedStop;
		System.out.println(label+" expected=["+expectedStart+","+expectedStop+"] staleStart="+
				staleStart+" actual="+rangeLabels(ranges)+" pass="+ok);
		return ok;
	}

	/** Confirms the fixed dictionary uses the build-independent gene_names_37.1.txt resource. */
	private static boolean checkFixedDictionary(){
		final Integer id=Data.geneNameToIdTable().get("FixedName");
		final String name=Data.geneIdToNameTable().get(Integer.valueOf(FIXED_ID));
		final boolean ok=Integer.valueOf(FIXED_ID).equals(id) && "FixedName".equals(name);
		System.out.println("fixed-dictionary id="+id+" name="+name+" pass="+ok);
		return ok;
	}

	/** Checks public range cache slots before and after unloadAll. */
	private static boolean checkPublicRangeSlots(String label, boolean loaded){
		boolean ok=true;
		ok&=checkRangeSlots(label+"-geneTxRangeMatrix", Data.geneTxRangeMatrix, loaded);
		ok&=checkRangeSlots(label+"-geneCodeRangeMatrix", Data.geneCodeRangeMatrix, loaded);
		ok&=checkRangeSlots(label+"-exonRangeMatrix", Data.exonRangeMatrix, loaded);
		ok&=checkRangeSlots(label+"-geneNearbyRangeMatrix", Data.geneNearbyRangeMatrix, loaded);
		ok&=checkRangeSlots(label+"-geneSetRangeMatrix", Data.geneSetRangeMatrix, loaded);
		return ok;
	}

	/** Checks whether all legacy chromosome range slots have the expected nullness. */
	private static boolean checkRangeSlots(String label, Range[][] matrix, boolean loaded){
		int loadedCount=0;
		for(int chrom=1; chrom<=25; chrom++){
			if(matrix[chrom]!=null){loadedCount++;}
		}
		final boolean ok=loaded ? loadedCount==25 : loadedCount==0;
		System.out.println(label+" expectedLoaded="+loaded+" loadedSlots="+loadedCount+
				"/25 pass="+ok);
		return ok;
	}

	/** Checks the active ID-table cache before and after unloadAll. */
	private static boolean checkActiveIdSlot(String label, boolean loaded){
		return checkSlot(label+"-geneIDTable", Data.geneIDTable, loaded);
	}

	/** Checks whether a sentinel slot has the expected nullness. */
	private static boolean checkSlot(String label, Object value, boolean loaded){
		final boolean ok=loaded ? value!=null : value==null;
		System.out.println(label+" expectedLoaded="+loaded+" actualLoaded="+(value!=null)+" pass="+ok);
		return ok;
	}

	/** Returns true iff the active ID table contains the expected singleton. */
	private static boolean containsId(HashMap<Integer, ArrayList<GeneSet>> table, int id){
		if(table==null){return false;}
		final ArrayList<GeneSet> list=table.get(Integer.valueOf(id));
		return list!=null && list.size()==1 && list.get(0).id==id;
	}

	/** Formats gene IDs compactly without using annotation formatting as an oracle. */
	private static String geneLabels(Gene[] genes){
		if(genes==null){return "null";}
		final StringBuilder sb=new StringBuilder();
		for(Gene gene : genes){
			if(sb.length()>0){sb.append(',');}
			sb.append(gene==null ? "null" : gene.id+"."+gene.symbol+"."+gene.mrnaAcc);
		}
		return sb.toString();
	}

	/** Formats gene-set IDs compactly without using annotation formatting as an oracle. */
	private static String geneSetLabels(GeneSet[] sets){
		if(sets==null){return "null";}
		final StringBuilder sb=new StringBuilder();
		for(GeneSet set : sets){
			if(sb.length()>0){sb.append(',');}
			sb.append(set==null ? "null" : set.id+"."+set.name);
		}
		return sb.toString();
	}

	/** Formats inclusive ranges compactly. */
	private static String rangeLabels(Range[] ranges){
		if(ranges==null){return "null";}
		final StringBuilder sb=new StringBuilder();
		for(Range range : ranges){
			if(sb.length()>0){sb.append(',');}
			sb.append(range==null ? "null" : "["+range.a+","+range.b+"]");
		}
		return sb.toString();
	}

	/** Whether the requested case name is supported. */
	private static boolean validMode(String mode){
		return "control".equals(mode) || "genes".equals(mode) || "sets".equals(mode) ||
				"names".equals(mode) || "transcripts".equals(mode) || "ids".equals(mode) ||
				"coding-range".equals(mode);
	}

	/** Keeps finite-oracle checks active under -da. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Build-specific gene ID. */
	private static int id(int build, int chrom){return build*1000+chrom;}
	/** Build-specific start coordinate. */
	private static int start(int build, int chrom){return 1000*build+10*chrom;}
	/** Build-specific stop coordinate. */
	private static int stop(int build, int chrom){return start(build, chrom)+8;}
	/** Build-specific gene symbol. */
	private static String symbol(int build, int chrom){return "B"+build+"_"+chrom;}
	/** Build-specific transcript accession. */
	private static String transcript(int build, int chrom){return "B"+build+"_TX_"+chrom;}
	/** The other generated build, used to detect stale retained caches. */
	private static int other(int build){return build==1 ? 2 : 1;}

	private static final int CHROM=1;
	private static final int FIXED_ID=3700;
	private static final String GENE_MAP="seqRefGene";
}
