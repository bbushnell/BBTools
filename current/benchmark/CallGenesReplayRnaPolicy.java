package benchmark;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.CompareGff;
import parse.LineParser1;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyRegistry;
import prok.ProkObject;

/** Research-only RNA startup and strict CDS views; changes no production defaults.
 * Each replay invocation owns a fresh JVM and calls the ordinary startup once.
 * Grading is serial, after all biological calling has finished.
 * Global flag mutation and count collection are not concurrent operations.
 * @author Sayu, Fischl
 */
public final class CallGenesReplayRnaPolicy{

	/*--------------------------------------------------------------*/
	/*----------------        Startup policy        ----------------*/
	/*--------------------------------------------------------------*/

	/** RNA startup choice; NORMAL delegates effective family/default selection to CallGenes. */
	public enum Mode{OFF, NORMAL}
	/** Fixed supplied PGM versus ordinary taxonomy-selected per-genome models. */
	public enum ModelPolicy{FIXED, PER_PHYLUM}
	/** Parse exactly fixed or perphylum; null and other spellings are rejected. */
	public static ModelPolicy parseModelPolicy(String value){
		if("fixed".equals(value)){return ModelPolicy.FIXED;}
		if("perphylum".equals(value)){return ModelPolicy.PER_PHYLUM;}
		throw new IllegalArgumentException("pgmpolicy must be fixed or perphylum");
	}

	/** Parse exactly off or normal; null and other spellings are rejected. */
	public static Mode parseMode(String value){
		if("off".equals(value)){return Mode.OFF;}
		if("normal".equals(value)){return Mode.NORMAL;}
		throw new IllegalArgumentException("rnacalling must be off or normal");
	}

	/** Build fixed-model control arguments. OFF preserves the old argument array.
	 * NORMAL keeps CDS and lets ordinary startup resolve its actual RNA defaults,
	 * including generic families. */
	public static String[] controlArgs(String fasta, String output, String pgm, Mode mode, Path fmnRoot){
		return controlArgs(fasta, output, pgm, mode, fmnRoot, ModelPolicy.FIXED);
	}

	/** Per-phylum mode leaves the ordinary fallback PGM default unresolved here.
	 * The benchmark separately binds the actual configured path before calling.
	 * Returns a new single-thread, one-pass argument array; this method neither starts
	 * CallGenes nor opens the named input, output, model or resource root. */
	public static String[] controlArgs(String fasta, String output, String pgm, Mode mode, Path fmnRoot, ModelPolicy policy){
		require(mode!=null, "RNA mode is required");
		require(policy!=null, "PGM policy is required");
		validateRoot(mode, fmnRoot);
		if(policy==ModelPolicy.PER_PHYLUM){
			require(mode==Mode.NORMAL, "Per-phylum comparison requires RNA ON");
			return new String[]{"in="+fasta, "out="+output, "taxonomy=t", "percontig=f", "passes=1",
				"callcds=t", "annotation=off", "ordered=t", "t=1", "overwrite=f", "fmnroot="+fmnRoot};
		}
		final String[] base=new String[]{"in="+fasta, "out="+output, "pgm="+pgm,
			"taxonomy=f", "passes=1", mode==Mode.OFF ? "cdsonly=t" : "callcds=t",
			"annotation=off", "ordered=t", "t=1", "overwrite=f"};
		if(mode==Mode.OFF){return base;}
		final String[] normal=Arrays.copyOf(base, base.length+1);
		normal[base.length]="fmnroot="+fmnRoot;
		return normal;
	}

	/** Require no root in OFF mode and a normalized absolute root in NORMAL mode.
	 * This is a path-shape check, not a resource existence or content check. */
	public static void validateRoot(Mode mode, Path fmnRoot){
		require(mode!=null, "RNA mode is required");
		if(mode==Mode.OFF){
			require(fmnRoot==null, "fmnroot is only valid with rnacalling=normal");
		}else{
			require(fmnRoot!=null && fmnRoot.isAbsolute() && fmnRoot.normalize().equals(fmnRoot),
				"rnacalling=normal requires a pinned normalized absolute fmnroot");
		}
	}

	/** Write actual startup state as a new UTF-8 table after ordinary initialization.
	 * The caller supplies the already-validated mode and keeps startup globals stable.
	 * Effective legacy flags can be false when a registered family replaces them. */
	public static void writeRuntime(Path path, Mode mode) throws IOException{
		require(ProkObject.callCDS, "RNA comparison must preserve CDS calling");
		require(!Files.exists(path), "RNA runtime receipt must be new");
		final StringBuilder text=new StringBuilder("kind\tname\tvalue\n");
		text.append("mode\trnacalling\t").append(mode==Mode.OFF ? "off" : "normal").append('\n');
		final boolean[] flags=flags();
		final String[] names={"CDS", "tRNA", "16S", "23S", "5S", "18S"};
		for(int i=0; i<names.length; i++){
			text.append("legacy_flag\t").append(names[i]).append('\t').append(flags[i]).append('\n');
		}
		text.append("generic_umbrella\tncrna\t").append(CallGenes.NCRNA_FAMILIES_ENABLED).append('\n');
		text.append("registered_families\tcount\t").append(GeneCaller.ncrnaFamilies.size()).append('\n');
		for(NcrnaFamily family:GeneCaller.ncrnaFamilies){
			require(family.name!=null && family.name.indexOf('\t')<0 && family.name.indexOf('\n')<0,
				"Invalid registered family name");
			text.append("registered_family\t").append(family.name).append('\t').append(family.outputType).append('\n');
			text.append("family_display_name\t").append(family.name).append('\t')
				.append(NcrnaFamilyRegistry.displayName(family.name)).append('\n');
		}
		Files.write(path, text.toString().getBytes(StandardCharsets.UTF_8));
	}

	/*--------------------------------------------------------------*/
	/*----------------    CDS views and grading     ----------------*/
	/*--------------------------------------------------------------*/

	/** Copy original CDS rows without reserializing coordinates/attributes. Non-CDS
	 * rows, including unknown generic RNA types, never enter a grading view.
	 * Preserve comments except embedded FASTA, which is rejected. Count all feature
	 * types and emitted family markers before filtering. Own and close the reader/writer;
	 * a failure can leave a partial output file.
	 * Omit blank lines: the retained GffLine reader dereferences line[0]. */
	public static Counts copyCds(Path input, Path output){
		require(!input.equals(output) && !Files.exists(output), "CDS view must be a fresh distinct file");
		final ByteFile reader=ByteFile.makeByteFile1(input.toString(), false);
		final ByteStreamWriter writer=new ByteStreamWriter(output.toString(), false, false, false);
		writer.start();
		final LineParser1 parser=new LineParser1('\t');
		final Counts counts=new Counts();
		for(NcrnaFamily family:GeneCaller.ncrnaFamilies){counts.registerFamily(family.name);}
		boolean failed=false;
		try{
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='#'){
					require(!new String(line, StandardCharsets.UTF_8).startsWith("##FASTA"),
						"Embedded FASTA is not a GFF feature view");
					writer.println(line);
					continue;
				}
				parser.set(line);
				require(parser.terms()==9, "Require nine-column GFF feature rows");
				require(parser.length(2)>0, "Empty GFF feature type");
				int category=-1;
				for(int i=0; i<Counts.TYPE_NAMES.length; i++){
					if(parser.termEquals(Counts.TYPE_NAMES[i], 2)){
						category=i;
						break;
					}
				}
				counts.total++;
				if(category>=0){
					counts.typeCounts[category]++;
				}else{
					Counts.increment(counts.otherTypes, parser.parseString(2));
				}
				if(category==0){
					counts.cds++;
					writer.println(line);
				}else{
					parser.setBounds(8);
					final int marker=familyMarker(line, parser.a(), parser.b());
					if(marker>=0){
						final int start=marker+8;
						int end=start;
						while(end<parser.b() && line[end]!=',' && line[end]!=';'){end++;}
						require(end>start, "Empty emitted RNA family");
						final String family=new String(line, start, end-start, StandardCharsets.UTF_8);
						Counts.increment(counts.families, family);
					}
				}
			}
		}catch(RuntimeException | Error problem){
			failed=true;
			throw problem;
		}finally{
			final boolean error=reader.close() | writer.poisonAndWait();
			if(error && !failed){throw new IllegalStateException("CDS view input/output failed");}
		}
		return counts;
	}

	/** Only call after all caller work. Restore the six flags even if grading fails.
	 * Inputs must be the explicit CDS views, not merely flag-filtered full GFFs. */
	public static void gradeFilteredCds(Path query, Path reference, Path report){
		final boolean[] prior=flags();
		try{
			ProkObject.callCDS=true;
			ProkObject.calltRNA=ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;
			CompareGff.main(new String[]{"in="+query, "ref="+reference, "outstream="+report});
		}finally{restore(prior);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Feature counts        ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable, thread-confined feature counts collected before CDS filtering.
	 * Known types have fixed row order; other types and display-family names are sorted. */
	public static final class Counts{
		/** All feature rows and the CDS subset; comments and blank lines are excluded. */
		public long total, cds;
		/** Fixed type order; index zero is the CDS category used by copyCds. */
		private static final String[] TYPE_NAMES={"CDS", "tRNA", "rRNA", "RNA"};
		/** Counts parallel to TYPE_NAMES, including zero-count categories. */
		private final long[] typeCounts=new long[TYPE_NAMES.length];
		/** Owned one-long counters, sorted by unknown type and emitted display-family name. */
		private final TreeMap<String, long[]> otherTypes=new TreeMap<>(), families=new TreeMap<>();
		/** Include a registry display name with a zero count unless already present. */
		public void registerFamily(String internalName){
			final String display=NcrnaFamilyRegistry.displayName(internalName);
			if(!families.containsKey(display)){families.put(display, new long[1]);}
		}
		/** Return the count for a display name, or zero when absent. */
		public long familyCount(String display){
			long[] count=families.get(display);
			return count==null ? 0 : count[0];
		}
		/** Return the count for an exact type name, or zero when absent. */
		public long typeCount(String type){
			for(int i=0; i<TYPE_NAMES.length; i++){
				if(TYPE_NAMES[i].equals(type)){return typeCounts[i];}
			}
			long[] count=otherTypes.get(type);
			return count==null ? 0 : count[0];
		}
		/** Increment the owned mutable counter, creating it on first observation. */
		private static void increment(TreeMap<String, long[]> map, String key){
			long[] count=map.get(key);
			if(count==null){
				count=new long[1];
				map.put(key, count);
			}
			count[0]++;
		}
		/** Return a fresh TSV fragment with totals, types and emitted-family counts.
		 * The caller supplies TSV-safe genome and arm labels; no header is added. */
		public String rows(String genome, String arm){
			final String prefix=genome+'\t'+arm+'\t';
			final StringBuilder out=new StringBuilder();
			out.append(prefix).append("total\tall\t").append(total).append('\n');
			out.append(prefix).append("total\tCDS\t").append(cds).append('\n');
			out.append(prefix).append("total\tnon_CDS\t").append(total-cds).append('\n');
			for(int i=0; i<TYPE_NAMES.length; i++){
				out.append(prefix).append("type\t").append(TYPE_NAMES[i]).append('\t').append(typeCounts[i]).append('\n');
			}
			for(Map.Entry<String, long[]> e:otherTypes.entrySet()){
				out.append(prefix).append("type\t").append(e.getKey()).append('\t').append(e.getValue()[0]).append('\n');
			}
			for(Map.Entry<String, long[]> e:families.entrySet()){
				out.append(prefix).append("emitted_family\t").append(e.getKey()).append('\t').append(e.getValue()[0])
					.append('\n');
			}
			return out.toString();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Find the emitted comma-family marker inside [start, end), or return -1.
	 * This recognizes the caller marker, not a general GFF attribute grammar. */
	private static int familyMarker(byte[] line, int start, int end){
		for(int i=start; i+8<=end; i++){
			if(line[i]==',' && line[i+1]=='f' && line[i+2]=='a' && line[i+3]=='m' && line[i+4]=='i'
				&& line[i+5]=='l' && line[i+6]=='y' && line[i+7]==':'){return i;}
		}
		return -1;
	}

	/** Snapshot six legacy flags in CDS, tRNA, 16S, 23S, 5S, 18S order. */
	private static boolean[] flags(){
		return new boolean[]{ProkObject.callCDS, ProkObject.calltRNA, ProkObject.call16S, ProkObject.call23S,
			ProkObject.call5S, ProkObject.call18S};
	}
	/** Restore a six-element snapshot produced by flags(), in the same order. */
	private static void restore(boolean[] f){
		ProkObject.callCDS=f[0];
		ProkObject.calltRNA=f[1];
		ProkObject.call16S=f[2];
		ProkObject.call23S=f[3];
		ProkObject.call5S=f[4];
		ProkObject.call18S=f[5];
	}
	/** Reject invalid policy or input independently of assertion enablement. */
	private static void require(boolean value, String reason){
		if(!value){throw new IllegalArgumentException(reason);}
	}
	/** Static startup/grading helpers; no instance state. */
	private CallGenesReplayRnaPolicy(){}
}
