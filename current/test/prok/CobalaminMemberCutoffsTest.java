package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CallGenes;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Optional member identities preserve defaults, explicit CLI precedence and
 * atomic group loading when a declared resource is broken.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminMemberCutoffsTest {

	public static void main(String[] args)throws Exception{
		check(args.length==2, "Supply the real B12 resource root and a fresh fixture directory");
		final Path source=Paths.get(args[0]), root=Paths.get(args[1]);Files.createDirectory(root);
		final ArrayList<NcrnaFamily> previous=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		final boolean enabled=CallGenes.COBALAMIN_ENABLED, umbrella=CallGenes.NCRNA_FAMILIES_ENABLED;
		final String oldRoot=CallGenes.COBALAMIN_RESOURCE_ROOT, filter=CallGenes.NCRNA_FAMILY_FILTER;
		final float oldPass=CallGenes.NCRNA_ID_PASS_OVERRIDE, oldBorder=CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=true;CallGenes.COBALAMIN_ENABLED=true;
			CallGenes.NCRNA_FAMILY_FILTER=null;
			CallGenes.NCRNA_ID_PASS_OVERRIDE=CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=Float.NaN;
			load(source);identities(.64f, .64f, .66f, .66f, .66f, .66f);
			modelProfiles(source, root);
			final String[] values={row(0, ".64", ".64"), row(1, ".68", ".68"), row(2, ".70", ".70")};
			final Path valid=group(source, root.resolve("valid"), values, -1);
			load(valid);identities(.64f, .64f, .68f, .68f, .70f, .70f);
			CallGenes.NCRNA_FAMILY_FILTER="cobalamin_rf00174";CallGenes.NCRNA_ID_PASS_OVERRIDE=.64f;
			load(valid);identities(.64f, .64f, .68f, .68f, .70f, .70f);
			CallGenes.NCRNA_FAMILY_FILTER="cobalamin_rf01689";CallGenes.NCRNA_ID_PASS_OVERRIDE=Float.NaN;
			CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=.68f;
			load(valid);identities(.64f, .64f, .68f, .68f, .70f, .68f);
			CallGenes.NCRNA_FAMILY_FILTER=null;CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=Float.NaN;
			// The generic model path must retain its scalar policy, not switch aligners.
			for(NcrnaFamily family:GeneCaller.ncrnaFamilies){
				check(!family.reuseConsensusAlignment && family.modelThresholds==null,
					"Member scalar files must not enable the euk5S single-alignment model policy");
			}
			CallGenes.NCRNA_FAMILY_FILTER="cobalamin_rf01482";
			CallGenes.NCRNA_ID_PASS_OVERRIDE=.72f;CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=.71f;
			load(valid);identities(.64f, .64f, .72f, .71f, .70f, .70f);
			CallGenes.NCRNA_FAMILY_FILTER="cobalamin";
			load(valid);identities(.72f, .71f, .72f, .71f, .72f, .71f);
			CallGenes.NCRNA_ID_PASS_OVERRIDE=CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=Float.NaN;
			CallGenes.NCRNA_FAMILY_FILTER=null;
			load(group(source, root.resolve("optional"), new String[]{null, values[1], null}, -1));
			identities(.64f, .64f, .68f, .68f, .66f, .66f);
			fails(group(source, root.resolve("missing"), values, 2), "A declared missing third-member cutoff must fail atomically");
			fails(group(source, root.resolve("wrong_family"), new String[]{values[0], values[1], values[1]}, -1),
				"A cutoff cannot silently bind another member's identity");
			fails(group(source, root.resolve("inverted"), new String[]{values[0], values[1], row(2, ".64", ".70")}, -1),
				"Borderline above pass must fail before any group is published");
			fails(group(source, root.resolve("duplicate"), new String[]{values[0], values[1], values[2]+values[2]}, -1),
				"A duplicate row must not silently select its first or last cutoff");
			fails(group(source, root.resolve("nonfinite"), new String[]{values[0], values[1], row(2, "NaN", ".70")}, -1),
				"Nonfinite member identities cannot enter alignment comparisons");
			CallGenes.NCRNA_FAMILY_FILTER="cobalamin_rf01482";CallGenes.NCRNA_ID_PASS_OVERRIDE=.60f;
			fails(valid, "Partial CLI overrides must be checked against the effective packaged borderline value");
			CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=.70f;
			boolean failed=false;
			try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "An explicitly inverted B12 pair must still fail before resource loading");
			CallGenes.NCRNA_FAMILY_FILTER="rnasep";CallGenes.NCRNA_ID_PASS_OVERRIDE=.50f;CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=Float.NaN;
			failed=false;try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Other families retain their early registry-default identity validation");
			CallGenes.COBALAMIN_ENABLED=false;
			CallGenes.COBALAMIN_RESOURCE_ROOT=root.resolve("does_not_exist").toString();
			GeneCaller.ncrnaFamilies.clear();CallGenes.loadCobalaminResources();
			check(GeneCaller.ncrnaFamilies.isEmpty(), "Disabled B12 must not open optional member files");
			failed=false;
			try{new FamilyResourceResolver().optional("cobalamin_rf00174", "consensus", null);}
			catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Optional resolution must never weaken a registered required role");
		}finally{
			GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.addAll(previous);
			CallGenes.COBALAMIN_ENABLED=enabled;CallGenes.NCRNA_FAMILIES_ENABLED=umbrella;
			CallGenes.COBALAMIN_RESOURCE_ROOT=oldRoot;CallGenes.NCRNA_FAMILY_FILTER=filter;
			CallGenes.NCRNA_ID_PASS_OVERRIDE=oldPass;CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE=oldBorder;
		}
		System.out.println("COBALAMIN_MEMBER_CUTOFFS_TEST_PASS");
	}
	private static void load(Path root){
		GeneCaller.ncrnaFamilies.clear();CallGenes.COBALAMIN_RESOURCE_ROOT=root.toString();
		CallGenes.validateNcrnaSweepOverrides();CallGenes.loadCobalaminResources();
		check(GeneCaller.ncrnaFamilies.size()==3, "Successful loading must publish exactly the three B12 members");
	}
	/** Real member loading must propagate every model identity without changing its verifier. */
	private static void modelProfiles(Path source, Path root){
		for(NcrnaFamily family:GeneCaller.ncrnaFamilies){
			final Path table=root.resolve(family.name+".models.tsv");
			final ByteStreamWriter writer=RrnaResourceIO.writer(table.toString());final ByteBuilder b=new ByteBuilder();
			try{
				writer.println("model\tmodel_length\tmin_length\tmax_length\tidentity");
				for(int i=0; i<family.library.length; i++){
					writer.print(b.clear().append(family.modelNames[i]).tab().append(family.library[i].length)
						.append("\t1\t2147483647\t").append(i%2==0 ? .6f : .8f, 2).nl());
				}
			}finally{RrnaResourceIO.finish(writer);}
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			options.resources.root=source.toString();options.resources.lengthBounds=table.toString();
			final NcrnaFamily loaded=NcrnaFamilyLoader.load(family.name, options, new FamilyResourceResolver());
			check(!loaded.reuseConsensusAlignment && loaded.modelThresholds!=null,
				"Per-model profiles must reach the existing generic member verifier");
			for(int i=0; i<loaded.library.length; i++){
				final float pass=i%2==0 ? .6f : .8f;
				check(loaded.modelThresholds.pass(i)==pass && loaded.modelThresholds.borderline(i)==Math.min(family.idBorderline, pass),
					"Model identity must override member defaults, with a reachable direct-pass interval");
			}
		}
	}
	private static void identities(float... values){
		check(values.length==6, "Each of three members needs an independently specified identity pair");
		for(int i=0; i<3; i++){
			final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(i);
			check(family.idPass==values[2*i] && family.idBorderline==values[2*i+1],
				"Wrong effective identity for "+family.name+": "+family.idPass+"/"+family.idBorderline);
		}
	}
	private static void fails(Path root, String why){
		boolean failed=false;try{load(root);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && GeneCaller.ncrnaFamilies.isEmpty(), why);
	}
	/** Fixture payloads are small resource copies only, never program/class trees. */
	private static Path group(Path source, Path destination, String[] rows, int missing)throws Exception{
		check(rows.length==3 && missing>=-1 && missing<3, "Fixture must describe each declared group member");
		Files.createDirectory(destination);Files.copy(source.resolve("seeds.fa.gz"), destination.resolve("seeds.fa.gz"));
		for(String member:MEMBERS){
			Files.createDirectory(destination.resolve(member));
			for(String name:new String[]{"consensus.fa", "models.hbmc.gz"}){
				Files.copy(source.resolve(member).resolve(name), destination.resolve(member).resolve(name));
			}
		}
		final ByteFile in=RrnaResourceIO.open(source.resolve("manifest.tsv").toString());
		final ByteStreamWriter out=RrnaResourceIO.writer(destination.resolve("manifest.tsv").toString());
		final LineParser1 parser=new LineParser1('\t');
		try{
			for(byte[] line; (line=in.nextLine())!=null; ){
				parser.set(line);
				// Each fixture supplies its own optional cutoffs, including deliberately absent rows.
				if(parser.terms()==8 && parser.termEquals("cutoffs", 1)){continue;}
				out.println(line);
			}
			for(int i=0; i<3; i++){
				if(rows[i]==null){continue;}
				final String text="family\tidpass\tidborderline\n"+rows[i];
				final byte[] bytes=text.getBytes(java.nio.charset.StandardCharsets.UTF_8);
				if(i!=missing){
					final ByteStreamWriter cutoff=RrnaResourceIO.writer(destination.resolve(MEMBERS[i]).resolve("cutoffs.tsv").toString());
					try{cutoff.print(bytes);}finally{RrnaResourceIO.finish(cutoff);}
				}
				// Runtime manifest validation checks length; fixture digest is explicitly synthetic.
				out.println(MEMBERS[i]+"\tcutoffs\t.\t.\tcutoffs.tsv\ttsv\t"+bytes.length+"\t00000000000000000000");
			}
		}finally{check(!in.close(), "Fixture manifest read failed");RrnaResourceIO.finish(out);}
		return destination;
	}
	private static String row(int member, String pass, String borderline){return "cobalamin_"+MEMBERS[member].toLowerCase(java.util.Locale.ROOT)+'\t'+pass+'\t'+borderline+'\n';}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
}
