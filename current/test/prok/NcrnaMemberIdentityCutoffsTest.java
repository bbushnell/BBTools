package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMemberIdentityCutoffs;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaFamily;
import prok.NcrnaModelThresholds;
import prok.TrnaConsensusBuilder;
import prok.RrnaResourceIO;

/** Saved assigned identities, percentile convention, clamps and singleton provenance.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaMemberIdentityCutoffsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String model="cobalamin_RF00174_original", single="cobalamin_RF00174_singleton";
		final String modelHeader=model+" n=21", singletonHeader=single+" n=1";
		write(root.resolve("ref.fa"), ">"+modelHeader+"\nACGTAC\n>"+singletonHeader+"\nATGCAT\n");
		final StringBuilder fasta=new StringBuilder(), quality=new StringBuilder(NcrnaMemberIdentityCutoffs.QUALITY_HEADER).append('\n');
		for(int i=0;i<22;i++){
			fasta.append(">m").append(i).append(" member_family=RF00174 build_group=RF00174\n").append(i==0?"ATGCAT":"ACGTAC").append('\n');
			quality.append("m").append(i).append("\t6\t");
			if(i==21){quality.append("orphan\tNA\tNA");}
			else{quality.append("kept\t").append(model).append('\t').append(i==0?"0.580000":i==1?"0.700000":"0.900000");}
			quality.append('\t').append(single).append("\t1.000000\ttrue\n");
		}
		write(root.resolve("in.fa"), fasta.toString());write(root.resolve("quality.tsv"), quality.toString());
		final String added=NcrnaMemberIdentityCutoffs.ADDITIONS_HEADER+'\n'+"m0\t"+single+"\tRF00174\tRF00174\t6\tEXACT_SINGLE_MEMBER_NO_REASSIGNMENT\n";
		write(root.resolve("added.tsv"), added);
		final String members=NcrnaMemberLengthBounds.HEADER+'\n'+row(single, 1)+row(model, 21);
		write(root.resolve("members.tsv"), members);run(root, "valid.");
		check(table(root.resolve("valid.A.tsv"), 1).equals(modelHeader+"\t6\t4\t8\t0.600000"), "A uses assigned minimum and lower clamp, not best identities");
		check(table(root.resolve("valid.B.tsv"), 1).equals(modelHeader+"\t6\t4\t8\t0.700000"), "p5 of21 values uses zero-based index1");
		check(table(root.resolve("valid.C.tsv"), 1).equals(modelHeader+"\t6\t4\t8\t0.680000"), "C subtracts .02 before clamp");
		for(String rule:new String[]{"A", "B", "C", "D"}){
			check(table(root.resolve("valid."+rule+".tsv"), 2).equals(singletonHeader+"\t6\t4\t8\t0.800000"), "All rules preserve complete caller names and clamp exact singleton to.80");
			final NcrnaFamily loaded=loadProfile(root.resolve("ref.fa"), root.resolve("valid."+rule+".tsv"));
			check(loaded.modelNames[0].equals(modelHeader) && loaded.modelThresholds.pass(1)==.8f, "Actual caller loading must consume every exported rule with full FASTA header metadata");
		}
		check(table(root.resolve("valid.models.tsv"), 1).contains("\t21\tfalse\t0.580000\t0.700000\t"), "Orphan best-model ties must not add construction members");
		write(root.resolve("quality.tsv"), quality.toString().replace("0.900000", "0.760000"));run(root, "floor64.", "floor=.64");
		check(table(root.resolve("floor64.A.tsv"), 1).endsWith("\t0.640000"), "Explicit floor64 replaces historical floor60");
		check(table(root.resolve("floor64.B.tsv"), 1).endsWith("\t0.700000") && table(root.resolve("floor64.C.tsv"), 1).endsWith("\t0.680000"), "Raising floor must not change already eligible p5 cutoffs");
		check(table(root.resolve("floor64.D.tsv"), 1).endsWith("\t0.710000"), "D is member median.76 minus.05, not p5 or a capped median");
		check(loadProfile(root.resolve("ref.fa"), root.resolve("floor64.D.tsv")).modelThresholds.pass(0)==.71f, "Actual loader must receive D's derived cutoff");
		write(root.resolve("quality.tsv"), quality.toString());
		int invalid=0;for(String value:new String[]{"NaN", "Infinity", "-.01", ".81"}){
			final String prefix="invalid_floor_"+(invalid++)+'.';boolean rejected=false;
			try{run(root, prefix, "floor="+value);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !Files.exists(root.resolve(prefix+"A.tsv")), "Invalid floor must fail before producing any cutoff table");
		}
		write(root.resolve("ref.fa"), ">"+modelHeader+"\nACGTAC\n>"+singletonHeader+"\nNNNNNN\n");
		write(root.resolve("in.fa"), fasta.toString().replace("ATGCAT", "NNNNNN"));run(root, "ambiguous.");
		check(table(root.resolve("ambiguous.models.tsv"), 2).contains("\t1\ttrue\t0.500000\t0.500000\t"), "Identical literal Ns receive half credit from the actual Scrabble/Tracer identity");
		for(String rule:new String[]{"A", "B", "C", "D"}){check(table(root.resolve("ambiguous."+rule+".tsv"), 2).endsWith("\t0.600000"), "All rules clamp ambiguous singleton identity at the historical lower floor");}
		run(root, "ambiguous64.", "floor=.64");
		for(String rule:new String[]{"A", "B", "C", "D"}){check(table(root.resolve("ambiguous64."+rule+".tsv"), 2).endsWith("\t0.640000"), "All four rules obey explicit floor64 even for ambiguous singleton contributors");}
		write(root.resolve("ref.fa"), ">"+modelHeader+"\nACGTAC\n>"+singletonHeader+"\nATGCAT\n");write(root.resolve("in.fa"), fasta.toString());
		write(root.resolve("quality.tsv"), quality.toString()+quality.substring(quality.lastIndexOf("m21\t")));fails(root, "duplicate.");
		write(root.resolve("quality.tsv"), quality.toString().replace("0.580000", "NaN"));fails(root, "nonfinite.");
		write(root.resolve("quality.tsv"), quality.toString());
		write(root.resolve("members.tsv"), NcrnaMemberLengthBounds.HEADER+'\n'+row(single, 1)+row(model, 20));fails(root, "count.");
		write(root.resolve("members.tsv"), members);
		write(root.resolve("in.fa"), fasta.toString().replace("ATGCAT", "ATGCAA"));fails(root, "changed_sequence.");
		write(root.resolve("in.fa"), fasta.toString());
		write(root.resolve("added.tsv"), added.replace("m0\t", "m21\t"));fails(root, "wrong_contributor.");
		write(root.resolve("added.tsv"), added+added.substring(added.indexOf('\n')+1));fails(root, "duplicate_singleton.");
		System.out.println("NCRNA_MEMBER_IDENTITY_CUTOFFS_TEST_PASS");
	}
	private static String row(String model, int count){return model+"\tRF00174\t6\t"+count+"\t6\t6\t6\t6\t6\t1\tPURE\t"+count+"\t1\t0\t0\t0\t0\t0\t0\n";}
	/** Exercise the actual caller's FASTA-name loader and profile consumer. */
	static NcrnaFamily loadProfile(Path refs, Path table){
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(refs.toString());
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", library, null, TrnaConsensusBuilder.lastLoadedNames, new map.LongHashSet(32), 17, 1, 1);
		NcrnaModelThresholds.loadLengths(table.toString(), family);return family;
	}
	private static void run(Path root, String prefix, String... extra){
		final String[] base={"in="+root.resolve("in.fa"), "ref="+root.resolve("ref.fa"), "quality="+root.resolve("quality.tsv"), "membership="+root.resolve("members.tsv"), "additions="+root.resolve("added.tsv"), "outprefix="+root.resolve(prefix), "family=RF00174", "expected=22", "models=2"};
		final String[] args=java.util.Arrays.copyOf(base, base.length+extra.length);System.arraycopy(extra, 0, args, base.length, extra.length);
		NcrnaMemberIdentityCutoffs.main(args);
	}
	private static void fails(Path root, String prefix){boolean failed=false;try{run(root, prefix);}catch(IllegalArgumentException expected){failed=true;}check(failed && !Files.exists(root.resolve(prefix+"A.tsv")), "Invalid memberships must fail before publishing thresholds");}
	private static String table(Path path, int index){final ByteFile in=RrnaResourceIO.open(path.toString());try{byte[] line=null;for(int i=0;i<=index;i++){line=in.nextLine();check(line!=null, "Missing fixture output row");}return new String(line);}finally{check(!in.close(), "Fixture read failed");}}
	/** Rewrites are deliberate fixture mutations inside the newly created test directory. */
	private static void write(Path file, String text){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(file.toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
