package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;

import prok.NcrnaLearningExclusions;
import prok.NcrnaTextIO.Input;
import shared.Shared;
import structures.ByteBuilder;

/** Literal exact-path panel exclusion oracle: canonical and original aliases,
 * repeated panel membership, a same-basename nonmatch, and ambiguous aliases.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningExclusionsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaLearningExclusionsTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown exclusion fixture argument: "+arg);
			}
			check(parser.out1!=null, "Supply out=fresh_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_LEARNING_EXCLUSIONS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private static void run(final File root){
		check(root.mkdir(), "Exclusion fixture output must be fresh");
		final String prefix=root.getPath()+File.separator;
		NcrnaLearningPopulationTest.write(prefix+"roster.tsv", NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(prefix+"loci.tsv", NcrnaLearningPopulationTest.loci());
		NcrnaLearningPopulationTest.write(prefix+"graph.tsv", NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(prefix+"original.fa", ">common\nAAAA\n>rare\nCCCC\n");
		final ByteBuilder manifest=new ByteBuilder("assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role\n");
		for(int i=0; i<=20; i++){
			manifest.append("assembly").append(i).tab().append(i).append("\ttest\t/canonical/assembly").append(i)
				.append(".fa\t/gff\t.\t.\t/source/assembly").append(i).append(".fa\t/gff\tfasta\tNCRNA\tTEST\n");
		}
		NcrnaLearningPopulationTest.write(prefix+"manifest.tsv", manifest.toString());
		NcrnaLearningPopulationTest.write(prefix+"first.list", "/canonical/assembly0.fa\n");
		NcrnaLearningPopulationTest.write(prefix+"second.list", "/source/assembly1.fa\n/source/assembly0.fa\n/elsewhere/assembly2.fa\n");
		invoke(prefix, "manifest.tsv", "result");
		try(Input input=new Input(prefix+"result/exclude.txt")){
			for(final String id:new String[]{"gene0", "gene1", "gene20"}){
				check(Arrays.equals(input.reader.nextLine(), id.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
					"Exclude every locus of exact assemblies zero/one, without duplicates or a basename guess");
			}
			check(input.reader.nextLine()==null, "Same-basename assembly two must stay in teaching data");
		}
		try(Input input=new Input(prefix+"result/summary.tsv")){
			check(input.reader.nextLine()!=null && input.reader.nextLine()!=null, "Held-out metadata and count schema remain visible");
			check(Arrays.equals(input.reader.nextLine(), "21\t3\t18\t20\t18\t2\t1".getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Twenty-one originals retain18 after3 exclusions;20eligible become18,2assemblies selected and1path unmatched");
			check(input.reader.nextLine()==null, "Summary is exactly one member-population row");
		}
		NcrnaLearningPopulationTest.write(prefix+"ambiguous.tsv", manifest.toString()
			.replace("/source/assembly1.fa", "/source/assembly0.fa"));
		boolean rejected=false;
		try{invoke(prefix, "ambiguous.tsv", "bad");}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"bad/exclude.txt").exists(), "An alias naming different assemblies must fail before emitting exclusions");
		NcrnaLearningPopulationTest.write(prefix+"first.list", "/elsewhere/assembly0.fa\n");
		rejected=false;
		try{invoke(prefix, "manifest.tsv", "bad_required");}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"bad_required/exclude.txt").exists(), "Every required-panel genome must map before exclusions are emitted");
	}
	private static void invoke(final String root, final String manifest, final String out){
		NcrnaLearningExclusions.main(new String[]{"in="+root+manifest, "panel="+root+"first.list", "optionalpanel="+root+"second.list",
			"roster="+root+"roster.tsv", "loci="+root+"loci.tsv", "graph="+root+"graph.tsv",
			"original="+root+"original.fa", "out="+root+out, "t=1"});
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
