package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionAssemblies;
import prok.NcrnaConstructionBindings;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal assembly identities and dimension-only candidates, not sequence proof.
 * Distinct short sequences of equal length deliberately remain comparison candidates.
 * @author Ganyu
 */
public final class NcrnaConstructionAssembliesTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssembliesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown assembly fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea and out=fresh_fixture_directory without sampling/overwrite overrides");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Assembly fixture directory must be fresh");
			run(root.getPath());
			pp.outstream.println("NCRNA_CONSTRUCTION_ASSEMBLIES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final Fixture good=new Fixture(root, "good");
		good.run("result");
		checkResult(good, "result");
		NcrnaConstructionAssemblies.main(new String[]{"genomes="+good.path("original.tsv"),
			"evalcensus="+good.path("eval.tsv"), "evalmembers="+good.path("eval.members.tsv"),
			"newcensus="+good.path("new.tsv"), "newmembers="+good.path("new.members.tsv"),
			"expectedgenomes=4", "expectedeval=2", "expectednew=4", "out="+good.path("cli"), "t=1"});
		checkResult(good, "cli");
		final Fixture unicode=new Fixture(root, "unicode_\u00e9");
		unicode.run("result");
		checkResult(unicode, "result");
		check(rows(unicode.path("result.inputs.tsv.gz")).get(1).contains(unicode.path("a.fa")+'\t'+unicode.path("a.fa")),
			"Both literal and canonical UTF-8 paths must round-trip without ByteBuilder String byte-casting");
		final Fixture foreign=new Fixture(root, "foreign_membership");
		foreign.run("result");
		check(rows(foreign.path("result.summary.tsv")).get(1).equals("2\t2\t1\t1\t0\t1\t0\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT")
			&& rows(foreign.path("result.inputs.tsv.gz")).get(6).endsWith("\tNA\t5\tfalse\tfalse\ttrue\tASSEMBLY_IDENTITY_UNRESOLVED"),
			"A foreign membership FASTA must not falsely attach the evaluation assembly to an otherwise unbound census input");
		final Fixture equalForeign=new Fixture(root, "foreign_equal_dimensions");
		equalForeign.run("result");
		check(rows(equalForeign.path("result.summary.tsv")).get(1).equals("2\t2\t2\t0\t0\t2\t0\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT")
			&& rows(equalForeign.path("result.inputs.tsv.gz")).get(6).endsWith("\tNA\t5\tfalse\ttrue\ttrue\tASSEMBLY_SEQUENCE_CHECK_PENDING"),
			"Even equal dimensions cannot bind a foreign membership path; they only nominate a sequence comparison");
		final Fixture unresolved=new Fixture(root, "unresolved");
		unresolved.run("result");
		check(rows(unresolved.path("result.summary.tsv")).get(1).equals("2\t2\t1\t1\t0\t1\t0\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT"),
			"Unsupported assembly claims must remain unresolved rather than quietly admitted");
		final Fixture unmeasured=new Fixture(root, "unmeasured_original");
		unmeasured.run("result");
		check(rows(unmeasured.path("result.summary.tsv")).get(1).equals("2\t2\t1\t0\t1\t1\t1\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT")
			&& rows(unmeasured.path("result.inputs.tsv.gz")).get(5).endsWith("\t4\tfalse\ttrue\tfalse\tASSEMBLY_SEQUENCE_CHECK_PENDING"),
			"The original extractor's zero-locus paired zeros are unmeasured, so the positive census still supplies candidate dimensions");
		for(final String fault:new String[]{"wrong_size", "wrong_dimensions", "conflicting_claim", "wrong_roster_count",
				"zero_with_loci", "one_zero_dimension", "zero_loci_with_dimensions"}){
			final Fixture bad=new Fixture(root, fault);
			boolean rejected=false;
			try{bad.run("result");}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Malformed assembly binding was accepted: "+fault);
		}
	}

	private static void checkResult(final Fixture f, final String prefix){
		check(rows(f.path(prefix+".summary.tsv")).get(1).equals("2\t2\t1\t0\t1\t1\t0\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT"),
			"Two evaluation inputs, two known aliases, one unresolved sequence candidate and one training candidate must conserve six inputs");
		final ArrayList<String> inputs=rows(f.path(prefix+".inputs.tsv.gz"));
		check(inputs.size()==7 && inputs.get(3).endsWith("\t0\ttrue\tfalse\tfalse\tEXCLUDE_KNOWN_EVALUATION_ASSEMBLY")
			&& inputs.get(4).endsWith("\t0\ttrue\tfalse\tfalse\tEXCLUDE_KNOWN_EVALUATION_ASSEMBLY"),
			"Both the explicit assembly-record version and canonical symlink alias must join the evaluation group");
		check(inputs.get(5).endsWith("\t4\tfalse\ttrue\tfalse\tASSEMBLY_SEQUENCE_CHECK_PENDING")
			&& inputs.get(6).endsWith("\t5\tfalse\tfalse\tfalse\tTRAIN_CANDIDATE_ASSEMBLY_IDENTITY"),
			"Equal dimensions do not prove identity and different dimensions do not erase known assembly provenance");
		final ArrayList<String> pairs=rows(f.path(prefix+".candidate_pairs.tsv.gz"));
		check(pairs.size()==2 && pairs.get(1).equals("1\t3\t"+f.path("a.fa")+'\t'+f.path("c.fa")
			+"\t3\t1\tDIMENSION_COMPATIBILITY_NOT_SEQUENCE_IDENTITY"),
			"Candidate comparison must preserve cohort-qualified inputs and exact paths without asserting these different sequences are aliases");
	}

	private static final class Fixture{
		Fixture(final String root, final String fault_){
			fault=fault_;
			dir=root+'/'+fault;
			final boolean foreign=fault.equals("foreign_membership") || fault.equals("foreign_equal_dimensions");
			check(new File(dir).mkdir(), "Every assembly fixture needs a fresh directory");
			write("a.fa", ">a\nAAA\n");
			write("b.fa", ">b\nAAAA\n");
			write("c.fa", ">c\nCCC\n");
			write("d.fa", fault.equals("foreign_equal_dimensions") ? ">d\nGGG\n" : ">d\nAAAAA\n");
			write("e.fa", ">e\nAAAAAAA\n");
			if(foreign){write("unselected.fa", ">u\nAAAAA\n");}
			try{Files.createSymbolicLink(new File(path("alias.fa")).toPath(), new File(path("a.fa")).toPath());}
			catch(IOException error){throw new IllegalStateException("The native Dori fixture needs a real symlink alias", error);}
			write("original.tsv", NcrnaConstructionBindings.GENOMES_HEADER+'\n'
				+original(0, "a", "GCF_100.1", 3)+original(1, "b", "GCF_100.2", 4)
				+original(2, "c", "GCF_300.1", 3)
				+original(3, foreign ? "unselected" : "d", "GCF_400.1", 5));
			write("eval.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+census(1, "a.fa", 3)+census(2, "e.fa", 7));
			write("new.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+census(1, "b.fa", 4)+census(2, "alias.fa", 3)
				+census(3, "c.fa", fault.equals("wrong_dimensions") ? 4 : 3)
				+census(4, "d.fa", fault.equals("foreign_equal_dimensions") ? 3 : 5));
			write("eval.members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+membership(1, "a.fa", "NA")+membership(2, "e.fa", "NA"));
			write("new.members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+membership(1, "b.fa", "NA")
				+membership(2, "alias.fa", "NA")+membership(3, "c.fa", "NA")
				+membership(4, foreign ? "a.fa" : "d.fa",
					foreign ? "GCF_100.1" : fault.equals("unresolved") ? "unknown-record"
					: fault.equals("conflicting_claim") ? "GCF_999.1" : "NA"));
		}
		void run(final String prefix){
			NcrnaConstructionAssemblies.process(path("original.tsv"),
				NcrnaEvidenceCatalog.load(path("eval.tsv"), path("eval.members.tsv")),
				NcrnaEvidenceCatalog.load(path("new.tsv"), path("new.members.tsv")),
				fault.equals("wrong_roster_count") ? 5 : 4, 2, 4, path(prefix));
		}
		String original(final int index, final String name, final String assembly, final int bases){
			final boolean target=name.equals("c");
			final boolean unmeasured=target && fault.equals("unmeasured_original");
			final int loci=unmeasured || target && fault.equals("zero_loci_with_dimensions") ? 0 : 1;
			final int observedBases=unmeasured || target && (fault.equals("zero_with_loci") || fault.equals("one_zero_dimension")) ? 0 : bases;
			final int contigs=unmeasured || target && fault.equals("zero_with_loci") ? 0 : 1;
			return index+"\t"+name+'\t'+assembly+"\t1\tbacteria\tP\t"+path(name+".gff")+'\t'+path(name+".fa")
				+'\t'+loci+'\t'+observedBases+'\t'+contigs+'\n';
		}
		String census(final int id, final String file, final int bases){
			return id+"\t"+path(file)+'\t'+(new File(path(file)).length()+(fault.equals("wrong_size") && file.equals("c.fa") ? 1 : 0))
				+'\t'+bases+"\t1\n";
		}
		String membership(final int id, final String file, final String assembly){
			return "fixture\t"+id+'\t'+id+'\t'+path(file)+"\tNA\t"+assembly+"\tNA\tprokaryotic\tNA\tNA\n";
		}
		void write(final String file, final String text){
			try(Output out=new Output(path(file))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		}
		String path(final String file){return dir+'/'+file;}
		final String dir, fault;
	}

	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){result.add(new String(line, StandardCharsets.UTF_8));}
		}
		check(!result.isEmpty(), "Assembly fixture output must include its header");
		return result;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
