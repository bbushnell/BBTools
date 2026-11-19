package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionBindings;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal old-member binding fixtures; no biological eligibility or sequence
 * equivalence is inferred from the deliberately synthetic catalog claims.
 * @author Ganyu
 */
public final class NcrnaConstructionBindingsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionBindingsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown binding fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use assertions and out=fresh_fixture_directory without input/sampling overrides");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root.getPath());
			pp.outstream.println("NCRNA_CONSTRUCTION_BINDINGS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final Fixture good=new Fixture(root, "good");
		final NcrnaConstructionBindings result=good.load();
		check(result.members()==3 && result.kept()==2, "The orphan must remain a seed input outside the two kept consensus members");
		result.write(good.dir+"/bound");
		checkSummary(good.dir+"/bound", 3);
		NcrnaConstructionBindings.main(new String[]{"genomes="+good.dir+"/genomes.tsv", "loci="+good.dir+"/loci.tsv",
			"assignments="+good.dir+"/assignments.tsv", "evalcensus="+good.dir+"/eval.census.tsv",
			"evalmembers="+good.dir+"/eval.members.tsv", "newcensus="+good.dir+"/new.census.tsv",
			"newmembers="+good.dir+"/new.members.tsv", "expectedgenomes=3", "expectedmembers=3", "out="+good.dir+"/cli", "t=1"});
		checkSummary(good.dir+"/cli", 3);
		try(Input in=new Input(good.dir+"/bound.members.tsv.gz")){
			check(in.reader.nextLine()!=null, "Member output requires a header");
			final String a=line(in), b=line(in), c=line(in);
			check(a.startsWith("true\ttrue\tEVALUATION_PATH\tUNVERIFIED\tNOT_EVALUATED\tRF00050_0_0\t"),
				"Evaluation path membership remains an exclusion identity, not a sequence proof");
			check(b.startsWith("true\tfalse\tNO_RECORDED_EVALUATION_IDENTITY\tUNVERIFIED\tNOT_EVALUATED\tRF00050_1_0\t"),
				"Discovery input1 must not collide with evaluation input1 or silently lose its seed-only member");
			check(c.startsWith("true\ttrue\tEVALUATION_ASSEMBLY_CLAIM\tUNVERIFIED\tNOT_EVALUATED\tRF00050_2_0\t"),
				"A different path with the evaluation assembly claim must remain marked");
			check(in.reader.nextLine()==null, "Exactly three original members must survive");
		}
		final Fixture multi=new Fixture(root, "multiple_bindings");
		multi.load().write(multi.dir+"/bound");
		checkSummary(multi.dir+"/bound", 4);
		final Fixture multipart=new Fixture(root, "multipart");
		multipart.load().write(multipart.dir+"/bound");
		try(Input in=new Input(multipart.dir+"/bound.members.tsv.gz")){
			line(in);
			line(in);
			check(line(in).contains("\t21-30,41-50\t20\t"), "A multipart source must preserve its gap instead of acquiring outer-span bases");
			line(in);
			check(in.reader.nextLine()==null, "Multipart preservation must not add or drop members");
		}
		final String[] faults={"missing_path", "missing_assignment", "duplicate_assignment", "changed_header",
			"blank_kept_model", "orphan_model", "wrong_count", "wrong_taxonomy", "bad_interval", "duplicate_locus",
			"bad_roster_index", "unknown_status"};
		final String[] reasons={"no literal cohort", "cover all original", "exactly one assignment", "exactly one assignment",
			"Only kept", "Only kept", "declared construction", "identity claims", "exceed contig bounds", "unique identity",
			"roster order", "Unknown original builder"};
		check(faults.length==reasons.length, "Every malformed join needs its expected diagnostic");
		for(int i=0; i<faults.length; i++){
			final Fixture bad=new Fixture(root, faults[i]);
			try{bad.load();}
			catch(IllegalArgumentException error){
				check(error.getMessage().contains(reasons[i]), "Wrong rejection for "+faults[i]+": "+error.getMessage());
				continue;
			}
			throw new AssertionError("Malformed original binding was accepted: "+faults[i]);
		}
	}

	private static void checkSummary(final String prefix, final int bindings){
		try(Input in=new Input(prefix+".summary.tsv")){
			check(in.reader.nextLine()!=null, "Summary header is required");
			check(line(in).equals("3\t3\t2\t1\t"+bindings+"\tCONSTRUCTION_IDENTITIES_BOUND_NOT_ELIGIBILITY_OR_SEQUENCE_VALIDATION"),
				"Original population and all unique cohort/input bindings must be conserved");
			check(in.reader.nextLine()==null, "Summary has one population row");
		}
	}

	private static String line(final Input in){
		final byte[] row=in.reader.nextLine();
		check(row!=null, "Expected fixture output row is missing");
		return new String(row, StandardCharsets.UTF_8);
	}

	private static final class Fixture{
		Fixture(final String root, final String fault){
			dir=root+'/'+fault;
			check(new File(dir).mkdir(), "Each fixture case needs a fresh directory");
			write("eval.census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/eval/current.fa\t100\t100\t1\n");
			write("eval.members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
				+membership("eval", 1, 1, "/old/eval.fa", "GCF_SHARED.1")
				+membership("eval", 2, 1, "/old/eval.fa", "GCF_SHARED.1"));
			write("new.census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/old/train.fa\t100\t100\t1\n2\t/old/alias.fa\t100\t100\t1\n");
			write("new.members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
				+membership("new", 1, 1, "/old/train.fa", "GCF_TRAIN.1")
				+membership("new", 2, 2, "/old/alias.fa", "GCF_SHARED.1")
				+(fault.equals("multiple_bindings") ? membership("new", 3, 2, "/old/train.fa", "GCF_SHARED.1") : ""));
			write("genomes.tsv", NcrnaConstructionBindings.GENOMES_HEADER+'\n'
				+genome(fault.equals("bad_roster_index") ? 4 : 0, "eval", "GCF_SHARED.1", "/old/eval.fa", 1)
				+genome(1, "train", "GCF_TRAIN.1", fault.equals("missing_path") ? "/missing.fa" : "/old/train.fa", fault.equals("wrong_count") ? 2 : 1)
				+genome(2, "alias", "GCF_SHARED.1", "/old/alias.fa", 1));
			final String a=locus(0, "eval", "GCF_SHARED.1", "1", "1-20");
			String trainLocus=locus(1, "train", "GCF_TRAIN.1", fault.equals("wrong_taxonomy") ? "2" : "1",
				fault.equals("bad_interval") ? "90-109" : fault.equals("multipart") ? "21-30,41-50" : "21-40");
			if(fault.equals("multipart")){trainLocus=trainLocus.replace("ANNOTATION_FULL_COPY", "MULTIPART_UNRESOLVED");}
			write("loci.tsv", NcrnaConstructionBindings.LOCI_HEADER+'\n'+a+trainLocus
				+locus(2, "alias", "GCF_SHARED.1", "1", "41-60")+(fault.equals("duplicate_locus") ? a : ""));
			final String assigned=assignment(0, "eval", fault.equals("changed_header") ? " changed" : "",
				fault.equals("blank_kept_model") ? "" : "FMN_model0", fault.equals("unknown_status") ? "unknown" : "kept");
			String trainAssignment=assignment(1, "train", "", fault.equals("orphan_model") ? "FMN_model1" : "", "orphan");
			if(fault.equals("multipart")){trainAssignment=trainAssignment.replace("ANNOTATION_FULL_COPY", "MULTIPART_UNRESOLVED");}
			write("assignments.tsv", NcrnaConstructionBindings.ASSIGNMENTS_HEADER+'\n'+assigned+trainAssignment
				+(fault.equals("missing_assignment") ? "" : assignment(2, "alias", "", "FMN_model0", "kept"))
				+(fault.equals("duplicate_assignment") ? assigned : ""));
		}
		NcrnaConstructionBindings load(){
			return NcrnaConstructionBindings.load(dir+"/genomes.tsv", dir+"/loci.tsv", dir+"/assignments.tsv",
				NcrnaEvidenceCatalog.load(dir+"/eval.census.tsv", dir+"/eval.members.tsv"),
				NcrnaEvidenceCatalog.load(dir+"/new.census.tsv", dir+"/new.members.tsv"), 3, 3);
		}
		void write(final String file, final String text){
			try(Output out=new Output(dir+'/'+file)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		}
		final String dir;
	}
	private static String membership(final String cohort, final int order, final int id, final String fasta, final String assembly){
		return cohort+'\t'+order+'\t'+id+'\t'+fasta+"\tclaim\t"+assembly+"\t1\tprokaryote\t/original.gff\tNA\n";
	}
	private static String genome(final int index, final String name, final String assembly, final String fasta, final int count){
		return index+"\t"+name+'\t'+assembly+"\t1\tbacteria\tP\t/"+name+".gff\t"+fasta+'\t'+count+"\t100\t1\n";
	}
	private static String locus(final int index, final String name, final String assembly, final String taxid, final String intervals){
		return "RF00050_"+index+"_0\t"+name+'\t'+assembly+'\t'+taxid+"\tbacteria\tP\tcontig\t+\t"+intervals
			+"\t20\t7\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t100\tfalse\t0\t0\t0\t7\n";
	}
	private static String assignment(final int index, final String name, final String suffix, final String model, final String status){
		return "RF00050_"+index+"_0 source="+name+" family=RF00050 status=ANNOTATION_FULL_COPY"+suffix
			+"\t1\tunknown\tRF00050_cluster\t"+model+'\t'+status+'\n';
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
