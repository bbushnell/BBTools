package test.prok;

import java.io.File;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlHoldout;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaCmTraining;
import shared.Shared;

/** Literal accepted-call support with shuffled source audits, alias-group
 * exclusions, unmatched calls and repeated winners from the same assembly.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaCmTrainingTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaCmTrainingTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Use out=fresh_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Fixture owns fresh output files");
			run(parser.out1);
			System.err.println("NCRNA_CM_TRAINING_INPUTS_FIXTURE_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final String root){
		try(Output output=new Output(root+"/sources.tsv")){
			output.writer.println(NcrnaDdlHoldout.PANEL_HEADER);
			for(int i=0; i<4; i++){output.writer.println(i+"\tasm"+i+"\tin"+i+"\t/literal/"+i+".fa\t/literal/"+i+".gff");}
		}
		for(int member=0; member<3; member++){
			try(Output output=new Output(root+"/audit"+member+".tsv")){
				output.writer.println("assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled");
				for(final int i:new int[]{3, 1, 0, 2}){
					output.writer.println("asm"+i+"\tin"+i+"\t/literal/"+i+".fa\t-1\t"+Math.min(i, 2)
						+'\t'+(member==1 && i==3)+"\tfalse\tfalse");
				}
			}
			try(Output output=new Output(root+"/ref"+member+".fa")){
				output.writer.println('>'+NAMES[member]+" n=10");
				output.writer.println("ACGT");
			}
		}
		try(Output output=new Output(root+"/exclude.tsv")){
			output.writer.println("# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result");
			output.writer.println("panel\trow1\tfasta\tassembly\tstatus");
			output.writer.println("/old_panel\t1\t/old_unknown.fa\t.\tUNMATCHED_MANIFEST_PATH_ALIAS");
		}
		try(Output output=new Output(root+"/prefixes.tsv")){
			output.writer.println("prefix");
			output.writer.println(root+"/ordinary");
		}
		try(Output output=new Output(root+"/ordinary.queries.tsv")){
			output.writer.println(PrepareNcrnaCmTraining.QUERY_HEADER);
			query(output, 0, 10, 0, "TP");
			query(output, 0, 30, 0, "TP");
			query(output, 1, 10, 0, "UNMATCHED");
			query(output, 1, 30, 1, "TP");
			query(output, 2, 10, 0, "UNMATCHED");
		}
		invoke(root, "sources.tsv", "pass");
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(root+"/pass.summary.tsv")){
			check(input.reader.nextLine()!=null, "Summary requires its declared header");
			p.set(input.reader.nextLine());
			final long[] expected={4, 2, 5, 1, 4, 3, 1, 1, 1, 4};
			check(p.terms()==expected.length, "Summary must retain every independent denominator");
			for(int i=0; i<expected.length; i++){check(p.parseLong(i)==expected[i], "Literal summary differs at column "+i);}
			check(input.reader.nextLine()==null, "Exactly one summary row expected");
		}
		try(Input input=new Input(root+"/pass.models.tsv")){
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			check(p.parseInt(3)==3 && p.parseInt(4)==2 && p.termEquals("true", 5), "Three wins from two assemblies pass, including the unmatched call");
			p.set(input.reader.nextLine());
			check(p.parseInt(3)==1 && p.parseInt(4)==1 && p.termEquals("false", 5), "The second member remains below support");
			p.set(input.reader.nextLine());
			check(p.parseInt(3)==0 && p.termEquals("false", 5), "Unsupported third model remains visible and unrefined");
		}
		try(Input input=new Input(root+"/pass.calls.tsv.gz")){
			input.reader.nextLine();
			int rows=0;
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.parseInt(0)==rows++ && !p.termEquals("asm2", 1), "Excluded alias group must never enter saved path inputs");
			}
			check(rows==4, "Every eligible accepted query survives, including unmatched calls");
		}
		// Two literal source IDs can denote one nonexcluded physical assembly.
		// Keep their calls, but do not manufacture a second assembly of support.
		for(int member=0; member<3; member++){
			try(Output output=new Output(root+"/grouped_audit"+member+".tsv")){
				output.writer.println("assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled");
				for(final int i:new int[]{3, 1, 0, 2}){
					output.writer.println("asm"+i+"\tin"+i+"\t/literal/"+i+".fa\t-1\t"+(i<2 ? 0 : 2)
						+'\t'+(member==1 && i==3)+"\tfalse\tfalse");
				}
			}
		}
		invoke(root, "sources.tsv", "grouped", "grouped_audit");
		try(Input input=new Input(root+"/grouped.models.tsv")){
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			check(p.parseInt(3)==3 && p.parseInt(4)==1 && p.termEquals("false", 5),
				"Three calls from two aliases of one assembly cannot satisfy the two-assembly support rule");
		}
		try(Output output=new Output(root+"/bad_sources.tsv")){
			output.writer.println(NcrnaDdlHoldout.PANEL_HEADER);
			for(int i=0; i<4; i++){output.writer.println(i+"\tasm"+i+"\t"+(i==1 ? "in0" : "in"+i)+"\t/literal/"+i+".fa\t/literal/"+i+".gff");}
		}
		boolean rejected=false;
		try{invoke(root, "bad_sources.tsv", "bad");}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Ambiguous source IDs cannot produce a training census");
	}

	private static void query(final Output output, final int source, final int start, final int member, final String outcome){
		assert(source>=0 && source<4 && member>=0 && member<3) : "Literal query must index the declared fixture source/model axes";
		output.writer.println("in"+source+"\tcontig\t"+start+'\t'+(start+10)+"\t+\t"+MEMBERS[member]+'\t'+NAMES[member]
			+" n=10\t"+outcome+"\ttrue\t"+(outcome.equals("TP") ? "reference"+source+'_'+start : "."));
	}

	private static void invoke(final String root, final String sources, final String stem){
		invoke(root, sources, stem, "audit");
	}
	private static void invoke(final String root, final String sources, final String stem, final String audit){
		assert(root!=null && sources!=null && stem!=null) : "Each fixture invocation owns distinct complete inputs and an output prefix";
		PrepareNcrnaCmTraining.main(new String[]{"in="+root+"/prefixes.tsv", "bindings="+root+'/'+sources,
			"audits="+root+'/'+audit+"0.tsv,"+root+'/'+audit+"1.tsv,"+root+'/'+audit+"2.tsv",
			"refs="+root+"/ref0.fa,"+root+"/ref1.fa,"+root+"/ref2.fa", "exclude="+root+"/exclude.tsv",
			"expectedinputs=4", "expectedmodels=3", "expectedprefixes=1", "mincalls=3", "minassemblies=2", "out="+root+'/'+stem});
	}

	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] NAMES={"cobalamin_RF00174_a", "cobalamin_RF01482_b", "cobalamin_RF01689_c"};
}
