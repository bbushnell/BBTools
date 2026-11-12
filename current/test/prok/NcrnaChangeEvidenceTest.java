package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.Parse;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaChangeEvidence;
import prok.NcrnaStageDiagLogSink;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Independent coordinate and original-input attribution fixtures; no calling.
 * The first input has no sidecars, so selecting it for the second input fails.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaChangeEvidenceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaChangeEvidenceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture parameter: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=fresh_directory");
			new NcrnaChangeEvidenceTest(new File(parser.out1)).run();
			pp.outstream.println("NCRNA_CHANGE_EVIDENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaChangeEvidenceTest(final File root_){
		root=root_.getAbsoluteFile();
		check(root.mkdir(), "Fixture directory must be fresh under an existing parent");
	}

	private void run(){
		check(NcrnaChangeEvidence.genomicStart(99, 109, 0, 1000)==100
			&& NcrnaChangeEvidence.genomicStop(99, 109, 0, 1000)==110, "Plus closed coordinates must be 100–110");
		check(NcrnaChangeEvidence.genomicStart(890, 900, 1, 1000)==100
			&& NcrnaChangeEvidence.genomicStop(890, 900, 1, 1000)==110, "Minus closed coordinates must be 100–110");
		write("policies.tsv", "policy\tmember\tmode\ttopn\tfloor\tminhits\ttopfrac\tqfrac\tpatience\n"
			+"shipped\tcobalamin_test\tfixed\t1\t0\t1\t0\t0\t-1\n"
			+"changed%name\tcobalamin_test\tadaptive\t1\t1\t1\t0\t0\t0\n");
		final String inputs="input_id\tpolicy\nfirst\tshipped\nfirst\tchanged%25name\nsecondβ%25\tshipped\nsecondβ%25\tchanged%25name\n";
		write("inputs.tsv", inputs);
		final String changes="input_id\tpolicy\tchange\treference_id\tseqid\tstart\tstop\tstrand\n"
			+"secondβ%25\tchanged%25name\tLOST\tparent\tseq%25α\t100\t110\t+\n";
		write("changes.tsv", changes);
		final String stage=NcrnaStageDiagLogSink.HEADER+"\n"
			+"SEED\tcobalamin_test\tseq%α description\t0\t1000\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t2\t99,109\n"
			+"WINDOW\tcobalamin_test\tseq%α description\t0\t1000\t1\t89\t119\tACCEPT_IDPASS\t2\t1\t0\tmodel\t0.7\t0\t0\t99\t109\ttrue\t.\t.\n"
			+"MODEL\tcobalamin_test\tseq%α description\t0\t1000\t1\t89\t119\tALIGNMENT\t.\t.\t0\tmodel\t0.7\t.\t.\t.\t.\t.\t.\t.\n"
			+"PATH\tcobalamin_test\tseq%α description\t0\t1000\t.\t.\t.\tDP_INPUT\t.\t.\t.\tmodel\t.\t.\t.\t99\t109\t.\t.\t.\n"
			+"PATH\tcobalamin_test\tseq%α description\t1\t1000\t.\t.\t.\tDP_RETAINED\t.\t.\t.\tother_model\t.\t.\t.\t890\t900\t.\t.\t.\n"
			+"WINDOW\tcobalamin_test\tseq%α description\t0\t1000\t1\t789\t819\tREJECT_ID\t2\t1\t0\tfar\t0.2\t0\t0\t-1\t-1\tfalse\t.\t.\n";
		final String gff="##gff-version 3\nseq%α\tBBTools\tCDS\t90\t105\t.\t+\t0\tID=near_cds\n"
			+"seq%α\tBBTools\tRNA\t100\t110\t.\t-\t.\tmodel:opposite\n"
			+"seq%α\tBBTools\trRNA\t300\t400\t.\t+\t.\tID=far_rrna\n";
		for(int policy=0; policy<2; policy++){
			final String marker=policy==0 ? "baseline_model" : "changed_model";
			write("sweep.input1.policy"+policy+".stages.tsv", stage.replace("\tmodel\t", "\t"+marker+"\t"));
			write("sweep.input1.policy"+policy+".gff", gff.replace("near_cds", policy==0 ? "baseline_near_cds" : "changed_near_cds"));
		}
		invoke("evidence.tsv", false);
		final String evidence=read("evidence.tsv");
		check(evidence.contains("\t90\t120\t100\t110\tSAVED_ROW\tWINDOW%09"), "Window and accepted spans must retain distinct converted coordinates");
		check(evidence.contains("\t100\t110\t100\t110\tSAVED_ROW\tPATH%09")
			&& evidence.contains("DP_RETAINED"), "Opposite-strand path evidence must survive");
		check(evidence.contains("near_cds") && evidence.contains("model:opposite")
			&& !evidence.contains("far_rrna") && !evidence.contains("REJECT_ID"), "Only overlapping closed contexts are selected at padding zero");
		check(!evidence.contains("ALIGNMENT") && evidence.contains("seq%25α"), "Default view omits model rows and preserves escaped Unicode contig identity");
		check(evidence.contains("\tsecondβ%25\tchanged%25name\tLOST\t"),
			"Original input and policy identities must decode for binding while raw metadata framing remains unchanged");
		invoke("models.tsv", true);
		check(read("models.tsv").contains("ALIGNMENT"), "Explicit models=t exposes observed model rows");
		attribution("evidence.tsv");
		write("empty_changes.tsv", changes+"secondβ%25\tchanged%25name\tLOST\tabsent_parent\tabsent_contig\t100\t110\t+\n");
		invoke("empty_context.tsv", false, "empty_changes.tsv", "inputs.tsv");
		emptyContexts("empty_context.tsv");
		write("missing_inputs.tsv", inputs+"third\tshipped\nthird\tchanged%25name\n");
		write("missing_changes.tsv", changes+"third\tchanged%25name\tLOST\tmissing_files\tseq%25α\t100\t110\t+\n");
		try{
			invoke("missing_sidecar.tsv", false, "missing_changes.tsv", "missing_inputs.tsv");
			throw new IllegalStateException("A missing later-batch sidecar must fail before output creation");
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("Cannot read native evidence files"), "Wrong missing-sidecar failure: "+expected);
			check(!new File(root, "missing_sidecar.tsv").exists(), "Preflight must finish every batch before any output opens");
		}
		write("invalid_inputs.tsv", inputs.replace("secondβ%25\tchanged%25name\n", ""));
		try{
			invoke("invalid.tsv", false, "changes.tsv", "invalid_inputs.tsv");
			throw new IllegalStateException("Truncated input-policy group must fail");
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("complete input/policy"), "Wrong negative fixture failure: "+expected);
			check(!new File(root, "invalid.tsv").exists(), "Malformed bindings must fail before output opens");
		}
	}

	/** Distinct model/GFF contents must remain under their actual source policy/file. */
	private void attribution(final String file){
		final LineParser1 p=new LineParser1('\t');
		int baseline=0, changed=0;
		try(Input in=new Input(path(file))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("change_row", 0)){continue;}
				check(p.terms()==20, "Witness schema must name both searched files and evidence kind");
				p.setBounds(19);
				final String original=new String(row, p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
				if(original.contains("baseline_model") || original.contains("baseline_near_cds")){
					check(p.termEquals("shipped", 9) && p.parseString(12).contains(".policy0."), "Baseline content was attributed to another policy");
					baseline++;
				}
				if(original.contains("changed_model") || original.contains("changed_near_cds")){
					check(p.termEquals("changed%25name", 9) && p.parseString(12).contains(".policy1."), "Changed-policy content was attributed to baseline");
					changed++;
				}
			}
		}
		check(baseline>0 && changed>0, "Both distinct source-policy payloads must actually be observed");
	}

	/** Two policies search existing files but find no rows for one declared contig. */
	private void emptyContexts(final String file){
		final LineParser1 p=new LineParser1('\t');
		int empty=0;
		try(Input in=new Input(path(file))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("change_row", 0) || !p.termEquals("NO_SAVED_EVIDENCE", 18)){continue;}
				check(p.terms()==20 && p.termEquals("absent_contig", 5) && p.parseInt(13)==0
					&& p.termEquals('.', 12) && p.termEquals('.', 19), "An empty witness must not invent a native path/line/row");
				check(p.parseString(10).endsWith(".stages.tsv") && p.parseString(11).endsWith(".gff"), "Empty evidence must identify both files searched");
				empty++;
			}
		}
		check(empty==2, "One no-row locus must have explicit baseline and changed-policy witnesses");
	}

	private void invoke(final String output, final boolean models){
		invoke(output, models, "changes.tsv", "inputs.tsv");
	}
	/** Each scenario owns fresh declarations; the production no-overwrite contract stays active. */
	private void invoke(final String output, final boolean models, final String changes, final String inputs){
		NcrnaChangeEvidence.main(new String[]{"in="+path(changes), "inputs="+path(inputs),
			"policies="+path("policies.tsv"), "calls="+path("sweep"), "family=cobalamin_",
			"models="+models, "out="+path(output)});
	}
	private void write(final String file, final String text){
		try(Output out=new Output(path(file))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String read(final String file){
		final ByteBuilder buffer=new ByteBuilder();
		try(Input in=new Input(path(file))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){buffer.append(line).nl();}
		}
		return new String(buffer.array, 0, buffer.length(), StandardCharsets.UTF_8);
	}
	private String path(final String file){return new File(root, file).toString();}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}
	private final File root;
}
