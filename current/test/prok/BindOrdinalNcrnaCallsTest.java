package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Locale;

import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BindOrdinalNcrnaCalls;
import prok.GradeCdsStartChanges;
import prok.GradeNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Independent binding oracle with repeated contigs, global IDs distinct from
 * ordinals/assembly names, empty middle/final inputs and unchanged feature bytes.
 * Nonexistent FASTA/reference paths prove this operation reads no sequence/truth.
 * @author Raiden
 */
public final class BindOrdinalNcrnaCallsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, BindOrdinalNcrnaCallsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture option: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null, "Supply out=fresh_fixture_directory");
			check(new File(parser.out1).mkdir(), "Fixture output must be fresh");
			new BindOrdinalNcrnaCallsTest(parser.out1).run();
			pp.outstream.println("ORDINAL_NCRNA_CALLS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private BindOrdinalNcrnaCallsTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		write("sources.tsv", "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n"
			+"0\tassembly_A\t10\t/unavailable/a.fa\t/unavailable/10.gff\n"
			+"1\tassembly_B\t11\t/unavailable/b.fa\t/unavailable/11.gff\n"
			+"2\tassembly_C\t12\t/unavailable/c.fa\t/unavailable/12.gff\n"
			+"3\tassembly_D\t13\t/unavailable/d.fa\t/unavailable/13.gff\n");
		final String first="same\tCallGenes\tRNA\t2\t9\t.\t+\t.\tmodel:cobalamin_RF00174_consensus_0,family:cobalamin_RF00174\n";
		final String second="same\tCallGenes\tRNA\t2\t9\t.\t-\t.\tmodel:cobalamin_RF01482_consensus_1,family:cobalamin_RF01482\n";
		final String cds="same\tCallGenes\tCDS\t15\t20\t.\t+\t0\tID:kept\n";
		final String raw="##gff-version 3\n##input_ordinal0\t0\n"+first+"##input_ordinal0\t1\n##input_ordinal0\t2\n"
			+second+cds+"##input_ordinal0\t3\n";
		write("calls.gff", raw);
		final String config="in=/unavailable/a.fa\nin=/unavailable/b.fa\nin=/unavailable/c.fa\nin=/unavailable/d.fa\nout="+path("calls.gff")+"\ngffinputordinal=t\n";
		write("caller.config", config);
		invoke("valid", "calls.gff", "caller.config");
		check(read("valid/complete.tsv").equals("#ORDINAL_NCRNA_CALL_BINDING_V2\ninputs\tfeatures\tsource_gff\tqualification\n4\t3\t"
			+path("calls.gff")+"\tDECLARED_CONFIG_AND_ORDINAL_BINDINGS_NOT_INDEPENDENT_PRODUCER_PROOF\n"),
			"Completion binds actual source path and counts without an unverified hash declaration");
		boolean obsoleteRejected=false;
		try{
			BindOrdinalNcrnaCalls.main(new String[]{"in="+path("calls.gff"), "bindings="+path("sources.tsv"), "callerconfig="+path("caller.config"),
				"inputsha80=0123456789abcdef0123", "expectedinputs=4", "out="+path("obsolete_option")});
		}catch(IllegalArgumentException expected){obsoleteRejected=expected.getMessage().contains("Unknown ordinal-binding option: inputsha80");}
		check(obsoleteRejected && !new File(path("obsolete_option")).exists(), "Removed provenance-only options must fail loudly before creating binding outputs");
		check(read("valid/queries/0.gff").equals("##gff-version 3\n"+first), "First input must retain exact feature bytes");
		check(read("valid/queries/2.gff").equals("##gff-version 3\n"+second+cds), "Repeated contig names and non-RNA rows keep their own source and order");
		check(read("valid/queries/1.gff").equals("##gff-version 3\n") && read("valid/queries/3.gff").equals("##gff-version 3\n"), "Empty middle and terminal inputs need explicit query files");
		check(read("valid/calls.tsv").equals("input_id\tfasta\tgff\n"
			+"10\t/unavailable/a.fa\t"+path("valid/queries/0.gff")+'\n'
			+"11\t/unavailable/b.fa\t"+path("valid/queries/1.gff")+'\n'
			+"12\t/unavailable/c.fa\t"+path("valid/queries/2.gff")+'\n'
			+"13\t/unavailable/d.fa\t"+path("valid/queries/3.gff")+'\n'), "Global reference IDs must not become ordinals or assembly names");
		check(read("valid/source_lines.tsv").equals("input_ordinal0\tinput_id\tquery_line\tsource_line\n0\t10\t2\t3\n2\t12\t2\t6\n2\t12\t3\t7\n"), "Every projected feature must retain its independent original line");
		write("skipped.gff", raw.replace("##input_ordinal0\t2", "##input_ordinal0\t3"));
		write("skipped.config", config.replace(path("calls.gff"), path("skipped.gff")));
		reject("skipped", "skipped.gff", "skipped.config", "contiguous");
		write("missing.gff", raw.replace("##input_ordinal0\t3\n", ""));
		write("missing.config", config.replace(path("calls.gff"), path("missing.gff")));
		reject("missing", "missing.gff", "missing.config", "Missing input boundary");
		write("early.gff", first+raw);
		write("early.config", config.replace(path("calls.gff"), path("early.gff")));
		reject("early", "early.gff", "early.config", "precedes");
		write("wrong.config", config.replace("in=/unavailable/a.fa", "in=/unavailable/c.fa"));
		reject("wrong", "calls.gff", "wrong.config", "input order");
		reject("output", "skipped.gff", "caller.config", "exact saved GFF");
		for(final String marker:new String[]{"-1", "4", "0"}){
			final String name="bad_marker_"+marker;
			write(name+".gff", raw.replace("##input_ordinal0\t2", "##input_ordinal0\t"+marker));
			write(name+".config", config.replace(path("calls.gff"), path(name+".gff")));
			reject(name, name+".gff", name+".config", "contiguous");
		}
		final String sources=read("sources.tsv");
		write("duplicate_id.tsv", sources.replace("2\tassembly_C\t12", "2\tassembly_C\t10"));
		rejectBindings("duplicate_id", "duplicate_id.tsv", "ordinal=2 input_id=10");
		write("missing_binding.tsv", sources.replace("2\tassembly_C\t12\t/unavailable/c.fa\t/unavailable/12.gff\n", ""));
		rejectBindings("missing_binding", "missing_binding.tsv", "ordinal=3 input_id=13 expected_ordinal=2");
		write("missing_last_binding.tsv", sources.substring(0, sources.indexOf("3\tassembly_D")));
		rejectBindings("missing_last_binding", "missing_last_binding.tsv", "next_ordinal=3 input_id=UNMAPPED");
		check(read("valid/references.tsv").equals("input_id\tfasta\tgff\n"
			+"10\t/unavailable/a.fa\t/unavailable/10.gff\n11\t/unavailable/b.fa\t/unavailable/11.gff\n"
			+"12\t/unavailable/c.fa\t/unavailable/12.gff\n13\t/unavailable/d.fa\t/unavailable/13.gff\n"),
			"Every global reference ID must be reached exactly once, including empty inputs");
		compressed();
	}

	/** Actual public consumers read compressed query views; no FASTA is created or read. */
	private void compressed() throws Exception{
		final String first="same\tBBTools\tRNA\t2\t9\t.\t+\t0\tRNA,model:A_c1,family:A\n";
		final String second="same\tBBTools\tRNA\t2\t9\t.\t-\t0\tRNA,model:B_c1,family:B\n";
		final String cds="same\tBBTools\tCDS\t15\t20\t.\t+\t0\tCDS,startScr:1\n";
		final String trna="same\tBBTools\ttRNA\t30\t40\t.\t+\t0\ttRNA,startScr:0\n";
		write("interop.gff", "##gff-version 3\n##input_ordinal0\t0\n"+first
			+"##input_ordinal0\t1\n##input_ordinal0\t2\n"+second+cds+trna+"##input_ordinal0\t3\n");
		write("interop.config", read("caller.config").replace(path("calls.gff"), path("interop.gff")));
		final boolean bgzip=ReadWrite.USE_BGZIP;
		for(final String arm:new String[]{"plain", "compressed"}){
			BindOrdinalNcrnaCalls.main(new String[]{"in="+path("interop.gff"), "bindings="+path("sources.tsv"),
				"callerconfig="+path("interop.config"), "expectedinputs=4", "out="+path(arm), "compressed="+arm.equals("compressed")});
			check(ReadWrite.USE_BGZIP==bgzip, "The binding tool must restore the caller's compressor preference after its output scope");
		}
		boolean failed=false;
		try{
			BindOrdinalNcrnaCalls.main(new String[]{"in="+path("skipped.gff"), "bindings="+path("sources.tsv"),
				"callerconfig="+path("skipped.config"), "expectedinputs=4", "out="+path("compressed_failure"), "compressed=t"});
		}catch(IllegalArgumentException expected){failed=expected.getMessage().contains("contiguous");}
		check(failed && ReadWrite.USE_BGZIP==bgzip && !new File(path("compressed_failure/complete.tsv")).exists(),
			"A malformed compressed input must fail without completion and restore the compressor preference after closing partial outputs");
		for(final String force:new String[]{"forcebgzip", "forcepigz"}){rejectCompressor(force);}
		check(readNative("compressed/queries/0.gff.gz").equals("##gff-version 3\n"+first)
			&& readNative("compressed/queries/2.gff.gz").equals("##gff-version 3\n"+second+cds+trna),
			"Compressed query bytes preserve repeated contigs, both strands and non-RNA features");
		check(readNative("compressed/queries/1.gff.gz").equals("##gff-version 3\n")
			&& readNative("compressed/queries/3.gff.gz").equals("##gff-version 3\n"),
			"Compressed empty middle and final inputs remain explicit gzip query files");
		check(!new File(path("compressed/source_lines.tsv")).exists()
			&& !new File(path("compressed/queries/0.gff")).exists(), "Compression must not leave expanded duplicate payloads");
		check(readNative("compressed/source_lines.tsv.gz").equals("input_ordinal0\tinput_id\tquery_line\tsource_line\n"
			+"0\t10\t2\t3\n2\t12\t2\t6\n2\t12\t3\t7\n2\t12\t4\t8\n"),
			"Compressed source-line ledger retains every hand-counted feature origin");
		check(read("plain/queries/0.gff").equals("##gff-version 3\n"+first)
			&& read("plain/source_lines.tsv").equals(readNative("compressed/source_lines.tsv.gz")),
			"Explicit compressed=false preserves the original plain output contract");
		write("interop.families.tsv", GradeNcrnaReference.RULE_HEADER+"\nA\tA\tA\ttrue\nB\tB\tB\ttrue\n");
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final String[] fasta={"/unavailable/a.fa", "/unavailable/b.fa", "/unavailable/c.fa", "/unavailable/d.fa"};
		final ByteBuilder refs=new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl();
		final ByteBuilder corpus=new ByteBuilder("assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role\n");
		for(int i=0; i<4; i++){
			final String id=""+(10+i), ref=path("reference_"+i+".gff"), stem="original_"+i;
			final Row row=i==0 || i==2 ? codec.create(id, new Locus("same", i==0 ? "A" : "B", 2, 9,
				i==0 ? '+' : '-'), false, null, null, "fixture:compressed_binding") : null;
			codec.write(ref, id, row==null ? Collections.<Row>emptyList() : Collections.singletonList(row));
			refs.append(id).tab().append(fasta[i]).tab().append(ref).nl();
			corpus.append("assembly_").append(i).tab().append(id).append("\tfixture\t").append(fasta[i])
				.append("\t/new.gff\t\t\t").append(path(stem+".fna.gz")).append("\t/new.gff\tfasta_gzip\tncrna\tcore\n");
			write(stem+".gff.gz", "##gff-version 3\n##sequence-region same 1 100\n"
				+"same\tNCBI\tregion\t1\t100\t.\t+\t.\tID=region;Is_circular=false\n"
				+(i==2 ? "same\tNCBI\tCDS\t15\t20\t.\t+\t0\tID=only\n" : ""));
			check(!new File(path(stem+".fna.gz")).exists(), "Consumer fixture must not need a sequence payload");
		}
		write("interop.references.tsv", refs.toString());
		write("interop.corpus.tsv", corpus.toString());
		GradeNcrnaReference.main(new String[]{"in="+path("compressed/calls.tsv"), "reference="+path("interop.references.tsv"),
			"familymap="+path("interop.families.tsv"), "out="+path("interop_grade"), "t=1"});
		final LineParser1 p=new LineParser1('\t');
		int families=0;
		try(Input input=new Input(path("interop_grade.summary.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check((p.termEquals("A", 0) || p.termEquals("B", 0)) && p.parseInt(1)==4
					&& p.parseInt(2)==1 && p.parseInt(3)==0 && p.parseInt(4)==0
					&& p.parseInt(10)==1 && p.parseInt(11)==1 && p.parseInt(12)==1,
					"Actual RNA grader must recover one exact TP per family from compressed queries with empty inputs");
				families++;
			}
		}
		check(families==2, "The two independent reference families must both be graded");
		GradeCdsStartChanges.main(new String[]{"in="+path("plain/calls.tsv"), "in2="+path("compressed/calls.tsv"),
			"corpus="+path("interop.corpus.tsv"), "expectedinputs=4", "out="+path("interop_cds"),
			"rnacallchanges=t", "familymap="+path("interop.families.tsv"), "t=1"});
		int inputRows=0;
		try(Input input=new Input(path("interop_cds.inputs.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.parseLong(2)==(p.termEquals("12", 0) ? 1 : 0) && p.parseLong(3)==0 && p.parseLong(4)==0,
					"Actual collateral reader must preserve the one CDS TP and zero FP/FN in both plain and compressed arms");
				inputRows++;
			}
		}
		check(inputRows==8, "All four before/after input pairs, including empty inputs, must survive compressed collateral grading");
		int rnaRows=0;
		try(Input input=new Input(path("interop_cds.rna_call_changes.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				final int expected=p.termEquals("12", 0) ? 1 : 0;
				check(p.parseInt(1)==expected && p.parseInt(2)==expected && p.parseInt(3)==0 && p.parseInt(4)==0,
					"The actual other-RNA comparison must retain the tRNA with no loss/gain from compression");
				rnaRows++;
			}
		}
		check(rnaRows==4, "Other-RNA comparison retains every declared input");
		System.err.println("ORDINAL_COMPRESSED_CONSUMERS_TEST_PASS");
	}

	/** Conflicting CLI flags fail before any output, even if the compressor is installed. */
	private void rejectCompressor(final String force){
		final boolean bgzip=ReadWrite.USE_BGZIP, pigz=ReadWrite.USE_PIGZ;
		final boolean forceBgzip=ReadWrite.FORCE_BGZIP, forcePigz=ReadWrite.FORCE_PIGZ;
		boolean failed=false;
		try{
			BindOrdinalNcrnaCalls.main(new String[]{"in="+path("interop.gff"), "bindings="+path("sources.tsv"),
				"callerconfig="+path("interop.config"), "expectedinputs=4", "out="+path(force), "compressed=t", force+"=t"});
		}catch(IllegalArgumentException expected){failed=expected.getMessage().contains("forcebgzip/forcepigz");}
		finally{
			ReadWrite.USE_BGZIP=bgzip;
			ReadWrite.USE_PIGZ=pigz;
			ReadWrite.FORCE_BGZIP=forceBgzip;
			ReadWrite.FORCE_PIGZ=forcePigz;
		}
		check(failed && !new File(path(force)).exists(), "Conflicting "+force+" must fail before creating compressed binding outputs");
	}

	private void invoke(final String output, final String gff, final String config){
		invoke(output, gff, config, "sources.tsv");
	}
	private void invoke(final String output, final String gff, final String config, final String bindings){
		BindOrdinalNcrnaCalls.main(new String[]{"in="+path(gff), "bindings="+path(bindings), "callerconfig="+path(config),
			"expectedinputs=4", "out="+path(output)});
	}
	private void rejectBindings(final String output, final String bindings, final String reason){
		boolean failed=false;
		try{invoke(output, "calls.gff", "caller.config", bindings);}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains(reason);}
		check(failed && !new File(path(output)).exists(), "Invalid panel must fail before attribution/output: "+output);
	}
	private void reject(final String output, final String gff, final String config, final String reason){
		boolean failed=false;
		try{invoke(output, gff, config);}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains(reason);}
		check(failed && !new File(path(output+"/complete.tsv")).exists(), "Malformed binding must fail without completion: "+output);
	}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String readNative(final String name){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path(name))){
			for(byte[] line; (line=input.reader.nextLine())!=null;){text.append(line).nl();}
		}
		return text.toString();
	}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean condition, final String reason){if(!condition){throw new AssertionError(reason);}}
	private final String root;
}
