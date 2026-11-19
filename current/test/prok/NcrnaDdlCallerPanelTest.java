package test.prok;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlCallerPanel;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal populations with repeated source slots, zero-target inputs, unrankable
 * parents, distinct sequence/source states and nonexistent genome paths.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlCallerPanelTest{

	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, NcrnaDdlCallerPanelTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown caller-panel fixture argument: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.in1==null && parser.in2==null
				&& parser.out2==null && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Caller-panel fixture needs assertions and out=fresh_directory without input/sampling");
			final Path root=Paths.get(parser.out1).toAbsolutePath();
			Files.createDirectory(root);
			run(root);
			pp.outstream.println("NCRNA_DDL_CALLER_PANEL_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final Path root) throws IOException{
		final String capture=root.resolve("capture.tsv").toString();
		final String combined=root.resolve("combined.tsv").toString();
		final String oldAudit=root.resolve("old.tsv").toString();
		final String newAudit=root.resolve("new.tsv").toString();
		final String unicode="/nonexistent/紫.fa";
		write(capture, NcrnaDdlCallerPanel.CAPTURE_HEADER,
			"0\t11\t/nonexistent/a.fa\t1\t300\t1\t2\t2\t0\t0\t0\t0\t0",
			"1\t12\t/nonexistent/zero.fa\t1\t300\t0\t0\t0\t0\t0\t0\t0\t0",
			"2\t13\t"+unicode+"\t1\t300\t0\t0\t1\t0\t0\t0\t0\t0");
		write(combined, NcrnaDdlCallerPanel.COMBINED_HEADER,
			"reused\t1\t12\t/nonexistent/zero.fa\t102\t0\t0\t0\t0\t0\t0",
			"reused\t0\t11\t/nonexistent/a.fa\t101\t2\t1\t0\t1\t0\t0",
			"new\t0\t21\t/nonexistent/new.fa\t201\t1\t0\t0\t0\t1\t0");
		final String seen=parent(0, 11, "/nonexistent/a.fa", "p11a", true, 101, true, "SEEN", "NOT_SEEN", "SEEN");
		final String notSeen=parent(0, 11, "/nonexistent/a.fa", "p11b", false, 101, true, "NOT_SEEN", "SEEN", "SEEN");
		final String unknown=parent(2, 13, unicode, "p13", false, 103, false, "UNKNOWN", "UNKNOWN", "UNKNOWN");
		write(oldAudit, NcrnaDdlCallerPanel.AUDIT_HEADER, seen, notSeen, unknown);
		write(newAudit, NcrnaDdlCallerPanel.AUDIT_HEADER,
			parent(0, 21, "/nonexistent/new.fa", "p21", false, 201, true, "NOT_SEEN", "NOT_SEEN", "NOT_SEEN"));

		final String original=root.resolve("original").toString();
		invoke("capture", capture, oldAudit, original, 3);
		equal(original+".config", "in=/nonexistent/a.fa\nin=/nonexistent/zero.fa\nin="+unicode+"\n");
		equal(original+".strata.tsv", NcrnaDdlCallerPanel.STRATA_HEADER+"\n"
			+"11\tp11a\tSEEN\tNOT_SEEN\ttrue\ttrue\n"
			+"11\tp11b\tNOT_SEEN\tSEEN\ttrue\tfalse\n"
			+"13\tp13\tUNKNOWN\tUNKNOWN\tfalse\tfalse\n");
		equal(original+".summary.tsv", "inputs\ttargets\texact_seen\texact_not_seen\texact_unknown\n3\t3\t1\t1\t1\n");
		final String expanded=root.resolve("expanded").toString();
		invoke("combined", combined, oldAudit+","+newAudit, expanded, 3);
		equal(expanded+".config", "in=/nonexistent/zero.fa\nin=/nonexistent/a.fa\nin=/nonexistent/new.fa\n");
		equal(expanded+".inputs.tsv", NcrnaDdlCallerPanel.INPUT_HEADER+"\n"
			+"0\t12\t/nonexistent/zero.fa\treused\t1\t102\t0\n"
			+"1\t11\t/nonexistent/a.fa\treused\t0\t101\t2\n"
			+"2\t21\t/nonexistent/new.fa\tnew\t0\t201\t1\n");
		equal(expanded+".summary.tsv", "inputs\ttargets\texact_seen\texact_not_seen\texact_unknown\n3\t3\t1\t2\t0\n");

		final String missing=root.resolve("missing.tsv").toString();
		write(missing, NcrnaDdlCallerPanel.AUDIT_HEADER, seen, unknown);
		expect("capture", capture, missing, root.resolve("missing_out").toString(), "Missing or duplicated exact-sequence labels");
		final String duplicate=root.resolve("duplicate.tsv").toString();
		write(duplicate, NcrnaDdlCallerPanel.AUDIT_HEADER, seen, notSeen, unknown, notSeen);
		expect("capture", capture, duplicate, root.resolve("duplicate_out").toString(), "exactly one label");
		final String wrong=root.resolve("wrong.tsv").toString();
		write(wrong, NcrnaDdlCallerPanel.AUDIT_HEADER,
			parent(1, 11, "/nonexistent/a.fa", "p11a", true, 101, true, "SEEN", "SEEN", "SEEN"), notSeen, unknown);
		expect("capture", capture, wrong, root.resolve("wrong_out").toString(), "Audit slots/paths");
		final String badGroup=root.resolve("bad_group.tsv").toString();
		write(badGroup, NcrnaDdlCallerPanel.COMBINED_HEADER,
			"reused\t0\t11\t/nonexistent/a.fa\t101\t2\t1\t0\t1\t0\t0",
			"reused\t1\t12\t/nonexistent/zero.fa\t101\t0\t0\t0\t0\t0\t0",
			"new\t0\t21\t/nonexistent/new.fa\t201\t1\t0\t0\t0\t1\t0");
		expect("combined", badGroup, oldAudit+","+newAudit, root.resolve("bad_group_out").toString(), "one representative per resolved assembly group");
	}

	private static void invoke(final String mode, final String input, final String audits, final String output, final int inputs){
		NcrnaDdlCallerPanel.main(new String[]{"mode="+mode, "in="+input, "audits="+audits,
			"out="+output, "expectedinputs="+inputs, "t=1"});
	}
	private static void expect(final String mode, final String input, final String audits, final String output, final String message){
		try{invoke(mode, input, audits, output, 3);}
		catch(RuntimeException | AssertionError failure){
			check(failure.getMessage()!=null && failure.getMessage().contains(message), "Fixture failed for the wrong reason: "+failure);
			return;
		}
		throw new AssertionError("Expected malformed caller-panel binding to fail: "+message);
	}
	private static String parent(final int slot, final int id, final String fasta, final String parent,
			final boolean common, final int group, final boolean complete, final String exact, final String source, final String composite){
		check(slot>=0 && id>=0 && group>=0, "Literal fixture parent needs valid original bindings");
		return slot+"\t"+id+"\t"+fasta+"\t"+parent+"\tcontig\t10\t50\t+\tfalse\ttrue\tfalse\tfalse\tRF00174\t.\t"
			+common+"\t"+group+"\t"+complete+"\t"+exact+"\t"+source+"\t"+composite+"\t.";
	}
	private static void write(final String path, final String... lines){
		check(lines.length>0, "Every fixture metadata file needs its schema header");
		try(Output out=new Output(path)){
			for(final String line:lines){out.writer.println(line.getBytes(StandardCharsets.UTF_8));}
		}
	}
	private static void equal(final String path, final String expected) throws IOException{
		check(Arrays.equals(Files.readAllBytes(Paths.get(path)), expected.getBytes(StandardCharsets.UTF_8)), "Literal caller-panel output differs: "+path);
	}
	private static void check(final boolean valid, final String message){
		if(!valid){throw new AssertionError(message);}
	}
}
