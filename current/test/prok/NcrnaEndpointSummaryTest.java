package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;

import parse.Parser;
import parse.Parse;
import parse.PreParser;
import prok.NcrnaEndpointSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Finite signed-offset, joint-radius and malformed-report controls.
 * Expected tables are literal independent arithmetic for six paired loci.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointSummaryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointSummaryTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(parser.out1==null){throw new IllegalArgumentException("Specify out=fresh_directory");}
			final File root=new File(parser.out1);
			check(root.mkdir(), "Test output directory must be fresh");
			run(root);
			System.out.println("NCRNA_ENDPOINT_SUMMARY_TEST_PASS signed_strands joint_radii quantiles unique_pairs malformed_offsets overlap_scope empty_family");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String header="#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\n";
		final String rows="T\ta\t+\t100\t120\t100\t120\t0\t0\n"
			+"T\tb\t+\t200\t220\t202\t217\t2\t-3\n"
			+"T\tc\t-\t300\t320\t305\t318\t2\t-5\n"
			+"T\td\t+\t400\t420\t396\t421\t-4\t1\n"
			+"T\te\t-\t500\t520\t500\t520\t0\t0\n"
			+"T\tf\t+\t600\t620\t600\t617\t0\t-3\n"
			+"OTHER\tx\t+\t700\t720\t700\t720\t0\t0\n";
		final File input=new File(root, "matches.tsv");
		write(input, header+rows);
		final String prefix=new File(root, "result").getPath();
		NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+prefix, "family=T", "overlap=reciprocal_half",
			"radii=5,0,3,2", "expected=6", "exact5=3", "exact3=2", "both=2"});
		equal(prefix+".summary.tsv", "family\toverlap\tmatched\texact5\texact3\tboth_exact\nT\treciprocal_half\t6\t3\t2\t2\n");
		equal(prefix+".histogram.tsv", "family\tend\tmeasure\toffset_bp\tloci\n"
			+"T\t5prime\tsigned\t-4\t1\nT\t5prime\tsigned\t0\t3\nT\t5prime\tsigned\t2\t2\n"
			+"T\t5prime\tabsolute\t0\t3\nT\t5prime\tabsolute\t2\t2\nT\t5prime\tabsolute\t4\t1\n"
			+"T\t3prime\tsigned\t-5\t1\nT\t3prime\tsigned\t-3\t2\nT\t3prime\tsigned\t0\t2\nT\t3prime\tsigned\t1\t1\n"
			+"T\t3prime\tabsolute\t0\t2\nT\t3prime\tabsolute\t1\t1\nT\t3prime\tabsolute\t3\t2\nT\t3prime\tabsolute\t5\t1\n");
		equal(prefix+".radii.tsv", "family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n"
			+"T\treciprocal_half\t0\t6\t3\t2\t2\nT\treciprocal_half\t2\t6\t5\t3\t2\n"
			+"T\treciprocal_half\t3\t6\t5\t5\t4\nT\treciprocal_half\t5\t6\t6\t6\t6\n");
		equal(prefix+".quantiles.tsv", "#NCRNA_ENDPOINT_QUANTILES_V2\nfamily\tend\tmeasure\tmatched\tminimum\tp50\tp90\tp95\tp99\tmaximum\n"
			+"T\t5prime\tsigned\t6\t-4\t0\t2\t2\t2\t2\nT\t5prime\tabsolute\t6\t0\t0\t4\t4\t4\t4\n"
			+"T\t3prime\tsigned\t6\t-5\t-3\t1\t1\t1\t1\nT\t3prime\tabsolute\t6\t0\t1\t5\t5\t5\t5\n");
		failure(root, "sign", header+rows.replace("318\t2\t-5", "318\t-2\t-5"), "Stored offsets disagree", "reciprocal_half");
		failure(root, "duplicate", header+rows+"T\tb\t+\t201\t220\t202\t217\t1\t-3\n", "consumed more than once", "reciprocal_half");
		failure(root, "overflow", header+rows.replace("200\t220", "4294967496\t220"), "Invalid decimal field", "reciprocal_half");
		final String contained="T\tz\t+\t100\t199\t100\t119\t0\t-80\n";
		failure(root, "overlap", header+contained, "declared overlap rule", "reciprocal_half");
		final File shortInput=new File(root, "shorter.tsv");
		write(shortInput, header+contained);
		NcrnaEndpointSummary.main(new String[]{"in="+shortInput, "out="+new File(root, "shorter"),
			"family=T", "overlap=half_shorter", "expected=1"});
		NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+new File(root, "empty"),
			"family=absent", "overlap=reciprocal_half", "expected=0"});
		equal(new File(root, "empty.summary.tsv").getPath(),
			"family\toverlap\tmatched\texact5\texact3\tboth_exact\nabsent\treciprocal_half\t0\t0\t0\t0\n");
		equal(new File(root, "empty.quantiles.tsv").getPath(),
			"#NCRNA_ENDPOINT_QUANTILES_V2\nfamily\tend\tmeasure\tmatched\tminimum\tp50\tp90\tp95\tp99\tmaximum\n");
		filtered(root, header);
	}

	/** Source IDs differ from ordinals; two sources reuse the same contig/interval. */
	private static void filtered(final File root, final String header){
		final String bindingText="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n"
			+"0\tA\t77\t/a.fa\t/a.gff\n1\tB\t9\t/b.fa\t/b.gff\n2\tC\tu/v:x\t/c.fa\t/c.gff\n";
		final String auditText="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result\n"
			+"panel\trow1\tfasta\tassembly\tstatus\n"
			+"/panels/new1000.txt\t1\t/a.fa\tA\tEXCLUDED_ASSEMBLY\n"
			+"/panels/old300.txt\t1\t/a.fa\tA\tEXCLUDED_ASSEMBLY\n";
		final String rows="T\t2:77/ctg\t+\t100\t120\t100\t120\t0\t0\n"
			+"T\t1:9/ctg\t+\t100\t120\t102\t117\t2\t-3\n"
			+"T\t5:u/v:x/ctg\t-\t300\t320\t299\t324\t-4\t1\n"
			+"clan\t2:77/ctg\t+\t100\t120\t100\t120\t0\t0\n";
		final File bindings=new File(root, "sources.tsv"), audit=new File(root, "exclude.tsv"), input=new File(root, "qualified.tsv");
		write(bindings, bindingText);
		write(audit, auditText);
		write(input, header+rows);
		final String prefix=new File(root, "filtered").getPath();
		NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+prefix, "family=T", "overlap=reciprocal_half",
			"bindings="+bindings, "exclude="+audit, "radii=0,2,3,4", "expectedall=3", "expected=2", "exact5=0", "exact3=0", "both=0"});
		equal(prefix+".summary.tsv", "family\toverlap\tmatched\texact5\texact3\tboth_exact\nT\treciprocal_half\t2\t0\t0\t0\n");
		equal(prefix+".population.tsv", "family\tselected\tretained\texcluded\tbinding_sources\texcluded_assemblies\taudit_panels\taudit_rows\tunmatched_alias_rows\n"
			+"T\t3\t2\t1\t3\t1\t2\t2\t0\n");
		equal(prefix+".radii.tsv", "family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n"
			+"T\treciprocal_half\t0\t2\t0\t0\t0\nT\treciprocal_half\t2\t2\t1\t1\t0\n"
			+"T\treciprocal_half\t3\t2\t1\t2\t1\nT\treciprocal_half\t4\t2\t2\t2\t2\n");
		equal(prefix+".histogram.tsv", "family\tend\tmeasure\toffset_bp\tloci\n"
			+"T\t5prime\tsigned\t-4\t1\nT\t5prime\tsigned\t2\t1\n"
			+"T\t5prime\tabsolute\t2\t1\nT\t5prime\tabsolute\t4\t1\n"
			+"T\t3prime\tsigned\t-3\t1\nT\t3prime\tsigned\t1\t1\n"
			+"T\t3prime\tabsolute\t1\t1\nT\t3prime\tabsolute\t3\t1\n");
		filteredFailure(root, "length", header+rows.replace("1:9/ctg", "9:9/ctg"), bindings, audit, "Declared source-ID length");
		filteredFailure(root, "unknown", header+rows.replace("1:9/ctg", "7:unknown/ctg"), bindings, audit, "outside the declared source binding");
		filteredFailure(root, "excluded_bad_sign", header+rows.replace("100\t120\t100\t120\t0\t0", "100\t120\t100\t120\t1\t0"),
			bindings, audit, "Stored offsets disagree");
		final File conflict=new File(root, "conflict.sources.tsv");
		write(conflict, bindingText.replace("1\tB\t9\t", "1\tB\t77\t"));
		filteredFailure(root, "conflict", header+rows, conflict, audit, "unambiguous");
		final File badAudit=new File(root, "bad.exclude.tsv");
		write(badAudit, auditText.replace("old300.txt\t1", "old300.txt\t2"));
		filteredFailure(root, "audit_index", header+rows, bindings, badAudit, "row indices");
		boolean rejected=false;
		try{
			NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+new File(root, "unpaired"), "family=T",
				"overlap=reciprocal_half", "bindings="+bindings});
		}catch(IllegalArgumentException expected){
			rejected=expected.getMessage().contains("together");
		}
		check(rejected, "One source-filter path alone must not silently disable exclusion");
	}

	/** Bad selected rows remain invalid even when their assembly would be excluded. */
	private static void filteredFailure(final File root, final String name, final String data,
			final File bindings, final File audit, final String reason){
		final File input=new File(root, name+".qualified.tsv"), prefix=new File(root, name+".qualified");
		write(input, data);
		boolean rejected=false;
		try{
			NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+prefix, "family=T", "overlap=reciprocal_half",
				"bindings="+bindings, "exclude="+audit});
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Unexpected source-filter failure: "+expected);
			rejected=true;
		}
		check(rejected && !new File(prefix+".summary.tsv").exists(), "Invalid source data must fail before creating a successful summary");
	}

	/** Rejection must name the intended defect; unrelated failures cannot satisfy it. */
	private static void failure(final File root, final String name, final String data,
			final String reason, final String overlap){
		final File input=new File(root, name+".input.tsv");
		write(input, data);
		boolean rejected=false;
		try{
			NcrnaEndpointSummary.main(new String[]{"in="+input, "out="+new File(root, name), "family=T", "overlap="+overlap});
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Unexpected failure: "+expected);
			rejected=true;
		}
		check(rejected, "Malformed matched report must fail: "+name);
	}

	private static void write(final File path, final String data){
		try(Output output=new Output(path.getPath())){output.writer.print(data);}
	}

	private static void equal(final String path, final String expected){
		final ByteBuilder actual=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){actual.append(line).nl();}
		}
		check(java.util.Arrays.equals(actual.toBytes(), expected.getBytes(StandardCharsets.UTF_8)), "Literal expected table differs: "+path);
	}

	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}
}
