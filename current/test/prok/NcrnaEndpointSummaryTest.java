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
		equal(prefix+".radii.tsv", "family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n"
			+"T\treciprocal_half\t0\t6\t3\t2\t2\nT\treciprocal_half\t2\t6\t5\t3\t2\n"
			+"T\treciprocal_half\t3\t6\t5\t5\t4\nT\treciprocal_half\t5\t6\t6\t6\t6\n");
		equal(prefix+".quantiles.tsv", "family\tend\tmeasure\tmatched\tminimum\tp50\tp90\tp95\tp99\tmaximum\n"
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
