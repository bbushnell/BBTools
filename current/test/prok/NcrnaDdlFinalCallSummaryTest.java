package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlFinalCallSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal saved-GFF checks: final rows, duplicate retention and absent metadata.
 * No caller, sequence, matching or timing measurement.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlFinalCallSummaryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlFinalCallSummaryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture flag: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			check(new File(parser.out1).mkdir(), "Fixture directory must be fresh");
			run(parser.out1);
			pp.outstream.println("NCRNA_DDL_FINAL_CALL_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final String target=feature("cobalamin_RF00174_first", "");
		final String ordinary=root+"/ordinary.gff";
		write(ordinary, "##gff-version 3\n\n"+feature("cobalamin_RF01482_other", "")+
			"chr\tBBTools\tCDS\t1\t20\t.\t+\t0\tCDS\n"+target+target+feature("cobalamin_RF001740_lookalike", ""));
		summarize(ordinary, root+"/ordinary", "public");
		final ArrayList<String> calls=rows(root+"/ordinary.calls.tsv");
		check(calls.size()==3 && calls.get(1).equals("5\tcobalamin_RF00174_first\tunclassified\tNA\tNA\tchr\t1\t20\t+")
			&& calls.get(2).equals("6\tcobalamin_RF00174_first\tunclassified\tNA\tNA\tchr\t1\t20\t+"),
			"Comments/blanks/non-target rows retain original physical-line keys; duplicate target predictions remain separate");
		check(rows(root+"/ordinary.summary.tsv").get(1).startsWith("unclassified\tNA\t2\t"),
			"Public final calls have unavailable provenance, not invented native support zero");

		final String hybrid=root+"/hybrid.gff";
		write(hybrid, feature("cobalamin_RF00174_a", ";ddl_stream=native;ddl_native_shared=1")+
			feature("cobalamin_RF00174_b", ";ddl_stream=native;ddl_native_shared=2")+
			feature("cobalamin_RF00174_c", ";ddl_stream=sketch;ddl_native_shared=0")+
			feature("cobalamin_RF00174_d", ";ddl_stream=sketch;ddl_native_shared=1")+
			feature("cobalamin_RF00174_e", ";ddl_native_shared=2;ddl_stream=sketch")+
			feature("cobalamin_RF00174_c", ";ddl_stream=sketch;ddl_native_shared=0"));
		summarize(hybrid, root+"/hybrid", "hybrid");
		final ArrayList<String> summary=rows(root+"/hybrid.summary.tsv");
		check(summary.size()==5 && summary.get(1).startsWith("native\tBELOW_SHIPPED_FLOOR_2\t1\t")
			&& summary.get(2).startsWith("native\tAT_LEAST_SHIPPED_FLOOR_2\t1\t")
			&& summary.get(3).startsWith("sketch\tBELOW_SHIPPED_FLOOR_2\t3\t")
			&& summary.get(4).startsWith("sketch\tAT_LEAST_SHIPPED_FLOOR_2\t1\t")
			&& rows(root+"/hybrid.calls.tsv").size()==7,
			"Actual final rows conserve native/sketch proposals on both sides of floor2, including the duplicate sketch-zero row");

		final String sketch=root+"/sketch.gff";
		write(sketch, feature("cobalamin_RF00174_s", ";partial=11;ddl_stream=sketch_only;ddl_native_shared=NA"));
		summarize(sketch, root+"/sketch", "sketch");
		check(rows(root+"/sketch.summary.tsv").get(1).startsWith("sketch_only\tNA\t1\t"),
			"Partial sketch-only calls keep unmeasured native support NA");
		final String empty=root+"/empty.gff";
		write(empty, "##gff-version 3\n"+feature("cobalamin_RF01482_other", ""));
		summarize(empty, root+"/empty", "hybrid");
		check(rows(root+"/empty.calls.tsv").size()==1 && rows(root+"/empty.summary.tsv").size()==1,
			"A real zero-target population writes headers without fabricating one call");

		final String[] modes={"public", "hybrid", "sketch", "hybrid", "sketch", "hybrid", "hybrid", "hybrid", "hybrid", "hybrid", "hybrid"};
		final String[] bad={
			";ddl_stream=native;ddl_native_shared=2", "", ";ddl_stream=sketch_only;ddl_native_shared=0",
			";ddl_stream=sketch_only;ddl_native_shared=NA", ";ddl_stream=native;ddl_native_shared=2",
			";ddl_stream=unknown;ddl_native_shared=2", ";ddl_stream=native;ddl_native_shared=NA",
			";ddl_stream=native;ddl_native_shared=65536", ";ddl_stream=native;ddl_native_shared=-1",
			";ddl_stream=native;ddl_stream=sketch;ddl_native_shared=2", ";model:second"};
		for(int i=0; i<bad.length; i++){
			final String path=root+"/bad_"+i+".gff";
			write(path, feature("cobalamin_RF00174_bad", bad[i]));
			reject(path, root+"/bad_"+i, modes[i]);
		}
		reject(ordinary, root+"/ordinary", "public");
	}

	private static String feature(final String model, final String extra){
		return "chr\tBBTools\tRNA\t1\t20\t.\t+\t0\tRNA,startScr:0,model:"+model+extra+'\n';
	}
	private static void summarize(final String input, final String prefix, final String mode){
		NcrnaDdlFinalCallSummary.main(new String[]{"in="+input, "out="+prefix, "mode="+mode, "t=1"});
	}
	private static void reject(final String input, final String prefix, final String mode){
		boolean rejected=false;
		try{summarize(input, prefix, mode);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Contradictory/missing proposal evidence or an output collision must fail loudly: "+input);
	}
	private static void write(final String path, final String text){
		try(Output out=new Output(path)){out.writer.print(text);}
	}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
