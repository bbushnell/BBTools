package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ExportNcrnaModelSubset;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec;
import shared.Shared;

/** Literal record oracle: retain C, add A, omit B, and preserve original A,C order.
 * Complete headers, bases and all five positional counts must survive dense and
 * compact encodings. Invalid selectors/pairs must fail before an export exists.
 * @author Raiden
 */
public final class ExportNcrnaModelSubsetTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ExportNcrnaModelSubsetTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown paired-export fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new ExportNcrnaModelSubsetTest(parser.out1).run();
			pp.outstream.println("NCRNA_MODEL_SUBSET_EXPORT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private ExportNcrnaModelSubsetTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		write("original.fa", ">A n=4\nAC\nGT\n>B n=9\nTTTT\n>C n=7\nGGCC\n");
		write("original.hbm", "#HBM\t1\n"+A+B+C);
		write("keep.fa", ">C n=7\nGG\nCC\n");
		write("names.txt", "A\n");
		invoke("original.hbm", "keep.fa", "names.txt", "dense");
		check(read("dense/consensus.fa").equals(">A n=4\nACGT\n>C n=7\nGGCC\n"), "Subset must preserve full original headers/bases and original A,C order, not request order");
		check(read("dense/models.hbm").equals("#HBM\t1\n"+A+C), "Every original header/base/count byte in the selected dense records must survive");
		check(read("dense/summary.tsv").contains("3\t1\t1\t2\tORIGINAL_PAIRED_RECORDS_VERIFIED_NOT_CALLER_VERIFIED"), "One original pair is added to one retained pair without caller verification");
		NucleotideHbmCodec.convert(path("original.hbm"), path("original.hbmc"));
		invoke("original.hbmc", "keep.fa", "names.txt", "compact");
		NucleotideHbmCodec.convert(path("compact/models.hbmc"), path("compact/decoded.hbm"));
		check(read("compact/decoded.hbm").equals("#HBM\t1\n"+A+C), "Compact export must preserve every original decoded record exactly");
		write("unknown.txt", "UNKNOWN\n");
		reject("original.hbm", "keep.fa", "unknown.txt", "unknown", "absent from original");
		write("duplicate.txt", "A\nA\n");
		reject("original.hbm", "keep.fa", "duplicate.txt", "duplicate", "Duplicate addition");
		write("retained.txt", "C\n");
		reject("original.hbm", "keep.fa", "retained.txt", "retained", "already-retained");
		write("bad_keep.fa", ">C n=8\nGGCC\n");
		reject("original.hbm", "bad_keep.fa", "names.txt", "changed", "preserve a unique original");
		write("bad.hbm", ("#HBM\t1\n"+A+B+C).replace(">B n=9", ">WRONG n=9"));
		reject("bad.hbm", "keep.fa", "names.txt", "unselected_bad_pair", "pair differs");
		write("empty.txt", "");
		invoke("original.hbm", "keep.fa", "empty.txt", "empty");
		check(read("empty/models.hbm").equals("#HBM\t1\n"+C), "A member with no additions must still export its original retained pair");
	}

	private void invoke(final String hbm, final String keep, final String names, final String out){
		ExportNcrnaModelSubset.main(new String[]{"in="+path("original.fa"), "hbm="+path(hbm),
			"keep="+path(keep), "names="+path(names), "out="+path(out)});
	}
	private void reject(final String hbm, final String keep, final String names, final String out, final String reason){
		boolean rejected=false;
		try{invoke(hbm, keep, names, out);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected && !new File(path(out)).exists(), "Invalid selection or original pair must fail before export: "+reason);
	}
	private void write(final String name, final String value){try(Output output=new Output(path(name))){output.writer.print(value);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final String root;
	private static final String A=">A n=4\tdepth=4\nACGT\n4\t1\t2\t3\t1\n2\t5\t1\t3\t2\n1\t2\t6\t3\t3\n1\t2\t3\t7\t4\n";
	private static final String B=">B n=9\tdepth=9\nTTTT\n1\t2\t3\t9\t1\n1\t2\t3\t9\t2\n1\t2\t3\t9\t3\n1\t2\t3\t9\t4\n";
	private static final String C=">C n=7\tdepth=7\nGGCC\n1\t2\t7\t3\t4\n2\t3\t8\t4\t5\n3\t9\t4\t5\t6\n4\t10\t5\t6\t7\n";
}
