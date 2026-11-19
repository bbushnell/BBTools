package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionAssemblyCompare;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;

/** Literal whole-assembly equality and conservative nonmatch fixtures.
 * No biological relationship is inferred from sequence inequality.
 * @author Ganyu
 */
public final class NcrnaConstructionAssemblyCompareTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssemblyCompareTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown assembly-comparison fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply out=fresh_directory without sampling, overwrite or extra inputs");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root.getPath());
			pp.outstream.println("NCRNA_CONSTRUCTION_ASSEMBLY_COMPARE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final String a=root+"/original_\u00e9.fa";
		final String b=root+"/renamed_reversed.fa.gz";
		final String c=root+"/multiplicity.fa";
		final String d=root+"/ambiguity.fa";
		write(a, ">a\nAACG\n>duplicate\nAACG\n>iupac\nARYKMBVDHN\n>palindrome\nAT\n");
		write(b, ">different names\nat\n>other\nndhbv\nkmryt\n>twice\ncgtt\n>once\ncgtt\n");
		write(c, ">a\nAACG\n>not_a_duplicate\nAAAA\n>iupac\nARYKMBVDHN\n>palindrome\nAT\n");
		write(d, ">a\nAACG\n>duplicate\nAACG\n>changed_ambiguity\nANYKMBVDHN\n>palindrome\nAT\n");
		write(root+"/circle.fa", ">circular\nACGAT\n");
		write(root+"/rotated.fa", ">rotated\nCGATA\n");
		write(root+"/u.fa", ">unsupported_DNA_U\nAU\n");
		final String same=pair("1", "1", a, b, 20, 4);
		final String nominated=same+pair("1", "2", a, c, 20, 4)+pair("1", "3", a, d, 20, 4)
			+pair("2", "4", root+"/circle.fa", root+"/rotated.fa", 5, 1);
		write(root+"/pairs.tsv", NcrnaConstructionAssemblyCompare.HEADER+'\n'+nominated);
		NcrnaConstructionAssemblyCompare.process(root+"/pairs.tsv", root+"/result", 4, 100);
		checkResult(root+"/result", same);
		NcrnaConstructionAssemblyCompare.main(new String[]{"in="+root+"/pairs.tsv", "out="+root+"/cli",
			"expectedpairs=4", "maxassemblybases=100", "t=1"});
		checkResult(root+"/cli", same);
		write(root+"/empty.tsv", NcrnaConstructionAssemblyCompare.HEADER+'\n');
		NcrnaConstructionAssemblyCompare.process(root+"/empty.tsv", root+"/empty", 0, 100);
		check(rows(root+"/empty.summary.tsv").get(1).equals("0\t0\t0\t0\t0\tPOSITIVE_IDENTITY_ONLY_NOT_FROZEN_SPLIT"),
			"An empty nomination list must not invent comparisons or read assembly payloads");

		fault(root, "duplicate_pair", same+same, 2, 100);
		fault(root, "wrong_count", same, 2, 100);
		fault(root, "wrong_bases", pair("1", "1", a, b, 21, 4), 1, 100);
		fault(root, "wrong_contigs", pair("1", "1", a, b, 20, 5), 1, 100);
		fault(root, "bounded_bases", same, 1, 19);
		fault(root, "changed_path", same+pair("1", "2", c, d, 20, 4), 2, 100);
		fault(root, "changed_dimensions", same+pair("1", "2", a, d, 21, 4), 2, 100);
		fault(root, "wrong_basis", same.replace(NcrnaConstructionAssemblyCompare.BASIS, "GUESS"), 1, 100);
		fault(root, "missing_file", pair("1", "1", a, root+"/absent.fa", 20, 4), 1, 100);
		fault(root, "unsupported_U", pair("1", "1", root+"/u.fa", root+"/u.fa", 2, 1), 1, 100);
		final boolean previous=Read.IUPAC_TO_N;
		try{
			Read.IUPAC_TO_N=true;
			fault(root, "lossy_conversion", same, 1, 100);
		}finally{Read.IUPAC_TO_N=previous;}
		boolean overwrite=false;
		try{NcrnaConstructionAssemblyCompare.process(root+"/pairs.tsv", root+"/result", 4, 100);}
		catch(IllegalArgumentException expected){overwrite=true;}
		check(overwrite, "A previous exact-identity result must never be overwritten");
	}

	private static void checkResult(final String prefix, final String first){
		final ArrayList<String> pairs=rows(prefix+".pairs.tsv.gz");
		check(pairs.size()==5, "Every nominated pair must be retained exactly once");
		check(pairs.get(1).equals(first.trim()+'\t'+NcrnaConstructionAssemblyCompare.MATCH+'\t'
			+NcrnaConstructionAssemblyCompare.CRITERION),
			"Names/order/case/wrapping/compression/orientation must normalize without changing UTF-8 provenance");
		for(int i=2; i<pairs.size(); i++){
			check(pairs.get(i).endsWith('\t'+NcrnaConstructionAssemblyCompare.UNRESOLVED+'\t'
				+NcrnaConstructionAssemblyCompare.CRITERION),
				"Different multiplicity, IUPAC symbols or circular serialization must stay biologically unresolved");
		}
		check(rows(prefix+".summary.tsv").get(1).equals("4\t1\t3\t6\t90\tPOSITIVE_IDENTITY_ONLY_NOT_FROZEN_SPLIT"),
			"Four nominated pairs conserve one positive and three unresolved outcomes; the left assembly is reused");
	}

	private static void fault(final String root, final String name, final String pairs, final int expected, final long limit){
		write(root+'/'+name+".tsv", NcrnaConstructionAssemblyCompare.HEADER+'\n'+pairs);
		boolean rejected=false;
		try{NcrnaConstructionAssemblyCompare.process(root+'/'+name+".tsv", root+'/'+name, expected, limit);}
		catch(IllegalArgumentException invalid){rejected=true;}
		check(rejected, "Malformed assembly comparison must reject: "+name);
	}

	private static String pair(final String eval, final String discovery, final String a, final String b, final long bases, final int contigs){
		return eval+'\t'+discovery+'\t'+a+'\t'+b+'\t'+bases+'\t'+contigs+'\t'+NcrnaConstructionAssemblyCompare.BASIS+'\n';
	}
	private static void write(final String path, final String text){
		try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
