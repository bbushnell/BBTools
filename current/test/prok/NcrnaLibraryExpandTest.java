package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaLibraryExpand;
import prok.NcrnaResourceReport;
import prok.NucleotideHbmCodec;
import prok.RrnaResourceIO;
import prok.TrnaConsensusBuilder;

/** Verifies additive single-member counts, immutable old models and full remeasurement.
 * @author Brian, Raiden
 */
public final class NcrnaLibraryExpandTest {

	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String a="AAAAAAAAAAAAAAAAAAAA", c="CCCCCCCCCCCCCCCCCCCC", g="GGGGGGGGGGGGGGGGGGGG";
		final String h0="q0 tid_1 member_family=RF01482 build_group=RF01482";
		final String h1="q1 tid_2 member_family=RF01482 build_group=RF01482";
		final String h2="q2 tid_3 member_family=B12clade build_group=RF01482";
		write(root.resolve("input.fa"), ">"+h0+'\n'+a+"\n>"+h1+'\n'+c+"\n>"+h2+'\n'+g+'\n');
		write(root.resolve("ref.fa"), ">cobalamin_RF01482_original n=2\n"+a+'\n');
		write(root.resolve("assignments.tsv"), "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"+
			h0+"\t1\tunknown\told\tcobalamin_RF01482_original n=2\tkept\n"+
			h1+"\t2\tunknown\told\tcobalamin_RF01482_original n=2\tkept\n"+h2+"\t3\tunknown\torphan\t\torphan\n");
		final ArrayList<BaseGraph> old=new ArrayList<BaseGraph>();
		old.add(new BaseGraph("cobalamin_RF01482_original n=2", a.getBytes(StandardCharsets.US_ASCII), null, 0, 0));
		old.get(0).ref[0].add((byte)'A', 3);
		NucleotideHbmCodec.writeModels(old, root.resolve("old.hbmc").toString(), false);
		report(root, "ref.fa", "old.hbmc", "before");
		expand(root, "input.fa", "ref.fa", "old.hbmc", "expanded", 3);
		final BaseGraph[] graphs=TrnaConsensusBuilder.loadModels(root.resolve("expanded.models.hbmc").toString());
		check(graphs.length==3 && graphs[0].name.equals(old.get(0).name), "Expansion must append two models without relabeling the original");
		for(int i=0; i<a.length(); i++){
			check(Arrays.equals(graphs[0].ref[i].acgtCount, old.get(0).ref[i].acgtCount)
				&& graphs[0].del[i].countSum==old.get(0).del[i].countSum, "Old HBM counts must be unchanged at every position");
		}
		for(int model=1; model<=2; model++){
			final byte[] bases=(model==1?c:g).getBytes(StandardCharsets.US_ASCII);
			final BaseGraph singleton=new BaseGraph("independent singleton", bases, null, 0, 0);
			check(Arrays.equals(graphs[model].original, bases), "Singleton FASTA/HBM must retain the exact admitted core");
			for(int i=0; i<bases.length; i++){
				check(Arrays.equals(graphs[model].ref[i].acgtCount, singleton.ref[i].acgtCount)
					&& graphs[model].del[i].countSum==0, "Added counts must match the builder's one-member BaseGraph semantics");
			}
		}
		report(root, "expanded.consensus.fa", "expanded.models.hbmc", "after");
		final ByteFile summary=RrnaResourceIO.open(root.resolve("after.summary.tsv").toString());final LineParser1 p=new LineParser1('\t');
		try{summary.nextLine();p.set(summary.nextLine());check(p.parseInt(0)==3 && p.parseInt(1)==3 && p.parseInt(6)==3 && p.parseInt(13)==9,
			"Expanded library must recover all three members in an independent exhaustive remeasurement");}
		finally{check(!summary.close(), "Summary read failed");}
		final ByteFile ledger=RrnaResourceIO.open(root.resolve("expanded.additions.tsv").toString());
		try{
			ledger.nextLine();p.set(ledger.nextLine());check(p.termEquals("q1", 0) && p.termEquals("RF01482", 3), "Own member addition must retain biological attribution");
			p.set(ledger.nextLine());check(p.termEquals("q2", 0) && p.termEquals("B12clade", 3) && ledger.nextLine()==null,
				"Operational placement must not relabel clan-only members or invent memberships");
		}finally{check(!ledger.close(), "Addition ledger read failed");}
		write(root.resolve("duplicate.fa"), ">"+h0+'\n'+a+"\n>"+h0+'\n'+a+"\n>"+h2+'\n'+g+'\n');
		reject(root, "duplicate.fa", "ref.fa", "old.hbmc", "duplicate", 3);
		write(root.resolve("missing.fa"), ">"+h0+'\n'+a+"\n>"+h1+'\n'+c+'\n');
		reject(root, "missing.fa", "ref.fa", "old.hbmc", "missing", 3);
		write(root.resolve("badref.fa"), ">cobalamin_RF01482_original n=2\n"+c+'\n');
		reject(root, "input.fa", "badref.fa", "old.hbmc", "mismatch", 3);
		reject(root, "input.fa", "ref.fa", "old.hbmc", "denominator", 4);
		System.out.println("NCRNA_LIBRARY_EXPAND_TEST_PASS");
	}
	private static void report(Path root, String ref, String hbm, String out){
		NcrnaResourceReport.main(new String[]{"in="+root.resolve("input.fa"), "ref="+root.resolve(ref), "hbm="+root.resolve(hbm),
			"assignments="+root.resolve("assignments.tsv"), "outprefix="+root.resolve(out), "expected=3", "reportcut=0.6"});
	}
	private static void expand(Path root, String input, String ref, String hbm, String out, int expected){
		NcrnaLibraryExpand.main(new String[]{"in="+root.resolve(input), "ref="+root.resolve(ref), "hbm="+root.resolve(hbm),
			"report="+root.resolve("before.members.tsv.gz"), "outprefix="+root.resolve(out), "family=RF01482", "prefix=cobalamin_RF01482_extra_", "expected="+expected});
	}
	private static void reject(Path root, String input, String ref, String hbm, String out, int expected){
		boolean failed=false;try{expand(root, input, ref, hbm, out, expected);}catch(IllegalArgumentException caught){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".consensus.fa")), "Invalid input must fail before publishing any expanded library: "+out);
	}
	private static void write(Path path, String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
