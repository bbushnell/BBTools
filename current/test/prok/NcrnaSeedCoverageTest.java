package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaSeedCoverage;
import prok.NcrnaClanSeedPrep;
import prok.RrnaResourceIO;

/** Repeated occurrences count once; every raw heldout copy remains in the denominator.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSeedCoverageTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root.resolve("seeds.fa"), ">a\nAAAAAAAAAAAAAAAAA\n>ac\nAAAAAAAAAAAAAAAAC\n>unused\nGGGGGGGGGGGGGGGGG\n");
		write(root.resolve("genes.fa"), ">one\nAAAAAAAAAAAAAAAAAAAA\n>two\nAAAAAAAAAAAAAAAAAC\n>zero\nTTTTTTTTTTTTTTTTTTTT\n");
		write(root.resolve("split.tsv"), "locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group\n"+
			"one\tg1\tg1\tbacteria/A\tone\tHELD_OUT\tfalse\n"+
			"two\tg2\tg2\tbacteria/A\ttwo\tHELD_OUT\tfalse\n"+
			"zero\tg3\tg3\tbacteria/B\tzero\tHELD_OUT\tfalse\n"+
			"training\tg4\tg4\tbacteria/B\ttraining\tTRAIN\tfalse\n");
		run(root, "good");
		run(root, "audit");
		for(String suffix:new String[]{"summary.tsv", "loci.tsv", "histogram.tsv", "missed.ids.txt"}){
			check(java.util.Arrays.equals(Files.readAllBytes(root.resolve("good."+suffix)), Files.readAllBytes(root.resolve("audit."+suffix))),
				"Optional per-key audit must preserve all pre-existing coverage outputs: "+suffix);
		}
		check(!Files.exists(root.resolve("good.keys.tsv")), "Disabled audit must not emit an extra table");
		final ByteFile keyRows=RrnaResourceIO.open(root.resolve("audit.keys.tsv").toString());final LineParser1 keys=new LineParser1('\t');int keyCount=0;
		try{keyRows.nextLine();for(byte[] row; (row=keyRows.nextLine())!=null; ){
			keys.set(row);keyCount++;final int count=keys.parseInt(1);
			check(keys.parseInt(2)==3 && keys.termEquals("HELD_OUT", 3) && keys.termEquals("false", 4), "Per-key support uses the same declared raw member population");
			if(keys.termEquals("AAAAAAAAAAAAAAAAA", 0)){check(count==2, "Repeated A17 occurrences in one member must count once; two different members give2");}
			else if(keys.termEquals("AAAAAAAAAAAAAAAAC", 0)){check(count==1, "A key in only one member must remain visible as singleton support");}
			else{check(keys.termEquals("GGGGGGGGGGGGGGGGG", 0) && count==0, "Unused selected keys must remain explicit zeros");}
		}}finally{check(!keyRows.close(), "Key-support fixture read failed");}
		check(keyCount==3, "Every loaded key must appear exactly once in the support audit");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.summary.tsv").toString());final LineParser1 p=new LineParser1('\t');
		try{
			in.nextLine();p.set(in.nextLine());
			check(p.termEquals("ALL", 0) && p.parseInt(1)==3 && p.parseInt(2)==2 && p.parseInt(3)==1 && p.parseInt(4)==1,
				"Three raw members: two have a seed, only one has two distinct keys; reverse complements are not forward matches");
		}finally{check(!in.close(), "Coverage summary read failed");}
		final ByteFile misses=RrnaResourceIO.open(root.resolve("good.missed.ids.txt").toString());
		try{check(java.util.Arrays.equals(misses.nextLine(), new byte[]{'z','e','r','o'}) && misses.nextLine()==null, "Missed-ID export must retain the exact absent-seed member");}
		finally{check(!misses.close(), "Misses read failed");}
		write(root.resolve("duplicate.fa"), ">one\nAAAAAAAAAAAAAAAAAAAA\n>one\nAAAAAAAAAAAAAAAAAAAA\n>zero\nTTTTTTTTTTTTTTTTTTTT\n");
		boolean failed=false;try{run(root, "duplicate");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Duplicate input IDs must fail even if raw FASTA count still matches");
		write(root.resolve("missing.fa"), ">one\nAAAAAAAAAAAAAAAAAAAA\n>zero\nTTTTTTTTTTTTTTTTTTTT\n");
		failed=false;try{run(root, "missing");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Missing heldout FASTA members must fail rather than reduce the denominator");
		folds(root);
		System.out.println("NCRNA_SEED_COVERAGE_TEST_PASS");
	}
	private static void folds(Path root)throws Exception{
		final String header="locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group\n";
		write(root.resolve("folds.tsv"), header+"one\tg1\tg1\tbacteria/A\tone\tHELD_OUT\tfalse\n"+
			"training\tg2\tg2\tbacteria/B\ttraining\tTRAIN\tfalse\n"+"training_alias\tg3\tg2\tbacteria/B\ttraining\tTRAIN\tfalse\n"+
			"lost\tg4\tg4\tbacteria/B\tlost\tTRAIN\tfalse\n");
		final String labeled=NcrnaClanSeedPrep.SEED_HEADER+'\n'+
			"one\tg1\tGCF_1.1\t1\tbacteria\tA\tc\t+\t1-20\tRF00174\tRF00174\tRF00174\n"+
			"training\tg2\tGCF_2.1\t2\tbacteria\tB\tc\t+\t1-20\tRF01482\tRF01482\tRF01482\n"+
			"training_alias\tg3\tGCF_3.1\t3\tbacteria\tB\tc\t+\t1-20\tRF01482\tRF01482\tRF01482\n"+
			"lost\tg4\tGCF_4.1\t4\tbacteria\tB\tc\t+\t1-20\tRF01689\tB12clade\tRF00174\n";
		write(root.resolve("labels.tsv"), labeled);
		write(root.resolve("train.fa"), ">training\nAAAAAAAAAAAAAAAAAAAA\n>training_alias\nAAAAAAAAAAAAAAAAAAAA\n>lost\nTTTTTTTTTTTTTTTTTTTT\n");
		write(root.resolve("all.fa"), ">one\nAAAAAAAAAAAAAAAAAAAA\n>training\nAAAAAAAAAAAAAAAAAAAA\n>training_alias\nAAAAAAAAAAAAAAAAAAAA\n>lost\nTTTTTTTTTTTTTTTTTTTT\n");
		write(root.resolve("representatives.fa"), ">one\nAAAAAAAAAAAAAAAAAAAA\n>training\nAAAAAAAAAAAAAAAAAAAA\n>lost\nTTTTTTTTTTTTTTTTTTTT\n");
		selected(root, "train", "TRAIN", false, 3, "labels.tsv");selected(root, "all", "ALL", false, 4, "labels.tsv");selected(root, "representatives", "ALL", true, 3, "labels.tsv");
		final ByteFile in=RrnaResourceIO.open(root.resolve("train.families.tsv").toString());final LineParser1 p=new LineParser1('\t');int source=0, biological=0, group=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);check(p.terms()==8, "Each label view needs counts plus coverage, with no denominator collapse");
			if(p.termEquals("SOURCE", 0)){
				source+=p.parseInt(2);if(p.termEquals("RF01482", 1)){check(p.parseInt(2)==2 && p.parseInt(3)==2, "Source view must retain both identical raw alias copies");}
				else{check(p.termEquals("RF01689", 1) && p.parseInt(2)==1 && p.parseInt(3)==0, "Separate cross-labeled miss must remain in its source denominator");}
			}
			else if(p.termEquals("BIOLOGICAL", 0)){biological+=p.parseInt(2);if(p.termEquals("B12clade", 1)){check(p.parseInt(2)==1 && p.parseInt(3)==0, "Cross-labeled member keeps its biological label and seed miss");}}
			else{check(p.termEquals("BUILD_GROUP", 0), "Unknown attribution view");group+=p.parseInt(2);}
		}}finally{check(!in.close(), "Family coverage fixture read failed");}
		check(source==3 && biological==3 && group==3, "Every attribution view must independently conserve raw TRAIN denominator");
		write(root.resolve("badlabels.tsv"), labeled.substring(0, labeled.indexOf("training_alias\t")));
		boolean failed=false;try{selected(root, "train_bad", "TRAIN", false, 3, "badlabels.tsv");}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("train_bad.loci.tsv")), "Missing family attribution must fail before publishing coverage");
	}
	private static void selected(Path root, String out, String fold, boolean representatives, int count, String metadata){
		final String fasta=out.equals("train_bad")?"train":out;
		NcrnaSeedCoverage.main(new String[]{"in="+root.resolve(fasta+".fa"), "seeds="+root.resolve("seeds.fa"), "split="+root.resolve("folds.tsv"),
			"metadata="+root.resolve(metadata), "outprefix="+root.resolve(out), "expected="+count, "fold="+fold, "representatives="+representatives});
	}
	private static void run(final Path root, final String name){
		NcrnaSeedCoverage.main(new String[]{"in="+root.resolve(name.equals("good") || name.equals("audit")?"genes.fa":name+".fa"), "seeds="+root.resolve("seeds.fa"), "split="+root.resolve("split.tsv"),
			"outprefix="+root.resolve(name), "expected=3", "keysupport="+name.equals("audit")});
	}
	private static void write(final Path path, final String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
