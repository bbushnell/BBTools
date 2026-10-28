package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaClanSeedPrep;
import prok.NcrnaSeedSplit;

/** Tests source-block separation, alias conservation and rare-phylum constraints.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSeedSplitTest {

	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final ArrayList<String> groups=new ArrayList<String>(), loci=new ArrayList<String>(), ids=new ArrayList<String>();
		for(int i=0; i<9; i++){
			final String genome="g"+i, group=i==1?"g0":i==5?"g4":genome, phylum=i<5?"A":i<8?"B":"SINGLETON";
			final String hash=String.format(Locale.ROOT, "%020x", i==1?1:i+1), id="RF00174_"+i;
			groups.add(genome+"\tGCF_"+i+".1\t"+(i+1)+'\t'+hash+'\t'+group+"\tfalse\n");
			loci.add(id+'\t'+genome+"\tGCF_"+i+".1\t"+(i+1)+"\tbacteria\t"+phylum+
				"\tcontig\t+\t101-290\t190\t1\tf\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t100\t200\t1\n");
			ids.add(id+'\n');
		}
		writeInputs(root, groups, loci, ids);run(root, "a", 9);
		final ObjectMap<String,String> folds=new ObjectMap<String,String>(String.class, String.class);
		final ByteFile in=ByteFile.makeByteFile(root.resolve("a.loci.tsv").toString(), false);
		final LineParser1 p=new LineParser1('\t');int count=0, aliases=0, trainRaw=0, heldRaw=0;
		try{
			check(in.nextLine()!=null, "Split locus header required");
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);check(p.terms()==7, "Split schema changed");count++;
				final String id=p.parseString(0), group=p.parseString(2), fold=p.parseString(5), prior=folds.get(group);
				check(prior==null || prior.equals(fold), "One source group cannot leak across folds");folds.put(group, fold);
				if(fold.equals("TRAIN")){trainRaw++;}else{check(fold.equals("HELD_OUT"), "Unknown fold");heldRaw++;}
				if(id.equals("RF00174_1")){check(p.termEquals("RF00174_0", 4), "Identical source bytes/coordinates need the lexical representative");aliases++;}
				else{check(p.termEquals(id, 4), "Distinct source loci must not be deduplicated by group alone");}
			}
		}finally{check(!in.close(), "Fixture ledger read failed");}
		check(count==9 && aliases==1 && trainRaw+heldRaw==9, "All raw loci and exactly one alias must survive");
		check(folds.get("g4").equals("TRAIN") && folds.get("g8").equals("TRAIN"), "Mixed-phyla and singleton groups must stay in TRAIN");
		check(heldRaw>0 && trainRaw>0, "Both folds must be populated");
		final int trainIds=Files.readAllLines(root.resolve("a.train.ids.txt"), StandardCharsets.US_ASCII).size();
		final int heldIds=Files.readAllLines(root.resolve("a.holdout.ids.txt"), StandardCharsets.US_ASCII).size();
		check(heldIds==heldRaw && trainIds==trainRaw-(folds.get("g0").equals("TRAIN")?1:0),
			"Evaluation keeps raw copies while construction collapses only proven aliases");
		final ByteFile phyla=ByteFile.makeByteFile(root.resolve("a.phyla.tsv").toString(), false);int strata=0;
		try{
			phyla.nextLine();for(byte[] row; (row=phyla.nextLine())!=null; ){
				p.set(row);strata++;
				if(p.termEquals("bacteria/SINGLETON", 0)){check(p.parseInt(3)==0 && p.termEquals("SINGLETON_TRAIN_ONLY", 8), "Singleton exception must be visible");}
				else{check(p.parseInt(2)>0 && p.parseInt(3)==1, "Separable phyla must appear in both folds");}
			}
		}finally{check(!phyla.close(), "Fixture phylum read failed");}
		check(strata==3, "Phylum labels must be preserved");
		Collections.reverse(groups);Collections.reverse(loci);Collections.reverse(ids);writeInputs(root, groups, loci, ids);run(root, "b", 9);
		for(String suffix:new String[]{".train.ids.txt", ".train.raw.ids.txt", ".holdout.ids.txt", ".loci.tsv", ".groups.tsv", ".phyla.tsv", ".summary.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("a"+suffix)), Files.readAllBytes(root.resolve("b"+suffix))), "Row order changed split: "+suffix);
		}
		check(Files.readAllLines(root.resolve("a.train.raw.ids.txt")).size()==trainRaw, "Raw TRAIN export must retain aliases for per-family coverage denominators");
		clan(root, loci, ids);
		boolean failed=false;try{run(root, "wrong_count", 8);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Changed full-population denominator must fail");
		loci.remove(0);writeInputs(root, groups, loci, ids);failed=false;
		try{run(root, "missing_locus", 9);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "An eligible ID without a source row must fail");
		System.out.println("NCRNA_SEED_SPLIT_TEST_PASS");
	}
	private static void clan(Path root, ArrayList<String> original, ArrayList<String> ids)throws Exception{
		final ArrayList<String> compact=new ArrayList<String>(), allIds=new ArrayList<String>(ids);
		for(String row:original){
			final String[] fields=row.trim().split("\t");final StringBuilder b=new StringBuilder();
			for(int i=0; i<9; i++){if(i>0){b.append('\t');}b.append(fields[i]);}compact.add(b+"\tRF00174\tRF00174\tRF00174\n");
		}
		compact.add("RF01482_0\tg0\tGCF_0.1\t1\tbacteria\tA\tcontig\t+\t101-290\tRF01482\tB12clade\tRF00174\n");
		compact.add("RF01689_0\tg0\tGCF_0.1\t1\tbacteria\tA\tcontig\t+\t101-290\tRF01689\tRF01689\tRF01689\n");
		allIds.add("RF01482_0\n");allIds.add("RF01689_0\n");
		write(root.resolve("clan.tsv"), NcrnaClanSeedPrep.SEED_HEADER+'\n', compact);write(root.resolve("clan.ids"), "", allIds);
		NcrnaSeedSplit.main(new String[]{"eligible="+root.resolve("clan.ids"), "loci="+root.resolve("clan.tsv"), "groups="+root.resolve("groups.tsv"),
			"outprefix="+root.resolve("clan"), "families=RF00174,RF01482,RF01689", "expected=11", "seed=20261025", "holdout=0.2"});
		final ByteFile in=ByteFile.makeByteFile(root.resolve("clan.loci.tsv").toString(), false);final LineParser1 p=new LineParser1('\t');
		String sharedFold=null;int count=0, representatives=0, shared=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);count++;final String id=p.parseString(0);if(p.termEquals(id, 4)){representatives++;}
			if(p.termEquals("g0", 2)){
				shared++;final String fold=p.parseString(5);check(sharedFold==null || sharedFold.equals(fold), "Coincident cross-family loci must still share a source-group fold");sharedFold=fold;
				check(p.termEquals(id.equals("RF00174_1")?"RF00174_0":id, 4), "Alias collapse may merge same-family copies but must preserve distinct biological/source labels");
			}
		}}finally{check(!in.close(), "Clan split fixture read failed");}
		check(count==11 && representatives==10 && shared==4, "Multifamily split must conserve every raw member and both coincident distinct-family contributors");
	}

	private static void run(final Path root, final String name, final int count){
		NcrnaSeedSplit.main(new String[]{"eligible="+root.resolve("eligible.txt"), "loci="+root.resolve("loci.tsv"), "groups="+root.resolve("groups.tsv"),
			"outprefix="+root.resolve(name), "family=RF00174", "expected="+count, "seed=20261025", "holdout=0.2"});
	}
	private static void writeInputs(final Path root, final ArrayList<String> groups, final ArrayList<String> loci, final ArrayList<String> ids)throws Exception{
		write(root.resolve("groups.tsv"), "genome\tassembly\ttaxid\tfasta_sha80\tsource_group\tpilot_exposed_group\n", groups);
		write(root.resolve("loci.tsv"), NcrnaCmLocusAudit.HEADER+'\n', loci);write(root.resolve("eligible.txt"), "", ids);
	}
	private static void write(final Path path, final String header, final ArrayList<String> rows)throws Exception{
		final StringBuilder b=new StringBuilder(header);for(String row:rows){b.append(row);}Files.write(path, b.toString().getBytes(StandardCharsets.US_ASCII));
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
