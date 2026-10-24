package test.prok;

import prok.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import aligner.CovarianceModel;
import aligner.CovarianceModelConsensusMap;
import map.LongHashSet;

/** Hand-constructed global maps test exact annotation without scoring/alignment.
 * @author Brian Bushnell, Raiden */
public final class SeedModelPositionTableTest {
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Pass one new fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		final java.lang.reflect.Method factory=test.aligner.CovarianceModelCykTest.class.getDeclaredMethod("pair");factory.setAccessible(true);
		final CovarianceModel cm=(CovarianceModel)factory.invoke(null);
		final byte[] a=new byte[37], g=new byte[17];Arrays.fill(a, (byte)'A');Arrays.fill(g, (byte)'G');
		final byte[] rfA=new byte[37], rfG=new byte[17];Arrays.fill(rfA, (byte)'.');Arrays.fill(rfG, (byte)'.');
		rfA[17]='x';rfA[19]='x';rfG[0]='x';rfG[16]='x';
		final CovarianceModelConsensusMap[] maps={new CovarianceModelConsensusMap(cm,a,a,rfA),new CovarianceModelConsensusMap(cm,g,g,rfG)};
		long gKey=0;for(int i=0; i<17; i++){gKey=(gKey<<2)|2;}
		final LongHashSet seeds=new LongHashSet(8);seeds.add(0);seeds.add(gKey);seeds.add(SeedModelPositionTable.MASK);
		final String[] names={"repeat_A","single_G"}, hashes={"00000000000000000001","00000000000000000002","00000000000000000003","00000000000000000004"};
		final int[] lengths={37,17}, sites={1};
		final SeedModelPositionTable t=SeedModelPositionTable.build(seeds,names,new byte[][]{a,g},maps,sites,hashes);
		check(t.seedCount()==3 && t.get(0).size()==21 && t.get(0).ambiguous(), "All overlapping exact positions must survive, rather than keeping one convenient anchor");
		check(t.get(gKey).size()==1 && !t.get(gKey).ambiguous() && t.get(SeedModelPositionTable.MASK).size()==0,
			"Forward-only annotation retains absent T17 as unmapped, not as reverse-complement support from A17");
		check(t.side(t.get(0).get(1),0)==-1 && t.side(t.get(0).get(2),0)==0 && t.side(t.get(0).get(19),0)==1,
			"The RF boundary maps to query interval18..19; only whole seeds beyond that interval support a side");
		check(t.get(0).get(0).rfStartLow==0 && t.get(0).get(0).rfEndHigh==0,
			"A seed entirely inside a terminal insertion maps to that boundary, never fabricated model residues");
		final String path=dir.resolve("positions.tsv").toString();t.write(path);
		final SeedModelPositionTable loaded=SeedModelPositionTable.load(path,seeds,names,lengths,2,sites,hashes);
		check(loaded.get(0).size()==21 && loaded.get(SeedModelPositionTable.MASK).size()==0 && loaded.side(loaded.get(0).get(19),0)==1,
			"The native loader preserves repeated positions, the unmapped sentinel and projected-site semantics");
		final String[] wrong=hashes.clone();wrong[1]="00000000000000000005";boolean rejected=false;
		try{SeedModelPositionTable.load(path,seeds,names,lengths,2,sites,wrong);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A resource from another consensus generation must fail before it can place seeds");
		final ArrayList<String> lines=new ArrayList<String>(Files.readAllLines(Paths.get(path),StandardCharsets.UTF_8));lines.remove(lines.size()-1);
		final Path missing=dir.resolve("missing.tsv");Files.write(missing,lines,StandardCharsets.UTF_8);rejected=false;
		try{SeedModelPositionTable.load(missing.toString(),seeds,names,lengths,2,sites,hashes);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Unmapped keys are still required resource members, including the last key");
		subset(dir, t, seeds, names, new byte[][]{a,g}, gKey);
		System.out.println("SEED_MODEL_POSITION_TEST_PASS repeat unmapped orientation interval_sides terminal_insertion provenance conservation subset alignments=0 cmScores=0");
	}
	private static void subset(final Path dir, final SeedModelPositionTable source, final LongHashSet seeds,
			final String[] names, final byte[][] bases, final long gKey) throws Exception{
		final String subsetHash="00000000000000000005";
		final SeedModelPositionTable onlyG=source.subset(new int[]{1}, subsetHash);
		check(onlyG.seedCount()==3 && onlyG.modelCount()==1 && onlyG.modelName(0).equals("single_G")
			&& onlyG.get(0).size()==0 && onlyG.get(gKey).size()==1 && onlyG.get(gKey).get(0).model==0,
			"Removing the first model remaps the second to zero while its removed-only seed remains explicitly unmapped");
		final SeedModelPositionTable onlyA=source.subset(new int[]{0}, subsetHash);
		check(onlyA.get(0).size()==21 && onlyA.get(0).ambiguous() && onlyA.get(gKey).size()==0
			&& onlyA.side(onlyA.get(0).get(1),0)==-1 && onlyA.side(onlyA.get(0).get(2),0)==0
			&& onlyA.side(onlyA.get(0).get(19),0)==1 && source.modelCount()==2 && source.get(gKey).get(0).model==1,
			"Subsetting preserves all repeat postings and uncertain site intervals without mutating the source");
		for(int[] bad:new int[][]{new int[0],{1,0},{0,0},{2},{-1}}){
			boolean rejected=false;try{source.subset(bad, subsetHash);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected,"Invalid subset indices must fail before constructing an incompatible model map: "+Arrays.toString(bad));
		}
		boolean rejected=false;try{source.subset(new int[]{0}, "invalid");}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected,"New consensus provenance must satisfy the table's sha80 contract");
		final Path envelope=dir.resolve("source.resource.tsv"), modelFile=dir.resolve("models.fa"), selectedFile=dir.resolve("selected.fa"), seedFile=dir.resolve("seeds.fa");
		SeedModelPositionResource.write(envelope.toString(), source, bases);
		Files.write(modelFile, Arrays.asList(">"+names[0],new String(bases[0],StandardCharsets.US_ASCII),">"+names[1],new String(bases[1],StandardCharsets.US_ASCII)), StandardCharsets.US_ASCII);
		Files.write(selectedFile, Arrays.asList(">single_G",new String(bases[1],StandardCharsets.US_ASCII)), StandardCharsets.US_ASCII);
		Files.write(seedFile, Arrays.asList(">A","AAAAAAAAAAAAAAAAA",">G","GGGGGGGGGGGGGGGGG",">T","TTTTTTTTTTTTTTTTT"), StandardCharsets.US_ASCII);
		final String out=dir.resolve("selected.resource.tsv").toString();
		SeedModelPositionResource.main(new String[]{"positions="+envelope,"models="+modelFile,"seeds="+seedFile,
			"subsetmodels="+selectedFile,"subsetconsensussha80="+subsetHash,"out="+out});
		final SeedModelPositionTable loaded=SeedModelPositionResource.load(out, seeds, new String[]{"single_G"}, new byte[][]{bases[1]});
		check(loaded.modelCount()==1 && loaded.get(gKey).get(0).model==0 && loaded.get(0).size()==0
			&& loaded.get(SeedModelPositionTable.MASK).size()==0,
			"The actual subset CLI must round-trip a remapped model and both preexisting and newly unmapped keys");
		final String text=new String(Files.readAllBytes(Paths.get(out)), StandardCharsets.US_ASCII);
		check(text.contains("#Binding\tconsensusSha80\t"+subsetHash)
			&& text.contains("#Binding\tseedFastaSha80\t00000000000000000001")
			&& text.contains("#Binding\tanchorsSha80\t00000000000000000003")
			&& text.contains("#Binding\tcmSha80\t00000000000000000004"),
			"Only the consensus binding changes; seed, anchor and CM provenance must survive subsetting");
		final Path changed=dir.resolve("changed.fa");Files.write(changed,Arrays.asList(">single_G","CGGGGGGGGGGGGGGGG"),StandardCharsets.US_ASCII);
		rejected=false;try{SeedModelPositionResource.main(new String[]{"positions="+envelope,"models="+modelFile,"seeds="+seedFile,
			"subsetmodels="+changed,"subsetconsensussha80="+subsetHash,"out="+dir.resolve("must_not_exist.tsv")});}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !Files.exists(dir.resolve("must_not_exist.tsv")),
			"A changed consensus is not a subset; reject it before writing inherited positions against different bases");
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
