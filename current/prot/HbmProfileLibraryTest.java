package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/** Independent composition/order/missing-family and native round-trip checks. */
public final class HbmProfileLibraryTest {
	public static void main(String[] args) throws Exception{
		final Path out=Paths.get(HmmComparisonData.required(HmmComparisonData.options(args),"out"));Files.createDirectory(out);
		final byte[] a=Blosum62.encode("ACDE".getBytes(StandardCharsets.US_ASCII),"a"),b=Blosum62.encode("WYP".getBytes(StandardCharsets.US_ASCII),"b");
		final AAGraph ga=AAGraphScorer.buildModel(a,new byte[][]{a,a}),gb=AAGraphScorer.buildModel(b,new byte[][]{b,b,b});
		final HbmBundleLoader.Loaded la=new HbmBundleLoader.Loaded(new String[]{"a"},new AAGraph[]{ga}),lb=new HbmBundleLoader.Loaded(new String[]{"b"},new AAGraph[]{gb});
		final List<String> roster=Arrays.asList("a","b");
		final HbmBundleLoader.Loaded combined=HbmBundleLoader.Loaded.combine(Arrays.asList(la,lb),roster);
		if(combined.familyCount()!=2 || !combined.repId(1).equals("b")){throw new AssertionError("Combined roster changed");}
		combined.assertStructuralMatch(0,ga);combined.assertStructuralMatch(1,gb);
		reject(()->HbmBundleLoader.Loaded.combine(Collections.singletonList(la),roster),"Missing family accepted");
		reject(()->HbmBundleLoader.Loaded.combine(Arrays.asList(la,la),Arrays.asList("a","a")),"Duplicate family accepted");
		reject(()->HbmBundleLoader.Loaded.combine(Arrays.asList(lb,la),roster),"Wrong order accepted");
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][32];for(byte[] row : provenance){Arrays.fill(row,(byte)1);}
		combined.writeBundle(out.resolve("combined.mqhb"),provenance);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(out.resolve("combined.mqhb"),roster,id->id.equals("a") ? a : b,provenance);
		loaded.assertStructuralEquivalent(combined);
		final double[] bg=combined.positionBackground();
		final HbmPositionModel[] ordinary=combined.positionModels("logodds",0.01,true,-4),explicit=combined.positionModels("logodds",0.01,true,-4,bg);
		for(int i=0; i<2; i++){for(byte[] q : new byte[][]{a,b}){
			final HbmPositionModel.Result x=ordinary[i].align(q,true),y=explicit[i].align(q,true);
			if(x.score!=y.score || x.start!=y.start || x.end!=y.end || !Arrays.equals(x.path,y.path)){throw new AssertionError("Explicit same background changed alignment");}
		}}
		Files.write(out.resolve("ranks.txt"),"1\n0\n".getBytes(StandardCharsets.US_ASCII));
		if(!Arrays.equals(HbmProfileLibrary.readRanks(out.resolve("ranks.txt").toString(),2),new int[]{1,0})){throw new AssertionError("Rank selection order changed");}
		Files.write(out.resolve("duplicate.txt"),"0\n0\n".getBytes(StandardCharsets.US_ASCII));
		reject(()->HbmProfileLibrary.readRanks(out.resolve("duplicate.txt").toString(),2),"Duplicate rank accepted");
		Files.write(out.resolve("outside.txt"),"2\n".getBytes(StandardCharsets.US_ASCII));
		reject(()->HbmProfileLibrary.readRanks(out.resolve("outside.txt").toString(),2),"Out-of-range rank accepted");
		System.err.println("PROFILE_LIBRARY_TEST_PASS composition=PASS missing_duplicate_order=PASS roundtrip=PASS frozen_background=PASS");
	}
	private static void reject(Runnable code,String reason){
		boolean failed=false;try{code.run();}catch(IllegalArgumentException expected){failed=true;}
		if(!failed){throw new AssertionError(reason);}
	}
}
