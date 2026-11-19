package test.prok;

import prok.Orf;
import prok.ProkObject;
import prok.ScoreTracker;

/** Positive-count score averages and GFF score components reject nonfinite numbers. */
public final class NcrnaFiniteOutputBoundaryTest {
	public static void main(String[] args){
		final ScoreTracker st=new ScoreTracker(ProkObject.CDS);
		check(st.toString().isEmpty() && st.toJson().omapSize()==0 && st.toJson().jmapSize()==0,
			"Empty ScoreTracker output omits undefined averages");
		final byte[] bases="ATGAAATAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final Orf bad=new Orf("bad", 0, 8, 0, 0, bases, true, ProkObject.CDS, true, true);
		bad.startScore=Float.NaN;
		bad.stopScore=2;
		bad.kmerScore=0;
		st.add(bad);
		fails(()->st.toString(), "Nonfinite text score average");
		fails(()->st.toJson(), "Nonfinite JSON score average");
		final Orf orf=new Orf("ctg", 0, 8, 0, 0, bases, true, ProkObject.CDS, true, true);
		orf.orfScore=1;
		orf.startScore=2;
		orf.stopScore=3;
		orf.kmerScore=0;
		check(orf.toGff().contains("startScr:2"), "Finite GFF score components must still serialize");
		orf.stopScore=Float.NaN;
		fails(()->orf.toGff(), "Nonfinite GFF stop score");
		orf.stopScore=3;
		orf.kmerScore=Float.POSITIVE_INFINITY;
		fails(()->orf.toGff(), "Nonfinite GFF inner score");
		System.out.println("PASS NcrnaFiniteOutputBoundaryTest scoretracker_gff_nonfinite_rejection");
	}
	private static void fails(Runnable r, String why){
		try{r.run();}
		catch(IllegalStateException expected){return;}
		throw new AssertionError("Expected "+why);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
