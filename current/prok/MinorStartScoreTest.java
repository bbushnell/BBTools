package prok;

import java.nio.charset.StandardCharsets;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import dna.AminoAcid;
import dna.GeneticCode;

/** Finite score-policy oracles; these do not claim biological accuracy.
 * @author Keqing
 */
public final class MinorStartScoreTest {
	public static void main(String[] args)throws Exception{
		final float configured=GeneCaller.parseMinorStartMultiplier(System.getProperty("bbtools.diagnostic.minorStartMultiplier", "1.0"));
		for(float factor:new float[]{1f, .9f, .8f}){
			for(int codon=0; codon<64; codon++){
				final String bases=AminoAcid.codonToString(codon);
				final boolean minor=bases.equals("CTG") || bases.equals("ATT") || bases.equals("ATC") || bases.equals("ATA");
				final float expected=minor ? 100f*factor : 100f;
				check(GeneCaller.minorStartScore(100f, codon, false, factor)==expected, "Wrong full-score penalty for "+bases+" factor="+factor);
				check(GeneCaller.minorStartScore(100f, codon, true, factor)==100f, "Synthetic5prime starts must not acquire a biological codon penalty");
				check(GeneCaller.minorStartScore(-7f, codon, false, factor)==-7f, "A penalty must not improve a negative score");
				check(Float.floatToRawIntBits(GeneCaller.minorStartScore(-0f, codon, false, factor))==Float.floatToRawIntBits(-0f), "Zero scores retain their exact bits");
			}
		}
		for(String start:new String[]{"CTG", "ATT", "ATC", "ATA", "ATG", "GTG", "TTG"}){
			final byte[] bases=(start+"AAACCCTAA").getBytes(StandardCharsets.US_ASCII);
			final boolean minor=start.equals("CTG") || start.equals("ATT") || start.equals("ATC") || start.equals("ATA");
			for(int strand=0; strand<2; strand++){
				final Orf complete=new Orf("fixture", 0, bases.length-1, strand, 0, bases, true, ProkObject.CDS, false, false);
				final Orf edge=new Orf("fixture", 0, bases.length-1, strand, 0, bases, true, ProkObject.CDS, true, false);
				final Orf stopEdge=new Orf("fixture", 0, bases.length-1, strand, 0, bases, true, ProkObject.CDS, false, true);
				check(complete.startCodon==GeneticCode.codon(start), "Coordinate flipping must retain the biological start codon");
				check(GeneCaller.minorStartScore(100f, complete.startCodon, complete.startTruncated, configured)==(minor ? 100f*configured : 100f), "Both strands use the configured whole-score factor");
				check(GeneCaller.minorStartScore(100f, edge.startCodon, edge.startTruncated, configured)==100f, "Both-strand artificial starts remain exempt");
				check(GeneCaller.minorStartScore(100f, stopEdge.startCodon, stopEdge.startTruncated, configured)==(minor ? 100f*configured : 100f), "An observed start remains eligible when only the stop is truncated");
			}
		}
		check(GeneCaller.minorStartScore(100f, -1, false, .8f)==100f, "Unavailable codons must not be relabeled as minor starts");
		for(String bad:new String[]{"0", "-0.1", "1.1", "NaN", "Infinity", "", "wrong"}){
			boolean rejected=false;
			try{GeneCaller.parseMinorStartMultiplier(bad);}catch(IllegalArgumentException e){rejected=e.getMessage().contains("minorStartMultiplier");}
			check(rejected, "Invalid diagnostic policy requires a named-property diagnostic: "+bad);
		}
		pathScores(configured);
		System.out.println("PASS MinorStartScoreTest factor="+configured+" all64codons positive_negative_zero strand edge_stop_independence invalid_policy actual_DP_overlap_nonoverlap_once_RNA_unchanged");
	}

	/** Exercises both real DP predecessor paths, including cached scores already scaled once. */
	private static void pathScores(float factor)throws Exception{
		final float[] saved={GeneCaller.p0,GeneCaller.p1,GeneCaller.q1,GeneCaller.q2,GeneCaller.q3,GeneCaller.q4,GeneCaller.q5};
		try{
			// Zero only transition bonuses, isolating the candidate contribution to a known500-point predecessor.
			GeneCaller.p0=GeneCaller.p1=GeneCaller.q1=GeneCaller.q2=GeneCaller.q3=GeneCaller.q4=GeneCaller.q5=0f;
			final GeneCaller caller=new GeneCaller(80,100,100,0,0,0,0,0,new GeneModel(false));
			for(int prevStrand=0; prevStrand<2; prevStrand++){
				final Method method=GeneCaller.class.getDeclaredMethod(prevStrand==0 ? "calcPathScorePlus" : "calcPathScoreMinus",
					Orf.class,ArrayList.class,int.class,int.class,int.class);method.setAccessible(true);
				for(int strand=0; strand<2; strand++){
					for(int overlap:new int[]{0,20}){
						for(int type:new int[]{ProkObject.CDS,ProkObject.RNA,ProkObject.r18S}){
							for(boolean edge:new boolean[]{false,true}){
								final Orf prev=orf(20,319,prevStrand,0,ProkObject.CDS,false);
								final int start=overlap==0 ? 400 : 300;
								final Orf next=orf(start,start+299,strand,1,type,edge);
								prev.pathScorePlus=prev.pathScoreMinus=500f;
								next.startScore=next.stopScore=1f;next.kmerScore=300f;
								final float raw=next.calcOrfScore();check(raw>0,"The DP fixture requires an actual positive scoring CDS");
								next.orfScore=type==ProkObject.CDS && !edge ? raw*factor : raw;
								next.pathScorePlus=next.pathScoreMinus=-999999f;
								float expected=overlap==0 || type==ProkObject.RNA ? next.orfScore : next.calcOrfScore(overlap);
								if(overlap>0 && type==ProkObject.CDS && !edge && expected>0){expected*=factor;}
								final ArrayList<Orf> list=new ArrayList<Orf>();list.add(prev);
								check(next.isValidPrev(prev,100),"The fixture must reach an actual valid DP edge");
								method.invoke(caller,next,list,prevStrand,0,0);
								check(next.prev(prevStrand)==prev,"The real DP method must select the supplied predecessor");
								check(Float.floatToIntBits(next.pathScore(prevStrand))==Float.floatToIntBits(500f+expected),
									"DP score must apply the factor once to CDS only: prevStrand="+prevStrand+" strand="+strand+" overlap="+overlap+" type="+type+" edge="+edge);
							}
						}
					}
				}
			}
		}finally{
			GeneCaller.p0=saved[0];GeneCaller.p1=saved[1];GeneCaller.q1=saved[2];GeneCaller.q2=saved[3];
			GeneCaller.q3=saved[4];GeneCaller.q4=saved[5];GeneCaller.q5=saved[6];
		}
	}
	private static Orf orf(int start,int stop,int strand,int frame,int type,boolean edge){
		final byte[] bases=new byte[1000];Arrays.fill(bases,(byte)'A');
		final int a=strand==0 ? start : bases.length-stop-1, z=strand==0 ? stop : bases.length-start-1;
		System.arraycopy("CTG".getBytes(StandardCharsets.US_ASCII),0,bases,a,3);
		return new Orf("fixture",a,z,strand,frame,bases,true,type,edge,false);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
