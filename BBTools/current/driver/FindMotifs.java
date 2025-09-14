package driver;
import java.util.ArrayList;
import java.util.HashSet;

import dna.ChromosomeArray;
import dna.Data;
import dna.GeneSet;
import dna.Motif;
import dna.MotifProbsN;
import shared.Shared;


public class FindMotifs {
	
	
	public static void main(String[] args){
		
		int chrom=1;
		if(args.length>0){
			chrom=Integer.parseInt(args[0]);
		}
		
		int maxChrom=22;
		
		
		float[][] grid={
				
//				{0.19540f, 0.26751f, 0.34873f, 0.18835f},
//				{0.18987f, 0.33930f, 0.28953f, 0.18131f},
//				{0.19421f, 0.32921f, 0.28259f, 0.19399f},
				
				{0.19519f, 0.23856f, 0.38961f, 0.17664f},
				{0.17382f, 0.33995f, 0.30720f, 0.17903f},
				{0.24452f, 0.38376f, 0.25710f, 0.11462f},
				
				{0.46075f, 0.09954f, 0.38018f, 0.05953f},
				{0.29028f, 0.38874f, 0.19941f, 0.12156f},
				{0.17610f, 0.46129f, 0.28953f, 0.07309f},
				
				{0.99859f, 0.00108f, 0.00011f, 0.00022f},
				{0.00001f, 0.00001f, 0.00076f, 0.99924f},
				{0.00001f, 0.00001f, 0.99924f, 0.00076f},
//
				{0.20993f, 0.14877f, 0.51269f, 0.12861f},
				{0.26903f, 0.39861f, 0.18337f, 0.14899f},
//				{0.14812f, 0.26448f, 0.39286f, 0.19453f},
//				{0.23476f, 0.24073f, 0.35697f, 0.16753f},
//				{0.24886f, 0.32043f, 0.22511f, 0.20560f},
//				{0.16504f, 0.31956f, 0.34288f, 0.17252f},
//				{0.23444f, 0.26838f, 0.32531f, 0.17187f},
		};
		
		
		float[][] gridATG={
				{1, 0, 0, 0},
				{0, 0, 0, 1},
				{0, 0, 1, 0},
		};
		
		
		float[][] grid2={
//		{0.03793f, 0.04397f, 0.06643f, 0.02087f, 0.06272f, 0.11378f, 0.06177f, 0.07713f, 0.06526f, 0.10277f, 0.09535f, 0.04248f, 0.02564f, 0.06590f, 0.05943f, 0.05859f},
//		{0.04534f, 0.04238f, 0.07427f, 0.02956f, 0.06897f, 0.11463f, 0.06579f, 0.07702f, 0.05594f, 0.09821f, 0.08677f, 0.04206f, 0.02363f, 0.05943f, 0.05721f, 0.05880f},
//		{0.04397f, 0.04524f, 0.07766f, 0.02702f, 0.06537f, 0.10372f, 0.07130f, 0.07427f, 0.05234f, 0.09609f, 0.09397f, 0.04164f, 0.02808f, 0.05954f, 0.06399f, 0.05583f},
//
//		{0.04990f, 0.04164f, 0.07342f, 0.02479f, 0.06039f, 0.11294f, 0.06113f, 0.07013f, 0.06929f, 0.10234f, 0.09556f, 0.03973f, 0.02998f, 0.06102f, 0.05329f, 0.05445f},
//		{0.05075f, 0.04725f, 0.08391f, 0.02765f, 0.06590f, 0.11664f, 0.06388f, 0.07151f, 0.06770f, 0.10012f, 0.07840f, 0.03719f, 0.02458f, 0.06113f, 0.05424f, 0.04916f},
//		{0.04831f, 0.04005f, 0.08931f, 0.03125f, 0.06685f, 0.09249f, 0.09546f, 0.07035f, 0.05583f, 0.08359f, 0.10234f, 0.03867f, 0.02278f, 0.05287f, 0.06293f, 0.04693f},
//
//		{0.04587f, 0.05117f, 0.07045f, 0.02627f, 0.05488f, 0.09450f, 0.05541f, 0.06420f, 0.06664f, 0.13328f, 0.10478f, 0.04534f, 0.02087f, 0.06208f, 0.05912f, 0.04513f},
		{0.04598f, 0.04015f, 0.07321f, 0.02892f, 0.06261f, 0.12575f, 0.07331f, 0.07935f, 0.06282f, 0.10637f, 0.08370f, 0.03687f, 0.02246f, 0.05721f, 0.05318f, 0.04810f},
		{0.03952f, 0.03189f, 0.09704f, 0.02543f, 0.07639f, 0.08666f, 0.10266f, 0.06378f, 0.05424f, 0.07850f, 0.11166f, 0.03899f, 0.02426f, 0.04270f, 0.07819f, 0.04810f},
		
		{0.04312f, 0.04291f, 0.08878f, 0.01960f, 0.03666f, 0.09800f, 0.04577f, 0.05933f, 0.07098f, 0.14207f, 0.12395f, 0.05255f, 0.02129f, 0.05795f, 0.04990f, 0.04714f},
		{0.04662f, 0.04672f, 0.06335f, 0.01536f, 0.09069f, 0.14980f, 0.05488f, 0.04556f, 0.07670f, 0.11654f, 0.09026f, 0.02490f, 0.03009f, 0.07162f, 0.04895f, 0.02797f},
		{0.09503f, 0.01695f, 0.11802f, 0.01409f, 0.20945f, 0.03856f, 0.11622f, 0.02045f, 0.11929f, 0.02998f, 0.09440f, 0.01377f, 0.03464f, 0.01504f, 0.05255f, 0.01155f},
		
		{0.16580f, 0.13751f, 0.11325f, 0.04185f, 0.02532f, 0.03687f, 0.01409f, 0.02426f, 0.08645f, 0.19578f, 0.05700f, 0.04195f, 0.01165f, 0.01992f, 0.01451f, 0.01377f},
		{0.07257f, 0.06929f, 0.13042f, 0.01695f, 0.04693f, 0.25532f, 0.05943f, 0.02839f, 0.04079f, 0.06918f, 0.07649f, 0.01240f, 0.01420f, 0.06876f, 0.02373f, 0.01515f},
		{0.17417f, 0.00021f, 0.00011f, 0.00001f, 0.46213f, 0.00021f, 0.00001f, 0.00021f, 0.28944f, 0.00064f, 0.00001f, 0.00001f, 0.07289f, 0.00001f, 0.00001f, 0.00001f},
		
/* */	{0.00001f, 0.00001f, 0.00074f, 0.99788f, 0.00001f, 0.00001f, 0.00001f, 0.00106f, 0.00001f, 0.00001f, 0.00001f, 0.00011f, 0.00001f, 0.00001f, 0.00001f, 0.00021f},
		{0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00074f, 0.00001f, 0.00001f, 0.99926f, 0.00001f},
		{0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.00001f, 0.20934f, 0.14769f, 0.51330f, 0.12893f, 0.00032f, 0.00001f, 0.00042f, 0.00001f},
		
		{0.07416f, 0.04132f, 0.06346f, 0.03072f, 0.03252f, 0.05265f, 0.01843f, 0.04407f, 0.15447f, 0.23382f, 0.07946f, 0.04598f, 0.00795f, 0.07194f, 0.02182f, 0.02723f},
//		{0.04428f, 0.05859f, 0.11855f, 0.04767f, 0.06155f, 0.11622f, 0.12681f, 0.09514f, 0.02829f, 0.05922f, 0.06621f, 0.02945f, 0.01335f, 0.03157f, 0.08168f, 0.02140f},
//		{0.04015f, 0.02691f, 0.06166f, 0.01875f, 0.07850f, 0.06876f, 0.06314f, 0.05520f, 0.08910f, 0.09916f, 0.15319f, 0.05181f, 0.02712f, 0.04534f, 0.08020f, 0.04100f},
//		{0.07522f, 0.05785f, 0.06282f, 0.03899f, 0.05244f, 0.07151f, 0.04428f, 0.07194f, 0.10965f, 0.11601f, 0.08306f, 0.04948f, 0.01165f, 0.07416f, 0.03644f, 0.04450f},
//		{0.05149f, 0.06187f, 0.09959f, 0.03602f, 0.05996f, 0.11749f, 0.07787f, 0.06420f, 0.03094f, 0.08083f, 0.07342f, 0.04142f, 0.02214f, 0.06018f, 0.09185f, 0.03072f},
//		{0.04195f, 0.03380f, 0.06717f, 0.02161f, 0.08656f, 0.09164f, 0.07342f, 0.06876f, 0.07660f, 0.09810f, 0.12003f, 0.04799f, 0.02871f, 0.04460f, 0.06568f, 0.03337f},
//		{0.07819f, 0.05318f, 0.05710f, 0.04534f, 0.05668f, 0.08232f, 0.04471f, 0.08444f, 0.09948f, 0.10520f, 0.06590f, 0.05573f, 0.01685f, 0.06378f, 0.03697f, 0.05414f},
		};
		

		Motif gstartMotif=new MotifProbsN("Gene Starts MP1", grid, 6, 1);
		Motif gstartATG=new MotifProbsN("ATG Gene Starts MP1", gridATG, 0, 1);
		
		Motif gstartMotif2=new MotifProbsN("Gene Starts MP2", grid2, 8, 2);
//
//		Motif estartMotif_ag=new MotifProbs(grid_ag, 9);
//		Motif estartMotif_ac=new MotifProbs(grid_ac, 9);
//		Motif estartMotif_atg=new MotifProbs(grid_atg, 9);
//		Motif estartMotif_nonagac=new MotifProbs(grid_nonagac, 9);
//
//		Motif estartMotif2_ag=new MotifProbsN(grid2_ag, 10);
//		Motif estartMotif2_ac=new MotifProbsN(grid2_ac, 9);
//		Motif estartMotif2_nonagac=new MotifProbsN(grid2_nonagac, 10);

//		Motif estartMotif_multi=new MotifMulti(estartMotif_ag, estartMotif_ac, estartMotif_nonagac);
//		Motif estartMotif2_multi=new MotifMulti(estartMotif2_ag, estartMotif2_ac, estartMotif2_nonagac);

		Motif m=gstartMotif2;
		
		
		ArrayList<Integer> firstBeaten=new ArrayList<Integer>();
		
		long count=0;
		for(chrom=1; chrom<=maxChrom; chrom++){
//			count+=analyzeChromosomeGStarts(chrom, m, locations);
//			count+=analyzeChromosomeGStartsStronger(chrom, m, locations, firstBeaten);
//			count+=analyzeChromosomeGStartsStrongerInFrame(chrom, m, locations, firstBeaten, true, Gene.PLUS);
			count+=analyzeChromosomeGStartsStrongerInFrame(chrom, m, locations, firstBeaten, true, Shared.MINUS);
			Data.unload(chrom, true);
		}
		
		Shared.sort(locations);

		int[] histogram=new int[CLEN+1];
		int[] histogramBeaten=new int[CLEN+1];
		for(Integer i : locations){
			histogram[i]++;
		}
		for(Integer i : firstBeaten){
			histogramBeaten[i]++;
		}
		
		System.out.println(count+" sites analyzed.  ATG occurances:");
		for(int i=0; i<histogram.length; i++){
			System.out.println(i+"\t"+histogram[i]+(firstBeaten.size()==0 ? "" : "\t"+histogramBeaten[i]));
		}
		
	}
	
	
	public static long analyzeChromosomeGStarts(int chrom, Motif m, ArrayList<Integer> list, byte strand){
		GeneSet[] genes=Data.getGeneSets(chrom);
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray ca=Data.getChromosome(chrom);
		
		HashSet<Integer> eset=new HashSet<Integer>();
		for(GeneSet g : genes){
			if(g.strand==strand){
				if(strand==Shared.PLUS){
					eset.add(g.minStart);
				}else{
					eset.add(ca.maxIndex-g.maxEnd);
				}
			}
		}
		
		ArrayList<Integer> list2=new ArrayList<Integer>(eset.size());
		list2.addAll(eset);
		Shared.sort(list2);
		
		for(Integer x : list2){
			
			for(int i=CLEN; i>=0; i--){
				int pos=x-i;
				float f=analyze(pos, m, ca);
				if(f>=THRESH){
					list.add(i);
				}
			}
		}
		return list2.size();
	}
	
	
	public static long analyzeChromosomeGStartsStronger(int chrom, Motif m, ArrayList<Integer> list, ArrayList<Integer> listBeat, byte strand){
		GeneSet[] genes=Data.getGeneSets(chrom);
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray ca=Data.getChromosome(chrom);
		
		HashSet<Integer> eset=new HashSet<Integer>();
		for(GeneSet g : genes){
			if(g.strand==strand){
				if(strand==Shared.PLUS){
					eset.add(g.minStart);
				}else{
					eset.add(ca.maxIndex-g.maxEnd);
				}
			}
		}
		
		ArrayList<Integer> list2=new ArrayList<Integer>(eset.size());
		list2.addAll(eset);
		Shared.sort(list2);
		
		for(Integer x : list2){
			
//			for(int i=CLEN; i>=0; i--){
//				int pos=x-i;
//				float f=analyze(pos, list, m, ca);
//				if(f>=THRESH){
//					list.add(i);
//				}
//			}
			
			int firstBeaten=CLEN+1;
			float basis=analyze(x, m, ca);
			for(int i=0; i<=CLEN; i++){
				int pos=x-i;
				float f=analyze(pos, m, ca);
				if(f>=basis){
					if(i>0 && i<firstBeaten){
						firstBeaten=i;
						listBeat.add(firstBeaten);
					}
					list.add(i);
				}
			}
		}
		return list2.size();
	}
	
	
	public static long analyzeChromosomeGStartsStrongerInFrame(int chrom, Motif m, ArrayList<Integer> list, ArrayList<Integer> listBeat, boolean in, byte strand){
		GeneSet[] genes=Data.getGeneSets(chrom);
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray ca=Data.getChromosome(chrom);
		
		HashSet<Integer> eset=new HashSet<Integer>();
		for(GeneSet g : genes){
			if(g.strand==strand){
				if(strand==Shared.PLUS){
					eset.add(g.minStart);
				}else{
					eset.add(ca.maxIndex-g.maxEnd);
				}
			}
		}
		
		ArrayList<Integer> list2=new ArrayList<Integer>(eset.size());
		list2.addAll(eset);
		Shared.sort(list2);
		
		for(Integer x : list2){
			
//			for(int i=CLEN; i>=0; i--){
//				int pos=x-i;
//				float f=analyze(pos, list, m, ca);
//				if(f>=THRESH){
//					list.add(i);
//				}
//			}
			
			int firstBeaten=CLEN+1;
			float basis=analyze(x, m, ca);
			for(int i=0; i<=CLEN; i++){
				int pos=x-i;
				

				if((in && i%3==0) || (!in && i%3==1)){
					float f=analyze(pos, m, ca);
					if(f>=basis){
						if(i>0 && i<firstBeaten){
							firstBeaten=i;
							listBeat.add(firstBeaten);
						}
						list.add(i);
					}
				}
				
				
			}
		}
		return list2.size();
	}
	
	
	public static float analyze(int point, Motif m, ChromosomeArray ca){
		
		float f=m.matchStrength(ca.array, point);
		f=m.normalize(f);
		return f;
	}
	
	
	private static String padFront(String s, int len){
		int spaces=len-s.length();
		for(int i=0; i<spaces; i++){
			s=" "+s;
		}
		return s;
	}
	
	public static void swap(long[] a, int x, int y){
		long temp=a[x];
		a[x]=a[y];
		a[y]=temp;
	}
	
	public static void swap(char[] a, int x, int y){
		char temp=a[x];
		a[x]=a[y];
		a[y]=temp;
	}
	
	
	public static float THRESH=.2f;
	
	
	public static long analyses=0;
	
	
	/** Exon start */
	public static final int ESTART=0;
	
	/** Exon stop */
	public static final int ESTOP=1;
	
	/** Gene (and exon) start */
	public static final int GSTART=2;
	
	/** Gene (and exon) stop */
	public static final int GSTOP=3;
	
	
	/** Exon start using AG */
	public static final int ESTARTAG=4;
	
	/** Exon stop using GT */
	public static final int ESTOPGT=5;
	
	
	/** Exon start without AG */
	public static final int ESTARTATG=6;
	
	/** Exon stop without GT */
	public static final int ESTOPNONGT=7;

	
	/** Exon start without AG, AC, or ATG */
	public static final int ESTARTNON=8;
	
	/** Exon start using AC */
	public static final int ESTARTAC=9;
	
	private static final int CLEN=200;

	public static ArrayList<Integer> locations=new ArrayList<Integer>();
	
	

	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
}
