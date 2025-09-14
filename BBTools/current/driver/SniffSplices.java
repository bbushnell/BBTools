package driver;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.Motif;
import dna.MotifProbsN;
import shared.Tools;

public class SniffSplices {
	
	public static void main(String[] args){
		
//		MotifProbsN mAG=MotifProbsN.makeMotif("AG Exon Starts MP2", 11, 13, 11, 2);
//		MotifProbsN mGT=MotifProbsN.makeMotif("GT Exon Stops MP2", 3, 10, 3, 2);
//
//		MotifProbsN eStarts2=MotifProbsN.makeMotif("Exon Starts MP2", 9, 11, 9, 2);
//		MotifProbsN eStops2=MotifProbsN.makeMotif("Exon Stops MP2", 3, 10, 3, 2);
//
//		MotifProbsN gStarts2=MotifProbsN.makeMotif("Gene Starts MP2", 9, 11, 9, 2);
//		MotifProbsN gStops2=MotifProbsN.makeMotif("Gene Stops MP2", 3, 10, 3, 2);
		

		Motif m=eStops2;
//		Motif m=eStarts2;
//		Motif m=eStarts2_15;
		
		
		ArrayList<String> list=new ArrayList<String>();
		
		boolean rcomp=false;
		if(args.length>0){
			for(String s1 : args){
				String s=s1.toLowerCase();
				if(s.equalsIgnoreCase("rcomp")){rcomp=true;}

				if(s.contains("estart_ac")){m=eStarts2_AC;}
				else if(s.contains("estart_15")){m=eStarts2_15;}
				else if(s.contains("estart")){m=eStarts2;}
				else if(s.contains("estop_gc")){m=eStops2_GC;}
				else if(s.contains("estop")){m=eStops2;}
				else if(s.contains("gstart")){m=gStarts2;}
				else if(s.contains("gstop")){m=gStops2;}
				else{list.add(s.toUpperCase());}
			}
		}
		
		
		System.out.println("Using motif "+m);
		
		int initialLoc=0;
		int increment=1; //1 for plus strand, -1 for minus strand
		
//		String s="NNNNNNNNAGCGGGAATCGGGGGGTCCTTCTGCTCCCCTGAGCGTCCTTCCTGTGTTCCCAGGC"+
//			"ACTATCGCCTACCTGTTTTTCACCAACCGCCACGAGGTGAGGAAGATGACCCTGGACCGAAGCGAATACACCAGCCTCAT"+
//			"CCCAAACTTGAAGAACGTGGTCGCCCTGGACACCGAGGTGGCCAGCAACAGAATATACTGGTCCGACCTGTCCCAAAGGA"+
//			"AGATCTACAGGTGAGCCTTGGAGCCACACCCAGCGCTCAACCCCCGGTGGCGCGGGGGCCCCTCTCACTGACGCTCTCCT"+
//			"TCCCCTGCTCCTCCCCCTCAGCACCCAAATCGACAGAGCCCCCGGCTTCTCCTCCTATGACACCGTCGTCAGCGAGGACC"+
//			"TCCAGGCCCCTGATGGGCTGGCGGTGGACTGGATCCACAGCAACATATACTGGACAGACTCCATCCTGGGCACCGTCTCC"+
//			"GTGGCCGACACCAAGGGCGTGAAGAGAAAGACGCTCTTCAAGGAGAAAGGCTCTAAGCCACGTGCCATCGTGGTGGATCC"+
//			"CGTTCACGGGTGGGTGCTGCTAAAGCCGAGGGCCACGGAAGGAANNNNNNNN";
		
		//		"AAGTACAGGAATTATATGCCCCCAGGTAA * AGTACAGGAATTATATGCCCCCAGGTAAC"
//		String[] array={
//				"GCCTACTTTGTATGATGACCCTGTCCT",
//				"AGCCCTGGCCGCCTACTTTGTATGATGACCCTGTCCTCCCTCACCCA",
//		};
//		String[] array={
//				"TGGCCGCCGCCGACCGTAAGTTTTGCGCGCAAACTCCC",
//				"TGGCCGCCGCCGACCGTTAAGTTTTGCGCGCAAACTCCC",
//		};
//		String[] array={
//				"CAACTGCCAAGGGAAGGGCACGGTTAGCGGCACCCTCATAGGTAAGTGATGGCCCCAGACGCTGGTCTCTCTCCATCTGGACCTGGCCTGGGAGGTGGCTTGG",
//				"CAACTGCCAAGGGAAGGGCACGGTTAGCGGCACCCTCATAGGTGAGTGATGGCCCCAGACGCTGGTCTCTCTCCATCTGGACCTGGCCTGGGAGGTGGCTTGG",
//		};
		
//		String[] array={
//				"GTCTTTCTCATGTGGTCCTTGTGTTCGTCGAGCAGGCCAGCAAGTGTGACAGTCATGGCACCCACCTGGCAGGGG",
//				"GTCTTTCTCATGTGGTCCTTGTGTTCGTTGAGCAGGCCAGCAAGTGTGACAGTCATGGCACCCACCTGGCAGGGG",
//		};
		
//		String[] array={
//				"GCAGGGTCATGGTCACCGACTTCGAGAATGTGCCCGAGGAGGACGGGACCCGCCTCCACAGACAGGTAAGCACAGCCGTCTGATGGGAGGGCTGCCTCTGCCCATATCCCCATCCTGGAG",
//				"GCAGGGTCATGGTCACCGACTTCGAGAATGTGCCCGAGGAGGACGGGACCCGCTTCCACAGACAGGTAAGCACGGCCGTCTGATGGGAGGGCTGCCTCTGCCCATATCCCCATCCTGGAG",
//		};

		
//		String[] array={
//				"RTGTTTTCACTCCAGCCACGGAGCTGGGTCTCTGGTCTCGGGGGCAGCTGTGTGACAGAGCGT" +
//				"GCCTCTCCCTACAGTGCTCTTCGTCTTCCTTTGCCTGGGGGTCTTCCTTCTATGGAAGAACTG",
//				"RTGTTTTCACTCCAGCCACGGAGCTGGGTCTCTGGTCTCGGGGGCAGCTGTGTGACAGAGCGT" +
//				"GCCTCTCCTTACAGTGCTCTTCGTCTTCCTTTGCCTGGGGGTCTTCCTTCTATGGAAGAACTG",
//		};
		
//		String[] array={
////				"CAGCGAAGATGCGAAGGTGATTCCCGGGTGGG",
////				"CAGCGAAGATGCGAAGGTGATTTCCGGGTGGG",
//				"GCGGCCGAAGCGGGCCATGGACGCGCTCAAGT",
//				"GCGGCCGGAGCGGGCCATGGACGCGCTCAAGT",
//		};
		
		
//		String[] array={
//				"AAGTATGTTTTTGCTTTTAGGAGGATTCTCT",
//				"AAGTATGTTTTTGTTTTTAGGAGGATTCTCT",
//		};
		
//		String[] array={
//				"TTAGGTTGCTGGTGTCTGTATAATGTGTGT"+
//				"A"+
//				"TCTTTGTTGCAGGTTTGTTTTTTATTCTGC",
//
//				"TTAGGTTGCTGGTGTCTGTATAATGTGTGT"+
//				"G"+
//				"TCTTTGTTGCAGGTTTGTTTTTTATTCTGC"
//		};
		
//		ATGTATTCTACTTTT[TCTTTT]AAGTATGTTTTTGTTTTTAGGAGGATTCTCTATGG
		
//		String[] array={
//				"CAGGTCCTCGAGATCCTGGGATACAGGAAA",
//				"CAGGTCCTCGAGATCCTGGGATATAGGAAA"
//		};
		
//		String[] array={
//				"TGTTTTTGCTTTTAGGAGGATTCTCTATG",
//				"TGTTTTTGTTTTTAGGAGGATTCTCTATG"
//		};

		
		
		for(String s : list){
			if(rcomp){s=AminoAcid.reverseComplementBases(s);}
			System.out.println("For string "+s+":");

			if(!s.startsWith("N") || !s.endsWith("N")){
				s="NNNN"+s+"NNNN";
			}
			byte[] code=s.getBytes();

			for(int i=0; i<s.length(); i++){

				float strength=m.matchStrength(code, i);
				float norm=m.normalize(strength);
				float percent=-1;
				try {
					percent=m.percentile(norm);
				} catch (Exception e) {
					// TODO Auto-generated catch block
//					e.printStackTrace();
				}

				System.out.print((initialLoc+i*increment)+"\t");

				System.out.print(s.charAt(i)+"  Strength = "+Tools.format("%.4f   ",norm));
				if(percent!=-1){System.out.print(Tools.format("->   %.4f   ",percent));}
				float norm2=norm;
				while(norm2>0.1f){
					norm2-=.1f;
					System.out.print("*");
				}

//				System.out.print("\t"+Tools.format("%.3f   ",m.percentile(norm)));

				System.out.println();
				
			}

		}
		
	}
	

	private static final int N_MOTIF=2;
	
//	private static final MotifProbsN eStarts2=MotifProbsN.makeMotif("Exon Starts MP"+N_MOTIF, 12, 9, 2);
////	private static final MotifProbsN eStops2=MotifProbsN.makeMotif("Exon Stops MP"+N_MOTIF, 3, 11, 3, 2);
//	private static final MotifProbsN eStops2=MotifProbsN.makeMotif("Exon Stops MP"+N_MOTIF, 12, 3, 2);
//
//	private static final MotifProbsN gStarts2=MotifProbsN.makeMotif("Gene Starts MP"+N_MOTIF, 13, 9, 2);
//	private static final MotifProbsN gStops2=MotifProbsN.makeMotif("Gene Stops MP"+N_MOTIF, 11, 3, 2);
//
//	private static final MotifProbsN trStarts2=MotifProbsN.makeMotif("Tr Starts MP"+N_MOTIF, 12, 7, 2);
//	private static final MotifProbsN trStops2=MotifProbsN.makeMotif("Tr Stops MP"+N_MOTIF, 11, 6, 2);

	private static final MotifProbsN eStarts2=MotifProbsN.makeMotif("Exon Starts MP"+N_MOTIF, 13, 9, 2);
	private static final MotifProbsN eStarts2_AC=MotifProbsN.makeMotif("AC Exon Starts MP"+N_MOTIF, 13, 9, 2);
	private static final MotifProbsN eStarts2_15=MotifProbsN.makeMotif("Exon Starts MP"+N_MOTIF, 19, 15, 2);
	private static final MotifProbsN eStops2=MotifProbsN.makeMotif("Exon Stops MP"+N_MOTIF, 13, 4, 2);
	private static final MotifProbsN eStops2_GC=MotifProbsN.makeMotif("GC Exon Stops MP"+N_MOTIF, 13, 4, 2);
	
	private static final MotifProbsN gStarts2=MotifProbsN.makeMotif("Gene Starts MP"+N_MOTIF, 13, 9, 2);
	private static final MotifProbsN gStops2=MotifProbsN.makeMotif("Gene Stops MP"+N_MOTIF, 13, 4, 2);
	
	private static final MotifProbsN trStarts2=MotifProbsN.makeMotif("Tr Starts MP"+N_MOTIF, 13, 7, 2);
	private static final MotifProbsN trStops2=MotifProbsN.makeMotif("Tr Stops MP"+N_MOTIF, 13, 7, 2);
	
	
}
