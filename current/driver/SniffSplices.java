package driver;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.Motif;
import dna.MotifProbsN;
import shared.Tools;

/**
 * Analyzes DNA sequences for splice site motifs using position weight matrices.
 * Provides command-line interface for testing motif strength at each position
 * in input sequences. Supports multiple splice site types including exon starts,
 * exon stops, gene starts, and gene stops with various motif variants.
 *
 * @author Brian Bushnell
 */
public class SniffSplices {
	
	/**
	 * Program entry point for splice site analysis.
	 * Processes command-line arguments to select motif type and input sequences,
	 * then calculates and displays motif strength at each sequence position.
	 * @param args Command-line arguments including motif type and sequences
	 */
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
	

	/** Motif version number used for constructing MotifProbsN instances */
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

	/** Position weight matrix for detecting exon start sites */
	private static final MotifProbsN eStarts2=MotifProbsN.makeMotif("Exon Starts MP"+N_MOTIF, 13, 9, 2);
	/** Position weight matrix for detecting AC-type exon start sites */
	private static final MotifProbsN eStarts2_AC=MotifProbsN.makeMotif("AC Exon Starts MP"+N_MOTIF, 13, 9, 2);
	/**
	 * Extended position weight matrix for detecting exon start sites using 15-bp window
	 */
	private static final MotifProbsN eStarts2_15=MotifProbsN.makeMotif("Exon Starts MP"+N_MOTIF, 19, 15, 2);
	/** Position weight matrix for detecting exon stop sites */
	private static final MotifProbsN eStops2=MotifProbsN.makeMotif("Exon Stops MP"+N_MOTIF, 13, 4, 2);
	/** Position weight matrix for detecting GC-type exon stop sites */
	private static final MotifProbsN eStops2_GC=MotifProbsN.makeMotif("GC Exon Stops MP"+N_MOTIF, 13, 4, 2);
	
	/** Position weight matrix for detecting gene start sites */
	private static final MotifProbsN gStarts2=MotifProbsN.makeMotif("Gene Starts MP"+N_MOTIF, 13, 9, 2);
	/** Position weight matrix for detecting gene stop sites */
	private static final MotifProbsN gStops2=MotifProbsN.makeMotif("Gene Stops MP"+N_MOTIF, 13, 4, 2);
	
	/** Position weight matrix for detecting transcription start sites */
	private static final MotifProbsN trStarts2=MotifProbsN.makeMotif("Tr Starts MP"+N_MOTIF, 13, 7, 2);
	/** Position weight matrix for detecting transcription stop sites */
	private static final MotifProbsN trStops2=MotifProbsN.makeMotif("Tr Stops MP"+N_MOTIF, 13, 7, 2);
	
	
}
