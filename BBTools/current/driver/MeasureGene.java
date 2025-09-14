package driver;

import java.util.ArrayList;
import java.util.HashSet;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Exon;
import dna.Gene;
import dna.MotifMulti;
import dna.MotifProbsN;
import shared.Shared;
import shared.Tools;

public class MeasureGene {
	
	
	public static void main(String[] args){
		
		byte minChrom=19;
		byte maxChrom=22;
		
		double sum=0;
		long count=0;
		
		
		for(byte chrom=minChrom; chrom<=maxChrom; chrom++){
			Data.getChromosome(chrom);
			Gene[] genes=Data.getGenes(chrom, Shared.PLUS);
			genes=toNormalGenes(genes);
			
			for(Gene g : genes){
//				ArrayList<Exon> exons=getExons(g);
				
				analyzeGene(g);
//				System.out.println("\nchr"+g.chromosome+"\t"+g.name+"\t"+g.nameTranscript);
//
//				for(int i=0; i<g.exons.length; i++){
//					Exon e=g.exons[i];
//
//					float f2=measureExonFrequency(e.a, e.b, e.chromosome, e.strand);
//					float f1, f3;
//
//					if(i==0 && g.exons.length==1){
//						f1=mGStart.matchStrength(ca.array, e.a);
//						f3=mGStop.matchStrength(ca.array, e.b);
//					}else if(i==0){
//						f1=mGStart.matchStrength(ca.array, e.a);
//						f3=mEStop.matchStrength(ca.array, e.b);
//					}else if(i==g.exons.length-1){
//						f1=mEStart.matchStrength(ca.array, e.a);
//						f3=mGStop.matchStrength(ca.array, e.b);
//					}else{
//						f1=mEStart.matchStrength(ca.array, e.a);
//						f3=mEStop.matchStrength(ca.array, e.b);
//					}
//
//					if(f2!=0){
//						sum+=f2;
//						count++;
//						System.out.println(Tools.format("%.3f, %.3f, %.5f", f1, f3, f2));
//					}
//				}
			}
			
			
		}

		System.out.println("Sum: "+sum);
		System.out.println("Count: "+count);
		System.out.println("Average: "+sum/count);
		
	}
	
	
	public static float analyzeGene(Gene g){
		assert(g.strand==Shared.PLUS) : "TODO";
		ChromosomeArray ca=Data.getChromosome(g.chromosome);

		System.out.println("\nchr"+g.chromosome+"\t"+g.symbol+"\t"+g.mrnaAcc);
		
		double sum=0;
		
		for(int i=0; i<g.exons.length; i++){
			Exon e=g.exons[i];

			float f2=measureExonFrequency(e.a, e.b, e.chromosome, e.strand);
			float f1, f3;

			if(i==0 && g.exons.length==1){
				f1=mGStart.matchStrength(ca.array, e.a);
				f3=mGStop.matchStrength(ca.array, e.b);
			}else if(i==0){
				f1=mGStart.matchStrength(ca.array, e.a);
				f3=mEStop.matchStrength(ca.array, e.b);
			}else if(i==g.exons.length-1){
				f1=mEStart.matchStrength(ca.array, e.a);
				f3=mGStop.matchStrength(ca.array, e.b);
			}else{
				f1=mEStart.matchStrength(ca.array, e.a);
				f3=mEStop.matchStrength(ca.array, e.b);
			}
			
			sum=sum+f1+f3;

//			if(f2!=0){System.out.println(Tools.format("%.3f, %.3f, %.5f", f1, f3, f2));}
		}
		
		float avg=(float)(sum/(2*g.exons.length));
		
		System.out.println(Tools.format("Average: %.3f", avg));
		return avg;
	}
	
	
	public static Gene[] toNormalGenes(Gene[] genes){
		ArrayList<Gene> normal=new ArrayList<Gene>(genes.length);
		for(Gene g : genes){
			if(g.isNormalGene()){normal.add(g);}
		}
		return normal.toArray(new Gene[normal.size()]);
	}
	
	
	public static ArrayList<Exon> getExons(Gene...genes){
		HashSet<Exon> exonTable=new HashSet<Exon>();
		for(Gene g : genes){
			for(int i=0; i<g.exons.length; i++){
				Exon e=g.exons[i];
				exonTable.add(e);
			}
		}
		ArrayList<Exon> exons=new ArrayList<Exon>(exonTable.size());
		exons.addAll(exonTable);
		exonTable=null;
		Shared.sort(exons);
		return exons;
	}
	
	
	public static float measureExonFrequency(int a, int b, byte chrom, byte strand){
//		assert e.strand==Gene.PLUS;
		
		int start=a;
		int stop=b-1;
		
		double sum=0;
		int count=0;
		
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray ca=Data.getChromosome(chrom);
		
		for(int i=start; i<stop; i++){
			int number=0;
			boolean invalid=false;
			for(int j=0; j<length; j++){
				int code=AminoAcid.baseToNumberACGTN[ca.get(i+j)];
				invalid=invalid || (code<0 || code>3);
				number=((number<<2)|code);
			}
			if(!invalid){
				count++;
				sum+=freqDif[number];
			}else{
				return 0;
			}
		}
		
		return count>0 ? (float)(sum/count) : 0;
	}
	
	

	

	private static final MotifProbsN mAG=MotifProbsN.makeMotif("AG Exon Starts MP2", 13, 11, 2);
	private static final MotifProbsN mAC=MotifProbsN.makeMotif("AC Exon Starts MP2", 13, 11, 2);
	private static final MotifProbsN mATG=MotifProbsN.makeMotif("ATG Exon Starts MP2", 13, 11, 2);

	private static final MotifProbsN mGT=MotifProbsN.makeMotif("GT Exon Stops MP2", 10, 3, 2);
	private static final MotifProbsN mGC=MotifProbsN.makeMotif("GC Exon Stops MP2", 10, 3, 2);

	private static final MotifProbsN mGStartATG=MotifProbsN.makeMotif("Gene Starts MP2", 13, 11, 2);
	
	private static final MotifProbsN mGStopTAA=MotifProbsN.makeMotif("TAA Gene Stops MP2", 13, 11, 2);
	private static final MotifProbsN mGStopTAG=MotifProbsN.makeMotif("TAG Gene Stops MP2", 13, 11, 2);
	private static final MotifProbsN mGStopTGA=MotifProbsN.makeMotif("TGA Gene Stops MP2", 13, 11, 2);

	private static final MotifMulti mGStart=new MotifMulti("Gene Starts MP2", mGStartATG);
	private static final MotifMulti mEStart=new MotifMulti("Exon Starts MP2", mAG, mAC);
	private static final MotifMulti mEStop=new MotifMulti("Exon Stops MP2", mGT, mGC);
	private static final MotifMulti mGStop=new MotifMulti("Gene Stops MP2", mGStopTAA, mGStopTAG, mGStopTGA);
	
	
	private static final int length=2;
	
	//Overall Frequency Exonic

	public static final float[] exonicFreq1={0.259195f, 0.260530f, 0.260441f, 0.219835f};

	//Overall Frequency Non-Exonic

	public static final float[] nonExonicFreq1={0.277111f, 0.204189f, 0.213443f, 0.305257f};
	
	//Overall Frequency Exonic

	public static final float[] exonicFreq2={0.071395f, 0.055355f, 0.077256f, 0.052618f,
		0.079593f, 0.077505f, 0.032685f, 0.071248f, 0.075189f, 0.070017f, 0.070666f,
		0.045554f, 0.032210f, 0.057977f, 0.079080f, 0.051651f};

	//Overall Frequency Non-Exonic

	public static final float[] nonExonicFreq2={0.086472f, 0.047310f, 0.070451f, 0.072291f,
		0.069003f, 0.055260f, 0.011722f, 0.071913f, 0.058469f, 0.045772f, 0.056984f,
		0.054175f, 0.062555f, 0.059560f, 0.076273f, 0.101790f};
	
	public static final float[] freqDif=(
			length==2 ? makeDif(exonicFreq2, nonExonicFreq2) :
				length==1 ? makeDif(exonicFreq1, nonExonicFreq1) :
					null);
	
	public static final float[] makeDif(float[] a, float[] b){
		float[] dif=new float[a.length];
		for(int i=0; i<a.length; i++){
			dif[i]=a[i]-b[i];
		}
		return dif;
	}
	
}
