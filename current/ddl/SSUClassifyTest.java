package ddl;

import java.util.ArrayList;
import java.util.Arrays;

import dna.Data;
import fileIO.FileFormat;
import idaligner.QuantumAligner;
import prok.ProkObject;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;

/**
 * Tests tiered SSU/ITS classification scheme:
 *   1) Align to universal 16S and 18S consensus.
 *   2) If ANI>0.64, assign whichever is higher.
 *   3) Align to all 16S and 18S subtypes.
 *   4) If ANI>0.64, assign whichever is higher.
 *   5) If ANI<0.56 to both, assign ITS.
 *   6) Otherwise assign unknown.
 *
 * @author Noire
 * @date May 24, 2026
 */
public class SSUClassifyTest {

	static final float CONFIDENT=0.64f;
	static final float ITS_CEILING=0.56f;

	static final int TYPE_16S=0, TYPE_18S=1, TYPE_ITS=2, TYPE_UNKNOWN=3;
	static final String[] TYPE_NAMES={"16S", "18S", "ITS", "Unknown"};

	public static void main(String[] args){
		int threads=128;
		for(String arg : args){
			String[] split=arg.split("=");
			if(split[0].equals("t") || split[0].equals("threads")){threads=Integer.parseInt(split[1]);}
		}

		// Load consensus variants
		Read[] con16Sall=ProkObject.loadConsensusSequenceType("16S", true, true);
		Read[] con18Sall=ProkObject.loadConsensusSequenceType("18S", true, false);

		// Find universal consensus (first in each array)
		byte[] universal16S=con16Sall[0].bases;
		byte[] universal18S=con18Sall[0].bases;

		byte[][] all16S=toByteArrays(con16Sall);
		byte[][] all18S=toByteArrays(con18Sall);

		System.err.println("16S consensus variants: "+con16Sall.length);
		for(Read r : con16Sall){System.err.println("  "+r.id+": "+r.length()+"bp");}
		System.err.println("18S consensus variants: "+con18Sall.length);
		for(Read r : con18Sall){System.err.println("  "+r.id+": "+r.length()+"bp");}

		// Load reference sequences
		String r16sPath=Data.findPath("?all_prok_16S_best_taxsorted.fa.gz");
		String r18sPath=Data.findPath("?all_euk_18S_best_taxsorted.fa.gz");
		String rITSPath=Data.findPath("?all_ITS_best_taxsorted.fa.gz");

		System.err.println("\nLoading sequences...");
		ArrayList<Read> reads16S=loadAll(r16sPath);
		System.err.println("16S: "+reads16S.size());
		ArrayList<Read> reads18S=loadAll(r18sPath);
		System.err.println("18S: "+reads18S.size());
		ArrayList<Read> readsITS=(rITSPath!=null ? loadAll(rITSPath) : new ArrayList<>());
		System.err.println("ITS: "+readsITS.size());

		System.err.println("\nClassifying with "+threads+" threads...\n");

		int[] class16S=classifyAll(reads16S, universal16S, universal18S, all16S, all18S, threads);
		int[] class18S=classifyAll(reads18S, universal16S, universal18S, all16S, all18S, threads);
		int[] classITS=classifyAll(readsITS, universal16S, universal18S, all16S, all18S, threads);

		// Count classifications
		int[][] counts=new int[3][4]; // [trueType][assignedType]
		for(int c : class16S){counts[0][c]++;}
		for(int c : class18S){counts[1][c]++;}
		for(int c : classITS){counts[2][c]++;}

		System.err.println("=== Classification Results (confident="+CONFIDENT+", its_ceiling="+ITS_CEILING+") ===");
		System.err.println();
		System.err.printf("%-12s %8s %8s %8s %8s %8s%n", "True_Type", "Total", "as_16S", "as_18S", "as_ITS", "Unknown");
		String[] labels={"16S", "18S", "ITS"};
		int[] totals={reads16S.size(), reads18S.size(), readsITS.size()};
		for(int t=0; t<3; t++){
			System.err.printf("%-12s %8d %8d %8d %8d %8d%n",
				TYPE_NAMES[t], totals[t], counts[t][0], counts[t][1], counts[t][2], counts[t][3]);
		}
		System.err.println();
		System.err.printf("%-12s %8s %8s %8s %8s %8s%n", "True_Type", "Total", "%16S", "%18S", "%ITS", "%Unk");
		for(int t=0; t<3; t++){
			System.err.printf("%-12s %8d %7.3f%% %7.3f%% %7.3f%% %7.3f%%%n",
				TYPE_NAMES[t], totals[t],
				counts[t][0]*100.0/totals[t], counts[t][1]*100.0/totals[t],
				counts[t][2]*100.0/totals[t], counts[t][3]*100.0/totals[t]);
		}
	}

	static int[] classifyAll(ArrayList<Read> reads, byte[] uni16S, byte[] uni18S,
			byte[][] all16S, byte[][] all18S, int threads){
		final int[] results=new int[reads.size()];
		final int size=reads.size();
		Thread[] workers=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			workers[t]=new Thread(){
				@Override
				public void run(){
					QuantumAligner aligner=new QuantumAligner();
					for(int i=tid; i<size; i+=threads){
						results[i]=classifyOne(reads.get(i).bases, uni16S, uni18S, all16S, all18S, aligner);
					}
				}
			};
			workers[t].start();
		}
		for(Thread w : workers){try{w.join();}catch(InterruptedException e){}}
		return results;
	}

	static int classifyOne(byte[] bases, byte[] uni16S, byte[] uni18S,
			byte[][] all16S, byte[][] all18S, QuantumAligner aligner){
		// Step 1-2: universal consensus
		float u16=aligner.align(bases, uni16S);
		float u18=aligner.align(bases, uni18S);
		float uBest=Tools.max(u16, u18);
		if(uBest>=CONFIDENT){
			return u16>=u18 ? TYPE_16S : TYPE_18S;
		}

		// Step 3-4: all subtypes
		float best16=u16;
		for(byte[] con : all16S){best16=Tools.max(best16, aligner.align(bases, con));}
		float best18=u18;
		for(byte[] con : all18S){best18=Tools.max(best18, aligner.align(bases, con));}
		float bestAll=Tools.max(best16, best18);
		if(bestAll>=CONFIDENT){
			return best16>=best18 ? TYPE_16S : TYPE_18S;
		}

		// Step 5: below ITS ceiling → ITS
		if(bestAll<ITS_CEILING){
			return TYPE_ITS;
		}

		// Step 6: unknown
		return TYPE_UNKNOWN;
	}

	private static byte[][] toByteArrays(Read[] reads){
		byte[][] result=new byte[reads.length][];
		for(int i=0; i<reads.length; i++){result[i]=reads[i].bases;}
		return result;
	}

	private static ArrayList<Read> loadAll(String path){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		return StreamerFactory.getReads(-1, false, ff, null, null, null);
	}
}
