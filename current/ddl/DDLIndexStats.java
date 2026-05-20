package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import cardinality.DynamicDemiLog;
import shared.Timer;

/**
 * Analyzes DDL index posting list lengths at various minhits thresholds.
 * Loads a DDL sketch file, builds an index, then for each record queries
 * the index (excluding self) and counts how many refs share at least H
 * keys for H=3,5,10,20,40,80.  Outputs a histogram of list lengths.
 *
 * Usage: DDLIndexStats <sketch.tsv> [k=13] [t=32]
 *
 * @author Noire, Brian Bushnell
 * @date May 20, 2026
 */
public class DDLIndexStats {

	public static void main(String[] args){
		String inFile=null;
		int k=13, threads=32;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(inFile==null){inFile=arg;}
		}
		if(inFile==null){
			System.err.println("Usage: DDLIndexStats <sketch.tsv> [k=13] [t=32]");
			System.exit(1);
		}

		Timer t=new Timer();
		System.err.println("Loading "+inFile+"...");
		ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(inFile, k, threads);
		System.err.println("Loaded "+records.size()+" records.");

		final int buckets=records.get(0).ddl.buckets;
		System.err.println("Building index ("+buckets+" buckets)...");
		DDLIndex index=new DDLIndex(buckets);
		index.addAll(records, threads);
		System.err.println("Index built.");

		System.err.println("Occupancy: "+index.populatedCells()+" / "+(buckets*65535L)
			+" = "+String.format("%.4f%%", index.populatedCells()*100.0/(buckets*65535L)));

		final int n=records.size();
		final int[] thresholds={3, 5, 7, 8, 10, 15, 20};
		final int maxBin=5;
		final int nThresh=thresholds.length;

		final AtomicInteger nextQuery=new AtomicInteger(0);
		final int[][] mergedHist=new int[nThresh][maxBin+1];

		System.err.println("Querying "+n+" records with "+threads+" threads...");
		Thread[] workers=new Thread[threads];
		final int[][][] threadHists=new int[threads][nThresh][maxBin+1];

		for(int wi=0; wi<threads; wi++){
			final int wid=wi;
			workers[wi]=new Thread("Stats-"+wi){
				@Override
				public void run(){
					int[][] localHist=threadHists[wid];
					while(true){
						int qi=nextQuery.getAndIncrement();
						if(qi>=n){break;}
						int[] counts=index.query(records.get(qi).ddl);
						int[] numAbove=new int[nThresh];
						for(int ri=0; ri<n; ri++){
							int c=counts[ri];
							if(c<=0 || ri==qi){continue;}
							for(int ti=0; ti<nThresh; ti++){
								if(c>=thresholds[ti]){numAbove[ti]++;}
							}
						}
						for(int ti=0; ti<nThresh; ti++){
							int bin=numAbove[ti]>maxBin ? maxBin : numAbove[ti];
							localHist[ti][bin]++;
						}
					}
				}
			};
			workers[wi].start();
		}
		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){}
		}

		for(int wi=0; wi<threads; wi++){
			for(int ti=0; ti<nThresh; ti++){
				for(int bin=0; bin<=maxBin; bin++){
					mergedHist[ti][bin]+=threadHists[wi][ti][bin];
				}
			}
		}

		t.stop();
		System.err.println("Done in "+t+".");

		StringBuilder sb=new StringBuilder();
		sb.append("NumRefs");
		for(int h : thresholds){sb.append('\t').append("H>=").append(h);}
		sb.append('\n');
		for(int bin=0; bin<=maxBin; bin++){
			sb.append(bin==maxBin ? bin+"+" : ""+bin);
			for(int ti=0; ti<nThresh; ti++){
				sb.append('\t').append(mergedHist[ti][bin]);
			}
			sb.append('\n');
		}
		System.out.print(sb);
	}
}
