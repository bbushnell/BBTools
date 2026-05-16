package ddl;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

import java.util.ArrayList;

/**
 * Compares two genomes using any CardinalityTracker type that supports
 * bucketValues(). Reports WKID, ANI, and collision statistics.
 *
 * Usage: SketchCompare genome1.fa genome2.fa [k=31] [buckets=2048] [loglogtype=ddl]
 *
 * Supported types: ddl (DynamicDemiLog, 16-bit), ll6 (LogLog6, 6-bit),
 * and any other type where bucketValues() is non-null.
 *
 * @author Noire
 * @date May 11, 2026
 */
public class SketchCompare {

	public static void main(String[] args){
		if(args.length<2){
			System.err.println("Usage: SketchCompare <genome1.fa> <genome2.fa> [k=31] [buckets=2048] [loglogtype=ddl]");
			System.exit(1);
		}

		String file1=null, file2=null;
		int k=31, buckets=2048;
		String type="ddl";

		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("loglogtype") || a.equals("type")){type=b;}
			else if(file1==null){file1=args[i];}
			else if(file2==null){file2=args[i];}
		}

		Timer t=new Timer();

		CardinalityTracker ctA=CardinalityTracker.makeTracker(type, buckets, k, 12345L, 0f);
		CardinalityTracker ctB=CardinalityTracker.makeTracker(type, buckets, k, 12345L, 0f);

		long basesA=hashFile(file1, ctA, k);
		long basesB=hashFile(file2, ctB, k);

		long cardA=ctA.cardinality();
		long cardB=ctB.cardinality();

		int[] valsA=ctA.bucketValues();
		int[] valsB=ctB.bucketValues();

		if(valsA==null || valsB==null){
			System.err.println("Type '"+type+"' does not support bucketValues(). Cannot compare.");
			System.exit(1);
		}

		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<valsA.length; i++){
			int a=valsA[i], b=valsB[i];
			if(a==0 && b==0){bothEmpty++; continue;}
			if(a<b){lower++;}
			else if(a>b){higher++;}
			else{equal++;}
		}

		int minDiv=Math.min(equal+lower, equal+higher);
		if(minDiv<6){minDiv=6;}
		float wkid=(equal<=0 ? 0 : Math.min(1f, (float)equal/minDiv));
		float ani=(wkid>0 ? (float)Math.exp(Math.log(wkid)/k) : 0);
		float contAB=(equal+higher>0 ? Math.min(1f, (float)equal/(equal+higher)) : 0);
		float contBA=(equal+lower>0 ? Math.min(1f, (float)equal/(equal+lower)) : 0);

		t.stop();

		System.out.println("Type:\t"+type);
		System.out.println("File1:\t"+file1);
		System.out.println("File2:\t"+file2);
		System.out.println("Bases1:\t"+basesA);
		System.out.println("Bases2:\t"+basesB);
		System.out.println("Cardinality1:\t"+cardA);
		System.out.println("Cardinality2:\t"+cardB);
		System.out.println("K:\t"+k);
		System.out.println("Buckets:\t"+buckets);
		System.out.println("Lower:\t"+lower);
		System.out.println("Equal:\t"+equal);
		System.out.println("Higher:\t"+higher);
		System.out.println("BothEmpty:\t"+bothEmpty);
		System.out.println("WKID:\t"+String.format("%.6f", wkid));
		System.out.println("ANI:\t"+String.format("%.6f", ani));
		System.out.println("Containment(1in2):\t"+String.format("%.6f", contAB));
		System.out.println("Containment(2in1):\t"+String.format("%.6f", contBA));
		System.out.println("Time:\t"+t);

		System.err.println(String.format(
			"[%s] ANI: %.2f%%  WKID: %.4f  Equal: %d/%d  BothEmpty: %d",
			type, ani*100, wkid, equal, lower+equal+higher, bothEmpty));
	}

	private static long hashFile(String path, CardinalityTracker ct, int k){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		long bases=0;
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		while(ln!=null && reads!=null && reads.size()>0){
			for(Read r : reads){
				ct.hash(r);
				bases+=r.length();
				if(r.mate!=null){bases+=r.mate.length();}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		return bases;
	}
}
