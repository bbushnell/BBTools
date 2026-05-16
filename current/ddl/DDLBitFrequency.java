package ddl;

import java.util.ArrayList;
import cardinality.DynamicDemiLog;
import shared.Shared;

/**
 * Analyzes the frequency of 1-bits at each of the 16 bit positions
 * across all non-empty registers in a DDL sketch database.
 *
 * Register format: (absNlz << 10) | invMantissa
 *   Bits 15-10: NLZ exponent (6 bits)
 *   Bits  9-0:  Inverted mantissa (10 bits)
 *
 * Usage: DDLBitFrequency <sketch_file.tsv.gz> [k=31] [t=4]
 *
 * @author Noire
 * @date May 15, 2026
 */
public class DDLBitFrequency {

	public static void main(String[] args){
		String path=null;
		int k=31;
		int threads=Shared.threads();

		for(String arg : args){
			if(arg.startsWith("k=")){k=Integer.parseInt(arg.substring(2));}
			else if(arg.startsWith("t=")){threads=Integer.parseInt(arg.substring(2));}
			else if(arg.startsWith("threads=")){threads=Integer.parseInt(arg.split("=")[1]);}
			else if(path==null){path=arg;}
		}

		if(path==null){
			System.err.println("Usage: DDLBitFrequency <sketch_file.tsv.gz> [k=31] [t=4]");
			System.exit(1);
		}

		System.err.println("Loading sketches from: "+path);
		ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(path, k, threads);
		System.err.println("Loaded "+records.size()+" sketches.");

		final long[] onesCount=new long[16];
		long totalRegisters=0;
		long emptyRegisters=0;

		for(DDLRecord rec : records){
			char[] arr=rec.ddl.toAbsoluteArray();
			for(char val : arr){
				if(val==0){emptyRegisters++; continue;}
				totalRegisters++;
				for(int bit=0; bit<16; bit++){
					if(((val>>bit)&1)==1){
						onesCount[bit]++;
					}
				}
			}
		}

		System.err.println("Total non-empty registers: "+totalRegisters);
		System.err.println("Empty registers: "+emptyRegisters);
		System.err.println("Sketches: "+records.size());
		System.err.println();

		System.out.println("Bit\tSection\tSubBit\tOnes\tTotal\tFrequency");
		for(int bit=15; bit>=0; bit--){
			String section=(bit>=10) ? "NLZ" : "Mantissa";
			int subBit=(bit>=10) ? (bit-10) : bit;
			double freq=(double)onesCount[bit]/totalRegisters;
			System.out.println(bit+"\t"+section+"\t"+subBit+"\t"
					+onesCount[bit]+"\t"+totalRegisters+"\t"
					+String.format("%.6f", freq));
		}
	}

}
