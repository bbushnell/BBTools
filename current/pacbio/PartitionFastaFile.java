package pacbio;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Timer;

/**
 * Partitions FASTA files into multiple smaller files based on sequence length.
 * Splits input FASTA by base count, creating new files when partition size is reached.
 * Ensures complete sequences are not split across partition boundaries.
 *
 * @author Brian Bushnell
 * @date Jul 10, 2012
 */
public class PartitionFastaFile {
	
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		String infile=args[0];
		String outfile=args[1];
		assert(!infile.equalsIgnoreCase(outfile));
		assert(outfile.contains("#"));
		long partition=Integer.parseInt(args[2]);
		//NOTE [pacbio/PartitionFastaFile#001] LOW/dev (no .sh, no callers): args[3] is SKIPPED — maxDataOut reads args[4], so a
		//value passed as args[3] is silently ignored (off-by-one gap in positional args). Also args[0..2] unguarded (AIOOBE<3),
		//parseInt/parseLong NFE. And the maxDataOut cutoff in split() (`dataOut<maxDataOut` checked at loop top) can stop MID-record
		//— header written, sequence partially written — which contradicts the javadoc's "complete sequences are not split"
		//(that guarantee holds only for the per-`partition` boundary, which IS checked at '>' headers). Dead dev tool → LOW.
		if(args.length>4){maxDataOut=Long.parseLong(args[4]);}
		
		if(ReadWrite.ZIPLEVEL<2){ReadWrite.ZIPLEVEL=2;}
		
		TextFile tf=new TextFile(infile, false);
		
		split(tf, outfile, partition);
		t.stop();
		System.out.println("Time:\t"+t);
		
	}
	
	public static void split(TextFile tf, String outfile, long partition) {
		long currentBases=0;
		int pnum=1;

		TextStreamWriter tsw=new TextStreamWriter(outfile.replace("#", ""+pnum), true, false, false);
		tsw.start();
		
		String s;
		for(s=tf.nextLine(); s!=null && dataOut<maxDataOut; s=tf.nextLine()){
			if(s.charAt(0)=='>'){
				if(currentBases>=partition){
					System.out.println("Ended partition "+pnum+" at "+currentBases);
					currentBases=0;
					pnum++;
					tsw.poison();
					tsw=new TextStreamWriter(outfile.replace("#", ""+pnum), true, false, false);
					tsw.start();
				}
			}else{
				int x=s.length();
				currentBases+=x;
				dataOut+=x;
			}
			tsw.println(s);
		}
		System.out.println("Ended partition "+pnum+" at "+currentBases);
		System.out.println("Total: "+dataOut);
		System.out.println("Avg:   "+(dataOut)/pnum);
//		System.out.println("\n"+s+"\n"+dataOut+"\n"+maxDataOut);
		
		//n [PartitionFastaFile] this synchronized(tsw){tsw.wait(100)} is a pointless 100ms sleep — nothing ever notifies tsw, so
		//n it just blocks 100ms before poison(). Harmless (poison drains and the non-daemon writer is joined at JVM exit) but a
		//n code smell / no-op. The per-partition split (currentBases>=partition, checked only at '>' headers) DOES keep whole records.
		try {
			synchronized(tsw){
				tsw.wait(100);
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		tsw.poison();
	}
	
	public static int MIN_CONTIG_TO_ADD=150; //Not currently used
	public static long MAX_OUTPUT_LEN=200000000000L;
	public static long maxDataOut=Long.MAX_VALUE;
	private static long dataOut=0;
	
}
