package pacbio;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Timer;

/**
 * @author Brian Bushnell
 * @date Jul 10, 2012
 *
 */
public class PartitionFastaFile {
	
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		String infile=args[0];
		String outfile=args[1];
		assert(!infile.equalsIgnoreCase(outfile));
		assert(outfile.contains("#"));
		long partition=Integer.parseInt(args[2]);
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
