package jgi;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * TODO
 * @author Brian Bushnell
 * @date Jan 14, 2014
 *
 */
public class CountUniqueness {

	
	public void process(){
		Timer t=new Timer();
		for(String s : in){
			process(s);
		}
		
		t.stop();

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(this.getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void process(Read r1, Read r2){
		if(r1==null || r2==null){return;}
		readsProcessed++;
		basesProcessed+=r1.length();
		readsProcessed++;
		basesProcessed+=r2.length();
		assert(false) : "TODO";
	}
	
	public void process(String fname){
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
			if(verbose){System.err.println("Starting cris");}
			cris.start(); //4567
		}
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					assert(false);
					process(r1, r2);
				}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris);
	
	}

	private static final int MAX=41;
	private static final int MAX2=MAX+1;
	private long[][][] goodMatrix=new long[MAX2][MAX2][MAX2];
	private long[][][] badMatrix=new long[MAX2][MAX2][MAX2];
	
	private PrintStream outstream=System.err;
	private boolean verbose=false;
	private long maxReads=-1;
	private String in[];
	private String out;
	private boolean overwrite=true;
	private boolean append=false;
	private long readsProcessed=0;
	private long basesProcessed=0;
	private boolean errorState=false;
	
	
}
