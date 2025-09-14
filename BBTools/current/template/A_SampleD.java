package template;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadInputStreamD;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class A_SampleD {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		A_SampleD x=new A_SampleD(args);
		assert(false) : "To support MPI, uncomment this.";
//		MPIWrapper.mpiInit(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
//		MPIWrapper.mpiFinalize();
	}
	
	public A_SampleD(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ConcurrentReadInputStreamD.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("A_SampleD: Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream crisD=ConcurrentReadInputStream.getReadInputStream(
				maxReads, false, ffin1, null, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
		crisD.start();
		final boolean paired=crisD.paired();
		
		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(paired && (in1==null || !in1.contains(".sam"))){
				outstream.println("A_SampleD: Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		{
			
			ListNum<Read> ln=crisD.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			if(verbose){outstream.println("A_SampleD: Initial A_SampleD list: "+(reads==null ? -1 : reads.size()));}
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==crisD.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("A_SampleD: Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//  *********  Process reads here  *********
					
					readsProcessed++;
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				crisD.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("A_SampleD: Returned a list.");}
				ln=crisD.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				crisD.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStreams(crisD, ros);
		if(verbose){outstream.println("A_SampleD: Finished.");}
		
		t.stop();
		outstream.println("A_SampleD: Time:                         \t"+t);
		outstream.println("A_SampleD: Rank "+ Shared.MPI_RANK + ": Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
