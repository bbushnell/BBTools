package template;

import java.io.File;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

/**
 * Accepts multiple input files.
 * Reads them each sequentially, and outputs everything to a single output file.
 * Generically, it can be used to concatenate files while recompressing them
 * and avoiding the use of stdio.
 * @author Brian Bushnell
 * @date May 10, 2024
 *
 */
public class A_SampleMultipleInput {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		A_SampleMultipleInput x=new A_SampleMultipleInput(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public A_SampleMultipleInput(String[] args){
		
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

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("in")){
				in.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){in.add(b3);}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
	}
	
	void process(Timer t){

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			for(String s : in) {
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		for(String s : in) {
			processInner(s, ros);
		}
		
		ReadWrite.closeStream(ros);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"
				+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	void processInner(String fname, ConcurrentReadOutputStream ros) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin==null || ffin.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					//  *********  Process reads here  *********
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStream(cris);
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String out1=null;
	
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long readsProcessed=0, basesProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
