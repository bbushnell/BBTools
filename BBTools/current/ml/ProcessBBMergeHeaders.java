package ml;

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
import stream.FASTQ;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class ProcessBBMergeHeaders {

	public static void main(String[] args){
		Timer t=new Timer();
		ProcessBBMergeHeaders x=new ProcessBBMergeHeaders(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ProcessBBMergeHeaders(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
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
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin1, null);
			cris.start();
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			ArrayList<Read> keep=new ArrayList<Read>(reads.size());
			keep.add(new Read(null, null, headerString(), 0));

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					Header h=makeHeader(r1.id);
					if(h!=null){
						r1.id=h.toString();
						keep.add(r1);
					}
					
					readsProcessed++;
				}
				
				if(ros!=null){ros.add(keep, ln.id);}
				keep=new ArrayList<Read>();

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private Header makeHeader(String line){
		if(!line.startsWith("insert=")){return null;}
		if(!line.contains(" mo=")){return null;}
		Header h=new Header(line);
		return h.valid ? h : null;
	}
	
	public String headerString(){
		return "#Correct\tminOverlap\tbestOverlap\tbestBadInt\tsecondBestOverlap\tsecondBestBadInt\t"
				+ "expectedErrors\tbestExpectedErrors\tbestRatio\tbestBad\tsecondBestRatio\tsecondBestBad\tprobability";
	}
	
	private class Header {
	
		//mo=14_r1ee=5.2728_r2ee=3.4856_bi=202_bo=98_bb=5.3063_br=0.0598_bbi=6_sbi=270_sbo=30_sbb=12.4775_sbr=0.4343_sbbi=14_be=6.5990_pr=0.0007
		
		Header(String line_){
			line=line_;
			String[] split=line.split(" ");
			trueInsert=Integer.parseInt(split[0].split("=")[1]);
			
			split=split[2].split("_");
			for(String s : split){
				String[] split2=s.split("=");
				String a=split2[0], b=split2[1];
				if(a.equals("mo")){
					minOverlap=Integer.parseInt(b);
				}else if(a.equals("bi")){
					bestInsert=Integer.parseInt(b);
				}else if(a.equals("bo")){
					bestOverlap=Integer.parseInt(b);
				}else if(a.equals("bbi")){
					bestBadInt=Integer.parseInt(b);
				}else if(a.equals("sbi")){
					secondBestInsert=Integer.parseInt(b);
				}else if(a.equals("sbo")){
					secondBestOverlap=Integer.parseInt(b);
				}else if(a.equals("sbbi")){
					secondBestBadInt=Integer.parseInt(b);
				}else if(a.equals("r1ee")){
					expectedErrors1=Float.parseFloat(b);
				}else if(a.equals("r2ee")){
					expectedErrors2=Float.parseFloat(b);
				}else if(a.equals("be")){
					bestExpected=Float.parseFloat(b);
				}else if(a.equals("pr")){
					probability=Float.parseFloat(b);
				}else if(a.equals("br")){
					bestRatio=Float.parseFloat(b);
				}else if(a.equals("bb")){
					bestBad=Float.parseFloat(b);
				}else if(a.equals("sbr")){
					secondBestRatio=Float.parseFloat(b);
				}else if(a.equals("sbb")){
					secondBestBad=Float.parseFloat(b);
				}else{
					throw new RuntimeException(s);
				}
			}
			
			correct=(bestInsert==trueInsert);
			valid=(split.length==15 && bestInsert>0 && secondBestInsert>0);
			assert(!valid || bestOverlap>0) : bestOverlap+", "+bestInsert;
		}
		
		@Override
		public String toString(){
			return Tools.format("%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.8f",
					correct ? 1 : 0, minOverlap, bestOverlap, bestBadInt, secondBestOverlap, secondBestBadInt,
							expectedErrors1+expectedErrors2, bestExpected, bestRatio, bestBad, secondBestRatio, secondBestBad, probability);
		}
		
		int trueInsert;
		
		int minOverlap;
		float expectedErrors1;
		float expectedErrors2;
		
		float bestExpected;
		float probability;
		
		int bestInsert;
		int bestOverlap;
		float bestRatio;
		float bestBad;
		int bestBadInt;
		
		int secondBestInsert;
		int secondBestOverlap;
		float secondBestRatio;
		float secondBestBad;
		int secondBestBadInt;

		boolean correct;
		boolean valid=false;
		
		String line;
		
	}
	
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
