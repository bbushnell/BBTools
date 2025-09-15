package jasper;

import java.util.ArrayList;
import java.util.HashSet;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;

/**
 * Read in file of high-throughput reads and report the positions at which kmers from reference start
 * @author Jasper Toscani Field
 * @date Jun 4, 2020
 *
 */
public class KmerPosition {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		KmerPosition x=new KmerPosition(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public KmerPosition(String[] args){
		
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
				
			}else if(a.equals("ref")){
				ref=b;
					
			}else if(a.equals("k")){
					k=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
	}
	
	void process(Timer t){
		HashSet<String> kr=kmerReturn();
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start();
		}
		boolean paired=cris.paired();
		//System.out.println(paired);
		
		long readsProcessed=0, basesProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					//  *********  Process reads here  *********
					//System.out.println(r1);
					processRead(r1, kr, counts1, totalEncounter1);
					if(r1.mate!=null) {
						processRead(r2, kr, counts2, totalEncounter2);
					}
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}
		
		//System.out.println(counts);
		outputResults(counts1, totalEncounter1, counts2, totalEncounter2);
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		assert(!errorState) : "An error was encountered.";
	}
	
	private void outputResults(LongList posCounts1, LongList readCounts1, LongList posCounts2, LongList readCounts2){
		if(ffout1==null) {return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		long[] readArray1 = readCounts1.toArray();
		long[] countArray1 = posCounts1.toArray();
		long[] readArray2 = readCounts2.toArray();
		long[] countArray2 = posCounts2.toArray();
		
		bsw.println("#pos\tread1_count\tread1_perc\tread2_count\tread2_perc");
		
		int maxLen = Tools.max(readArray1.length, readArray2.length);
		
		for(int i=0; i<maxLen; i++) {
			
			bsw.print(i);
			bsw.print('\t');
			bsw.print(countArray1.length>i ? countArray1[i] : 0);
			bsw.print('\t');
			bsw.print(countArray1.length>i ? (countArray1[i] / (float) readArray1[i]) * 100 : 0, 3);
			bsw.print('\t');
			bsw.print(countArray2.length>i ? countArray2[i] : 0);
			bsw.print('\t');
			bsw.print(countArray2.length>i ? (countArray2[i] / (float) readArray2[i]) * 100 : 0, 3);
			bsw.println();
		}
		//Write stuff to the bsw
		//bsw.println("Stuff");

		errorState=bsw.poisonAndWait() | errorState;
	}
	
	private HashSet<String> kmerReturn(){
		HashSet<String> hs=new HashSet<String>();
		ArrayList<Read> readArray=ConcurrentReadInputStream.getReads(maxReads, false, ffref, null, null, null);
		for(Read r : readArray) {
			addToSet(hs, r);
		}
		
		//System.out.println(hs);
		
		return hs;
	}
	
	private int addToSet(HashSet<String> hs, Read r) {
		int countRead=0;
		for(int i=0, j=k; j<=r.length(); i++, j++) {
			//String(byte[] bytes, int offset, int length)
			String s=new String(r.bases, i, k);
			hs.add(s);
			countRead++;
		}
		return countRead;
	}
	
	private LongList processRead(Read r, HashSet<String> hs, LongList count, LongList readCount) {
		for(int i=0, j=k; j<=r.length(); i++, j++) {
			//String(byte[] bytes, int offset, int length)
			
			int x=5;
			String s=new String(r.bases, i, k);
			readCount.increment(i);
			if(hs.contains(s)) {
				count.increment(i);
			}
			//System.out.println(totalEncounter);
//			System.out.println(counts);
			
		}
		return count;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	private String out1=null;
	private String ref=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	private final FileFormat ffref;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean errorState=false;
	private int k=6;
	private LongList counts1=new LongList();
	private LongList totalEncounter1=new LongList();
	private LongList counts2=new LongList();
	private LongList totalEncounter2=new LongList();
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
