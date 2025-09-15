/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jgi;

import static jgi.Dedupe.maxNmer;
import static jgi.Dedupe.nmerIndex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 *
 * @author syao
 * Last updated : 01102018
 */
public class TetramerFrequencies_old {
	public static void main(String[] args){

		System.out.println("Start Tetramer Frequencies analysis ...");
		
		final int oldThreads=Shared.threads();
		Shared.capThreads(16);

		Timer t=new Timer();

		//Create an instance of this class
		TetramerFrequencies_old x=new TetramerFrequencies_old(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		Shared.setThreads(oldThreads);
	}

	public TetramerFrequencies_old(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		int k_=4;
		Parser parser=new Parser();
		for (String arg : args) {
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if (a.equals("s") || a.equals("step")){
				step = Integer.parseInt(b);
			} else if (a.equals("w") || a.equals("window")){
				winSize = Integer.parseInt(b);
			} else if (a.equals("out") || a.equals("freq")){
				out1 = b;
			} else if (a.equals("dropshort")){
				keepShort=!Parse.parseBoolean(b);
			} else if (a.equals("keepshort") || a.equals("short")){
				keepShort=Parse.parseBoolean(b);
			} else if (a.equalsIgnoreCase("gc")){
				printGC=Parse.parseBoolean(b);
			} else if (a.equals("k")){
				k_ = Integer.parseInt(b);
			} else if(parser.parse(arg, a, b)){
				//do nothing
			} else {
				throw new RuntimeException("Unknown argument "+arg);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			in1=parser.in1;
		}

		k=k_;
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		assert(ffin1!=null) : "No input file.";
		assert(ffin1.exists() && ffin1.canRead()) : "Cannot read input file "+in1+".";


		// initialize the class member textstream writer here, so no need of keep results in memory
		if (out1==null || out1.equals("")){
			out1 = "stdout";
		}

		threads=Tools.max(1, Shared.threads());
		inq=new ArrayBlockingQueue<Line>(threads+1);
		
		FileFormat ff=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, true, false, true);
		bsw = new ByteStreamWriter(ff);
		bsw.start();

		//create and write the output header
		List<String> head = new ArrayList<String>();
		head.add("scaffold");
		head.add("range");
		if(printGC) {head.add("GC");}

		// the unit tetramer strings in header
		tetramerGen2(head);
		bsw.add(new ByteBuilder(String.join("\t", head)).append('\n'), nextID);
		nextID++;
	}

	void process(Timer t){
		
		ArrayList<PrintThread> alpt=spawnThreads();

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		long readsProcessed=0;
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				for (Read r1 : reads){
					windowedTetramerProfile(r1.bases, r1.id);
					readsProcessed++;
				}

				cris.returnList(ln);
				if(verbose){
					outstream.println("Returned a list.");
				}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		putLine(POISON_LINE);
		waitForFinish(alpt);

		// close the output file
		bsw.poisonAndWait();

		ReadWrite.closeStream(cris);
		if(verbose){outstream.println("Finished.");}

		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}

	// TODO: optimize the algorithm ...
	private void windowedTetramerProfileOpt(byte[] bases, String header){
		int sidx = 0;
		int eidx = winSize;

		final int[] counts = new int[maxNmer+1];

		int leadlen = 0;
		int leadtmer = -1;
		int endlen = 0;
		int endtmer = -1;

		int tmer = 0;
		int len = 0;

		for (int i=0; i<= eidx; i++){
			int x = AminoAcid.baseToNumber[bases[i]];
			if (x < 0){
				tmer = 0;
				len = 0;
			} else {
				tmer = (leadtmer<<2)|x ;   // this can produce -1 vaue if any base in tetramer is N!
				len++;
				if (len >= k){
					int idx = tmerIndex(tmer);
					counts[idx]++; 

					if( leadtmer==-1){
						leadtmer = tmer;   
					} else if(i == eidx){

					}


					len--;
				}   
			}
		}



	}

	private void windowedTetramerProfile(byte[] bases, String header){
		int sidx = 0;
		int eidx = winSize<1 ? bases.length : Tools.min(bases.length, winSize);
		
		while (eidx <= bases.length){
			Line line=new Line(header, bases, sidx, eidx, nextID);
			putLine(line);
			sidx+=windowsPerLine*step;
			eidx+=windowsPerLine*step;
			nextID++;
		}
	}
	
	//Slow version
	void append(Line line, int[] counts, StringBuilder sb){
		sb.append(line.header);
		sb.append('\t');
		sb.append(line.sidx+1);
		sb.append('-');
		sb.append(line.eidx);
		float mult=1f/(line.length()-k+1);
		for (int cnt: counts){
			sb.append('\t');
			sb.append(Tools.format("%.4f", cnt*mult));
		}
		sb.append('\n');
	}
	
	void append(Line line, int[] counts, float gc, ByteBuilder bb){
		bb.append(line.header);
		bb.tab();
		bb.append(line.sidx+1);
		bb.append('-');
		bb.append(line.eidx);
		if(printGC) {
			bb.tab().append(gc, 4);
		}
		for (int cnt: counts){
			bb.tab();
			bb.append(cnt);
		}
		bb.nl();
	}
	
//	// factor this out so we can work on reads
//	public int[] tetramerCounter(byte[] bases, int startidx, int endidx, int[] counts){
//		if (counts == null){
//			counts = new int[maxNmer+1];
//		}
//		
//		for (int i=startidx; i<=endidx-k; i++){ 
//			int kmer = 0;   // the binary representation of the tetramer in the sliding window
//			for (int j=0; j<k; j++){
//				int x = AminoAcid.baseToNumber[bases[i+j]];
//				kmer= (kmer<<2)|x ;   // this can produce -1 vaue if any base in tetramer is N!
//			}
//
//			// ignore tetramers containing "N" base
//			if (kmer > -1){
//				int idx = tmerIndex(kmer);
//				counts[idx]++; 
//			}
//		}
//
//		return counts;
//	}
	
	// factor this out so we can work on reads
	public int[] tetramerCounter(byte[] bases, int startidx, int endidx, int[] counts, float[] gc){
		if (counts == null){
			counts = new int[maxNmer+1];
		}
		final int mask=(1<<k)-1;
		final int[] acgtn=new int[5];
		
		int len=0, kmer=0;
		for (int i=startidx; i<endidx; i++){
			int x = AminoAcid.baseToNumberACGTother[bases[i]];
			acgtn[x]++;
			kmer=((kmer<<2)|x)&mask;   // this can produce -1 vaue if any base in tetramer is N!
			if(x<4) {
				len++;
				// ignore tetramers containing "N" base
				if (len>=k){
					int idx = tmerIndex(kmer);
					counts[idx]++; 
				}
			}else {
				len=0;
				kmer=0;
			}
		}

		int gcCount=acgtn[1]+acgtn[2];
		int atCount=acgtn[0]+acgtn[3];
		gc[0]=gcCount/Math.max(1f, (gcCount+atCount));
		return counts;
	}

	private int tmerIndex(int tmer){
		int rtmer = AminoAcid.reverseComplementBinaryFast(tmer, k); 
		int mtmer = Tools.min(tmer, rtmer);
		int idx = -1;
		try{
			idx = nmerIndex[rtmer];
		} catch (ArrayIndexOutOfBoundsException e){
			System.out.printf("ArrayIndexOutOfBoundsException tmer=%d: rtmer=%d; \n", tmer, rtmer);
			System.exit(-1);
		}
		return idx;
	}


	public List<String> tetramerGen(List<String> tlist){
		if (tlist == null){
			tlist = new ArrayList<String>();
		}

		int[] bases = {0, 1, 2, 3}; // A C G T

		for (int b1 : bases){
			for (int b2 : bases){
				for (int b3 : bases){
					for (int b4 : bases){
						int tcode = b1;
						tcode = (tcode<<2)|b2;
						tcode = (tcode<<2)|b3;
						tcode = (tcode<<2)|b4;
						String tmer = AminoAcid.kmerToString(tcode, k).toLowerCase();
						String rtmer = AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(tcode, k), k).toLowerCase(); 

						if (!tlist.contains(rtmer)){
							tlist.add(tmer);
						} 
					}
				}
			} 
		}
		// total of 136 unique tetramers
		return tlist;
	}

	public List<String> tetramerGen2(List<String> tlist){
		if (tlist == null){
			tlist = new ArrayList<String>();
		}

		final int max=(1<<(2*k))-1;

		for(int i=0; i<=max; i++){
			final int a=i, b=AminoAcid.reverseComplementBinaryFast(i, k);
			int min=Tools.min(a, b);
			if(min==a){
				tlist.add(AminoAcid.kmerToString(a, k).toLowerCase());
			}
		}
		return tlist;
	}

	public void printTetramerFromCode(long code){
		System.out.println(AminoAcid.kmerToString(code, k));
	}

//	public void printTetramers(byte[] bases, String header){
//		int sidx = 0;
//		int eidx = winSize;
//
//		while (true){
//			if (eidx > bases.length){
//				break;
//			}
//			System.out.printf("%s %d-%d\n", header, sidx+1, eidx);
//
//			// count tetramer here
//			for (int i=sidx; i<eidx-k; i++){
//				char[] tetramer = new char[k];
//				for (int j=0; j<k; j++){
//					tetramer[j] = (char)bases[i+j];
//				}
//				System.out.println(new String(tetramer));
//			}
//
//			sidx+=step;
//			eidx+=step;
//		}
//	}

	public static void printHelp(){
		List<String> helplist = new ArrayList<String>();
		helplist.add("Program Name : TetramerFrequencies v1.1");
		helplist.add("Usage : ");
		helplist.add(" -h : this page");
		helplist.add(" -s : step [500]");
		helplist.add(" -w : window size [2000]. If set to 0 the whole sequence is processed");
		System.out.println(String.join("\n", helplist));
	}

	/*--------------------------------------------------------------*/
	
	final Line takeLine(){
		Line line=null;
		while(line==null){
			try {
				line=inq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		System.err.println("takeLine("+line.id+")");
		return line;
	}
	
	final void putLine(Line line){
//		System.err.println("putLine("+line.id+")");
		while(line!=null){
			try {
				inq.put(line);
				line=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Spawn process threads */
	private ArrayList<PrintThread> spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Fill a list with PrintThreads
		ArrayList<PrintThread> alpt=new ArrayList<PrintThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new PrintThread());
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(PrintThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	private void waitForFinish(ArrayList<PrintThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(PrintThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
	}
	
	private class PrintThread extends Thread {
		
		PrintThread(){}
		
		@Override
		public void run(){
			Line line=takeLine();
			while(line!=null && line!=POISON_LINE){
				processLine(line);
				line=takeLine();
			}
			putLine(POISON_LINE);
		}
		
		private void processLine(Line line){
			ByteBuilder bb=new ByteBuilder(512);
			for(int i=0; i<windowsPerLine && line.eidx<=line.bases.length; i++){
				Arrays.fill(counts, 0);
				tetramerCounter(line.bases, line.sidx, line.eidx, counts, gc);
				if(keepShort || line.length()>=winSize){
					append(line, counts, gc[0], bb);
				}
				line.sidx+=step;
				line.eidx+=step;
			}
			bsw.add(bb, line.id);
		}
		
		private final int[] counts=new int[maxNmer+1];
		private final float[] gc=new float[1];
	}
	
	private class Line {
		
		Line(String header_, byte[] bases_, int sidx_, int eidx_, long id_){
			header=header_;
			bases=bases_;
			sidx=sidx_;
			eidx=eidx_;
			id=id_;
		}
		
		public int length() {
			return eidx-sidx+1;
		}

		final String header;
		final byte[] bases;
		int sidx;
		int eidx;
		final long id;
		
	}
	
	/*--------------------------------------------------------------*/
	
	final Line POISON_LINE=new Line(null, null, -1, -1, -1);
	private final ArrayBlockingQueue<Line> inq;
	private final int threads;
	private long nextID=0;
	
	/*--------------------------------------------------------------*/

	private String in1 = null;
	private String out1 = null;
	private ByteStreamWriter bsw = null; // for output

	private final FileFormat ffin1;

	private java.io.PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	int step = 500;
	private int winSize = 2000;
	private final int k;
	private boolean keepShort=false;
	private boolean printGC=false;

	/*--------------------------------------------------------------*/

	private final static int windowsPerLine=8;
	public static boolean verbose=false;
}
