package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.Aligner;
import aligner.MultiStateAligner9PacBioAdapter2;
import aligner.SingleStateAlignerFlat;
import aligner.SingleStateAlignerFlat2;
import aligner.SingleStateAlignerFlat2Amino;
import aligner.SingleStateAlignerFlat2_1D;
import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SiteScore;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class FindPrimers implements Accumulator<FindPrimers.ProcessThread> {

	public static void main(String[] args){
		Timer t=new Timer();
		FindPrimers x=new FindPrimers(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FindPrimers(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBufferLen(8);
		
		float cutoff_=0;
		String literal_=null;
		String ref_=null;
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("usemsa2") || a.equals("msa2")){
				useMSA2=Parse.parseBoolean(b);
			}else if(a.equals("usessa2") || a.equals("ssa2")){
				useSSA2=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("usessa1d") || a.equalsIgnoreCase("ssa1d")){
				useSSA1D=Parse.parseBoolean(b);
			}else if(a.equals("swap")){
				swapQuery=Parse.parseBoolean(b);
			}else if(a.equals("printzeros")){
				printZeros=Parse.parseBoolean(b);
			}else if(a.equals("twocolumn") || a.equals("2column")){
				oneColumn=!Parse.parseBoolean(b);
			}else if(a.equals("onecolumn")){
				oneColumn=Parse.parseBoolean(b);
			}else if(a.equals("addr")){
				addR=Parse.parseBoolean(b);
			}else if(a.equals("replicate") || a.equals("expand")){
				replicateAmbiguous=Parse.parseBoolean(b);
			}else if(a.equals("literal")){
				literal_=b;
			}else if(a.equals("cutoff") || a.equals("minid")){
				cutoff_=Float.parseFloat(b);
				if(cutoff_>1){cutoff_=cutoff_/100;}
				assert(cutoff_>=0 && cutoff_<=1) : "Cutoff should range from 0 to 1";
			}else if(a.equals("primer") || a.equals("query") || a.equals("ref")){
				if(new File(b).exists()){ref_=b;}
				else{literal_=b;}
			}else if(a.equals("idhist")){
				outIdHist=b;
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
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
		cutoff=cutoff_;
		
//		ArrayList<byte[]> sharedHeader=new ArrayList<byte[]>();
		ByteBuilder sharedHeader=new ByteBuilder();
		sharedHeader.append("@HD\tVN:1.4\tSO:unsorted\n");
		if(ref_!=null){
			ArrayList<Read> list=FastaReadInputStream.toReads(ref_, FileFormat.FASTA, -1);
			int max=0;
			queries=new ArrayList<Read>();
//			assert(false) : list;
			for(int i=0; i<list.size(); i++){
				Read r=list.get(i);
				max=Tools.max(max, r.length());
				if(swapQuery){
					sharedHeader.append("@SQ\tSN:"+r.name()+"\tLN:"+r.length()).nl();
				}
				queries.add(r);
//				if(rcomp){
//					r=r.copy();
//					r.reverseComplement();
//					if(addR){r.id="r_"+r.id;}
//					r.setStrand(1);
//					queries.add(r);
//				}
			}
			maxqlen=max;
		}else if(literal_!=null){
			int max=0;
			String[] s2=literal_.split(",");
			queries=new ArrayList<Read>();
			for(int i=0; i<s2.length; i++){
				Read r=new Read(s2[i].getBytes(), null, "query"+i, i);
				max=Tools.max(max, r.length());
				queries.add(r);
				if(swapQuery){
					sharedHeader.append("@SQ\tSN:"+r.name()+"\tLN:"+r.length()).nl();
				}
			}
			maxqlen=max;
		}else{
			queries=null;
			maxqlen=0;
		}
		
		if(sharedHeader.length()>0) {
			ReadStreamWriter.HEADER=sharedHeader;
		}
		
		if(replicateAmbiguous){
			queries=Tools.replicateAmbiguous(queries, 1);
		}
		if(Shared.AMINO_IN){rcomp=false;}
		if(rcomp){
			for(int i=0, max=queries.size(); i<max; i++){
				Read r=queries.get(i).copy();
				r.reverseComplement();
				if(addR){r.id="r_"+r.id;}
				r.setStrand(1);
				queries.add(r);
			}
		}
		
//		assert(false) : queries;
		
		ffout1=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, "attachment", true, true, false, ordered);
//		assert(ffout1.type()==FileFormat.ATTACHMENT) : ffout1;
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();

		final ConcurrentReadOutputStream ros=makeCros();
		
//		final ByteStreamWriter bsw;
//		if(out1!=null){
//
//			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
//			
//			bsw=new ByteStreamWriter(ffout1);
//			bsw.start();
//		}else{bsw=null;}
		
		spawnThreads(cris, ros);
		
//		if(bsw!=null){bsw.poisonAndWait();}
		ReadWrite.closeStreams(cris, ros);
		
		int minid=100;
		{
			ByteBuilder bb=new ByteBuilder();
			bb.append(oneColumn ? ReadWrite.stripToCore(outIdHist) : "ID\tCount").nl();
			for(int i=0; i<idHist.length(); i++){
				int count=idHist.get(i);
				if(count>0 || printZeros) {
					if(!oneColumn){bb.append(i).tab();}
					bb.append(count).nl();
					if(count>0){minid=Tools.min(minid, i);}
				}
			}
			if(outIdHist!=null){ReadWrite.writeString(bb, outIdHist);}
		}
		
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Average Identity:   "+Tools.format("%.3f%%", (identitySum*100/identityCount)));
		outstream.println("Min Identity:       "+minid);
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, true);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		identitySum+=pt.identitySumT;
		identityCount+=pt.identityCountT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_){
			cris=cris_;
			ros=ros_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				assert(r.mate==null);
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r.length();

				//Increment counters
				readsProcessedT+=r.pairCount();
				basesProcessedT+=initialLength1;
				
				{
					//Reads are processed in this block.
					boolean keep=processRead(r);
					
					if(!keep){reads.set(idx, null);}
					else{
						readsOutT+=r.pairCount();
						basesOutT+=r.pairLength();
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processRead(final Read r){
			
			final int a=0, b=(swapQuery ? maxqlen-1 : r.length()-1);
			int[] max;
			
			
			SiteScore bestSite=null;
			Read bestQuery=null;
			for(int qnum=0; qnum<queries.size(); qnum++){
				final Read query=queries.get(qnum);
				
//				System.err.println("msa.maxColumns()="+msa.maxColumns()+", msa.maxRows()="+msa.maxRows()+
//						", maxqlen="+maxqlen+", rlen="+r.length()+", qlen="+query.length());
				
				
				//{rows, maxCol, maxState, maxScore, maxStart}
				if(swapQuery){
					max=msa.fillUnlimited(r.bases, query.bases, a, b, -999999);
				}else{
					max=msa.fillUnlimited(query.bases, r.bases, a, b, -999999);
				}
//				assert(false) : new String(r.bases)+"\n"+new String(queries.get(0).bases);
				
//				System.err.println("\nAligned: "+new String(r.bases).substring(0, 10)+"..."+new String(r.bases).substring(r.bases.length-11)+
//						", "+new String(query.bases).substring(0, 10)+"..."+new String(query.bases).substring(query.bases.length-11));
				
				if(max!=null){
					
					final int rows=max[0];
					final int maxCol=max[1];
					final int maxState=max[2];
//					final int maxScore=max[3];
//					final int maxStart=max[4];
					
					//{score, bestRefStart, bestRefStop} 
					//byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int maxRow, int maxCol, int maxState
//					int[] score=msa.score(query.bases, r.bases, a, b, max[0], max[1], max[2]);
					int[] score;
					if(swapQuery){
						score=msa.score(r.bases, query.bases, a, b, rows, maxCol, maxState);
					}else{
						//returns {score, bestRefStart, bestRefStop} 
						score=msa.score(query.bases, r.bases, a, b, rows, maxCol, maxState);
					}
					SiteScore ss=new SiteScore(1, query.strand(), score[1], score[2], 1, score[0]);
//					System.err.println("Score: "+ss.quickScore);
					if(bestSite==null || ss.quickScore>bestSite.quickScore){
						bestQuery=query;
						ss.setSlowScore(ss.quickScore);
						ss.score=ss.quickScore;
//						ss.match=msa.traceback(query.bases, r.bases, a, b, score[3], score[4], score[5], false);
						
						//(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state)
						//int refStartLoc, int refEndLoc, int row, int col, int state
						if(swapQuery){
							ss.match=msa.traceback(r.bases, query.bases, a, b, rows, maxCol, maxState);
						}else{
							ss.match=msa.traceback(query.bases, r.bases, a, b, rows, maxCol, maxState);
						}
//						System.err.println(new String(ss.match));
//						System.err.println(Read.identity(ss.match));
						bestSite=ss;
					}
				}
			}
			
			ByteBuilder bb2=new ByteBuilder(toBytes(r, bestQuery, bestSite));
			r.obj=bb2;
			
			return true;
		}
		

		/*--------------------------------------------------------------*/
		
		private ByteBuilder toBytes(Read r, Read query, SiteScore ss){
			bb.clear();
			if(ss==null){return toBytes(r, query);}
			float identity=Read.identity(ss.match);
			idHist.incrementAndGet((int)(100*identity));
			if(identity<cutoff){return toBytes(r, query);}
			if(swapQuery){
				bb.append(r.id).append('\t');
			}else{
				bb.append(query.id).append('\t');
			}
			bb.append(makeFlag(ss)).append('\t');
			if(swapQuery){
				bb.append(query.id).append('\t');
			}else{
				bb.append(r.id.replace('\t', '_')).append('\t');
			}
			bb.append(Tools.max(0, ss.start)+1).append('\t');
			bb.append(Tools.max(ss.score/query.length(), 4)).append('\t');
			String cigar;
			if(swapQuery){
				cigar=SamLine.toCigar14(ss.match, ss.start, ss.stop, query.length(), r.bases);
//				assert(SamLine.calcCigarReadLength(cigar, true, false)==r.length()) : 
//					"cigar="+SamLine.calcCigarReadLength(cigar, true, false)+", r="+r.length()+", q="+query.length()+", swap="+swapQuery+", rid="+r.id;
			}else{
				cigar=SamLine.toCigar14(ss.match, ss.start, ss.stop, r.length(), query.bases);
//				assert(SamLine.calcCigarReadLength(cigar, true, false)==query.length()) : 
//					"cigar="+SamLine.calcCigarReadLength(cigar, true, false)+", r="+r.length()+", q="+query.length()+", swap="+swapQuery+", rid="+r.id;
			}
			if(cigar==null){bb.append('*').append('\t');}else{bb.append(cigar).append('\t');}
			bb.append('*').append('\t');
			bb.append('0').append('\t');
			bb.append('0').append('\t');
			
			if(swapQuery){
				if(ss.strand()==Shared.MINUS){
					for(int i=r.bases.length-1; i>=0; i--){bb.append(AminoAcid.baseToComplementExtended[r.bases[i]]);}
					bb.tab();
				}else{
					bb.append(r.bases).append('\t');
				}
			}else{
				bb.append(query.bases).append('\t');
			}
//			assert(ss.strand()!=Shared.MINUS) : bb+"\n"+new String(r.bases);
			bb.append('*').append('\t');
			
			identitySumT+=identity;
			identityCountT++;
			bb.append("YI:f:").append(100*identity, 2);
			
			return bb;
		}
		
		private ByteBuilder toBytes(Read r, Read query){
			bb.clear();
//			if(cutoff>0){return bb;}
			if(swapQuery){
				bb.append(r.id).append('\t');
			}else{
				bb.append(query.id).append('\t');
			}
			bb.append(4).append('\t');
			bb.append('*').append('\t');
			bb.append(0).append('\t');
			bb.append(0).append('\t');
			bb.append('*').append('\t');
			bb.append('*').append('\t');
			bb.append('0').append('\t');
			bb.append('0').append('\t');
			
			if(swapQuery){
				bb.append(r.bases).append('\t');
			}else{
				bb.append(query.bases).append('\t');
			}
			bb.append('*').append('\t');
			
//			identitySumT+=f;
//			identityCountT++;
//			bb.append("YI:f:").append(100*f, 2);
			
			return bb;
		}
		
		//This seems to be up to 30% faster if the raw class rather than interface is used.
		Aligner msa=Shared.AMINO_IN ? new SingleStateAlignerFlat2Amino() : useMSA2 ? new MultiStateAligner9PacBioAdapter2() : 
			useSSA1D ? new SingleStateAlignerFlat2_1D() : useSSA2 ? new SingleStateAlignerFlat2() : new SingleStateAlignerFlat();
//		SingleStateAlignerFlat msa=new SingleStateAlignerFlat(maxqlen+5, maxqlen+5);
		
		SingleStateAlignerFlat xxx0; //Never fastest
		MultiStateAligner9PacBioAdapter2 xxx1; //Best for indels, but slowest
		SingleStateAlignerFlat2 xxx2; //Fastest for 16S
		SingleStateAlignerFlat2_1D xxx3; //Fastest for 23S
		Aligner xxx9;

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		double identitySumT=0;
		long identityCountT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		private ByteBuilder bb=new ByteBuilder(1000);
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
	}
	
	public static int makeFlag(SiteScore ss){
		int flag=0;
		if(ss.strand()==Shared.MINUS){flag|=0x10;}
		return flag;
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String outIdHist=null;
	
	private AtomicIntegerArray idHist=new AtomicIntegerArray(101);
	
	private final float cutoff;
	private boolean rcomp=true;
	private boolean replicateAmbiguous=false;
	private boolean swapQuery=false;
	private boolean addR=false;
	private boolean printZeros=true;
	private boolean oneColumn=false;
	private boolean useMSA2=false;
	private boolean useSSA2=true;
	private boolean useSSA1D=false;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private ArrayList<Read> queries;
	private final int maxqlen;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	protected long readsOut=0;
	protected long basesOut=0;

	private long maxReads=-1;
	
	double identitySum=0;
	long identityCount=0;
	
	boolean ordered=true;
	boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
