package jgi;

import java.util.ArrayList;

import align2.BandedAligner;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentCollectionReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;

/**
 * Calculates an all-to-all identity matrix.
 * @author Brian Bushnell
 * @date Nov 23, 2014
 *
 */
public class IdentityMatrix {

	/**
	 * Program entry point.
	 * Creates an IdentityMatrix instance and processes sequences to generate the matrix.
	 * @param args Command-line arguments for input/output files and parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		IdentityMatrix x=new IdentityMatrix(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs IdentityMatrix and parses command-line arguments.
	 * Sets up file formats, alignment parameters, and output formatting options.
	 * @param args Command-line arguments including input files, output paths,
	 * alignment parameters (edits, width), and formatting options
	 */
	public IdentityMatrix(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		FileFormat.PRINT_WARNING=false;
		int maxEdits_=-1;
		int maxWidth_=-1;
		
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
			}else if(a.equals("edits") || a.equals("maxedits")){
				maxEdits_=Integer.parseInt(b);
			}else if(a.equals("width") || a.equals("maxwidth")){
				maxWidth_=Integer.parseInt(b);
			}else if(a.equals("percent")){
				percent=Parse.parseBoolean(b);
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
		FASTQ.FORCE_INTERLEAVED=false;
		FASTQ.TEST_INTERLEAVED=false;
		
		maxEdits=maxEdits_==-1 ? BandedAligner.big : maxEdits_;
		maxWidth=maxWidth_==-1 ? (int)(Tools.min(maxEdits, BandedAligner.big)*2L+1) : maxWidth_;
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	/**
	 * Main processing method that generates the all-to-all identity matrix.
	 * Loads all sequences, creates worker threads to perform pairwise alignments,
	 * computes similarity scores, and outputs the final matrix.
	 * @param t Timer for tracking execution time and reporting performance metrics
	 */
	void process(Timer t){
		
		allReads=load();
		Shared.setBufferLen(4);
		ConcurrentCollectionReadInputStream cris=new ConcurrentCollectionReadInputStream(allReads, null, -1);
		cris.start(); //4567
		
		
		ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>();
		final int tmax=Tools.max(Shared.threads(), 1);
		for(int i=0; i<tmax; i++){
			threads.add(new ProcessThread(cris));
		}
		for(ProcessThread pt : threads){pt.start();}
		for(ProcessThread pt : threads){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		ReadWrite.closeStreams(cris);
		
		final int numReads=allReads.size();
		for(int i=1; i<numReads; i++){
			Read r1=allReads.get(i);
			assert(r1.numericID==i);
			for(int j=0; j<i; j++){
				Read r2=allReads.get(j);
				assert(r2.numericID==j);
				((float[])r2.obj)[i]=((float[])r1.obj)[j];
			}
		}
		
		if(ffout1!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffout1);
			tsw.start();
			for(Read r : allReads){
				float[] obj=(float[])r.obj;
				tsw.print(r.id);
				if(percent){
					for(float f : obj){
						tsw.print(Tools.format("\t%.2f", f));
					}
				}else{
					for(float f : obj){
						tsw.print(Tools.format("\t%.4f", f));
					}
				}
				tsw.print("\n");
				r.obj=null;
			}
			tsw.poisonAndWait();
		}
		
		t.stop();
		outstream.println("Total Time:                   \t"+t);
		outstream.println("Reads Processed:    "+allReads.size()+" \t"+Tools.format("%.2fk alignments/sec", (allReads.size()*(long)(allReads.size())/(double)(t.elapsed))*1000000));
		outstream.println("Min Similarity:     "+Tools.format("%.5f", minID));
		outstream.println("Max Similarity:     "+Tools.format("%.5f", maxID));
		outstream.println("Avg Similarity:     "+Tools.format("%.5f", avgID));
	}
	
	/**
	 * Loads all sequences from the input file into memory.
	 * Reads sequences using ConcurrentReadInputStream and stores them in an ArrayList.
	 * Tracks the longest sequence length for alignment parameter optimization.
	 * @return ArrayList containing all loaded sequences with numeric IDs assigned
	 */
	private ArrayList<Read> load(){
		Timer t=new Timer();
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		assert(!paired) : "This program is not designed for paired reads.";
		
		long readsProcessed=0;
		int maxLen=0;
		ArrayList<Read> bigList=new ArrayList<Read>();
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
					
					bigList.add(r1);
					maxLen=Tools.max(maxLen, r1.length());
					
					readsProcessed++;
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
		ReadWrite.closeStreams(cris);
		if(verbose){outstream.println("Finished loading "+readsProcessed+" sequences.");}
		
		longestSequence=maxLen;
		
		t.stop();
		outstream.println("Load Time:                    \t"+t);
		
		return bigList;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread that performs pairwise alignments between sequences.
	 * Each thread processes a subset of sequences and computes similarity scores
	 * using BandedAligner for efficient alignment with edit distance calculation.
	 */
	private class ProcessThread extends Thread {
		
		/**
		 * Constructs a processing thread with optimized alignment parameters.
		 * Creates a BandedAligner instance with width calculated from maximum edits
		 * and longest sequence length for efficient alignment processing.
		 * @param cris_ Input stream for reading sequences to process
		 */
		ProcessThread(ConcurrentReadInputStream cris_){
			cris=cris_;
			maxEdits2=Tools.min(maxEdits, longestSequence);
			int width=Tools.min(maxEdits2*2+1, maxWidth);
			bandy=BandedAligner.makeBandedAligner(width);
		}
		
		@Override
		public void run(){
			final int numReads=allReads.size();
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			double sum=0;
			long compares=0;
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					float[] obj=new float[numReads];
					r1.obj=obj;
					for(Read r2 : allReads){
						if(r2.numericID>r1.numericID){break;}
//						int edits=bandy.alignQuadruple(r1.bases, r2.bases, maxEdits2, false);
						int edits=bandy.alignQuadrupleProgressive(r1.bases, r2.bases, 10, maxEdits2, false);
						System.err.println(r1.id+"->"+r2.id+": Edits="+edits);
						float editRate=edits/(float)Tools.max(r1.length(), r2.length());
						float similarity=1-editRate;
						if(r1!=r2){
							compares++;
							sum+=similarity;
							minID=Tools.min(minID, similarity);
							maxID=Tools.max(maxID, similarity);
						}
						if(percent){
							float id=100*similarity;
							obj[(int)r2.numericID]=id;
						}else{
							obj[(int)r2.numericID]=similarity;
						}
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
			
			avgID=sum/compares;
		}
		
		/** Input stream for reading sequences in this processing thread */
		private final ConcurrentReadInputStream cris;
		/**
		 * Alignment engine optimized for calculating edit distances between sequences
		 */
		private final BandedAligner bandy;
		/**
		 * Thread-local maximum edit distance limit, constrained by longest sequence length
		 */
		private final int maxEdits2;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Input file path containing sequences for matrix calculation */
	private String in1=null;
	/** Output file path for writing the identity matrix */
	private String out1=null;
	
	/** File format specification for input sequence file */
	private final FileFormat ffin1;
	/** File format specification for output matrix file */
	private final FileFormat ffout1;
	/**
	 * Whether to output similarity values as percentages (true) or fractions (false)
	 */
	private boolean percent=false;
	
	/** Container holding all loaded sequences for matrix calculation */
	private ArrayList<Read> allReads;
	
	/*--------------------------------------------------------------*/
	
	/** Maximum number of reads to process, or -1 for unlimited */
	private long maxReads=-1;
	/** Maximum edit distance allowed in alignments before giving up */
	private final int maxEdits;
	/** Maximum alignment band width for BandedAligner optimization */
	private final int maxWidth;
	/** Length of the longest sequence in the input set */
	private int longestSequence;
	
	/** Average similarity value across all pairwise comparisons */
	/** Maximum similarity value observed across all pairwise comparisons */
	/** Minimum similarity value observed across all pairwise comparisons */
	private double minID=1, maxID=0, avgID=0;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for logging and status messages */
	private java.io.PrintStream outstream=System.err;
	/** Enable verbose logging of processing steps and debugging information */
	public static boolean verbose=false;
	
}
