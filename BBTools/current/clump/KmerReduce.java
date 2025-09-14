package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import assemble.AbstractRemoveThread;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import kmer.KmerTableSet;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Reduces reads to their feature kmer.
 * @author Brian Bushnell
 * @date Nov 10, 2015
 *
 */
public class KmerReduce {

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		Timer t=new Timer();
		KmerReduce x=new KmerReduce(args);
		x.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * @param fname0 Input filename of reads
	 * @param k Kmer length
	 * @param cutoff Minimum count to retain
	 * @return Set of pivot kmers
	 */
	public static KmerTableSet getValidKmersFromReads(final String fname0, int k, int cutoff){
		final String fname=fname0+"_"+(new Random().nextLong()>>>1)+".fa.gz";
		assert(!new File(fname).exists());

		ArrayList<String> arglist=new ArrayList<String>();
		arglist.add("in="+fname0);
		arglist.add("out="+fname);
		arglist.add("k="+k);
		String[] args=arglist.toArray(new String[0]);

		main(args);

		assert(false) : fname+", "+k+", "+cutoff;
		KmerTableSet set=getValidKmers(fname, k, cutoff);
		File f=new File(fname);
		if(f.exists()){f.delete();}
		
		return set;
	}
	
	/**
	 * @param fname Input filename of pivot kmers
	 * @param k Kmer length
	 * @param cutoff Minimum count to retain
	 * @return Set of pivot kmers
	 */
	public static KmerTableSet getValidKmers(final String fname, int k, int cutoff){
		ArrayList<String> arglist=new ArrayList<String>();
		arglist.add("in="+fname);
		arglist.add("k="+k);
		if(cutoff>1 && prefilter){
			arglist.add("prefilter="+(cutoff-1));
		}
		
		String[] args=arglist.toArray(new String[0]);
		KmerTableSet set=new KmerTableSet(args, 12);
		
		Timer t=new Timer();
		
		set.process(t);
//		errorState|=set.errorState;
		assert(!set.errorState);
		t.stop();
		
		set.prefilterArray=null;
		AbstractRemoveThread.process(Shared.threads(), cutoff, Integer.MAX_VALUE, set, true);
		
		return set;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerReduce(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("comparisons") || a.equals("c")){
				//do nothing
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("rename") || a.equals("addname")){
				//do nothing
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//do nothing
			}else if(a.equals("condense") || a.equals("consensus")){
				//do nothing
			}else if(a.equals("correct") || a.equals("ecc")){
				//do nothing
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				//do nothing
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		
		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=Tools.max(4, Shared.threads());

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Manage threads */
	public void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, false, false);
		
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, ros, kc));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class HashThread extends Thread{

		HashThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream ros_, KmerComparator kc_){
			cris=cris_;
			ros=ros_;
			kc=kc_;
		}

		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				ArrayList<Read> out=new ArrayList<Read>(reads.size());
				for(Read r : reads){
					if(ecco && r.mate!=null){
						if(r.mate!=null){BBMerge.findOverlapStrict(r, r.mate, true);}
					}
					final long kmer=kc.hash(r, null, 0, false);
					readsProcessedT++;
					basesProcessedT+=r.length();
					if(kmer>=0){
						Read temp=new Read(toBytes(kmer), null, header, r.numericID);
						out.add(temp);
					}
				}
				if(ros!=null){ros.add(out, ln.id);}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		final ConcurrentReadInputStream cris;
		final ConcurrentReadOutputStream ros;
		final KmerComparator kc;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		private static final String header="1";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public byte[] toBytes(final long kmer){
		byte[] dest=KillSwitch.allocByte1D(k);
		fill(kmer, dest, 0);
		return dest;
	}
	
	public void fill(final long kmer, final byte[] dest, int pos){
		for(int i=k-1; i>=0; i--, pos++){
			int x=(int)((kmer>>(2*i))&3);
			dest[pos]=AminoAcid.numberToBase[x];
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int k=31;
	static boolean prefilter=true;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	protected boolean ecco=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
