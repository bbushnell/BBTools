package stream;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import structures.ListNum;

/**
 * Abstract superclass of all ConcurrentReadStreamInterface implementations.
 * ConcurrentReadInputStreams allow paired reads from twin files to be treated as a single stream.
 * @author Brian Bushnell
 * @date Nov 26, 2014
 *
 */
public abstract class ConcurrentReadInputStream implements ConcurrentReadStreamInterface {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected ConcurrentReadInputStream(String fname_){fname=fname_;}
	
	/**
	 * Special method for testing, handles some parsing.
	 * Used by MultiCros.
	 * Not really necessary; safe to remove.
	 */
	protected static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, boolean allowSubprocess, String...args){
		assert(args.length>0) : Arrays.toString(args);
		for(int i=0; i<args.length; i++){
			if("null".equalsIgnoreCase(args[i])){args[i]=null;}
		}
		assert(args[0]!=null) : Arrays.toString(args);
		
		assert(args.length<2 || !args[0].equalsIgnoreCase(args[1]));
		String in1=args[0], in2=null, qf1=null, qf2=null;
		if(args.length>1){in2=args[1];}
		if(args.length>2){qf1=args[2];}
		if(args.length>3){qf2=args[3];}

		final FileFormat ff1=FileFormat.testInput(in1, null, allowSubprocess);
		final FileFormat ff2=FileFormat.testInput(in2, null, allowSubprocess);
		
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2);
	}
	
	/** @See primary method */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/** @See primary method */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2,
			final boolean mpi, final boolean keepAll){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null, mpi, keepAll);
	}
	
//	/** @See primary method */
//	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, String qf1){
//		return getReadInputStream(maxReads, keepSamHeader, ff1, null, qf1, null, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
//	}
	
	/** @See primary method */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/**
	 * @param maxReads Quit producing after this many reads (or pairs)
	 * @param keepSamHeader If the input is sam, store the header in the static shared header object 
	 * @param ff1 Read 1 file (required)
	 * @param ff2 Read 2 file (optional)
	 * @param qf1 Qual file 1 (optional)
	 * @param qf2 Qual file 2 (optional)
	 * @param mpi True if MPI will be used
	 * @param keepAll In MPI mode, tells this stream to keep all reads instead of just a fraction
	 * @return
	 */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2, final boolean mpi, final boolean keepAll){
		if(mpi){
			final int rank=Shared.MPI_RANK;
			final ConcurrentReadInputStream cris0;
			if(rank==0){
				cris0=getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, false, true);
				cris0.start();
			}else{
				cris0=null;
			}
			final ConcurrentReadInputStream crisD;
			if(Shared.USE_CRISMPI){
				assert(false) : "To support MPI, uncomment this.";
//				crisD=new ConcurrentReadInputStreamMPI(cris0, rank==0, keepAll);
				crisD=null;
			}else{
				crisD=new ConcurrentReadInputStreamD(cris0, rank==0, keepAll);
			}
			return crisD;
		}
		
		assert(ff1!=null);
		assert(ff2==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(ff2.name())) : ff1.name()+", "+ff2.name();
		assert(qf1==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(qf2));
		assert(qf1==null || qf2==null || qf1.equalsIgnoreCase(qf2));
		
		final ConcurrentReadInputStream cris;
		
		if(ff1.fastq()){
			
			ReadInputStream ris1, ris2;
			
			ris1=new FastqReadInputStream(ff1);
			try {
				ris2=(ff2==null ? null : new FastqReadInputStream(ff2));
			} catch (AssertionError e) {//Handles problems with quality score autodetection
				ris1.close();
				throw e;
			}
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.oneline()){
			
			ReadInputStream ris1=new OnelineReadInputStream(ff1);
			ReadInputStream ris2=(ff2==null ? null : new OnelineReadInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

		}else if(ff1.fasta()){
			
			ReadInputStream ris1;
			ReadInputStream ris2;
			if(ff1.preferShreds()){
				ris1=new FastaShredInputStream(ff1, Shared.AMINO_IN, ff2==null ? Shared.bufferData() : -1);
				ris2=(ff2==null ? null : new FastaShredInputStream(ff2, Shared.AMINO_IN, -1));
			}else{
				ris1=(qf1==null ? new FastaReadInputStream(ff1, (FASTQ.FORCE_INTERLEAVED && ff2==null), Shared.AMINO_IN, ff2==null ? Shared.bufferData() : -1)
						: new FastaQualReadInputStream(ff1, qf1));
				ris2=(ff2==null ? null : qf2==null ? new FastaReadInputStream(ff2, false, Shared.AMINO_IN, -1) : new FastaQualReadInputStream(ff2, qf2));
			}
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
//			cris.start();
//			ListNum<Read> ln=cris.nextList();
//			System.out.println(ln);
//			
//			assert(false) : ff1+", "+ff2;
		}else if(ff1.scarf()){
			
			ReadInputStream ris1=new ScarfReadInputStream(ff1);
			ReadInputStream ris2=(ff2==null ? null : new ScarfReadInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.samOrBam()){
			
			ReadInputStream ris1=new SamReadInputStream(ff1, keepSamHeader, FASTQ.FORCE_INTERLEAVED);
			ReadInputStream ris2=(ff2==null ? null : new SamReadInputStream(ff2, false, false));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.bread()){
//			assert(false) : ff1;
			RTextInputStream rtis=new RTextInputStream(ff1, ff2, maxReads);
			cris=new ConcurrentLegacyReadInputStream(rtis, maxReads); //TODO: Change to generic
			
		}else if(ff1.header()){
			
			HeaderInputStream ris1=new HeaderInputStream(ff1);
			HeaderInputStream ris2=(ff2==null ? null : new HeaderInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.sequential()){
			
			SequentialReadInputStream ris=new SequentialReadInputStream(maxReads, 200, 50, 0, false);
			cris=new ConcurrentLegacyReadInputStream(ris, maxReads);
			
		}else if(ff1.csfasta()){
			
			throw new RuntimeException("csfasta is no longer supported.");
			
		}else if(ff1.random()){
			
			RandomReadInputStream3 ris=new RandomReadInputStream3(maxReads, FASTQ.FORCE_INTERLEAVED);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);
			
		}else if(ff1.embl()){
			
			EmblReadInputStream ris=new EmblReadInputStream(ff1);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);
			
		}else if(ff1.gbk()){
			
			GbkReadInputStream ris=new GbkReadInputStream(ff1);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);
			
		}else{
			cris=null;
			throw new RuntimeException(""+ff1);
		}
		
		return cris;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ArrayList<Read> getReads(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		ConcurrentReadInputStream cris=getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
		cris.start();
		return cris.getReads();
	}
	
	public ArrayList<Read> getReads(){
		
		ListNum<Read> ln=nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		ArrayList<Read> out=new ArrayList<Read>();
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			out.addAll(reads);
			returnList(ln.id, ln.list.isEmpty());
			ln=nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		boolean error=ReadWrite.closeStream(this);
		if(error){
			System.err.println("Warning - an error was encountered during read input.");
		}
		return out;
	}
	
	@Override
	public void start(){
//		System.err.println("Starting "+this);
		new Thread(this).start(); //Prevents a strange deadlock in ConcurrentCollectionReadInputStream
		started=true;
	}
	
	public final boolean started(){return started;}

	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public abstract ListNum<Read> nextList();
	
	@Override
	public final void returnList(ListNum<Read> ln){
		if(ln!=null){returnList(ln.id, ln.isEmpty());}
	}
	
	@Override
	public abstract void returnList(long listNum, boolean poison);
	
	@Override
	public abstract void run();
	
	@Override
	public abstract void shutdown();
	
	@Override
	public abstract void restart();
	
	@Override
	public abstract void close();

	@Override
	public abstract boolean paired();

	@Override
	public String fname(){return fname;}
	
	@Override
	public abstract Object[] producers();
	
	@Override
	public abstract boolean errorState();
	
	@Override
	public abstract void setSampleRate(float rate, long seed);
	
	@Override
	public abstract long basesIn();
	
	@Override
	public abstract long readsIn();
	
	@Override
	public abstract boolean verbose();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final int BUF_LEN=Shared.bufferLen();;
	final int NUM_BUFFS=Shared.numBuffers();
	final long MAX_DATA=Shared.bufferData();
	public final String fname;
	public boolean ALLOW_UNEQUAL_LENGTHS=false;
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean SHOW_PROGRESS=false;
	public static boolean SHOW_PROGRESS2=false; //Indicate time in seconds between dots.
	public static long PROGRESS_INCR=1000000;
	public static boolean REMOVE_DISCARDED_READS=false;
	
}
