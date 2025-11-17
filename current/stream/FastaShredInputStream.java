package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Feb 13, 2013
 *
 */
public class FastaShredInputStream extends ReadInputStream {
	
	/**
	 * Program entry point for testing FASTA shredding functionality.
	 * Processes command-line arguments for read limits and shredding parameters.
	 * @param args Command-line arguments [filename, start_index, end_index, min_length, target_length]
	 */
	public static void main(String[] args){
		
		int a=20, b=Integer.MAX_VALUE;
		if(args.length>1){a=Integer.parseInt(args[1]);}
		if(args.length>2){b=Integer.parseInt(args[2]);}
		if(args.length>3){MIN_READ_LEN=Integer.parseInt(args[3]);}
		if(args.length>4){TARGET_READ_LEN=Integer.parseInt(args[4]);}
		
		Timer t=new Timer();
		
		FastaShredInputStream fris=new FastaShredInputStream(args[0], false, false, Shared.bufferData());
		Read r=fris.nextList().get(0);
		int i=0;
		
		while(r!=null){
			if(i<a){System.out.println(r.toText(false));}
			r=fris.next();
			if(++i>=a){break;}
		}
		while(r!=null && i++<b){r=fris.next();}
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	/**
	 * Creates a FASTA shredding stream from filename with specified parameters.
	 *
	 * @param fname Input FASTA filename
	 * @param amino_ Whether sequences are amino acid (true) or nucleotide (false)
	 * @param allowSubprocess_ Whether to allow subprocess for compressed files
	 * @param maxdata Maximum data to buffer in memory
	 */
	public FastaShredInputStream(String fname, boolean amino_, boolean allowSubprocess_, long maxdata){
		this(FileFormat.testInput(fname, FileFormat.FASTA, FileFormat.FASTA, 0, allowSubprocess_, false, false), amino_, maxdata);
	}
	
	/**
	 * Creates a FASTA shredding stream from FileFormat with validation.
	 * Sets up buffering parameters and validates file format expectations.
	 *
	 * @param ff FileFormat object specifying input file details
	 * @param amino_ Whether sequences are amino acid (true) or nucleotide (false)
	 * @param maxData_ Maximum data to buffer in memory
	 */
	public FastaShredInputStream(FileFormat ff, boolean amino_, long maxData_){
		name=ff.name();
		amino=amino_;
		flag=(amino ? Read.AAMASK : 0);
		
		if(!fileIO.FileFormat.hasFastaExtension(name) && !name.startsWith("stdin")){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+name);
		}
		
		allowSubprocess=ff.allowSubprocess();
		minLen=MIN_READ_LEN;
		maxLen=TARGET_READ_LEN;
		maxData=maxData_>0 ? maxData_ : Shared.bufferData();
		
		bf=open();
		
		assert(settingsOK());
	}
	
	/**
	 * Single read access is not supported.
	 * Use nextList() for batch processing instead.
	 * @return Never returns; always throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unsupported
	 */
	public Read next() {
		return nextList().get(0);
	}
	
	@Override
	public ArrayList<Read> nextList() {
		if(currentList==null){
			boolean b=fillList();
		}
		ArrayList<Read> list=currentList;
		currentList=null;
		if(list==null || list.isEmpty()){
			list=null;
		}else{
			consumed+=list.size();
		}
		return list;
	}
	
	@Override
	public boolean hasMore() {
		if(currentList==null || currentList.size()==0){
			if(open){
				fillList();
			}else{
//				assert(generated>0) : "Was the file empty?";
			}
		}
		return (currentList!=null && currentList.size()>0);
	}
	
	@Override
	public void restart() {
		if(bf!=null){close();}
		assert(bf==null);
//		generated=0;
		consumed=0;
		nextReadID=0;
		currentList=null;
		
		if(bf==null){
			bf=open();
		}else{
			assert(false) : "bf should be null";
		}
	}
	
	@Override
	public final boolean close(){
		synchronized(this){
			if(!open){return false;}
			open=false;
			assert(bf!=null);
			errorState|=bf.close();
			bf=null;
		}
		return false;
	}
	
	@Override
	public boolean paired() {return false;}
	
	/**
	 * Fills the current read list with shredded reads up to buffer limits.
	 * Generates reads until buffer length or data limits are reached.
	 * @return true if any reads were generated, false if end of file
	 */
	private final boolean fillList(){
//		assert(open);
		if(!open){
			currentList=null;
			return false;
		}
		assert(currentList==null);
		currentList=new ArrayList<Read>(BUF_LEN);
		
		long len=0;
		for(int i=0; i<BUF_LEN && len<maxData; i++){
			Read r=generateRead();
			if(r==null){
				close();
				break;
			}
			currentList.add(r);
			len+=r.length();
			nextReadID++;
			if(verbose){System.err.println("Generated a read; i="+i+", BUF_LEN="+BUF_LEN);}
//			if(i==1){assert(false) : r.numericID+", "+r.mate.numericID;}
		}
		
		return currentList.size()>0;
	}
	
	/**
	 * Generates a single shredded read from the internal buffer.
	 * Reads lines from file until sufficient sequence is accumulated,
	 * then extracts a read of target length with specified overlap.
	 * @return Read object of target length, or null if end of file reached
	 */
	private final Read generateRead(){
		Read r=null;
		boolean eof=false;
		while(r==null && !eof){
			while(buffer.length()<maxLen){
				byte[] line=bf.nextLine();
				if(line==null){
					eof=true;
					break;
				}
				if(line.length>0 && line[0]==carrot){break;}
				buffer.append(line);
			}
			if(buffer.length>=minLen){
				byte[] bases=buffer.expelAndShift(Tools.min(maxLen, buffer.length()), TARGET_READ_OVERLAP);
				r=new Read(bases, null, Long.toString(nextReadID), nextReadID, flag);
				nextReadID++;
				if(verbose){System.err.println("Made read:\t"+(r.length()>1000 ? r.id : r.toString()));}
				if(bases.length<maxLen){buffer.clear();}
			}else{
				buffer.clear();
			}
		}
		return r;
	}
	
	/**
	 * Opens the input FASTA file for reading.
	 * @return ByteFile handle for the opened file
	 * @throws RuntimeException if file is already open
	 */
	private final ByteFile open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened fasta file "+name);
		}
		open=true;
		ByteFile bf=ByteFile.makeByteFile(name, allowSubprocess);
		return bf;
	}
	
	/** Returns whether the input stream is currently open */
	public boolean isOpen(){return open;}
	
	/**
	 * Validates that shredding parameters are within acceptable ranges.
	 * Checks minimum and target read lengths for consistency and bounds.
	 * @return true if all settings are valid, false otherwise
	 * @throws RuntimeException if settings are invalid with detailed error message
	 */
	public static final boolean settingsOK(){
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1){
			throw new RuntimeException("Minimum FASTA read length is too long: "+MIN_READ_LEN);
		}
		if(MIN_READ_LEN<1){
			throw new RuntimeException("Minimum FASTA read length is too short: "+MIN_READ_LEN);
		}
		if(TARGET_READ_LEN<1){
			throw new RuntimeException("Target FASTA read length is too short: "+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>TARGET_READ_LEN){
			throw new RuntimeException("Minimum FASTA read length is longer than maximum read length: "+MIN_READ_LEN+">"+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1 || MIN_READ_LEN<1){return false;}
		if(TARGET_READ_LEN<1 || MIN_READ_LEN>TARGET_READ_LEN){return false;}
		return true;
	}
	
	@Override
	public String fname(){return bf.name();}
	
	/** Name of the input FASTA file */
	public final String name;
	
	/** Current batch of shredded reads ready for consumption */
	private ArrayList<Read> currentList=null;

	/** Whether the input stream is currently open for reading */
	private boolean open=false;
	/** Buffer for accumulating sequence data before shredding into reads */
	private ByteBuilder buffer=new ByteBuilder();
	/** ByteFile handle for the input FASTA file */
	public ByteFile bf;
	
	/** Counter for assigning unique IDs to generated reads */
	private long nextReadID=0;
	/** Number of reads that have been consumed from the stream */
	private long consumed=0;

	/** Whether to allow subprocess execution for compressed file handling */
	public final boolean allowSubprocess;
	/** Whether sequences are amino acid (true) or nucleotide (false) */
	public final boolean amino;
	/** Read flag indicating amino acid sequences (AAMASK) or nucleotide (0) */
	public final int flag;
	/** Maximum number of reads to buffer in a single batch */
	private final int BUF_LEN=Shared.bufferLen();;
	/** Maximum amount of sequence data to buffer in memory */
	private final long maxData;
	/** Maximum length for generated reads */
	private final int maxLen, minLen;
	
	
	/** Enable verbose debugging output during read generation */
	public static boolean verbose=false;
	/** ASCII character '>' used to identify FASTA header lines */
	private static final byte carrot='>';
	
	/** Target length for shredded reads in bases */
	public static int TARGET_READ_LEN=800;
	/** Number of bases of overlap between consecutive shredded reads */
	public static int TARGET_READ_OVERLAP=31;
	/** Minimum read length required to generate a read */
	public static int MIN_READ_LEN=31;
	/** Whether to generate fake quality scores for FASTA reads */
	public static boolean FAKE_QUALITY=false;
	/** Whether to warn if FASTA entries contain no sequence data */
	public static boolean WARN_IF_NO_SEQUENCE=true;
	/** Whether to limit no-sequence warnings to first occurrence only */
	public static boolean WARN_FIRST_TIME_ONLY=true;
	/** Whether to abort processing if FASTA entries contain no sequence data */
	public static boolean ABORT_IF_NO_SEQUENCE=false;
	
}
