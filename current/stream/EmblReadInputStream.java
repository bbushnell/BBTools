package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Input stream for reading biological sequences from EMBL format files.
 * Parses ID lines for sequence identifiers and SQ sections for sequence data.
 * Converts sequence characters to uppercase and filters non-letter characters.
 * @author Brian Bushnell
 */
public class EmblReadInputStream extends ReadInputStream {
	
	/** Test method that reads the first sequence from an EMBL file.
	 * @param args Command-line arguments where args[0] is the EMBL file path */
	public static void main(String[] args){
		
		EmblReadInputStream fris=new EmblReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Constructs an EMBL read input stream from a filename.
	 * @param fname Path to the EMBL format file
	 * @param allowSubprocess_ Whether to allow subprocess execution for compressed files
	 */
	public EmblReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.EMBL, null, allowSubprocess_, false));
	}
	
	/**
	 * Constructs an EMBL read input stream from a FileFormat object.
	 * Initializes the ByteFile and sets sequence flags for amino acid detection.
	 * @param ff FileFormat object specifying the input file and options
	 */
	public EmblReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.embl()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		bf=ByteFile.makeByteFile(ff);
//		assert(false) : interleaved;
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(bf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}
	
	/**
	 * Fills the internal buffer with reads from the EMBL file.
	 * Closes the file if fewer reads than expected are returned.
	 * Updates generated count and sets error state if buffer is null.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=toReadList(bf, BUF_LEN, nextReadID, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){bf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}
	}
	

	
	public static ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn, long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		
		int added=0;
		
		String idLine=null;
		ByteBuilder bb=new ByteBuilder();
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
			if(Tools.startsWith(s, "ID")){
				idLine=new String(s, 2, s.length-2).trim();
//				System.err.println(idLine);
			}else if(Tools.startsWith(s, "SQ")){
//				System.err.println(new String(s));
				byte[] line=null;
				for(line=bf.nextLine(); line!=null && line[0]==' '; line=bf.nextLine()){
					for(byte b : line){
						if(Tools.isLetter(b)){
							bb.append(Tools.toUpperCase(b));
						}
					}
				}
				assert(line==null || Tools.startsWith(line, "//")) : new String(line);

				Read r=new Read(bb.toBytes(), null, idLine==null ? ""+numericID : idLine, numericID, flag);
				list.add(r);
				added++;
				numericID++;
				
				bb.clear();
				idLine=null;
				
				if(added>=maxReadsToReturn){break;}
			}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		errorState|=bf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		bf.reset();
	}

	@Override
	public boolean paired() {return false;}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}
	
	@Override
	public String fname(){return bf.name();}

	/** Internal buffer for storing reads before they are consumed */
	private ArrayList<Read> buffer=null;
	/** Index of the next read to return from the buffer */
	private int next=0;
	
	/** ByteFile object for reading the input file */
	private final ByteFile bf;
	/** Bit flags for read properties such as amino acid detection */
	private final int flag;

	/** Maximum number of reads to store in buffer at once */
	private final int BUF_LEN=Shared.bufferLen();;
	/** Maximum amount of data to buffer (unused in current implementation) */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	/** Total number of reads generated from the input file */
	public long generated=0;
	/** Total number of reads consumed by the caller */
	public long consumed=0;
	/** Numeric ID to assign to the next read */
	private long nextReadID=0;
	
	/** True if reading from standard input, false if reading from file */
	public final boolean stdin;
	/** Controls verbose output for debugging purposes */
	public static boolean verbose=false;

}
