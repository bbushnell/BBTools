package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Input stream for reading sequences from GenBank format (.gbk) files.
 * Parses GenBank files to extract sequence data from ORIGIN sections.
 * Extends ReadInputStream to provide buffered reading capabilities for GenBank format.
 * @author Brian Bushnell
 */
public class GbkReadInputStream extends ReadInputStream {
	
	/**
	 * Test method for GbkReadInputStream functionality.
	 * Creates a stream from the first command-line argument and prints the first read.
	 * @param args Command-line arguments where args[0] is the input filename
	 */
	public static void main(String[] args){
		
		GbkReadInputStream fris=new GbkReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Creates a GbkReadInputStream from a filename.
	 * @param fname Input filename for GenBank format file
	 * @param allowSubprocess_ Whether to allow subprocess for decompression
	 */
	public GbkReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.GBK, null, allowSubprocess_, false));
	}
	
	/**
	 * Creates a GbkReadInputStream from a FileFormat object.
	 * Sets up buffering, amino acid flags, and validates the file format.
	 * Issues warning if file extension doesn't match GenBank format.
	 * @param ff FileFormat object containing file information and settings
	 */
	public GbkReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.gbk()){
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
	 * Fills the internal buffer with reads from the GenBank file.
	 * Parses GenBank format to extract sequences and creates Read objects.
	 * Closes the file when fewer reads than buffer size are returned.
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
	

	
	/**
	 * Parses GenBank format file and converts to Read objects.
	 * Looks for ORIGIN sections and extracts sequence data, filtering out
	 * non-letter characters and converting to uppercase.
	 *
	 * @param bf ByteFile to read from
	 * @param maxReadsToReturn Maximum number of reads to return
	 * @param numericID Starting numeric ID for reads
	 * @param flag Read flags (e.g., amino acid mask)
	 * @return ArrayList of Read objects parsed from GenBank format
	 */
	public static ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn, long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		
		int added=0;
		
		String idLine=null;
		ByteBuilder bb=new ByteBuilder();
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
//			if(Tools.startsWith(s, "ID")){
//				idLine=new String(s, 2, s.length-2).trim();
//			}else 
			if(Tools.startsWith(s, "ORIGIN")){
//				System.err.println(new String(s));
				byte[] line=null;
				for(line=bf.nextLine(); line!=null && line[0]!='/'; line=bf.nextLine()){
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
	
	@Override
	public String fname(){return bf.name();}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}

	/** Internal buffer holding parsed Read objects */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;
	
	/** ByteFile for reading from the GenBank input file */
	private final ByteFile bf;
	/** Read flags such as amino acid mask */
	private final int flag;

	/** Buffer length for reads from shared buffer settings */
	private final int BUF_LEN=Shared.bufferLen();;
	/** Maximum data size from shared buffer settings */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	/** Total number of reads generated from the file */
	public long generated=0;
	/** Total number of reads consumed by caller */
	public long consumed=0;
	/** Numeric ID for the next read to be created */
	private long nextReadID=0;
	
	/** Whether input is from standard input */
	public final boolean stdin;
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;

}
