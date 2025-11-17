package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Input stream for reading FASTA sequences with separate quality files.
 * Combines sequence data from FASTA files with quality scores from
 * corresponding quality files, creating complete Read objects.
 * Supports both numeric and ASCII-encoded quality scores.
 *
 * @author Brian Bushnell
 */
public class FastaQualReadInputStream extends ReadInputStream {
	
	/** Program entry point for testing.
	 * @param args Command-line arguments [sequence_file, quality_file] */
	public static void main(String[] args){
		
		FastaQualReadInputStream fris=new FastaQualReadInputStream(args[0], args[1], true);
		
		Read r=fris.nextList().get(0);
		int i=0;
		while(r!=null){
			System.out.println(r.toText(false));
			r=fris.nextList().get(0);
			if(i++>3){break;}
		}
		
	}
	
	/**
	 * Constructs input stream from sequence and quality file names.
	 * @param fname FASTA sequence file name
	 * @param qfname Quality scores file name
	 * @param allowSubprocess_ Whether to allow subprocess reading
	 */
	public FastaQualReadInputStream(String fname, String qfname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false), qfname);
	}
	
	/**
	 * Constructs input stream from FileFormat and quality file name.
	 * Validates that the FileFormat represents a FASTA file.
	 * @param ff FileFormat for the sequence file
	 * @param qfname Quality scores file name
	 */
	public FastaQualReadInputStream(FileFormat ff, String qfname){
		
		if(!ff.fasta() && !ff.stdio()){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+ff.name());
		}
		
		btf=ByteFile.makeByteFile(ff);
		qtf=ByteFile.makeByteFile(FileFormat.testInput(qfname, FileFormat.QUAL, null, ff.allowSubprocess(), false));
		interleaved=false;
	}
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(btf.isOpen()){
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
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return r;
	}
	
	/**
	 * Fills the internal read buffer by parsing sequences and quality scores.
	 * Closes input streams when fewer reads than buffer size are read.
	 * Updates generation statistics and error state.
	 */
	private synchronized void fillBuffer(){
		if(builder==null){builder=new ByteBuilder(2000);}
		if(verbose){System.err.println("Filling buffer.  buffer="+(buffer==null ? null : buffer.size()));}
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReads(BUF_LEN, nextReadID, interleaved);
		final int count=(buffer==null ? 0 : buffer.size());

		if(verbose){System.err.println("Filled buffer.  size="+count);}
		
		nextReadID+=count;
		if(count<BUF_LEN){
			if(verbose){System.err.println("Closing tf");}
			errorState|=close();
		}
		
		generated+=count;
		if(verbose){System.err.println("generated="+generated);}
	}
	
	/**
	 * Wrapper method for toReadList with assertion checks.
	 *
	 * @param maxReadsToReturn Maximum number of reads to return
	 * @param numericID Starting numeric ID for reads
	 * @param interleaved Whether reads are paired/interleaved
	 * @return ArrayList of Read objects
	 */
	private ArrayList<Read> toReads(int maxReadsToReturn, long numericID, boolean interleaved){
		ArrayList<Read> list=toReadList(maxReadsToReturn, numericID, interleaved);
		if(list==null){assert(finished);}
		else{assert(list.size()<=maxReadsToReturn);}
		return list;
	}
	
	/**
	 * Parses sequence and quality files to create Read objects.
	 * Handles both single and paired/interleaved read formats.
	 * Validates that sequence and quality headers match.
	 *
	 * @param maxReadsToReturn Maximum number of reads to return
	 * @param numericID Starting numeric ID for reads
	 * @param interleaved Whether reads are paired/interleaved
	 * @return ArrayList of Read objects
	 */
	private ArrayList<Read> toReadList(int maxReadsToReturn, long numericID, boolean interleaved){
		if(finished){return null;}
		if(verbose){System.err.println("FastaQualRIS fetching a list.");}
		
		if(currentHeader==null && numericID==0){
//			assert(numericID==0) : numericID;
			nextBases(btf, builder);
			nextQualities(qtf, builder);
			if(nextHeaderB==null){
				finish();
				return null;
			}
			assert(Tools.equals(nextHeaderB, nextHeaderQ)) : "Quality and Base headers differ for read "+numericID;
			currentHeader=nextHeaderB;
			nextHeaderB=nextHeaderQ=null;
			if(currentHeader==null){
				finish();
				return null;
			}
		}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(1000, maxReadsToReturn));
		
		int added=0;
		
		Read prev=null;
		
		while(added<maxReadsToReturn){
			Read r=makeRead(numericID);
			if(verbose){System.err.println("Made "+r);}
			if(r==null){
				finish();
				if(verbose){System.err.println("makeRead returned null.");}
				break;
			}
			if(interleaved){
				if(prev==null){prev=r;}
				else{
					prev.mate=r;
					r.mate=prev;
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}
		}
		
		assert(list.size()<=maxReadsToReturn);
		if(verbose){System.err.println("FastaQualRIS returning a list.  Size="+list.size());}
		return list;
	}
	
	/**
	 * Reads sequence bases from FASTA file until next header is encountered.
	 * Accumulates sequence lines into a single byte array.
	 *
	 * @param btf ByteFile for reading sequence data
	 * @param bb ByteBuilder for accumulating sequence bytes
	 * @return Byte array containing sequence bases
	 */
	private final byte[] nextBases(ByteFile btf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=btf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			bb.append(line);
			line=btf.nextLine();
		}
		
		
		if(line==null){//end of file
			//do nothing
		}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderB=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/**
	 * Reads quality scores from quality file until next header is encountered.
	 * Handles both numeric (space-separated) and ASCII-encoded quality formats.
	 * Converts ASCII qualities by subtracting FASTQ offset.
	 *
	 * @param qtf ByteFile for reading quality data
	 * @param bb ByteBuilder for accumulating quality bytes
	 * @return Byte array containing quality scores
	 */
	private final byte[] nextQualities(ByteFile qtf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=qtf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			if(NUMERIC_QUAL && line.length>0){
				int x=0;
				for(int i=0; i<line.length; i++){
					byte b=line[i];
					if(b==space){
						assert(i>0);
						bb.append((byte)x);
						x=0;
					}else{
						x=10*x+(b-zero);
					}
				}
				bb.append((byte)x);
			}else{
				for(byte b : line){bb.append((byte)(b-FASTQ.ASCII_OFFSET));}
			}
			line=qtf.nextLine();
		}
//		assert(bb.length()<1) : "'"+Arrays.toString(bb.toBytes())+"'";
		
		if(line==null){//end of file
			//do nothing
		}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderQ=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/**
	 * Creates a Read object from current sequence and quality data.
	 * Validates that sequence and quality arrays have equal length.
	 * Converts sequence bases to uppercase and removes header prefix.
	 *
	 * @param numericID Numeric identifier for the read
	 * @return Read object or null if no more data available
	 */
	private Read makeRead(long numericID){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		if(currentHeader==null){return null;}
		assert(nextHeaderB==null);
		assert(nextHeaderQ==null);
		
		final byte[] bases=nextBases(btf, builder);
		final byte[] quals=nextQualities(qtf, builder);
		final byte[] header=currentHeader;
		
		currentHeader=nextHeaderB;
		nextHeaderB=nextHeaderQ=null;
		
		if(bases==null){
			if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
			return null;
		}
		
		assert(bases.length==quals.length) :
			"\nFor sequence "+numericID+", name "+new String(header)+":\n" +
					"The bases and quality scores are different lengths, "+bases.length+" and "+quals.length;
		
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Tools.toUpperCase(bases[i]);
		}
//		for(int i=0; i<quals.length; i++){
//			quals[i]=(byte)(quals[i]-FASTQ.ASCII_OFFSET);
//		}
		
		assert(bases[0]!=carrot) : new String(bases)+"\n"+numericID+"\n"+header[0];
		String hd=new String(header, 1, header.length-1);
		Read r=new Read(bases, quals, hd, numericID);
		return r;
	}
	
	@Override
	public synchronized boolean close(){
		if(closed){return errorState;}
		if(verbose){System.err.println("FastaQualRIS closing.");}
//		if(verbose){new Exception().printStackTrace(System.err);}
		builder=null;
		finish();
		boolean a=btf.close();
		boolean b=qtf.close();
		closed=true;
		return a|b;
	}

	@Override
	public synchronized void restart() {
		if(verbose){System.err.println("FastaQualRIS restarting.");}
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		
		finished=false;
		closed=false;

		buffer=null;
		nextHeaderB=null;
		nextHeaderQ=null;
		currentHeader=null;
		builder=null;

		btf.reset();
		qtf.reset();
	}

	@Override
	public boolean paired() {
		return interleaved;
	}
	
	/** Marks the input stream as finished reading.
	 * Sets internal finished flag to prevent further read attempts. */
	private synchronized void finish(){
		if(verbose){System.err.println("FastaQualRIS setting finished "+finished+" -> "+true);}
		finished=true;
	}
	
	@Override
	public String fname(){return btf.name()+","+qtf.name();}

	/** Internal buffer for storing parsed reads */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;

	/** ByteFile for reading sequence data */
	private final ByteFile btf;
	/** ByteFile for reading quality score data */
	private final ByteFile qtf;
	/** Whether reads are paired and interleaved */
	private final boolean interleaved;

	/** Buffer length for batch read operations */
	private final int BUF_LEN=Shared.bufferLen();;

	/** Total number of reads generated from input files */
	public long generated=0;
	/** Total number of reads consumed by client code */
	public long consumed=0;
	/** Numeric ID to assign to the next read */
	private long nextReadID=0;
	
	/** Whether quality scores are numeric (space-separated) rather than ASCII */
	public static boolean NUMERIC_QUAL=true;
	
	/** Enable verbose debug output */
	public static boolean verbose=false;

	/** Next header from sequence file */
	private byte[] nextHeaderB=null;
	/** Next header from quality file */
	private byte[] nextHeaderQ=null;
	
	/** Current header being processed */
	private byte[] currentHeader=null;
	
	/** ByteBuilder for accumulating sequence and quality data */
	private ByteBuilder builder=null;
	
	/** Whether input stream has finished reading all data */
	private boolean finished=false, closed=false;
	/** FASTA header prefix character '>' */
	private final byte carrot='>', space=' ', zero='0';
	
}
