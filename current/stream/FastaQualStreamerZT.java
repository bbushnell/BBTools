package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Zero-threaded (Host-driven) Fasta+Qual streamer.
 * Adapted from the legacy FastaQualReadInputStream to implement the Streamer interface.
 * Uses raw byte scanning for maximum performance on space-delimited integers.
 * @author Brian Bushnell
 * @contributor Collei
 * @date November 22, 2025
 */
public class FastaQualStreamerZT implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public FastaQualStreamerZT(FileFormat ffFa, String qf, int pairnum_, long maxReads_){
		fname=ffFa.name();
		qfname=qf;
		FileFormat ffQual=FileFormat.testInput(qfname, FileFormat.QUAL, null, ffFa.allowSubprocess(), false);
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		btf=ByteFile.makeByteFile(ffFa, 2);
		qtf=ByteFile.makeByteFile(ffQual, 2);
		
		// Legacy logic used a static flag; we make it instance-level here
		numericQual=true; 
		
		if(verbose){outstream.println("Made FastaQualStreamerZT for "+fname);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start() {
		// No threads to start
	}

	@Override
	public synchronized ListNum<Read> nextList() {
		if(finished) {return null;}
		
		// Generate the list on the calling thread (Host-Driven)
		ArrayList<Read> list=toReadList(TARGET_LIST_SIZE);
		
		if(list==null || list.isEmpty()){
			finished=true;
			close();
			return null;
		}
		
		return new ListNum<Read>(list, listID++);
	}
	
	@Override
	public synchronized void close() {
		if(closed){return;}
		btf.close();
		qtf.close();
		closed=true;
	}

	@Override
	public synchronized boolean hasMore() {
		return !finished;
	}

	@Override
	public synchronized boolean errorState() {
		return errorState; // Helper field, though strictly not used much in ZT
	}

	@Override
	public boolean paired() {
		return false;
	}

	@Override
	public int pairnum() {
		return pairnum;
	}

	@Override
	public synchronized long readsProcessed() {
		return generated;
	}

	@Override
	public synchronized long basesProcessed() {
		return consumedBases; // Approximate tracking
	}

	@Override
	public synchronized void setSampleRate(float rate, long seed) {
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	@Override
	public ListNum<SamLine> nextLines() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String fname() {return "("+fname+","+qfname+")";}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses sequence and quality files to create Read objects.
	 * @param maxReadsToReturn Maximum number of reads to return
	 * @return ArrayList of Read objects
	 */
	private synchronized ArrayList<Read> toReadList(int maxReadsToReturn){
		if(finished){return null;}
		if(builder==null){builder=new ByteBuilder(2000);}
		
		if(currentHeader==null && generated==0){
			nextBases(btf, builder);
			nextQualities(qtf, builder);
			if(nextHeaderB==null){
				finish();
				return null;
			}
			// Legacy validation
			if(!Tools.equals(nextHeaderB, nextHeaderQ)) {
				// Relaxed check: just check IDs? Legacy code was strict.
				// We will assume strict for now to match legacy behavior.
				String hb=new String(nextHeaderB);
				String hq=new String(nextHeaderQ);
				if(!hb.equals(hq)){
					// Attempt to recover if it's just a description mismatch?
					// For now, crash to be safe, or implement loose check.
					throw new RuntimeException("Quality and Base headers differ:\n"+hb+"\n"+hq);
				}
			}
			currentHeader=nextHeaderB;
			nextHeaderB=nextHeaderQ=null;
			if(currentHeader==null){
				finish();
				return null;
			}
		}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(1000, maxReadsToReturn));
		int added=0;
		
		while(added<maxReadsToReturn && generated<maxReads){
			Read r=makeRead(generated);
			if(r==null){
				finish();
				break;
			}
			
			if(samplerate>=1f || randy.nextFloat()<samplerate){
				r.setPairnum(pairnum);
				list.add(r);
				added++;
			}
			generated++;
		}
		
		return list;
	}
	
	/**
	 * Creates a Read object from current sequence and quality data.
	 */
	private Read makeRead(long numericID){
		if(currentHeader==null){return null;}
		
		final byte[] bases=nextBases(btf, builder);
		final byte[] quals=nextQualities(qtf, builder);
		final byte[] header=currentHeader;
		
		currentHeader=nextHeaderB;
		nextHeaderB=nextHeaderQ=null;
		
		if(bases==null){return null;}
		
		if(bases.length != quals.length){
			throw new RuntimeException("\nFor sequence "+numericID+", name "+new String(header)+":\n" +
					"The bases and quality scores are different lengths, "+bases.length+" and "+quals.length);
		}
		
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Tools.toUpperCase(bases[i]);
		}
		
		String hd=new String(header, 1, header.length-1); //Strip '>'
		Read r=new Read(bases, quals, hd, numericID);
		consumedBases+=r.length();
		return r;
	}
	
	/**
	 * Reads sequence bases from FASTA file until next header is encountered.
	 * Matches FastaQualReadInputStream logic.
	 */
	private final byte[] nextBases(ByteFile btf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=btf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			bb.append(line);
			line=btf.nextLine();
		}
		
		if(line!=null){
			nextHeaderB=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/**
	 * Reads quality scores from quality file until next header is encountered.
	 * Includes the optimized manual integer parsing loop.
	 */
	private final byte[] nextQualities(ByteFile qtf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=qtf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			if(numericQual && line.length>0){
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
				// ASCII encoded support
				for(byte b : line){bb.append((byte)(b-FASTQ.ASCII_OFFSET));}
			}
			line=qtf.nextLine();
		}
		
		if(line!=null){
			nextHeaderQ=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	private void finish(){
		finished=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final String fname;
	public final String qfname;
	final int pairnum;
	final long maxReads;
	
	private final ByteFile btf;
	private final ByteFile qtf;
	
	private long generated=0;
	private long consumedBases=0;
	private long listID=0;
	
	private ByteBuilder builder;
	private byte[] currentHeader=null;
	private byte[] nextHeaderB=null;
	private byte[] nextHeaderQ=null;
	
	private boolean finished=false;
	private boolean closed=false;
	public boolean errorState=false;
	
	private float samplerate=1f;
	private shared.Random randy=null;
	
	public boolean numericQual=true; //Defaults to true, supports legacy ASCII check if adapted
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	private final byte carrot='>', space=' ', zero='0';
	
	protected PrintStream outstream=System.err;
	public static final boolean verbose=false;

}