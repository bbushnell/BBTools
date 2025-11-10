package stream.bam;

import java.util.ArrayList;

import fileIO.FileFormat;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Reads BAM files and outputs Read objects.
 * Drop-in replacement for SamReadInputStream when reading BAM files.
 *
 * @author Brian Bushnell
 * @date October 2025
 */
public class BamReadInputStream extends ReadInputStream {

	public static void main(String[] args){
		BamReadInputStream bris=new BamReadInputStream(args[0], false, false, false, true);
		bris.start();

		Read r=bris.nextList().get(0);
		System.out.println(r.toText(false));
		System.out.println();
		if(r.samline!=null){
			System.out.println(r.samline.toText());
			System.out.println();
		}

		bris.close();
	}

	public BamReadInputStream(String fname, boolean loadHeader_, boolean ordered_, boolean interleaved_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, ordered_, interleaved_, -1);
	}

	public BamReadInputStream(String fname, boolean loadHeader_, boolean ordered_, boolean interleaved_, boolean allowSubprocess_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, ordered_, interleaved_, maxReads_);
	}

	public BamReadInputStream(FileFormat ff, boolean loadHeader_, boolean ordered_, boolean interleaved_, long maxReads_){
		loadHeader=loadHeader_;
		ordered=ordered_;
		interleaved=interleaved_;
		maxReads=maxReads_;

		stdin=ff.stdio();
		if(!ff.bam()){
			System.err.println("Warning: Did not find expected bam file extension for filename "+ff.name());
		}

		fname=ff.name();
		header=new ArrayList<byte[]>();
		
		bls=StreamerFactory.makeSamOrBamStreamer(ff, -1, loadHeader, ordered_, maxReads, true);
	}
	
	public void start() {
		if(bls!=null){bls.start();}
	}
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(!finished){
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

	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());

		buffer=null;
		next=0;

		// Get next batch of SamLines from Streamer
		ListNum<SamLine> samLines=bls.nextLines();
		if(samLines==null || samLines.size()==0){
			finished=true;
			buffer=new ArrayList<Read>(0);
			return;
		}

		// Convert SamLine objects to Read objects
		buffer=toReadList(samLines.list, nextReadID);
		nextReadID+=buffer.size();
		generated+=buffer.size();
	}

	private final ArrayList<Read> toReadList(ArrayList<SamLine> samLines, long nextReadID2) {
		ArrayList<Read> list=new ArrayList<Read>(samLines.size());
		for(int idx=0; idx<samLines.size(); idx++){
			SamLine sl1=samLines.get(idx);

			Read r1=sl1.toRead(FASTQ.PARSE_CUSTOM);
			r1.samline=sl1;
			r1.numericID=nextReadID2;
			list.add(r1);

			// Handle paired reads if interleaved
			if(interleaved && (sl1.flag&0x1)!=0 && (sl1.flag&0x40)!=0){
				// This is first of pair, next should be second
				if(idx+1<samLines.size()){
					SamLine sl2=samLines.get(idx+1);
					if((sl2.flag&0x1)!=0 && (sl2.flag&0x80)!=0){
						Read r2=sl2.toRead(FASTQ.PARSE_CUSTOM);
						r2.samline=sl2;
						r2.numericID=nextReadID2;
						r1.mate=r2;
						r2.mate=r1;
						idx++; // Skip the mate since we processed it
					}
				}
			}

			nextReadID2++;
		}
		return list;
	}

	@Override
	public boolean close(){
		finished=true;
		return true;
	}

	@Override
	public synchronized void restart() {
		throw new RuntimeException("BamReadInputStream.restart() not supported - BAM streams cannot be reset");
	}

	@Override
	public String fname(){return fname;}

	@Override
	public boolean paired() {return interleaved;}

	private ArrayList<Read> buffer=null;
	private ArrayList<byte[]> header=null;
	private int next=0;

	private final Streamer bls;
	private final boolean interleaved;
	private final boolean loadHeader;
	private final boolean ordered;
	private final String fname;
	private final long maxReads;
	private boolean finished=false;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;

	public final boolean stdin;

}
