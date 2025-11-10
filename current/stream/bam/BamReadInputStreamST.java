package stream.bam;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.SamReadInputStream;
import structures.ByteBuilder;

/**
 * Single-threaded BAM reader that reads BAM files directly.
 * Drop-in replacement for SamReadInputStream when reading BAM files.
 * Unlike BamReadInputStream, this version has no background threads.
 *
 * @author Brian Bushnell
 * @date October 2025
 */
public class BamReadInputStreamST extends ReadInputStream {

	public static void main(String[] args){
		BamReadInputStreamST bris=new BamReadInputStreamST(args[0], false, false, true);

		Read r=bris.nextList().get(0);
		System.out.println(r.toText(false));
		System.out.println();
		if(r.samline!=null){
			System.out.println(r.samline.toText());
			System.out.println();
		}

		bris.close();
	}

	public BamReadInputStreamST(String fname, boolean loadHeader_, boolean interleaved_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, interleaved_);
	}

	public BamReadInputStreamST(FileFormat ff, boolean loadHeader_, boolean interleaved_){
		loadHeader=loadHeader_;
		interleaved=interleaved_;

		stdin=ff.stdio();
		if(!ff.bam()){
			System.err.println("Warning: Did not find expected bam file extension for filename "+ff.name());
		}

		fname=ff.name();
		header=new ArrayList<byte[]>();

		try {
            fis=new FileInputStream(fname);
			if(BgzfSettings.USE_MULTITHREADED_BGZF) {
				int threads = Math.max(1, BgzfSettings.READ_THREADS);
				bgzf=new BgzfInputStreamMT(fis, threads);
			} else {
				bgzf=new BgzfInputStream(fis);
			}
            reader=new BamReader(bgzf);

			// Read BAM magic
			byte[] magic=reader.readBytes(4);
			if(!Arrays.equals(magic, new byte[]{'B', 'A', 'M', 1})){
				throw new RuntimeException("Not a BAM file: "+fname);
			}

			// Read header text
			long l_text=reader.readUint32();
			byte[] text=reader.readBytes((int)l_text);

			// Parse header if requested
			if(loadHeader){
				int start=0;
				for(int i=0; i<text.length; i++){
					if(text[i]=='\n'){
						if(i>start){
							byte[] line=Arrays.copyOfRange(text, start, i);
							if(Shared.TRIM_RNAME){line=SamReadInputStream.trimHeaderSQ(line);}
							header.add(line);
						}
						start=i+1;
					}
				}
				// Add last line if not ending with newline
				if(start<text.length){
					byte[] line=Arrays.copyOfRange(text, start, text.length);
					header.add(line);
				}
				SamReadInputStream.setSharedHeader(header);
			}

			// Read reference sequence dictionary
			int n_ref=reader.readInt32();
			refNames=new String[n_ref];
			for(int i=0; i<n_ref; i++){
				long l_name=reader.readUint32();
				refNames[i]=reader.readString((int)l_name-1); // Exclude NUL
				reader.readUint8(); // Skip NUL terminator
				long l_ref=reader.readUint32(); // Reference length (unused here)
			}

			converter=new BamToSamConverter(refNames);

		} catch(IOException e){
			throw new RuntimeException("Error opening BAM file: "+fname, e);
		}
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

		int BUF_LEN=Shared.bufferLen();
		buffer=new ArrayList<Read>(BUF_LEN);
		final ByteBuilder cigar=new ByteBuilder(1024);

		// Read alignment records until buffer full or EOF
		try {
			while(buffer.size()<BUF_LEN){
				long block_size=reader.readUint32();
				byte[] bamRecord=reader.readBytes((int)block_size);
				SamLine sl=converter.toSamLine(bamRecord, cigar);
				Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
				r.samline=sl;
				r.numericID=nextReadID;
				buffer.add(r);

				nextReadID++;
			}
		} catch(EOFException e){
			// Normal end of file
			finished=true;
		} catch(IOException e){
			throw new RuntimeException("Error reading BAM file: "+fname, e);
		}

		generated+=buffer.size();
	}

	@Override
	public boolean close(){
		finished=true;
		try {
			if(bgzf!=null){bgzf.close();}
			if(fis!=null){fis.close();}
		} catch(IOException e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	@Override
	public synchronized void restart() {
		throw new RuntimeException("BamReadInputStreamST.restart() not supported - BAM streams cannot be reset");
	}

	@Override
	public String fname(){return fname;}

	@Override
	public boolean paired() {return interleaved;}

	private ArrayList<Read> buffer=null;
	private ArrayList<byte[]> header=null;
	private int next=0;

    private final FileInputStream fis;
	private final InputStream bgzf;
	private final BamReader reader;
	private final BamToSamConverter converter;
	private final String[] refNames;

	private final boolean interleaved;
	private final boolean loadHeader;
	private final String fname;
	private boolean finished=false;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;

	public final boolean stdin;

}
