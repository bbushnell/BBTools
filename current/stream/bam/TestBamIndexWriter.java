package stream.bam;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;

/**
 * Generates an index with BamIndexWriter and applies selected structural checks.
 * The checks require indexed content, a pseudo-bin and a trailing no-coordinate count;
 * they are not a general BAI validator. Known parsing limitations are annotated below.
 *
 * Usage: {@code java -ea stream.bam.TestBamIndexWriter [bamFile] [indexFile]}.
 * With no BAM argument the legacy relative path {@code ../Chloe/mapped.bam} is used.
 * An omitted index path is temporary; cleanup is attempted after generation and validation.
 */
public final class TestBamIndexWriter{

	/**
	 * Generates an index, checks it and prints success only after validation returns.
	 * An explicitly supplied index path is retained; later arguments are ignored.
	 * @param args Optional BAM path followed by an optional index path
	 * @throws Exception If temporary-path preparation, generation or validation fails
	 */
	public static void main(String[] args) throws Exception{
		String bamPath=args.length>0 ? args[0] : "../Chloe/mapped.bam";
		String indexPath=args.length>1 ? args[1] : createTemporaryIndexPath();

		try{
			BamIndexWriter.writeIndex(bamPath, indexPath);
			validateIndex(indexPath);
			System.out.println("BamIndexWriter test passed for "+bamPath);
		}finally{
			if(args.length<=1){deleteQuietly(indexPath);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates and deletes a temporary file to obtain an available index pathname.
	 * The returned path is not reserved after deletion.
	 * @return Absolute path of the deleted temporary file
	 * @throws IOException If creating or deleting the temporary file fails
	 */
	private static String createTemporaryIndexPath() throws IOException{
		File temp=File.createTempFile("bam-index-", ".bai");
		String path=temp.getAbsolutePath();
		if(!temp.delete()){throw new IOException("Failed to prepare temporary index file: "+path);}
		return path;
	}

	/**
	 * Deletes the supplied path if present, ignoring an IOException from cleanup.
	 * @param path File to remove
	 */
	private static void deleteQuietly(String path){
		try{
			Files.deleteIfExists(new File(path).toPath());
		}catch(IOException ignored){
			// Best effort clean up
		}
	}

	/**
	 * Reads the complete generated index and checks its magic, counts and offset values.
	 * Requires at least one chunk and pseudo-bin, an eight-byte no-coordinate count,
	 * and no trailing bytes. This routine retains the pseudo-bin parsing issue below.
	 * @param indexPath Index file to read completely into memory
	 * @throws IOException If reading fails or an explicit structural check fails
	 * @throws java.nio.BufferUnderflowException If a requested field exceeds the buffer
	 */
	private static void validateIndex(String indexPath) throws IOException{
		byte[] data=Files.readAllBytes(new File(indexPath).toPath());
		if(data.length<16){throw new IOException("BAI file too short: "+indexPath);}

		ByteBuffer bb=ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);

		int magic=bb.getInt();
		int expectedMagic=('B')|('A'<<8)|('I'<<16)|(1<<24);
		if(magic!=expectedMagic){throw new IOException("Invalid BAI magic number");}

		int nRef=bb.getInt();
		if(nRef<0){throw new IOException("Negative reference count in index");}

		boolean sawIndexedContent=false;
		boolean sawPseudoBin=false;

		for(int ref=0; ref<nRef; ref++){
			long nBin=Integer.toUnsignedLong(bb.getInt());
			for(long b=0; b<nBin; b++){
				int bin=bb.getInt();
				long nChunk=Integer.toUnsignedLong(bb.getInt());
				if(nChunk>0){sawIndexedContent=true;}
				//TODO: Probable bug [stream/bam/TestBamIndexWriter#113] - BamIndexWriter
				//writes pseudo-bin n_chunk=2 followed by four longs: start, end, mapped,
				//unmapped. This loop consumes all four as offset pairs (also comparing the
				//counts as offsets), then the pseudo-bin branch reads two additional longs.
				//The consumer therefore disagrees with this writer's metadata layout.
				for(long c=0; c<nChunk; c++){
					long beg=bb.getLong();
					long end=bb.getLong();
					if(beg<0 || end<beg){throw new IOException("Invalid chunk offsets for bin "+bin);}
				}
				if(bin==37450){
					sawPseudoBin=true;
					long mapped=bb.getLong();
					long unmapped=bb.getLong();
					if(mapped<0 || unmapped<0){throw new IOException("Negative read totals in pseudo-bin");}
				}
			}

			long nIntv=Integer.toUnsignedLong(bb.getInt());
			for(long i=0; i<nIntv; i++){
				long offset=bb.getLong();
				if(offset<0){throw new IOException("Negative linear index offset encountered");}
			}
		}

		//TODO: Possible bug [stream/bam/TestBamIndexWriter#112] - historical review
		//noted an optional n_no_coor field. This harness always reads eight bytes;
		//absent or partial data raises BufferUnderflowException. BamIndexWriter writes
		//the count unconditionally, so support for other producers is separate scope.
		long nNoCoord=bb.getLong();
		if(nNoCoord<0){throw new IOException("Negative n_no_coor value in index");}

		if(bb.hasRemaining()){throw new IOException("Unexpected trailing bytes in index");}

		if(!sawIndexedContent){throw new IOException("Generated index contains no chunks");}
		if(!sawPseudoBin){throw new IOException("Generated index missing pseudo-bin 37450");}
	}
}
