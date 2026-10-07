package stream;

import java.io.File;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Locale;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.bam.BamOutputStream;
import stream.bam.BgzfOutputStream;

/**
 * Opens output backends that perform compression on the calling thread.
 * Does not consult global MT compressor preferences or spawn subprocesses.
 * FileFormat selects compression; .bgz/.bgzip select native BGZF within GZIP.
 * This internal factory serves ST/ZT Writers, whose formatting policy is separate.
 * @author Shinobu
 * @date October 1, 2026
 */
final class LightweightOutputStream{

	/** Utility class; all operations use per-output state. */
	private LightweightOutputStream(){}

	/** Checks backend support and overwrite/append options without opening the file. */
	static void validate(final FileFormat ff){
		if(ff==null || !ff.write()){throw new IllegalArgumentException("An output descriptor is required: "+ff);}
		final int compression=ff.compression();
		if(ff.bam()){
			if(compression!=FileFormat.RAW){throw new IllegalArgumentException("BAM already supplies BGZF compression: "+ff);}
		}else if(compression!=FileFormat.RAW && compression!=FileFormat.GZIP && compression!=FileFormat.ZIP){
			throw new IllegalArgumentException("No lightweight compression backend for "+ff.compressionString()+": "+ff.name());
		}
		if(compression==FileFormat.ZIP && ff.append()){
			throw new IllegalArgumentException("ZIP outputs cannot be reopened for append: "+ff.name());
		}
		if(ff.file() && new File(ff.name()).exists() && !ff.append() && !ff.overwrite()){
			throw new IllegalArgumentException("Output exists and overwrite is disabled: "+ff.name());
		}
	}

	/** Opens a checked descriptor, without background compression or process resources. */
	static OutputStream open(final FileFormat ff){
		validate(ff);
		final int level=Math.max(-1, Math.min(9, ReadWrite.ZIPLEVEL));
		if(ff.bam() && ff.file()){
			final File parent=new File(ff.name()).getAbsoluteFile().getParentFile();
			if(!parent.isDirectory() && !parent.mkdirs() && !parent.isDirectory()){
				throw new IllegalArgumentException("Cannot create output directory: "+parent);
			}
			try{return new BamOutputStream(ff.name(), Math.min(level, 6), 1, ff.append());}
			catch(IOException e){throw new RuntimeException("Cannot open BAM output: "+ff.name(), e);}
		}
		final OutputStream raw=ReadWrite.getRawOutputStream(ff.name(), ff.append(), ff.compression()==FileFormat.RAW && !ff.bam());
		final OutputStream out=(raw==System.out || raw==System.err ? new BorrowedOutput(raw) : raw);
		try{
			if(ff.bam()){return new BamOutputStream(out, true, Math.min(level, 6), 1);}
			if(ff.compression()==FileFormat.GZIP){
				final String name=ff.name().toLowerCase(Locale.ROOT);
				if(name.endsWith(".bgz") || name.endsWith(".bgzip")){return new BgzfOutputStream(out, level);}
				return new GZIPOutputStream(out, 8192){
					{def.setLevel(level);}
				};
			}else if(ff.compression()==FileFormat.ZIP){
				final ZipOutputStream zip=new ZipOutputStream(out);
				zip.setLevel(level);
				zip.putNextEntry(new ZipEntry(ReadWrite.basename(ff.name())));
				return zip;
			}
			return out;
		}catch(IOException|RuntimeException e){
			try{out.close();}catch(IOException closeError){e.addSuppressed(closeError);}
			throw new RuntimeException("Cannot open lightweight output: "+ff.name(), e);
		}
	}

	/** Lets compressors finalize without closing a process-wide standard stream. */
	private static final class BorrowedOutput extends FilterOutputStream{
		/** Borrows a standard stream. */
		BorrowedOutput(final OutputStream out){super(out);}
		/** Writes an entire byte range without FilterOutputStream's per-byte loop. */
		@Override
		public void write(final byte[] bytes, final int offset, final int length) throws IOException{
			out.write(bytes, offset, length);
		}
		/** Flushes borrowed output; ownership remains with the process. */
		@Override
		public void close() throws IOException{flush();}
	}
}
