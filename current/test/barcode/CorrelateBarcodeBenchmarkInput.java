package test.barcode;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Generates and independently verifies a bounded literal FASTQ for barcode timing.
 * Records have sequential IDs, a fixed four-base Q40 barcode, 150 defined bases,
 * and uniform Q40 body qualities. Generation and verification are not timed workload.
 * @author Brian Bushnell, Shinobu
 */
public final class CorrelateBarcodeBenchmarkInput{

	/** Requires mode=generate|verify|selftest, path=FILE, reads=1..8000000 and length=150. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, CorrelateBarcodeBenchmarkInput.class, false);
		try{
			String mode=null, path=null;
			long reads=-1;
			int length=150;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("path")){path=value;}
				else if(key.equals("reads")){reads=Long.parseLong(value);}
				else if(key.equals("length")){length=Integer.parseInt(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			require(path!=null && reads>0 && reads<=8000000 && length==150,
					"Set path, reads in [1,8000000], and length=150; larger workloads require a new cost plan.");
			boolean assertions=false;
			assert(assertions=true);
			require(assertions, "Benchmark input qualification requires -ea.");
			final long bytes;
			if("generate".equals(mode)){bytes=generate(path, reads, length);}
			else if("verify".equals(mode)){bytes=verify(path, reads, length);}
			else if("selftest".equals(mode)){
				require(reads==2, "Framing self-test requires reads=2.");
				bytes=selftest(path, length);
			}else{throw new IllegalArgumentException("Set mode=generate, mode=verify or mode=selftest.");}
			require(Files.size(Paths.get(path))==bytes, "File size differs from the exact FASTQ byte count.");
			System.out.println("BENCHMARK_INPUT mode="+mode+" reads="+reads+" length="+length+" bytes="+bytes);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Writes a new file with one reusable record builder and propagates writer failure. */
	private static long generate(final String path, final long reads, final int length) throws IOException{
		require(!Files.exists(Paths.get(path)), "Refusing to overwrite benchmark input: "+path);
		final byte[] bases=new byte[length], qualities=new byte[length];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[i&3];}
		Arrays.fill(qualities, (byte)'I');
		final ByteBuilder record=new ByteBuilder(512);
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		long bytes=0;
		try{
			for(long i=0; i<reads; i++){
				record.clear();
				record.append("@ACGT_IIII_").append(i).nl().append(bases).nl().append('+').nl().append(qualities).nl();
				bytes+=record.length();
				writer.print(record);
			}
		}finally{
			if(writer.poisonAndWait()){throw new IOException("Benchmark input writer failed: "+path);}
		}
		return bytes;
	}

	/** Checks exact field content, increasing identifiers, EOF and all four LF terminators per record. */
	private static long verify(final String path, final long reads, final int length) throws IOException{
		final ByteBuilder expectedHeader=new ByteBuilder(32);
		final String alphabet="ACGT";
		long bytes=0;
		try(final RawInput input=new RawInput(path)){
			for(long i=0; i<reads; i++){
				expectedHeader.clear();
				expectedHeader.append("@ACGT_IIII_").append(i);
				for(int j=0; j<expectedHeader.length(); j++){input.expect(expectedHeader.array[j], i);}
				input.expect('\n', i);
				for(int j=0; j<length; j++){input.expect(alphabet.charAt(j&3), i);}
				input.expect('\n', i);
				input.expect('+', i);
				input.expect('\n', i);
				for(int j=0; j<length; j++){input.expect('I', i);}
				input.expect('\n', i);
				bytes+=expectedHeader.length()+2L*length+5;
			}
			require(input.next()==-1, "Unexpected bytes after requested record count.");
		}
		return bytes;
	}

	/** Rejects CRLF plus a removed final LF even though the two changes preserve file size. */
	private static long selftest(final String path, final int length) throws IOException{
		final long bytes=generate(path, 2, length);
		require(verify(path, 2, length)==bytes, "Valid two-record control failed.");
		try(final RawInput input=new RawInput(path)){
			while(input.next()!=-1){}
			require(input.next()==-1, "EOF must remain EOF on repeated reads.");
		}
		final byte[] original=Files.readAllBytes(Paths.get(path)), malformed=new byte[original.length];
		int firstLF=0;
		while(firstLF<original.length && original[firstLF]!='\n'){firstLF++;}
		require(firstLF<original.length && original[original.length-1]=='\n', "Generated control lacks required LF framing.");
		System.arraycopy(original, 0, malformed, 0, firstLF);
		malformed[firstLF]='\r';
		System.arraycopy(original, firstLF, malformed, firstLF+1, original.length-firstLF-1);
		final String broken=path+".bad-framing";
		Files.write(Paths.get(broken), malformed, StandardOpenOption.CREATE_NEW);
		boolean rejected=false;
		try{verify(broken, 2, length);}
		catch(final IOException expected){rejected=true;}
		require(rejected && Files.size(Paths.get(broken))==bytes, "Equal-size malformed framing was not rejected.");
		final String truncated=path+".missing-lf";
		Files.write(Paths.get(truncated), Arrays.copyOf(original, original.length-1), StandardOpenOption.CREATE_NEW);
		rejected=false;
		try{verify(truncated, 2, length);}
		catch(final IOException expected){rejected=true;}
		require(rejected, "A missing final LF was not rejected.");
		return bytes;
	}

	/** Buffered raw-byte cursor: no CR stripping or implicit acceptance of a missing terminator. */
	private static final class RawInput implements AutoCloseable{
		RawInput(final String path) throws IOException{input=Files.newInputStream(Paths.get(path));}
		int next() throws IOException{
			if(limit<0){return -1;}
			if(offset==limit){
				do{limit=input.read(buffer);}while(limit==0);
				offset=0;
				if(limit<0){return -1;}
			}
			return buffer[offset++]&255;
		}
		void expect(final int expected, final long record) throws IOException{
			final int observed=next();
			if(observed!=expected){
				throw new IOException("Unexpected byte at record "+record+": expected="+expected+", observed="+observed);
			}
		}
		@Override
		public void close() throws IOException{input.close();}
		final InputStream input;
		final byte[] buffer=new byte[65536];
		int offset=0, limit=0;
	}

	/** Enforces generator/verification contracts independently of assertions in downstream code. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
	private CorrelateBarcodeBenchmarkInput(){}
}
