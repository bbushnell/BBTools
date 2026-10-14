package fileIO;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Finite synchronous writeString checks with literal bytes and controlled standard streams.
 * Run one case per JVM with -ea; this fixture temporarily replaces System.out or
 * System.err and must not share a JVM with concurrent users of those global streams.
 * Both streams are restored before reporting a standard-stream case.
 *
 * Required arguments are mode=baseline|fixed, name=stdout|stdout.txt|stderr|stderr.txt|file,
 * size=small|large|unicode|trailing and sink=normal|error|append|overwrite.
 * mode selects the expected ownership behavior; it does not configure ReadWrite.
 * Raw stream cases use normal or error; file cases use append or overwrite and
 * create or replace output.txt in the current directory with a literal starting prefix.
 * Run small/large cases with charset=US-ASCII; they contain 4/65536 A bytes.
 * File cases require those sizes.
 * Unicode/trailing stream cases use charset=US-ASCII|UTF-8|UTF-16LE (default US-ASCII)
 * and target the JDK18+ PrintWriter constructor charset contract. Their literal byte
 * oracles cover non-ASCII/supplementary characters and an incomplete terminal surrogate.
 *
 * A successful invocation prints one TSV row: name, size, sink, whether the expected
 * RuntimeException occurred, backend close count (NA for files), and retained bytes.
 * Fixed normal-stream byte counts include the subsequent caller's marker. Expected
 * errors are caught and verified; an unexpected observation fails the invocation.
 * File controls check contents, not operating-system descriptor closure.
 * Compression, asynchronous writing and exceptional cleanup are outside this fixture.
 * @author Shinobu
 */
public final class WriteStringOwnershipTest{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one case and checks literal output, reported errors and backend closure. */
	public static void main(String[] args) throws Exception{
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run ownership checks with -ea.");}
		String mode=null, name=null, size=null, sink=null, charset="US-ASCII";
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equals("mode")){
				mode=value;
			}else if(key.equals("name")){
				name=value;
			}else if(key.equals("size")){
				size=value;
			}else if(key.equals("sink")){
				sink=value;
			}else if(key.equals("charset")){
				charset=value;
			}else{
				throw new IllegalArgumentException("Unknown argument: "+arg);
			}
		}
		if(!"baseline".equals(mode) && !"fixed".equals(mode)){
			throw new IllegalArgumentException("mode=baseline|fixed required");
		}
		if(!"US-ASCII".equals(charset) && !"UTF-8".equals(charset) && !"UTF-16LE".equals(charset)){
			throw new IllegalArgumentException("Unsupported fixture charset: "+charset);
		}
		final byte[] payload;
		final String text;
		if("small".equals(size) || "large".equals(size)){
			payload=new byte[size.equals("small") ? 4 : 65536];
			Arrays.fill(payload, (byte)'A');
			text=new String(payload, StandardCharsets.US_ASCII);
		}else if("unicode".equals(size)){
			text="A\u00e9\u20ac\ud83d\ude00";
			payload=charset.equals("US-ASCII") ? new byte[] {65, 63, 63, 63} :
				charset.equals("UTF-8") ? new byte[] {65, (byte)0xc3, (byte)0xa9, (byte)0xe2, (byte)0x82,
					(byte)0xac, (byte)0xf0, (byte)0x9f, (byte)0x98, (byte)0x80} :
				new byte[] {65, 0, (byte)0xe9, 0, (byte)0xac, 0x20, 0x3d, (byte)0xd8, 0, (byte)0xde};
		}else if("trailing".equals(size)){
			text="A\ud800";
			payload=charset.equals("UTF-16LE") ? new byte[] {65, 0, (byte)0xfd, (byte)0xff} : new byte[] {65, 63};
		}else{throw new IllegalArgumentException("Unsupported fixture size: "+size);}
		final byte[] after=charset.equals("UTF-16LE") ?
			new byte[] {124, 0, 97, 0, 102, 0, 116, 0, 101, 0, 114, 0} : AFTER;
		if("file".equals(name)){
			if(!"small".equals(size) && !"large".equals(size)){
				throw new IllegalArgumentException("File controls require small or large ASCII payloads");
			}
			if(!"append".equals(sink) && !"overwrite".equals(sink)){
				throw new IllegalArgumentException("File sink must be append or overwrite");
			}
			final Path path=Paths.get("output.txt");
			Files.write(path, PREFIX);
			ReadWrite.writeString(text, path.toString(), sink.equals("append"));
			final byte[] expected=sink.equals("append") ? join(PREFIX, payload) : payload;
			assert(Arrays.equals(Files.readAllBytes(path), expected)) :
				"Owned-file control must preserve literal overwrite/append bytes; sink="+sink;
			System.out.println(name+"\t"+size+"\t"+sink+"\tfalse\tNA\t"+expected.length);
			return;
		}
		if(!"stdout".equals(name) && !"stdout.txt".equals(name) &&
			!"stderr".equals(name) && !"stderr.txt".equals(name)){
			throw new IllegalArgumentException("Unrecognized raw standard-stream name: "+name);
		}
		if(!"normal".equals(sink) && !"error".equals(sink)){
			throw new IllegalArgumentException("Standard-stream sink must be normal or error");
		}
		final boolean fail=sink.equals("error"), fixed=mode.equals("fixed");
		final ControlledOutput backend=new ControlledOutput(fail);
		final PrintStream controlled=new PrintStream(backend, false, charset);
		final PrintStream savedOut=System.out, savedErr=System.err;
		RuntimeException failure=null;
		try{
			if(name.startsWith("stdout")){System.setOut(controlled);}else{System.setErr(controlled);}
			try{
				ReadWrite.writeString(text, name, true, false);
			}catch(RuntimeException e){failure=e;}
			assert((failure!=null)==fail) :
				"FIO036 requires retained output errors to throw; name="+name+" sink="+sink+" failure="+failure;
			assert(failure==null || failure.getMessage().equals("Error writing to "+name)) :
				"Expected writeString retained-error report, not an unrelated failure: "+failure;
			assert(backend.closes==(fixed ? 0 : 1)) :
				"Borrowed raw stream must remain open after repair; mode="+mode+" closes="+backend.closes;
			assert(Arrays.equals(backend.bytes.toByteArray(), fail ? new byte[0] : payload)) :
				"writeString must flush the literal payload; name="+name+" size="+size;
			if(!fail){
				controlled.print("|after");
				controlled.flush();
				assert(Arrays.equals(backend.bytes.toByteArray(), fixed ? join(payload, after) : payload)) :
					"A later caller must still be able to use the borrowed stream; mode="+mode;
			}
		}finally{
			System.setOut(savedOut);
			System.setErr(savedErr);
		}
		System.out.println(name+"\t"+size+"\t"+sink+"\t"+(failure!=null)+"\t"+backend.closes+"\t"+backend.bytes.size());
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Concatenates independent literal expectations; never reads implementation output. */
	private static byte[] join(byte[] first, byte[] second){
		final byte[] result=Arrays.copyOf(first, first.length+second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Controlled Sink       ----------------*/
	/*--------------------------------------------------------------*/

	/** A finite sink that retains bytes, counts closure and can reject every write. */
	private static final class ControlledOutput extends OutputStream{
		/** Selects whether every write throws before retaining bytes. */
		ControlledOutput(boolean fail_){fail=fail_;}
		/** Retains one byte, or throws the controlled write error. */
		@Override
		public void write(int value) throws IOException{
			if(fail){throw new IOException("Controlled write failure");}
			bytes.write(value);
		}
		/** Retains a block, or throws before retaining any part of it. */
		@Override
		public void write(byte[] data, int offset, int length) throws IOException{
			if(fail){throw new IOException("Controlled write failure");}
			bytes.write(data, offset, length);
		}
		/** Counts closure without discarding the bytes needed by the oracle. */
		@Override
		public void close(){closes++;}
		/** Whether writes fail before retaining any bytes. */
		final boolean fail;
		/** Bytes actually received, including any successful later caller marker. */
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		/** Number of delegated backend close calls. */
		int closes;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Existing file bytes that append must retain and overwrite must replace. */
	private static final byte[] PREFIX="prefix\n".getBytes(StandardCharsets.US_ASCII);
	/** Literal subsequent caller marker in ASCII and UTF-8. */
	private static final byte[] AFTER="|after".getBytes(StandardCharsets.US_ASCII);
}
