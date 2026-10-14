package fileIO;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Arrays;

/** Finite injected-input checks for OpenFile's read, close and error contracts.
 * @author Shinobu
 */
public class OpenFileInputTest {

	/** Runs one case; an observed target exception is rethrown after checking its identity. */
	public static void main(String[] args) throws Throwable{
		if(args.length!=2){throw new IllegalArgumentException("Expected baseline|fixed and a case name");}
		final boolean fixed=args[0].equals("fixed");
		if(!fixed && !args[0].equals("baseline")){throw new IllegalArgumentException(args[0]);}
		if(new File("stdin").exists()){throw new IllegalStateException("Fixture requires an absent stdin filename");}
		final String name=args[1];
		final boolean readError=name.equals("readError") || name.equals("bothError");
		final boolean closeError=name.equals("closeError") || name.equals("bothError");
		byte[] payload="abc".getBytes("US-ASCII");
		int chunk=100;
		if(name.equals("empty")){payload=new byte[0];}
		else if(name.equals("shortRead")){payload="abcdef".getBytes("US-ASCII"); chunk=2;}
		else if(name.equals("long")){payload=new byte[150]; Arrays.fill(payload, (byte)'A');}
		else if(!name.equals("ordinary") && !readError && !closeError){
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		final ControlledInput input=new ControlledInput(payload, chunk, readError, closeError);
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final ByteArrayOutputStream captured=new ByteArrayOutputStream();
		final PrintStream captureStream=new PrintStream(captured);
		Throwable failure=null;
		try{
			System.setIn(input);
			System.setErr(captureStream);
			try{OpenFile.main(new String[]{"stdin"});}
			catch(Throwable t){failure=t;}
		}finally{
			captureStream.flush();
			System.setIn(originalIn);
			System.setErr(originalErr);
		}
		originalErr.write(captured.toByteArray());
		final boolean expectFailure=fixed && (readError || closeError);
		assert((failure!=null)==expectFailure) : "OpenFile must propagate input failures after repair; case="+name+", failure="+failure;
		assert(input.reads==1) : "OpenFile's diagnostic performs one bulk read, not a fill loop; reads="+input.reads;
		assert(input.closes==1) : "OpenFile must attempt close once, including after a read failure; closes="+input.closes;
		if(expectFailure){
			assert(!readError || contains(failure, input.readFailure)) : "Propagation must retain the original read error; case="+name;
			assert(!closeError || contains(failure, input.closeFailure)) : "Propagation must retain the original close error, including combined failure; case="+name;
		}
		if(!readError && !closeError){
			final int size=Math.min(payload.length, Math.min(chunk, 100));
			final String expected=size==0 ? "" : "'"+new String(payload, 0, size)+"'"+System.lineSeparator();
			assert(captured.toString().equals(expected)) : "Successful output must preserve one-read/default-charset rendering; case="+name;
		}else if(!readError){
			assert(captured.toString().startsWith("'abc'"+System.lineSeparator())) : "Close failure must not erase already printed bytes";
		}
		System.out.println(name+"\t"+(failure!=null)+"\t"+input.reads+"\t"+input.closes);
		if(failure!=null){throw failure;}
	}

	/** Searches only the finite exception chain produced by this fixture and the target. */
	private static boolean contains(Throwable failure, final Throwable expected){
		if(failure==null){return false;}
		if(failure==expected){return true;}
		for(Throwable suppressed : failure.getSuppressed()){
			if(contains(suppressed, expected)){return true;}
		}
		return contains(failure.getCause(), expected);
	}

	/** Returns one bounded chunk, or the named read/close errors without filesystem faults. */
	private static class ControlledInput extends InputStream {
		ControlledInput(byte[] payload_, int chunk_, boolean readError_, boolean closeError_){
			payload=payload_; chunk=chunk_; readError=readError_; closeError=closeError_;
		}
		@Override
		public int read(){throw new AssertionError("OpenFile should use the bulk InputStream read overload");}
		@Override
		public int read(byte[] dest, int offset, int length) throws IOException{
			reads++;
			if(readError){throw readFailure;}
			if(payload.length==0){return -1;}
			final int size=Math.min(length, Math.min(chunk, payload.length));
			System.arraycopy(payload, 0, dest, offset, size);
			return size;
		}
		@Override
		public void close() throws IOException{
			closes++;
			if(closeError){throw closeFailure;}
		}
		final byte[] payload;
		final int chunk;
		final boolean readError, closeError;
		final IOException readFailure=new IOException("controlled-read-failure");
		final IOException closeFailure=new IOException("controlled-close-failure");
		int reads, closes;
	}
}
