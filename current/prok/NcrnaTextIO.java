package prok;

import java.nio.charset.StandardCharsets;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import stream.Streamer;

/** Shared ownership of native text streams used by the ncRNA model tools.
 * Use these owners in try-with-resources so cleanup failures remain suppressed
 * behind an earlier body failure. Record parsing and buffering belong to callers.
 * Each owner is confined to one thread and belongs to one try-with-resources scope.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaTextIO{

	/*--------------------------------------------------------------*/
	/*----------------          Input Owner         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns a ByteFile; callers retain their separate format-specific parsing state. */
	public static final class Input implements AutoCloseable{
		/** Opens a native byte reader; compression follows the input filename. */
		public Input(final String path_){
			path=path_;
			reader=ByteFile.makeByteFile(path, false);
		}

		/** Preserves a caller's explicit native input/compression options. */
		public Input(final FileFormat format){
			path=format.name();
			reader=ByteFile.makeByteFile(format);
		}

		/** Closes the reader and reports its native error state. */
		@Override
		public final void close(){
			if(reader.close()){throw new IllegalStateException("Could not finish reading "+path);}
		}

		/** Borrowed for record reads; the owner alone closes it. */
		public final ByteFile reader;
		/** Input path retained for failure diagnostics. */
		private final String path;
	}

	/** Owns native sequence-reader closure after consumers reach terminal input.
	 * The caller starts it and joins consumers before closing this scope. Native
	 * startup/cancellation limitations remain those of the selected implementation.
	 */
	public static final class SequenceInput implements AutoCloseable{
		/** Takes ownership of the unstarted reader returned by StreamerFactory. */
		public SequenceInput(final Streamer stream_){
			if(stream_==null){throw new IllegalArgumentException("A sequence input owner requires a reader");}
			stream=stream_;
		}

		/** Reports native close/error state; resource scopes suppress it under prior failures. */
		@Override
		public void close(){
			if(ReadWrite.closeStream(stream)){throw new IllegalStateException("Sequence input did not finish cleanly: "+stream.fname());}
		}

		/** Borrowed for start/read operations; this owner alone closes it. */
		public final Streamer stream;
	}

	/** Owns a reader for small key=value configuration files.
	 * Skips blank lines and whitespace-prefixed comments, trims key/value edges,
	 * and splits only at the first equals sign. Key case and empty values are
	 * preserved for the consumer's own policy checks; no flags or globals change.
	 */
	public static final class ConfigInput implements AutoCloseable{
		/** Opens a configuration through the native text owner. */
		public ConfigInput(final String path_){
			path=path_;
			input=new Input(path);
		}

		/** Advances to one framed record; key and value are valid only after true. */
		public boolean next(){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				int start=0;
				while(start<line.length && line[start]<=32){start++;}
				if(start==line.length || line[start]=='#'){continue;}
				parser.set(line, 0);
				final int separator=parser.b();
				if(separator>=line.length){throw new IllegalArgumentException("Expected configuration key=value in "+path);}
				key=parser.parseString(0).trim();
				value=new String(line, separator+1, line.length-separator-1, StandardCharsets.UTF_8).trim();
				return true;
			}
			key=null;
			value=null;
			return false;
		}

		/** Closes the owned reader; try-with-resources preserves earlier parse errors. */
		@Override
		public void close(){input.close();}

		/** Current trimmed key, retaining its supplied case. */
		public String key;
		/** Current trimmed value, retaining further equals signs and allowing empty text. */
		public String value;
		private final String path;
		private final Input input;
		private final LineParser1 parser=new LineParser1('=');
	}

	/*--------------------------------------------------------------*/
	/*----------------          Output Owner        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns one fresh native writer, including cleanup if worker startup fails. */
	public static final class Output implements AutoCloseable{
		/** Opens and starts a fresh output; existing files are not overwritten. */
		public Output(final String path){this(new ByteStreamWriter(path, false, false, false));}

		/** Preserves the caller's native format, compression and output policies. */
		public Output(final FileFormat format){this(new ByteStreamWriter(format));}

		/** Takes ownership of a newly constructed writer before attempting startup.
		 * A constructor that fails before returning a writer remains the native backend's responsibility.
		 */
		public Output(final ByteStreamWriter writer_){
			if(writer_==null){throw new IllegalArgumentException("An output owner requires its native writer");}
			writer=writer_;
			path=writer.fname;
			try{
				if(writer.getState()!=Thread.State.NEW){throw new IllegalArgumentException("Ownership must precede writer startup: "+path);}
				writer.start();
			}catch(RuntimeException | Error failure){
				try{
					close();
				}catch(RuntimeException | Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
				throw failure;
			}
		}

		/** Finishes and joins a started writer, or closes constructor-opened output
		 * after failed startup. ByteStreamWriter.closeBeforeStart requires state NEW.
		 */
		@Override
		public final void close(){
			final boolean error=writer.getState()==Thread.State.NEW ? writer.closeBeforeStart() : writer.poisonAndWait();
			if(error){throw new IllegalStateException("Could not finish writing "+path);}
		}

		/** Borrowed for writes; the owner alone requests shutdown and joins it. */
		public final ByteStreamWriter writer;
		/** Output path retained for failure diagnostics. */
		private final String path;
	}

	/** A thread-confined collection of output owners, closed in reverse acquisition order.
	 * Workers borrow writers only after acquisition and finish before this scope closes.
	 * Closing attempts every owner, even when an earlier close reports an error.
	 */
	public static final class Outputs implements AutoCloseable{
		/** Reserves one ownership slot per declared output, before any file is opened. */
		public Outputs(final int capacity){
			if(capacity<1){throw new IllegalArgumentException("An output collection needs a positive file count");}
			owners=new Output[capacity];
		}

		/** Opens the next output without altering the requested native backend options. */
		public ByteStreamWriter open(final FileFormat format){
			if(size==owners.length){throw new IllegalStateException("All declared output slots are already owned");}
			final Output owner=new Output(format);
			owners[size++]=owner;
			return owner.writer;
		}

		/** Preserves every cleanup failure; an enclosing resource scope preserves its body failure. */
		@Override
		public void close(){
			Throwable failure=null;
			while(size>0){
				final Output owner=owners[--size];
				owners[size]=null;
				assert(owner!=null) : "Every successful open occupies its slot until reverse cleanup; index="+size;
				try{
					owner.close();
				}catch(RuntimeException | Error cleanup){
					if(failure==null){
						failure=cleanup;
					}else if(failure!=cleanup){
						failure.addSuppressed(cleanup);
					}
				}
			}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			if(failure instanceof Error){throw (Error)failure;}
		}

		private final Output[] owners;
		private int size;
	}

	private NcrnaTextIO(){}
}
