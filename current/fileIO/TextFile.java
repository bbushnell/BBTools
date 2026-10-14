package fileIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Line-oriented reader using the default charset and a 32,768-character buffer.
 * ReadWrite selects the backend from the input name and subprocess setting.
 * Logical-line methods skip the blanks described by readLine(boolean), retain
 * other lines without trimming, and leave input open at EOF. Instances and their
 * public stream fields are mutable and unsynchronized; callers serialize access.
 *
 * Public close/reset and static collectors can close stdin through reader wrappers.
 * The diagnostic main separately preserves borrowed stdin and stdout.
 * @author Brian Bushnell, Shinobu
 */
public class TextFile{
	/**
	 * Reads a logical-line range, or times reading when the second argument is speedtest.
	 * Owns regular input cleanup after successful construction; preserves borrowed stdin
	 * and stdout. Print mode checks recorded stdout errors after successful reading.
	 * Earlier failures stay primary if owned-input cleanup also fails.
	 * @param args Optional filename, start index or speedtest, and exclusive stop index
	 */
	public static void main(String[] args){
		final String fname=(args.length>0 ? args[0] : "stdin").replace('\\', '/');
		final String lower=fname.toLowerCase();
		// Match ReadWrite.getRawInputStream: an existing stdin-named file is owned.
		final boolean borrowedInput=!new File(fname).exists() &&
			(lower.equals("stdin") || lower.startsWith("stdin."));
		final TextFile tf=new TextFile(fname, false);
		Throwable failure=null;
		try{
			int first=0;
			long last=100;
			boolean speedtest=false;
			if(args.length>1){
				if(args[1].equalsIgnoreCase("speedtest")){
					speedtest=true;
					first=0;
					last=Long.MAX_VALUE;
				}else{
					first=Integer.parseInt(args[1]);
					last=first+100;
				}
			}
			if(args.length>2){
				last=Integer.parseInt(args[2]);
			}
			speedtest(tf, first, last, !speedtest);
			// FIO-031: PrintStream records write failures instead of throwing them.
			if(!speedtest && System.out.checkError()){
				throw new RuntimeException("TextFile: error writing output");
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			// FIO-030: close owned input on parse/read exits, preserving the primary failure.
			// Reader wrappers would close borrowed stdin despite ReadWrite's raw-stream guard.
			if(!borrowedInput){
				try{
					if(tf.close()){throw new RuntimeException("Error closing "+tf.name);}
				}catch(RuntimeException | Error t){
					if(failure==null){throw t;}
					if(failure!=t){failure.addSuppressed(t);}
				}
			}
		}
	}

	/**
	 * Skips first logical lines, then reads until last or EOF. Does not close the reader.
	 * The diagnostic count labeled Bytes sums String.length(), not encoded bytes or terminators.
	 * @param tf Reader borrowed from main
	 * @param first Number of logical lines to skip for a nonnegative start
	 * @param last Exclusive logical stop position
	 * @param reprint Print retained lines and totals; otherwise print timing and totals
	 */
	private static void speedtest(TextFile tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				String s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length();
				System.out.println(s);
			}
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				String s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length();
			}
		}
		t.stop();
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}

	/**
	 * Opens name immediately with subprocess decompression disabled.
	 * @param name Input name interpreted by ReadWrite
	 */
	public TextFile(String name){this(name, false);}

	/**
	 * Opens ff.name() using ff.allowSubprocess(); ReadWrite selects the actual backend.
	 * @param ff Nonnull descriptor supplying the name and subprocess permission
	 */
	public TextFile(FileFormat ff){
		file=new File(ff.name());
		allowSubprocess=ff.allowSubprocess();
		name=ff.name();
		br=open();
	}

	/**
	 * Normalizes backslashes and immediately opens a default-charset buffered reader.
	 * @param fname Nonnull input name interpreted by ReadWrite
	 * @param allowSubprocess_ Whether ReadWrite may use a subprocess backend
	 */
	public TextFile(String fname, boolean allowSubprocess_){
		fname=fname.replace('\\', '/');
		file=new File(fname);
		allowSubprocess=allowSubprocess_;
		name=fname;
		br=open();
	}

	/**
	 * Collects remaining logical lines from a newly opened reader, then attempts closure.
	 * Uses the bounded overload with an unlimited logical-line count.
	 * @param ff Descriptor supplying the input name and subprocess permission
	 * @return New array of retained lines
	 */
	public static final String[] toStringLines(FileFormat ff){return toStringLines(ff, -1);}

	/**
	 * Collects logical lines from a newly opened reader with subprocesses disabled.
	 * Uses the bounded overload with an unlimited logical-line count and attempts closure.
	 * @param fname Input name
	 * @return New array of retained lines
	 */
	public static final String[] toStringLines(String fname){return toStringLines(fname, -1);}

	/**
	 * Collects logical lines and attempts closure after successful reader construction.
	 * Cleanup runs on normal and exceptional collection exits, including for stdin.
	 * Reported close errors reject the result; distinct cleanup failures are suppressed
	 * onto an existing primary collection failure.
	 * @param ff Descriptor supplying the input name and subprocess permission
	 * @param max Maximum retained lines; negative means unlimited and zero consumes none
	 * @return New array of retained lines
	 */
	public static final String[] toStringLines(FileFormat ff, int max){
		return collectAndClose(new TextFile(ff), max);
	}

	/**
	 * Collects logical lines with subprocesses disabled, then attempts reader closure.
	 * Cleanup begins after successful construction and also applies to stdin. Reported
	 * close errors reject the result; an existing collection failure remains primary.
	 * @param fname Input name
	 * @param max Maximum retained lines; negative means unlimited and zero consumes none
	 * @return New array of retained lines
	 */
	public static final String[] toStringLines(String fname, int max){
		return collectAndClose(new TextFile(fname), max);
	}

	/**
	 * Owns collection cleanup after construction, preserving a primary collection failure.
	 * @param tf Successfully constructed reader whose closure is attempted
	 * @param max Maximum logical lines, or negative for unlimited
	 * @return Collected lines when collection and closure succeed
	 */
	private static String[] collectAndClose(final TextFile tf, final int max){
		Throwable failure=null;
		try{
			return tf.toStringLines(max);
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			//FIO-029: always attempt closure; cleanup must not replace a collection failure.
			try{
				if(tf.close()){throw new RuntimeException("Error closing "+tf.name);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
	}

	/**
	 * Collects remaining logical lines to EOF without closing or resetting this reader.
	 * @return New array of retained lines
	 */
	public final String[] toStringLines(){
		return toStringLines(-1);
	}

	/**
	 * Collects at most max remaining logical lines without consuming the following line.
	 * A negative limit reads to EOF; zero leaves the logical cursor unchanged. The
	 * underlying readers may prefetch bytes. Does not close or reset this reader.
	 * @param max Maximum retained lines, or negative for unlimited
	 * @return New array of retained lines
	 */
	public final String[] toStringLines(long max){
		if(max<0){max=Long.MAX_VALUE;}
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		//FIO-028: test the collection limit before consuming another logical line.
		while(list.size()<max && (s=nextLine())!=null){
			list.add(s);
		}
		return list.toArray(new String[list.size()]);
	}

	/**
	 * Counts remaining logical lines, then closes and reopens the input name.
	 * Reset occurs only after successful traversal; it is not exceptional cleanup.
	 * Reopening is not a rewind guarantee for stdin or another nonreopenable source.
	 * @return Number of remaining retained lines encountered before EOF
	 */
	public final long countLines(){
		String s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		return count;
	}

	/**
	 * Splits each supplied string on tabs, preserving trailing empty fields.
	 * Allocates new row arrays without changing the input array or its strings.
	 * @param lines Nonnull array of nonnull strings
	 * @param trim Apply String.trim() before splitting each string
	 * @return New mutable array of field arrays
	 */
	public static String[][] doublesplitTab(String[] lines, boolean trim){
		String[][] lines2=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			if(trim){
				lines2[i]=lines[i].trim().split("\t", -1);
			}else{
				lines2[i]=lines[i].split("\t", -1);
			}
		}
		return lines2;
	}

	/**
	 * Splits on runs of Java whitespace using String.split's default empty-field rules.
	 * Trailing empty fields are discarded. Does not change the input array or strings.
	 * @param lines Nonnull array of nonnull strings
	 * @param trim Apply String.trim() before splitting; its definition differs from Java whitespace
	 * @return New mutable array of field arrays
	 */
	public static String[][] doublesplitWhitespace(String[] lines, boolean trim){
		String[][] lines2=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			if(trim){
				lines2[i]=lines[i].trim().split("\\p{javaWhitespace}+");
			}else{
				lines2[i]=lines[i].split("\\p{javaWhitespace}+");
			}
		}
		return lines2;
	}

	/**
	 * Closes and reopens the input name, ignoring the close return value.
	 * Sticky errorState is retained. Ordinary files restart; stdin is not seeked and
	 * may have been closed through the reader wrappers.
	 */
	public final void reset(){
		close();
		br=open();
	}

	/**
	 * Applies a name/File.exists heuristic without opening or validating accessibility.
	 * The case-sensitive stdin, stdin.* and jar: shortcuts return true directly;
	 * these spelling rules differ from ReadWrite's raw-input alias selection.
	 * @return Whether a shortcut matches or the stored File exists
	 */
	public boolean exists(){
		return name.equals("stdin") || name.startsWith("stdin.") || name.startsWith("jar:") || file.exists(); //TODO Ugly and unsafe hack for files in jars
	}

	/**
	 * Attempts reader and input closure once, folding reported closure errors into errorState.
	 * Reader wrappers may close stdin even when the raw-input closer exempts it.
	 * A normally completed open-reader close clears references and sets lineNum=-1.
	 * A thrown cleanup exception can leave references intact after open becomes false.
	 * @return Sticky closure error flag, also returned on repeated closed calls
	 */
	public final boolean close(){
		//FIO-027: retain the close-status contract used by CallPeaks and other boolean-return callers.
		if(!open){return errorState;}
		open=false;
		assert(br!=null);
		errorState|=ReadWrite.finishReading(is, name, allowSubprocess, br, isr);
		br=null;
		is=null;
		isr=null;
		lineNum=-1;
		return errorState;
	}

	/**
	 * Reads the next retained logical line using readLine(true).
	 * @return Untrimmed line, or null at EOF or when closed
	 */
	public String nextLine(){return readLine(true);}

	/**
	 * Reads the next retained logical line using readLine(true).
	 * @return Untrimmed line, or null at EOF or when closed
	 */
	public final String readLine(){return readLine(true);}

	/**
	 * Reads through the buffered reader, optionally skipping its exact blank predicate.
	 * A skipped line is empty, or has Character.isWhitespace endpoints and an empty
	 * String.trim() result. Retained lines are not trimmed; terminators are removed
	 * by BufferedReader. This predicate is not a general Unicode-whitespace test.
	 * Each attempted buffered read increments lineNum, including skipped lines, EOF
	 * and failed attempts. EOF leaves the reader open. Caught read exceptions print
	 * diagnostics and are wrapped in RuntimeException; they do not latch errorState.
	 * @param skipBlank Whether to skip lines matching the predicate above
	 * @return Retained line, or null at EOF or when closed/the buffered reader is null
	 */
	public final String readLine(boolean skipBlank){
		String currentLine=null;
		// No readiness probe: BufferedReader.readLine may block until input is available.
		if(!open || br==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name);}
			return null;
		}
		while(true){
			try{
				lineNum++;
				currentLine=br.readLine();
			}catch(Exception e){
				System.err.println("Oops! Bad read in file "+name+" at line "+lineNum);
				System.err.println(""+open+", "+(br==null));
				try{
					File f=new File(name);
					System.err.println("path and length: \t"+f.getAbsolutePath()+"\t"+f.length());
				}catch(Exception e1){
					// Preserve the original read failure if diagnostic file inspection also fails.
				}
				throw new RuntimeException(e);
			}
			if(currentLine==null){return null;}
			// Test endpoints before trim to avoid that work for ordinary retained lines.
			if(skipBlank && (currentLine.length()==0 ||
					(Character.isWhitespace(currentLine.charAt(0)) &&
							(Character.isWhitespace(currentLine.charAt(currentLine.length()-1)))) &&
							currentLine.trim().length()==0)){
				continue; //Skips blank lines
			}
			return currentLine;
		}
	}

	/**
	 * Opens a ReadWrite backend and wraps it in default-charset character readers.
	 * Marks open before backend construction; constructor failures are not cleaned up here.
	 * Does not reset lineNum or errorState.
	 * @return BufferedReader with capacity 32,768 characters
	 */
	private final BufferedReader open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name);
		}
		open=true;
		is=ReadWrite.getInputStream(name, true, allowSubprocess);
		isr=new InputStreamReader(is);
		BufferedReader b=new BufferedReader(isr, 32768);
		return b;
	}

	/** Returns the local open flag; this is not an EOF or backend-health check. */
	public boolean isOpen(){return open;}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Local open flag; set before backend construction and cleared before closure. */
	private boolean open=false;
	/** Sticky reported closure errors; read exceptions are thrown rather than latched here. */
	public boolean errorState=false;

	/** Input name used to open/reopen and to report diagnostics. */
	public final String name;
	/** File used by exists(); replacing it does not change the input name or current reader. */
	public File file;
	/** Subprocess permission passed to ReadWrite for opening and finishing input. */
	private final boolean allowSubprocess;

	/** Underlying backend stream; cleared after a normally completed close. */
	public InputStream is;
	/** Default-charset decoder wrapping is; cleared after a normally completed close. */
	public InputStreamReader isr;
	/** Buffered character reader with capacity 32,768; cleared after normal close. */
	public BufferedReader br;

	/** Read-attempt index starting at -1; includes skipped lines, EOF and failed attempts. */
	public long lineNum=-1;

	/** Legacy verbosity field; this class does not currently consult it. */
	public static boolean verbose=false;
}
