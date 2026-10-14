package fileIO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.zip.ZipOutputStream;

import dna.Data;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * Copies bytes through ReadWrite's input and output backends, including recompression.
 * Matching filename extensions do not bypass the stream factories; codec selection still
 * depends on suffixes and global ReadWrite settings, including RAWMODE for input.
 * File-identity preflight rejects existing aliases, except stream/sink destinations. Copies are not
 * atomic: an error can leave partial output, and concurrent path changes are not prevented.
 * The synchronized copy methods serialize calls on this class; backend settings are global.
 *
 * @author Brian Bushnell, Shinobu
 * @date Jan 23, 2013
 */
public class CopyFile{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses filenames and common/compression options, then copies and reports timing.
	 * Supports in/out, overwrite/ow (default true), append/app (default false), or
	 * positional input/output in the first two argument slots. Remaining unrecognized
	 * arguments are ignored. Append also updates ReadStats.append.
	 * Speed uses File.length() on the supplied input name before opening, divided by
	 * elapsed copy time in decimal MB/s; it does not count decoded bytes and can be zero
	 * for stream names or paths resolved by the backend elsewhere. Diagnostics use stderr.
	 * @param args Command-line arguments, including optional config/help preprocessing
	 */
	public static void main(String[] args){
		{//Preparse help and config files.
			PreParser pp=new PreParser(args, new Object(){}.getClass().getEnclosingClass(), false);
			args=pp.args;
		}
		String in=null, out=null;
		boolean overwrite=true;
		boolean append=false;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(Parser.parseCommonStatic(arg, a, b)){
				//Handled by the shared parser.
			}else if(Parser.parseZip(arg, a, b)){
				//Handled by the compression parser.
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(in==null && i==0 && !args[i].contains("=")){
				in=args[i];
			}else if(out==null && i==1 && !args[i].contains("=")){
				out=args[i];
			}
		}
		assert(in!=null && out!=null);
		long bytes=new File(in).length();
		Timer t=new Timer();
		copyFile(in, out, false, overwrite, append);
		t.stop();
		double mbps1=bytes*1000d/t.elapsed;
		System.err.println("Time:  \t"+t);
		System.err.println(Tools.format("Speed: \t%.2f MB/s", mbps1));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies without appending, using a 16 KiB buffer and the configured ReadWrite backends.
	 * Completion is attempted for both acquired streams, including their registered codec
	 * processes. Direct System.in/out/err objects are borrowed and remain open; output
	 * standard streams are flushed. Wrappers follow the backend's normal close behavior.
	 * A copy failure remains primary, with distinct cleanup failures suppressed on it.
	 * Reported completion errors fail the copy even if all bytes were transferred.
	 * @param source Nonnull input path or supported input stream name
	 * @param dest Nonnull output path or supported output stream name
	 * @param createPathIfNeeded Whether to attempt parent-directory creation before acquisition;
	 * the ReadWrite output factory also attempts parent creation
	 * @param overwrite Whether an existing normalized destination is allowed
	 * @throws IllegalArgumentException If existing input and output paths identify the same file
	 * @throws AssertionError With assertions enabled, if overwrite is false and the destination exists
	 * @throws RuntimeException If path comparison, copying or completion reports an I/O failure
	 */
	public static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded, boolean overwrite){
		copyFile(source, dest, createPathIfNeeded, overwrite, false);
	}

	/** Copies with the public helper's completion policy and an explicit append setting.
	 * Append permits existing output without overwrite; individual codec backends may reject
	 * append. The alias preflight applies to both append and replacement before acquisition.
	 * @param source Nonnull input path or stream name
	 * @param dest Nonnull output path or stream name
	 * @param createPathIfNeeded Whether to attempt parent creation before the output factory
	 * @param overwrite Whether to replace existing output when not appending
	 * @param append Whether to append through the selected output backend
	 */
	private static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded,
			boolean overwrite, boolean append){
		checkCopyTarget(source, dest);
		//Match the raw output backend's filename normalization before checking overwrite.
		assert(append || overwrite || !new File(dest.replace('\\', '/').replace('|', '_')).exists()) : "Destination file already exists: "+dest;
		if(createPathIfNeeded){
			File parent=new File(dest).getParentFile();
			if(parent!=null && !parent.exists()){parent.mkdirs();}
		}
		InputStream in=null;
		OutputStream out=null;
		Throwable failure=null;
		try{
			in=ReadWrite.getInputStream(source, false, true, true);
			out=ReadWrite.getOutputStream(dest, append, false, true);
			final byte[] buffer=new byte[16384];
			int len;
			while((len=in.read(buffer))>0){out.write(buffer, 0, len);}
			out.flush();
			if(out.getClass()==ZipOutputStream.class){
				ZipOutputStream zos=(ZipOutputStream)out;
				zos.closeEntry();
				zos.finish();
			}
			//The historical XZ-specific finish branch remains inactive.
		}catch(FileNotFoundException e){
			failure=new RuntimeException(e);
		}catch(IOException e){
			failure=new RuntimeException(e);
		}catch(RuntimeException | Error e){
			failure=e;
		}finally{
			//Attempt both completions, preserving a copy failure ahead of cleanup failures.
			boolean completionError=false;
			try{
				if(in!=null){completionError|=ReadWrite.finishReading(in, source, true);}
			}catch(RuntimeException | Error e){
				if(failure==null){failure=e;}
				else if(failure!=e){failure.addSuppressed(e);}
			}
			try{
				if(out!=null){completionError|=ReadWrite.finishWriting(null, out, dest, true);}
			}catch(RuntimeException | Error e){
				if(failure==null){failure=e;}
				else if(failure!=e){failure.addSuppressed(e);}
			}
			if(completionError){
				final RuntimeException e=new RuntimeException("Error completing copy from "+source+" to "+dest);
				if(failure==null){failure=e;}
				else{failure.addSuppressed(e);}
			}
		}
		if(failure instanceof Error){throw (Error)failure;}
		if(failure!=null){throw (RuntimeException)failure;}
	}
	
	/** Rejects existing file aliases before output can truncate or extend its input.
	 * Mirrors raw filename normalization and the bare-input fallback; standard output
	 * aliases and /dev/null are exempt from this file-identity check. Concurrent path changes
	 * after this check are not prevented.
	 * @param source Nonnull input filename or stream name
	 * @param dest Nonnull output filename or stream name
	 * @throws IllegalArgumentException If both existing paths identify the same file
	 * @throws RuntimeException If Files.isSameFile reports an I/O failure
	 */
	private static void checkCopyTarget(String source, String dest){
		final String outputName=dest.replace('\\', '/');
		if(outputName.equals("stdout") || outputName.startsWith("stdout.") ||
				outputName.equals("stderr") || outputName.startsWith("stderr.") || outputName.equals("/dev/null")){return;}
		final File output=new File(outputName.replace('|', '_'));
		if(!output.exists()){return;}
		final String inputName=source.replace('\\', '/');
		File input=new File(inputName);
		if(!input.exists() && inputName.indexOf('/')<0 &&
				!inputName.equalsIgnoreCase("stdin") && !inputName.regionMatches(true, 0, "stdin.", 0, 6)){
			input=new File(Data.ROOT_CURRENT+"/"+inputName);
		}
		if(!input.exists()){return;}
		try{
			if(Files.isSameFile(input.toPath(), output.toPath())){
				throw new IllegalArgumentException("Input file and output file refer to the same file: "+source+" and "+dest);
			}
		}catch(IOException e){
			throw new RuntimeException("Could not compare input and output before copying: "+source+" and "+dest, e);
		}
	}
}
