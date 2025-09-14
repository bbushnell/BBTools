package fileIO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.ZipOutputStream;

import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;


/**
 * Unlike ReadWrite's version, this one forces compression and decompression even with same extensions.
 * Mainly for benchmarking.
 * @author Brian Bushnell
 * @date Jan 23, 2013
 *
 */
public class CopyFile {
	
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
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
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
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
		copyFile(in, out, false, overwrite);
		t.stop();
		double mbps1=bytes*1000d/t.elapsed;
		System.err.println("Time:  \t"+t);
		System.err.println(Tools.format("Speed: \t%.2f MB/s", mbps1));
	}
	
	
	public static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded, boolean overwrite){

		assert(overwrite || !new File(dest).exists()) : "Destination file already exists: "+dest;
		if(createPathIfNeeded){
			File parent=new File(dest).getParentFile();
			if(parent!=null && !parent.exists()){
				parent.mkdirs();
			}
		}

		try{
			InputStream in=ReadWrite.getInputStream(source, false, true);
			OutputStream out=ReadWrite.getOutputStream(dest, false, false, true);

			final byte[] buffer=new byte[16384];
			int len;

			while((len = in.read(buffer)) > 0){
				out.write(buffer, 0, len);
			}

			in.close();
			out.flush();
			if(out.getClass()==ZipOutputStream.class){
				ZipOutputStream zos=(ZipOutputStream)out;
				zos.closeEntry();
				zos.finish();
			}
			//			else if(PROCESS_XZ && out.getClass()==org.tukaani.xz.XZOutputStream.class){
			//				org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)out;
			//				zos.finish();
			//			}
			out.close();

		}catch(FileNotFoundException e){
			throw new RuntimeException(e);
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
}
