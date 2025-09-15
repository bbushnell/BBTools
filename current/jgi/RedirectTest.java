package jgi;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import dna.Data;
import fileIO.PipeThread;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * @author Brian Bushnell
 * @date Jan 22, 2013
 *
 */
public class RedirectTest {
	
	public static void main(String[] args) throws IOException{
		
		String fin=args[0];
//		String fout=args[1];
		
		System.out.println("fin="+fin);
		
		InputStream in=null;
		final OutputStream os=System.out;
		InputStream es=null;
		Process p=null;

		System.out.println("Samtools="+Data.SAMTOOLS());
		System.out.println("Gzip="+Data.GZIP());
		System.out.println("Pigz="+Data.PIGZ());
		System.out.println("Gunzip="+Data.GUNZIP());
		
		if(Shared.WINDOWS){
			System.out.println("WINDOWS");
			in=ReadWrite.getInputStream(fin, false, false);
		}else{
			System.out.println("LINUX");
			p=Runtime.getRuntime().exec("gunzip -c -d "+fin);
			in=p.getInputStream();
			es=p.getErrorStream();
			assert(es!=null);
			PipeThread et=new PipeThread(es, System.err);
			et.start();
			System.out.println(p);
		}
		
		final byte[] buf=new byte[4096];
		for(int len=in.read(buf); len>0; len=in.read(buf)){
			os.write(buf, 0, len);
		}
		
		in.close();
		if(es!=null){es.close();}
		ReadWrite.close(os);
		
	}
	
	public static void main_0(String[] args) throws IOException{
		
		String fin=args[0];
		String fout=args[1];
		
		InputStream in=ReadWrite.getInputStream(fin, false, false);
		
		OutputStream os=System.out;
		
		byte[] buf=new byte[4096];
		
		for(int len=in.read(buf); len>0; len=in.read(buf)){
			os.write(buf, 0, len);
		}
		
	}
	
}
