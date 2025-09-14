package driver;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.zip.ZipOutputStream;

import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Timer;

public class TestCompressionSpeed {
	
	
	public static void main(String[] args){
		
		TextFile tf=new TextFile(args[0], false);
		String[] lines=tf.toStringLines();
		tf.close();

		Timer t=new Timer();
		
		for(int i=0; i<=9; i++){
			t.start();
			String fname=args[1].replaceFirst("#", ""+i);
			compress(lines, fname, i);
			t.stop();

			System.out.println("Level "+i+" compress:   "+t+"  \tsize:       "+new File(fname).length());
		}
		
		for(int i=0; i<=9; i++){
			t.start();
			String fname=args[1].replaceFirst("#", ""+i);
			String[] lines2=read(fname);
			assert(lines2.length>=lines.length);
			t.stop();
			
			System.out.println("Level "+i+" decompress: "+t);
		}
		
	}
	
	
	public static void compress(String[] text, String fname, int level){
		ReadWrite.ZIPLEVEL=level;
		OutputStream os=ReadWrite.getOutputStream(fname, false, true, true);
		PrintWriter writer=new PrintWriter(os);

		for(String s : text){writer.println(s);}
		for(String s : text){writer.println(s);}
		for(String s : text){writer.println(s);}
		for(String s : text){writer.println(s);}
		
		try {
			writer.flush();
			if(os.getClass()==ZipOutputStream.class){
				ZipOutputStream zos=(ZipOutputStream)os;
				zos.closeEntry();
				zos.finish();
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public static String[] read(String fname){
		TextFile tf=new TextFile(fname, false);
		String[] s=tf.toStringLines();
		tf.close();
		return s;
	}
	
}
