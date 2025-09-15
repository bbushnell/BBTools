package fileIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import shared.Shared;
import shared.Timer;
import shared.Tools;


public class TextFile {
	
	
	public static void main(String[] args){
		TextFile tf=new TextFile(args.length>0 ? args[0] : "stdin", false);
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
		
//		long lines=0;
//		long bytes=0;
//		if(args.length>1){
//			first=Integer.parseInt(args[1]);
//			last=first+100;
//		}
//		if(args.length>2){
//			last=Integer.parseInt(args[2]);
//		}
//		
//		for(int i=0; i<first; i++){tf.readLine();}
//		for(int i=first; i<last; i++){
//			String s=tf.readLine();
//			if(s==null){break;}
//
//			lines++;
//			bytes+=s.length();
//			System.out.println(s);
////			System.out.println(Arrays.toString(s.getBytes()));
//		}
//		
//		System.err.println("\n");
//		System.err.println("Lines: "+lines);
//		System.err.println("Bytes: "+bytes);
//		tf.close();
//		tf.reset();
//		tf.close();
//		
////		for(int i=first; i<last; i++){
////			String s=tf.readLine();
////			if(s==null){break;}
////
////			lines++;
////			bytes+=s.length();
////			System.out.println(s);
////		}
	}
	
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

	public TextFile(String name){this(name, false);}
	
	public TextFile(FileFormat ff){
		file=new File(ff.name());
		allowSubprocess=ff.allowSubprocess();
		name=ff.name();
		
		br=open();
	}
	
	public TextFile(String fname, boolean allowSubprocess_){
		fname=fname.replace('\\', '/');
		file=new File(fname);
		allowSubprocess=allowSubprocess_;
		name=fname;
		
		br=open();
	}
	
	public static final String[] toStringLines(FileFormat ff){
		TextFile tf=new TextFile(ff);
		String[] lines=tf.toStringLines();
		tf.close();
		return lines;
	}
	
	public static final String[] toStringLines(String fname){
		TextFile tf=new TextFile(fname);
		String[] lines=tf.toStringLines();
		tf.close();
		return lines;
	}
	
	/** Generate an array of the lines in this TextFile */
	public final String[] toStringLines(){
		
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list.toArray(new String[list.size()]);
		
	}
	
	public final long countLines(){
		
		String s=null;
		long count=0;
		
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		
		reset();
		
		return count;
		
	}
	
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
	
	public final void reset(){
		close();
		br=open();
	}
	
	public boolean exists(){
		return name.equals("stdin") || name.startsWith("stdin.") || name.startsWith("jar:") || file.exists(); //TODO Ugly and unsafe hack for files in jars
	}
	
	public final boolean close(){
		if(!open){return false;}
		open=false;
		assert(br!=null);
		
		errorState|=ReadWrite.finishReading(is, name, allowSubprocess, br, isr);
		
		br=null;
		is=null;
		isr=null;
		lineNum=-1;
		return false;
	}
	
	public String nextLine(){
		return readLine(true);
	}
	
	public final String readLine(){
		return readLine(true);
	}
	
	public final String readLine(boolean skipBlank){
		String currentLine=null;
		
		
		//Note:  Disabling this block seems to speed things up maybe 5%.
//		boolean ready=false;
//		try {
//			ready=br.ready();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		if(!ready){return null;}
		
		if(!open || br==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name);}
			return null;
		}
		try{
			lineNum++;
			currentLine=br.readLine();
//			System.out.println(lineNum+":\t"+currentLine);
		}catch(Exception e){
			System.err.println("Oops! Bad read in file "+name+" at line "+lineNum);
			System.err.println(""+open+", "+(br==null));
			try {
				File f=new File(name);
				System.err.println("path and length: \t"+f.getAbsolutePath()+"\t"+f.length());
			} catch (Exception e1) {
				//e1.printStackTrace();
			}
			throw new RuntimeException(e);
		}
		if(currentLine==null){return null;}
//		System.out.println("Read "+line);
		
//		currentLine=currentLine.trim();
		
		//Note! This may generate a new String for every line and thus be slow.
//		if(currentLine.trim().length()==0){return readLine();} //Skips blank lines
		if(skipBlank && (currentLine.length()==0 ||
				(Character.isWhitespace(currentLine.charAt(0)) &&
						(Character.isWhitespace(currentLine.charAt(currentLine.length()-1)))) &&
						currentLine.trim().length()==0)){
			return readLine(skipBlank); //Skips blank lines
		}
		
		return currentLine;
	}
	
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
	
	public boolean isOpen(){return open;}

	private boolean open=false;
	public boolean errorState=false;
	
	public final String name;
	public File file;
	private final boolean allowSubprocess;
	
	public InputStream is;
	public InputStreamReader isr;
	public BufferedReader br;
	
	public long lineNum=-1;

	public static boolean verbose=false;
	
}
