package fileIO;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;

import shared.Shared;
import structures.ListNum;


public abstract class ByteFile {
	
//	public static final ByteFile makeByteFile(String fname){
//		return makeByteFile(fname, false, true);
//	}
	
	public static final ByteFile makeByteFile1(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return new ByteFile1(ff);
	}
	
	public static final ByteFile makeByteFile(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return makeByteFile(ff);
	}
	
	public static final ByteFile makeByteFile(FileFormat ff){
		return makeByteFile(ff, 0);
	}
	
	public static final ByteFile makeByteFile(FileFormat ff, int type){
		type=pickType(type);
		if(type==4){return new ByteFile4(ff);}
		if(type==3){return new ByteFile3(ff);}
		if(type==2){return new ByteFile2(ff);}
		return new ByteFile1(ff);
	}
	
	protected ByteFile(FileFormat ff_){
		ff=ff_;
		assert(ff.read()) : ff;
	}
	
	public final ArrayList<byte[]> toByteLines(){
		
		byte[] s=null;
		ArrayList<byte[]> list=new ArrayList<byte[]>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list;
	}
	
	public static final ArrayList<byte[]> toLines(FileFormat ff){
		ByteFile bf=makeByteFile(ff);
		ArrayList<byte[]> lines=bf.toByteLines();
		bf.close();
		return lines;
	}
	
	public static final ArrayList<byte[]> toLines(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, false);
		return toLines(ff);
	}
	
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		
		return count;
	}

	public abstract void reset();
	final void superReset(){
		nextID=0;
	}
	
	public synchronized ListNum<byte[]> nextList(){
		byte[] line=nextLine();
		if(line==null){return null;}
		final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
		ArrayList<byte[]> list=new ArrayList<byte[]>(slimit);
		list.add(line);
		int bytes=line.length;
		
		for(int i=1; i<slimit && bytes<blimit; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
			bytes+=line.length;
		}
		ListNum<byte[]> ln=new ListNum<byte[]>(list, nextID);
		nextID++;
		return ln;
	}
	
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists(); //TODO Ugly and unsafe hack for files in jars
	}

	public abstract InputStream is();
	public abstract long lineNum();
	
	/** Returns true if there was an error */
	public abstract boolean close();
	
	public abstract byte[] nextLine();
	
//	public final void pushBack(byte[] line){
//		assert(pushBack==null);
//		pushBack=line;
//	}
	
	public abstract void pushBack(byte[] line);
	
	public abstract boolean isOpen();
	
	public final String name(){return ff.name();}
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	private static final int pickType(int type) {
		if(type==4 && !ALLOW_BF4) {type=0;}
		else if(type==3 && !ALLOW_BF3) {type=0;}
		else if(type==2 && !ALLOW_BF2) {type=0;}
		else if(type==1 && !ALLOW_BF1) {type=0;}
		if(type>0) {return type;}
		assert(type==0) : type;
		
		final int threads=Shared.threads();
		if(FORCE_MODE_BF1) {return 1;}
		if(FORCE_MODE_BF4) {return 4;}
		if(FORCE_MODE_BF3) {return 3;}
		if(FORCE_MODE_BF2) {return 2;}
		
		if(Shared.LOW_MEMORY || threads<12) {return 1;}
		if(ALLOW_BF4) {return 4;}
		if(ALLOW_BF3) {return 3;}
		if(ALLOW_BF2) {return 2;}
		return 1;
	}
	
	public final FileFormat ff;
	
	/** Force usage of ByteFile1 */
	public static boolean FORCE_MODE_BF1=false;
	public static boolean FORCE_MODE_BF2=false;
	public static boolean FORCE_MODE_BF3=false;
	public static boolean FORCE_MODE_BF4=false;

	public static boolean ALLOW_BF1=true;
	public static boolean ALLOW_BF2=true;
	public static boolean ALLOW_BF3=false;//Hung on fasta input
	public static boolean ALLOW_BF4=true;//Hangs if there are limited reads and it is never closed. 
	
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';
	protected static final byte[] plusLine=new byte[] {plus};
	
	public static int TARGET_LIST_SIZE=800;
	public static int TARGET_LIST_BYTES=262144;
	
//	byte[] pushBack=null;
	protected long nextID=0;
	
}
