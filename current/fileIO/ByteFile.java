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
		if(type==1){return new ByteFile1(ff);}
		if(type==2){return new ByteFile2(ff);}
		if(!Shared.LOW_MEMORY && (FORCE_MODE_BF2 || (!FORCE_MODE_BF1 && Shared.threads()>4/* && (ReadWrite.isCompressed(fname) || ReadWrite.isSam(fname))*/))){
//			if(allowSubprocess && ((ReadWrite.USE_UNPIGZ || ReadWrite.USE_GUNZIP) && (fname.endsWith(".gz") || fname.endsWith(".gzip")))){}
			return new ByteFile2(ff);
		}
//		if(FORCE_MODE_BF3){return new QuickFile(ff);}
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
	
	public synchronized final ListNum<byte[]> nextList(){
		byte[] line=nextLine();
		if(line==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>(200);
		list.add(line);
		for(int i=1; i<200; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
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
	
	public final FileFormat ff;
	
	/** Force usage of ByteFile1 */
	public static boolean FORCE_MODE_BF1=false;//!(Data.GENEPOOL || Data.DENOVO || Data.CORI || Shared.WINDOWS);
	
	/** Force usage of ByteFile2 */
	public static boolean FORCE_MODE_BF2=false;
	
	/** Unused */
	@Deprecated
	public static boolean FORCE_MODE_BF3=false;
	
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';
	
//	byte[] pushBack=null;
	private long nextID=0;
	
}
