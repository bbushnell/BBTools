package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Feb 13, 2013
 *
 */
public class FastaShredInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		int a=20, b=Integer.MAX_VALUE;
		if(args.length>1){a=Integer.parseInt(args[1]);}
		if(args.length>2){b=Integer.parseInt(args[2]);}
		if(args.length>3){MIN_READ_LEN=Integer.parseInt(args[3]);}
		if(args.length>4){TARGET_READ_LEN=Integer.parseInt(args[4]);}
		
		Timer t=new Timer();
		
		FastaShredInputStream fris=new FastaShredInputStream(args[0], false, false, Shared.bufferData());
		Read r=fris.next();
		int i=0;
		
		while(r!=null){
			if(i<a){System.out.println(r.toText(false));}
			r=fris.next();
			if(++i>=a){break;}
		}
		while(r!=null && i++<b){r=fris.next();}
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	public FastaShredInputStream(String fname, boolean amino_, boolean allowSubprocess_, long maxdata){
		this(FileFormat.testInput(fname, FileFormat.FASTA, FileFormat.FASTA, 0, allowSubprocess_, false, false), amino_, maxdata);
	}
	
	public FastaShredInputStream(FileFormat ff, boolean amino_, long maxData_){
		name=ff.name();
		amino=amino_;
		flag=(amino ? Read.AAMASK : 0);
		
		if(!fileIO.FileFormat.hasFastaExtension(name) && !name.startsWith("stdin")){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+name);
		}
		
		allowSubprocess=ff.allowSubprocess();
		minLen=MIN_READ_LEN;
		maxLen=TARGET_READ_LEN;
		maxData=maxData_>0 ? maxData_ : Shared.bufferData();
		
		bf=open();
		
		assert(settingsOK());
	}
	
	@Override
	public Read next() {
		throw new RuntimeException("Unsupported");
	}
	
	@Override
	public ArrayList<Read> nextList() {
		if(currentList==null){
			boolean b=fillList();
		}
		ArrayList<Read> list=currentList;
		currentList=null;
		if(list==null || list.isEmpty()){
			list=null;
		}else{
			consumed+=list.size();
		}
		return list;
	}
	
	@Override
	public boolean hasMore() {
		if(currentList==null || currentList.size()==0){
			if(open){
				fillList();
			}else{
//				assert(generated>0) : "Was the file empty?";
			}
		}
		return (currentList!=null && currentList.size()>0);
	}
	
	@Override
	public void restart() {
		if(bf!=null){close();}
		assert(bf==null);
//		generated=0;
		consumed=0;
		nextReadID=0;
		currentList=null;
		
		if(bf==null){
			bf=open();
		}else{
			assert(false) : "bf should be null";
		}
	}
	
	@Override
	public final boolean close(){
		synchronized(this){
			if(!open){return false;}
			open=false;
			assert(bf!=null);
			errorState|=bf.close();
			bf=null;
		}
		return false;
	}
	
	@Override
	public boolean paired() {return false;}
	
	@Override
	public void start() {}
	
	private final boolean fillList(){
//		assert(open);
		if(!open){
			currentList=null;
			return false;
		}
		assert(currentList==null);
		currentList=new ArrayList<Read>(BUF_LEN);
		
		long len=0;
		for(int i=0; i<BUF_LEN && len<maxData; i++){
			Read r=generateRead();
			if(r==null){
				close();
				break;
			}
			currentList.add(r);
			len+=r.length();
			nextReadID++;
			if(verbose){System.err.println("Generated a read; i="+i+", BUF_LEN="+BUF_LEN);}
//			if(i==1){assert(false) : r.numericID+", "+r.mate.numericID;}
		}
		
		return currentList.size()>0;
	}
	
	private final Read generateRead(){
		Read r=null;
		boolean eof=false;
		while(r==null && !eof){
			while(buffer.length()<maxLen){
				byte[] line=bf.nextLine();
				if(line==null){
					eof=true;
					break;
				}
				if(line.length>0 && line[0]==carrot){break;}
				buffer.append(line);
			}
			if(buffer.length>=minLen){
				byte[] bases=buffer.expelAndShift(Tools.min(maxLen, buffer.length()), TARGET_READ_OVERLAP);
				r=new Read(bases, null, Long.toString(nextReadID), nextReadID, flag);
				nextReadID++;
				if(verbose){System.err.println("Made read:\t"+(r.length()>1000 ? r.id : r.toString()));}
				if(bases.length<maxLen){buffer.clear();}
			}else{
				buffer.clear();
			}
		}
		return r;
	}
	
	private final ByteFile open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened fasta file "+name);
		}
		open=true;
		ByteFile bf=ByteFile.makeByteFile(name, allowSubprocess);
		return bf;
	}
	
	public boolean isOpen(){return open;}
	
	public static final boolean settingsOK(){
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1){
			throw new RuntimeException("Minimum FASTA read length is too long: "+MIN_READ_LEN);
		}
		if(MIN_READ_LEN<1){
			throw new RuntimeException("Minimum FASTA read length is too short: "+MIN_READ_LEN);
		}
		if(TARGET_READ_LEN<1){
			throw new RuntimeException("Target FASTA read length is too short: "+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>TARGET_READ_LEN){
			throw new RuntimeException("Minimum FASTA read length is longer than maximum read length: "+MIN_READ_LEN+">"+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1 || MIN_READ_LEN<1){return false;}
		if(TARGET_READ_LEN<1 || MIN_READ_LEN>TARGET_READ_LEN){return false;}
		return true;
	}
	
	@Override
	public String fname(){return bf.name();}
	
	public final String name;
	
	private ArrayList<Read> currentList=null;

	private boolean open=false;
	private ByteBuilder buffer=new ByteBuilder();
	public ByteFile bf;
	
	private long nextReadID=0;
	private long consumed=0;

	public final boolean allowSubprocess;
	public final boolean amino;
	public final int flag;
	private final int BUF_LEN=Shared.bufferLen();;
	private final long maxData;
	private final int maxLen, minLen;
	
	
	public static boolean verbose=false;
	private static final byte carrot='>';
	
	public static int TARGET_READ_LEN=800;
	public static int TARGET_READ_OVERLAP=31;
	public static int MIN_READ_LEN=31;
	public static boolean FAKE_QUALITY=false;
	public static boolean WARN_IF_NO_SEQUENCE=true;
	public static boolean WARN_FIRST_TIME_ONLY=true;
	public static boolean ABORT_IF_NO_SEQUENCE=false;
	
}
