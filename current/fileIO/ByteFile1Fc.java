package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.bam.BgzfSettings;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/**
 * ByteFile variant that splits on FASTA record boundaries (\n>) instead of newlines.
 * Each "line" returned is a block of complete FASTA records (header + sequence),
 * with \r stripped and only 2 \n per record.  Specifically, the blocks 
 * start with > and end with \n which signifies the end of a record.
 * Uses SIMD for boundary detection.
 * 
 * @author Brian Bushnell
 * @date November 11, 2025
 */
public final class ByteFile1Fc extends ByteFile {

	public static void main(String[] args){
		long first=0, last=100;
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
			if(args[2].equalsIgnoreCase("simd")){
				Shared.SIMD=true;
			}else {
				last=Integer.parseInt(args[2]);
			}
		}
		if(args.length>3){
			if(args[3].equalsIgnoreCase("native")){
				ReadWrite.ALLOW_NATIVE_BGZF=BgzfSettings.USE_MULTITHREADED_BGZF=true;
			}
		}
		ByteFile1Fc tf=new ByteFile1Fc(args.length>0 ? args[0] : "stdin", true);
		speedtest(tf, first, last, !speedtest);

		tf.close();
		tf.reset();
		tf.close();
	}

	private static void speedtest(ByteFile1Fc bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long records=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}
				records++;
				bytes+=s.length+1;
				System.out.print("*");
				System.out.print(new String(s));
				System.out.println("*");
			}
			System.err.println("\n");
			System.err.println("Records: "+records);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}
				records++;
				bytes+=s.length+1;
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, records, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
	}

	public ByteFile1Fc(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false));
	}

	public ByteFile1Fc(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1F("+ff+")");}
		is=open();
	}

	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
		firstRecord=true;
		listPos=0;
		positions.clear();
	}

	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		errorState|=ReadWrite.finishReading(is, name(), (allowSubprocess() || FileFormat.isBamFile(name())));

		is=null;
		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public final byte[] nextLine(){
		return nextLine(new IntList());
	}
	
	public final byte[] nextLine(IntList newlines){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		int lastRecordLoc=fillBuffer();
		byte[] records=condense(buffer, positions, newlines, lastRecordLoc);
		return records;
	}
	
	@Override
	public final ListNum<byte[]> nextList(){
		return nextList(new IntList());
	}
	
	public final ListNum<byte[]> nextList(IntList newlines){
		byte[] records=nextLine();
		if(records==null) {return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>(1);
		list.add(records);
		return new ListNum<byte[]>(list, nextID++);
	}

	/** Fill buffer and find all FASTA boundaries */
	private int fillBuffer(){
		// Shift remaining data to start
		if(bstart>0 && bstart<bstop){
			int extra=bstop-bstart;
			System.arraycopy(buffer, bstart, buffer, 0, extra);
			bstop=extra;
			bstart=0;
		}else if(bstart>=bstop){
			bstart=0;
			bstop=0;
		}

		// Clear position list
		listPos=0;
		positions.clear();

		// Read data
		int lastRecordLoc=-1;
		while(lastRecordLoc<0){
			if(bstop==buffer.length){
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			
			int r=0;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}catch(NullPointerException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}
			
//			final int from=Math.max(0, bstop-1);//Not safe because buffer can come in with residual newlines
			bstop+=r;
			
			if(firstRecord){
			    assert(bstop<1 || buffer[0]==carrot) : "File does not start with '>' - "+ff.name()+" is not a valid FASTA file";
			    firstRecord=false;
			}

			// Find all \n> boundaries in the buffer
			Vector.findSymbols(buffer, 0, bstop, slashn, positions);
			for(int i=positions.size-1; i>=0 && lastRecordLoc<0; i--) {
				int loc=positions.array[i];
				if(loc+1<bstop && buffer[loc+1]==carrot) {lastRecordLoc=loc;}
			}
			
			if(r<=0){// EOF
				if(bstop>0 && buffer[bstop-1]!='\n') {positions.add(bstop);}
				break;
			}
		}
		return lastRecordLoc;
	}
	
	private byte[] condense(byte[] buffer, IntList newlines1, IntList newlines2, int lastRecordLoc) {
		int lastByte=lastRecordLoc>0 ? lastRecordLoc : bstop;
		newlines2.clear();
		if(lastByte<1) {
			close();
			return null;
		}
		ByteBuilder bb=new ByteBuilder(lastByte+1);
		int prev=-1;
		if(verbose) {System.err.println("Entered condense: lastRecordLoc="+
			lastRecordLoc+", bstop="+bstop+", newlines="+newlines1);}
		for(int i=0; i<newlines1.size && prev<lastByte; i++) {
			int loc=newlines1.get(i);
			int len=loc-prev-1;
			System.arraycopy(buffer, prev+1, bb.array, bb.length, len);//Copy excluding newline
			if(verbose) {System.err.println("Added "+len+" bytes.");}
			bb.length+=len;
			if(bb.endsWith('\r')){bb.length--;}//Trim \r
			if(buffer[prev+1]==carrot){//Record start
				newlines2.add(bb.length);
				bb.nl();
				if(verbose) {System.err.println("Added header start newline.");}
			}
			if(loc>=lastByte || buffer[loc+1]==carrot) {//Record end
				newlines2.add(bb.length);
				bb.nl();
				if(verbose) {System.err.println("Added record stop newline.");}
			}
			prev=loc;
		}
		//Shift residual
		int residual=bstop-lastRecordLoc-1;
		if(lastRecordLoc>0 && residual>0) {
			System.arraycopy(buffer, lastRecordLoc+1, buffer, 0, residual);
			bstop=residual;
		}else {
			bstop=0;
		}
		//assert(false) : "*"+bb.toString()+"*";
		return bb.array;
	}

	@Override
	public void pushBack(byte[] record){
		throw new UnsupportedOperationException("pushBack not supported for ByteFile1F");
	}

	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile1F "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=0;
		bstop=0;
		firstRecord=true;
		listPos=0;
		positions.clear();
		return is;
	}

	@Override
	public boolean isOpen(){return open;}

	@Override
	public final InputStream is(){return is;}

	@Override
	public final long lineNum(){return lineNum;}

	private boolean open=false;
	private byte[] buffer=new byte[bufferlen];
	private static final byte[] blankLine=new byte[0];
	private int bstart=0, bstop=0;
	public InputStream is;
	public long lineNum=-1;
	private IntList positions=new IntList();
	private int listPos=0;
	private boolean firstRecord=true;

	public static final boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=262144;

	private boolean errorState=false;
}