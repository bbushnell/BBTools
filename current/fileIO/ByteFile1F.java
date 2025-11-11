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
import structures.IntList;
import structures.ListNum;

/**
 * ByteFile variant that splits on FASTA record boundaries (\n>) instead of newlines.
 * Each "line" returned is a complete FASTA record (header + sequence).
 * Uses SIMD for boundary detection.
 * 
 * @author Brian Bushnell
 * @date January 2026
 */
public final class ByteFile1F extends ByteFile {

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
		ByteFile1F tf=new ByteFile1F(args.length>0 ? args[0] : "stdin", true);
		speedtest(tf, first, last, !speedtest);

		tf.close();
		tf.reset();
		tf.close();
	}

	private static void speedtest(ByteFile1F bf, long first, long last, boolean reprint){
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
				System.out.println(new String(s));
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

	public ByteFile1F(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false));
	}

	public ByteFile1F(FileFormat ff){
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
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		// For the first record, find the first '>'
		if(firstRecord){
			if(bstop==0){fillBuffer();}
			while(bstart<bstop && buffer[bstart]!=carrot){bstart++;}
			if(bstart>=bstop){
				close();
				return null;
			}
			firstRecord=false;
		}

		// Need more positions?
		if(listPos>=positions.size()){
			fillBuffer();
			if(listPos>=positions.size() && bstart>=bstop){
				// EOF
				close();
				return null;
			}
		}

		lineNum++;

		// Determine record boundaries
		final int limit;
		if(listPos<positions.size()){
			// Next boundary is at positions[listPos] (points to \n before >)
			limit=positions.get(listPos++);
		}else{
			// No more boundaries, use rest of buffer
			limit=bstop;
		}

		if(bstart>=limit){
			// Empty record or at end
			if(bstart<bstop){
				bstart++; // Skip the >
			}
			return blankLine;
		}

		byte[] record=KillSwitch.copyOfRange(buffer, bstart, limit);

		// Move past this record (\n> or end)
		if(limit<bstop){
			bstart=limit+2; // Skip \n>
		}else{
			bstart=bstop;
		}

		return record;
	}
	
	@Override
	public final ListNum<byte[]> nextList(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);
		for(int i=0; i<listSize; i++){
			byte[] record=nextLine();
			if(record==null){
				break;
			}
			list.add(record);
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Fill buffer and find all FASTA boundaries */
	private void fillBuffer(){
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
		while(true){
			if(bstop==buffer.length){
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			
			int r=-1;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}catch(NullPointerException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}
			
			if(r<=0){break;} // EOF
			
			bstop+=r;
			
			// Check if we have enough data or hit a natural boundary
			if(bstop>=65536){break;}
		}

		// Find all \n> boundaries in the buffer
		Vector.findFastaHeaders(buffer, 0, bstop, positions);
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

	public static boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=65536;

	private boolean errorState=false;
}