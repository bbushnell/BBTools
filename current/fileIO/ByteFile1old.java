package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.bam.BgzfSettings;
import structures.IntList;
import structures.ListNum;


/**
 * @author Brian Bushnell
 *
 */
public final class ByteFile1old extends ByteFile {


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
		ByteFile1old tf=new ByteFile1old(args.length>0 ? args[0] : "stdin", true);
		speedtest(tf, first, last, !speedtest);

		tf.close();
		tf.reset();
		tf.close();
	}

	private static void speedtest(ByteFile1old bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length+1;
				System.out.println(new String(s));
			}

			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			boolean nextList=false;
			if(nextList) {
				for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
					for(byte[] line : ln.list) {
						lines++;
						bytes+=line.length+1;
					}
				}
			}else if(!Shared.SIMD) {
				for(long i=first; i<last; i++){
					byte[] s=bf.nextLine();
					if(s==null){break;}
					lines++;
					bytes+=s.length+1;
				}
			}else {
				for(long i=first; i<last; i++){
					byte[] s=bf.nextLineList();
					if(s==null){break;}
					lines++;
					bytes+=s.length+1;
				}
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
	}

	public ByteFile1old(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	public ByteFile1old(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1("+ff+")");}
		is=open();
	}

	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
	}

	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		//		assert(false) : name()+","+allowSubprocess();
		errorState|=ReadWrite.finishReading(is, name(), (allowSubprocess() || FileFormat.isBamFile(name())));

		is=null;
		lineNum=-1;
		//		pushBack=null;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public final byte[] nextLine(){
		if(Shared.SIMD) {return nextLineList();}
		if(verbose){System.err.println("Reading line "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}

		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		int nlpos=bstart;

		while(nlpos<bstop && buffer[nlpos]!=slashn){nlpos++;}

		if(nlpos>=bstop){//At this point we are at the last character which may or may not be a newline
			nlpos=fillBuffer();
		}

		// Check if buffer is empty after fill attempt
		if(bstop<1){
			close();
			return null;
		}

		lineNum++;

		// Determine the limit for copying
		// If nlpos >= bstop, we have data but no newline (EOF case) - use bstop
		// Otherwise, nlpos points to the newline
		final int limit;
		if(nlpos >= bstop){
			// No newline found, but we have data - use everything remaining
			limit = bstop;
		} else {
			// Found newline - exclude it (and any preceding \r)
			limit = (nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		}

		if(bstart==limit){//Empty line.
			bstart = (nlpos < bstop) ? nlpos+1 : bstop;
			return blankLine;
		}

		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);

		assert(line.length>0) : bstart+", "+nlpos+", "+limit;

		// Advance bstart past the newline (if there was one)
		bstart = (nlpos < bstop) ? nlpos+1 : bstop;

		return line;
	}

	public final byte[] nextLineList(){

		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		if(listPos>=positions.size()) {
			listPos=0;
			positions.clear();
			fillBuffer();
			Vector.findSymbols(buffer, 0, bstop, slashn, positions);
		}

		final int nlpos;
		if(listPos>=positions.size()) {//No newlines
			nlpos=bstop;
		}else {
			nlpos=positions.get(listPos++);
		}

		// Check if buffer is empty after fill attempt
		if(bstop<1){
			close();
			return null;
		}

		lineNum++;

		// Determine the limit for copying
		// If nlpos >= bstop, we have data but no newline (EOF case) - use bstop
		// Otherwise, nlpos points to the newline
		final int limit;
		if(nlpos >= bstop){
			// No newline found, but we have data - use everything remaining
			limit = bstop;
		} else {
			// Found newline - exclude it (and any preceding \r)
			limit = (nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		}

		if(bstart==limit){//Empty line.
			bstart = (nlpos < bstop) ? nlpos+1 : bstop;
			return blankLine;
		}

		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);

		assert(line.length>0) : bstart+", "+nlpos+", "+limit;

		// Advance bstart past the newline (if there was one)
		bstart = (nlpos < bstop) ? nlpos+1 : bstop;

		return line;
	}
	
	private final ListNum<byte[]> nextListScalar(){
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);
		for(int i=0; i<listSize; i++) {
			byte[] line=nextLine();
			if(line==null) {
				close();
				break;
			}
			list.add(line);
			lineNum++;
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	@Override
	public final ListNum<byte[]> nextList(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		if(!Shared.SIMD) {return nextListScalar();}
		
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);

		while(list.size()<listSize && open) {
			if(listPos>=positions.size()) {
				listPos=0;
				positions.clear();
				fillBuffer();
				Vector.findSymbols(buffer, 0, bstop, slashn, positions);

				// If still no positions after fill, we're done
				if(positions.size()==0){
					break;
				}
			}

			// Process positions until we have enough lines OR run out of positions
			final int iters=Math.min(positions.size-listPos, listSize-list.size());
			for(int i=0; i<iters; i++) {
				final int nlpos = positions.get(listPos++);

				// Check if buffer is empty
				if(bstop<1){
					close();
					break;
				}
				lineNum++;

				// Determine the limit for copying
				final int limit;
				if(nlpos >= bstop){
					limit = bstop;
				} else {
					limit = (nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
				}

				if(bstart==limit){//Empty line
					bstart = (nlpos < bstop) ? nlpos+1 : bstop;
					list.add(blankLine);
				}else {
					byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);
					assert(line.length>0) : bstart+", "+nlpos+", "+limit;
					bstart = (nlpos < bstop) ? nlpos+1 : bstop;
					list.add(line);
				}
			}
		}

		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	private final void printBuffer(){
		for(int i=0; i<bstop; i++){
			char c=(char)buffer[i];
			if(c=='\n'){
				System.err.println("\\n");
			}else if(c==slashr){
				System.err.print("\\r");
			}else{
				System.err.print(c);
			}
		}
	}

	private int fillBuffer(){
		if(bstart<bstop){ //Shift end bytes to beginning
			//			System.err.println("Shift: "+bstart+", "+bstop);
			assert(bstart>0);
			//			assert(bstop==buffer.length);
			int extra=bstop-bstart;
			for(int i=0; i<extra; i++, bstart++){
				//				System.err.print((char)buffer[bstart]);
				//System.err.print('.');
				buffer[i]=buffer[bstart];
				//				assert(buffer[i]>=slasher || buffer[i]==tab);
				assert(buffer[i]!=slashn);
			}
			bstop=extra;
			//			System.err.println();

			//			{//for debugging only
			//				buffer=new byte[bufferlen];
			//				bstop=0;
			//				bstart=0;
			//			}
		}else{
			bstop=0;
		}

		bstart=0;
		int len=bstop;
		int r=-1;
		while(len==bstop){//hit end of input without encountering a newline
			if(bstop==buffer.length){
				//				assert(false) : len+", "+bstop;
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			try {
				r=is.read(buffer, bstop, buffer.length-bstop);
				//				byte[] x=new byte[buffer.length-bstop];
				//				r=is.read(x);
				//				if(r>0){
				//					for(int i=0, j=bstop; i<r; i++, j++){
				//						buffer[j]=x[i];
				//					}
				//				}
			} catch (IOException e) {//java.io.IOException: Stream Closed
				//TODO: This should be avoided rather than caught.  It happens when a stream is shut down with e.g. "reads=100".
				if(!Shared.anomaly){
					e.printStackTrace();
					System.err.println("open="+open);
				}
			} catch (NullPointerException e) {//Can be thrown by java.util.zip.Inflater.ensureOpen(Inflater.java:389)
				//TODO: This should be avoided rather than caught.  It happens when a stream is shut down with e.g. "reads=100".
				if(!Shared.anomaly){
					e.printStackTrace();
					System.err.println("open="+open);
				}
			}
			if(r>0){
				bstop=bstop+r;
				//				//while(len<bstop && (buffer[len]>slasher || buffer[len]==tab)){len++;}//Obsolete; handled old-style Mac convention

				while(len<bstop && buffer[len]!=slashn){len++;}
				//				len=Vector.find(buffer, slashn, len, bstop); //x0.85 speed in simd mode for short sequences; x1.00 for 150bp fastq
			}else{
				len=bstop;
				break;
			}
		}

		//		System.err.println("After Fill: ");
		//		printBuffer();
		//		System.err.println();

		//		System.out.println("Filled buffer; r="+r+", returning "+len);
		assert(r==-1 || buffer[len]==slashn);

		//		System.err.println("lasteol="+(lasteol=='\n' ? "\\n" : lasteol==slashr ? "\\r" : ""+(int)lasteol));
		//		System.err.println("First="+(int)buffer[0]+"\nLastEOL="+(int)lasteol);

		return len;
	}

	@Override
	public void pushBack(byte[] line) {
		if(bstart>line.length){
			bstart--;
			buffer[bstart]='\n';
			for(int i=0, j=bstart-line.length; i<line.length; i++, j++){
				buffer[j]=line[i];
			}
			bstart=bstart-line.length;
			return;
		}

		int bLen=bstop-bstart;
		int newLen=bLen+line.length+1;
		int rShift=line.length+1-bstart;
		assert(rShift>0) : bstop+", "+bstart+", "+line.length;
		while(newLen>buffer.length){
			//This could get big if pushback is used often,
			//unless special steps are taken to prevent it, like leaving extra space for pushbacks.
			buffer=Arrays.copyOf(buffer, buffer.length*2);
		}

		Tools.shiftRight(buffer, rShift);

		for(int i=0; i<line.length; i++){
			buffer[i]=line[i];
		}
		buffer[line.length]='\n';
		bstart=0;
		bstop=newLen;
	}

	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=-1;
		bstop=-1;
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

	public static boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=65536;

	private boolean errorState=false;


}