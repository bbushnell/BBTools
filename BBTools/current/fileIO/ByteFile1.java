package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * @author Brian Bushnell
 *
 */
public final class ByteFile1 extends ByteFile {
	
	
	public static void main(String[] args){
		ByteFile1 tf=new ByteFile1(args.length>0 ? args[0] : "stdin", true);
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
		speedtest(tf, first, last, !speedtest);
		
		tf.close();
		tf.reset();
		tf.close();
		
//		tf.reset();
//		byte[] b=tf.nextLine();
//		System.err.println(new String(b));
//		tf.pushBack(b);
//		b=tf.nextLine();
//		System.err.println(new String(b));
//		tf.bstart+=2;
//		tf.pushBack(b);
//		b=tf.nextLine();
//		System.err.println(new String(b));
//		b=tf.nextLine();
//		System.err.println(new String(b));
//		b=tf.nextLine();
//		System.err.println(new String(b));
//		b=tf.nextLine();
//		System.err.println(new String(b));
	}
	
	private static void speedtest(ByteFile1 tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length;
				System.out.println(new String(s));
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
			}
		}
		t.stop();
		
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}
	
	public ByteFile1(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}
	
	public ByteFile1(FileFormat ff){
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
//		if(pushBack!=null){//This was making the massive post-Meltdown-patch performance problem more common
//			byte[] temp=pushBack;
//			pushBack=null;
//			return temp;
//		}
		
		if(verbose){System.err.println("Reading line "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

//		System.out.println("\nCalled nextLine() for line "+lineNum);
//		System.out.println("A: bstart="+bstart+", bstop="+bstop);
		
		//if(bstart<bstop && lasteol==slasher && buffer[bstart]==slashn){bstart++;}
//		assert(bstart>=bstop || (buffer[bstart]!=slashn)/*buffer[bstart]>slasher || buffer[bstart]==slashn*/);
		int nlpos=bstart;
		
//		System.out.println("B: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
//		while(nlpos<bstop && (buffer[nlpos]>slasher || buffer[nlpos]==tab)){nlpos++;}
		while(nlpos<bstop && buffer[nlpos]!=slashn){nlpos++;}//Another vectorization place
//		System.out.println("C: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		if(nlpos>=bstop){
			nlpos=fillBuffer();
//			System.out.println("Filled buffer.");
		}
//		System.out.println("D: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		
		if(nlpos<0 || bstop<1){
			close();
			return null;
		}

		lineNum++;
		//Limit is the position after the last position to copy.
		//Limit equals nlpos unless there was a \r before the \n.
		final int limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		if(bstart==limit){//Empty line.
			bstart=nlpos+1;
//			System.out.println("E: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos+", returning='"+printNL(blankLine)+"'");
			return blankLine;
		}
		
		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);
//		byte[] line=Arrays.copyOfRange(buffer, bstart, limit);
//		byte[] line=new String(buffer, bstart, limit-bstart).getBytes();
//		byte[] line=new byte[limit-bstart];
//		byte[] line=dummy;
		
		assert(line.length>0) : bstart+", "+nlpos+", "+limit;
		bstart=nlpos+1;
//		System.out.println("F: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos+", returning='"+printNL(line)+"'");
		return line;
	}
	
	final byte[] dummy=new byte[100];
	
	private static final String printNL(byte[] b){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<b.length; i++){
			char c=(char)b[i];
			if(c=='\n'){
				sb.append("\\n");
			}else if(c==slashr){
				sb.append("\\r");
			}else{
				sb.append(c);
			}
		}
		return sb.toString();
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
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess());
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
	
	public static boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=16384;

	private boolean errorState=false;
	
	
}
