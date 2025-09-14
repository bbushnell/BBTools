package driver;

import java.io.PrintStream;
import java.util.Random;

import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date April 9, 2020
 *
 */
public class Life {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Life x=new Life(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Life(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		{//Parse the arguments
			final Parser parser=parse(args);
			parser.out1="stdout.txt";
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		
		current=new byte[rows][columns];
		prev=new byte[rows][columns];
		next=prev;
		
		xMax=columns-1;
		yMax=rows-1;
		bsw=makeBSW(ffout);
		buffer=new byte[columns+1];
		buffer[columns]='\n';
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equalsIgnoreCase("x") || a.equalsIgnoreCase("width") || a.equalsIgnoreCase("columns")){
				columns=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("y") || a.equalsIgnoreCase("height") || a.equalsIgnoreCase("rows")){
				rows=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("cycles") || a.equalsIgnoreCase("rounds")){
				maxCycles=Parse.parseKMG(b);
			}else if(a.equalsIgnoreCase("prob") || a.equalsIgnoreCase("load")){
				prob=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("display")){
				display=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("interval") || a.equalsIgnoreCase("delay")){
				delay=Parse.parseIntKMG(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		init();
		
		processInner();
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesOut, bytesOut, 8));
		outstream.println();
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void init(){
		Random randy=Shared.threadLocalRandom();
		for(int y=0; y<rows; y++){
			for(int x=0; x<columns; x++){
				prev[y][x]=(byte)(randy.nextFloat()<=prob ? 1 : 0);
			}
		}
	}
	
	private void processInner(){

		if(display){
			System.out.print("\033[2J\033[1;1H");
//			bsw.print("\033[2J\033[1;1H");
		}
		for(long cycle=0; cycle<maxCycles; cycle++){
			runCycle(cycle);
		}
	}
	
	private boolean runCycle(long cycle){

		if(display){
			System.out.print("\033[1;1H");
			System.out.println("cycle="+cycle);
			//		bsw.println("\033[1;1H");
			//		bsw.println(""+cycle);
		}
		
		int changes=advance();

//		printMatrix(prev); //For testing
//		printMatrix(current); //For testing
		
		prev=current;
		current=next;
		next=prev;
		
		return changes>0;
	}
	
	private int advance(){
		int changes=0;
		for(int y=0; y<rows; y++){
			changes+=fillRow(y);
			if(display){
				printRow2(current[y]);
				//			printRow(current[y]);
			}
			linesOut++;
			bytesOut+=buffer.length;
		}
		delay(delay);
		return changes;
	}
	
//	//For testing
//	private void printMatrix(byte[][] matrix){
//		System.out.println("-----------------");
//		for(int i=0; i<matrix.length; i++){
//			System.out.println(Arrays.toString(matrix[i]));
//		}
//		System.out.println("-----------------");
//	}
	
	private int fillRow(int y){
		final byte[] dest=current[y];
		final byte[] a=getRow(y-1, prev);
		final byte[] b=getRow(y, prev);
		final byte[] c=getRow(y+1, prev);
		int sum;
		int changed=0;
		
//		printMatrix(new byte[][] {a, b, c});
		
		//Calculate first cell
		sum=a[xMax]+b[xMax]+c[xMax];
		sum+=a[0]+b[0]+c[0];
		sum+=a[1]+b[1]+c[1];
		dest[0]=stateMap[b[0]][sum];
		sum-=(a[xMax]+b[xMax]+c[xMax]);
		
		//Calculate middle cells
		for(int x=1 ; x<xMax; x++){
			final int left=x-1, right=x+1;
			sum+=(a[right]+b[right]+c[right]);
			final byte prevState=b[x];
			final byte nextState=stateMap[prevState][sum];
			dest[x]=nextState;
			sum-=(a[left]+b[left]+c[left]);
//			System.out.println(prevState+" + "+sum+" -> "+nextState);
//			changed+=(ns^ps);//1 if something changed
		}
//		assert(false);

		//Calculate last cell
		sum+=a[0]+b[0]+c[0];
		final byte prevState=b[0];
		final byte nextState=stateMap[prevState][sum];
		dest[xMax]=nextState;

//		changed+=(ns^ps);//1 if something changed
		return changed;
	}
	
	//Using bsw
	private void printRow(byte[] row){
		for(int i=0; i<row.length; i++){
			byte state=row[i];
			buffer[i]=charMap[state];
		}
		bsw.print(new String(buffer));
		linesOut++;
		bytesOut+=row.length;
	}
	
	//Using System.out
	private void printRow2(byte[] row){
		for(int i=0; i<row.length; i++){
			byte state=row[i];
			buffer[i]=charMap[state];
		}
		System.out.print(new String(buffer));
		linesOut++;
		bytesOut+=row.length;
	}
	
	private byte[] getRow(int y, byte[][] matrix){
		return y<0 ? matrix[yMax] : y>yMax ? matrix[0] : matrix[y];
	}
	
	private static TextStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		TextStreamWriter bsw=new TextStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private void delay(int millis){
		if(millis<1){return;}
		
		long until=System.currentTimeMillis()+millis;
		while(System.currentTimeMillis()<until){
			try {
				Thread.sleep(millis);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String out=null;
	private long maxCycles=100;
	
	private int columns=50;
	private int rows=20;
	private float prob=0.25f;
	private boolean display=true;
	private int delay=0;
	
	private final int xMax, yMax;
	private final byte[] buffer;
	
	private final TextStreamWriter bsw;
	
	private byte[][] current;
	private byte[][] prev;
	private byte[][] next;
	
	private final byte[][] stateMap=new byte[][] {
		{0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, //Dead
		{0, 0, 0, 1, 1, 0, 0, 0, 0, 0} //Live (Add 1)
	};
	
	private final byte[] charMap=new byte[] {' ', '@'};
//	private final byte[] charMap=new byte[] {' ', block}; //Eclipse can't handle block character
	
	/*--------------------------------------------------------------*/
	
	private long linesOut=0;
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
