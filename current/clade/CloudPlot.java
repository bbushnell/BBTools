package clade;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintStream;

import javax.imageio.ImageIO;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.FloatList;

/**
 * Visualizes 3D compositional metrics (GC, HH, CAGA) as 2D scatter plots with color encoding.
 * Supports TSV input with future expansion for FASTA format via Scalars integration.
 * Generates PNG images with configurable scaling and point sizes.
 *
 * @author Brian Bushnell
 * @contributor G11
 * @date October 6, 2025
 */
public class CloudPlot {

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
		CloudPlot x=new CloudPlot(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CloudPlot(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			in1=parser.in1;
			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written

		ffout1=FileFormat.testOutput(out1, FileFormat.UNKNOWN, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("xmin")){
				xmin=Float.parseFloat(b);
			}else if(a.equals("xmax")){
				xmax=Float.parseFloat(b);
			}else if(a.equals("ymin")){
				ymin=Float.parseFloat(b);
			}else if(a.equals("ymax")){
				ymax=Float.parseFloat(b);
			}else if(a.equals("zmin")){
				zmin=Float.parseFloat(b);
			}else if(a.equals("zmax")){
				zmax=Float.parseFloat(b);
			}else if(a.equals("scale")){
				scale=Integer.parseInt(b);
			}else if(a.equals("pointsize")){
				pointsize=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		return parser;
	}

	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
		if(scale<1){
			throw new RuntimeException("scale must be >= 1");
		}
		if(pointsize<1){
			throw new RuntimeException("pointsize must be >= 1");
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){

		// Read input data
		readData();

		// Apply autoscaling if needed
		autoscale();

		// Render the plot
		BufferedImage img=renderPlot();

		// Write output
		try{
			writeOutput(img);
		}catch(Exception e){
			throw new RuntimeException("Error writing output file: "+out1, e);
		}

		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, pointsProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Points plotted:    \t"+pointsProcessed);

		//Throw an exception if there was an error
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Read data from input file (TSV or FASTA) */
	private void readData(){
		if(ffin1.isSequence()){
			// TODO: Phase 2 - FASTA input support
			// Will use: data = Scalars.toIntervals(in1);
			throw new RuntimeException("FASTA input not yet supported. Use TSV format: <GC>\\t<HH>\\t<CAGA>");
		}else{
			// TSV input
			data=readTSV();
		}
	}

	/** Read TSV format input */
	private FloatList[] readTSV(){
		FloatList gc=new FloatList();
		FloatList hh=new FloatList();
		FloatList caga=new FloatList();

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		LineParser1 parser=new LineParser1('\t');

		byte[] line;
		while((line=bf.nextLine())!=null){
			if(line.length>0 && line[0]!='#'){
				bytesProcessed+=(line.length+1);
				parser.set(line);
				if(parser.terms()>=3){
					gc.add(parser.parseFloat(0));
					hh.add(parser.parseFloat(1));
					caga.add(parser.parseFloat(2));
					pointsProcessed++;
				}
			}
		}
		bf.close();

		return new FloatList[]{gc, hh, caga};
	}

	/** Apply autoscaling to any axis with negative min/max values */
	private void autoscale(){
		if(data==null || data[0].size()<1){
			throw new RuntimeException("No data points to plot");
		}

		// X-axis (GC)
		if(xmin<0){xmin=findMin(data[0]);}
		if(xmax<0){xmax=data[0].max();}

		// Y-axis (HH)
		if(ymin<0){ymin=findMin(data[1]);}
		if(ymax<0){ymax=data[1].max();}

		// Z-axis/Color (CAGA)
		if(zmin<0){zmin=findMin(data[2]);}
		if(zmax<0){zmax=data[2].max();}
	}

	/** Find minimum value in FloatList (FloatList has max() but not min()) */
	private float findMin(FloatList list){
		if(list.size()<1){return 0;}
		float min=list.array[0];
		for(int i=1; i<list.size(); i++){
			min=Math.min(min, list.array[i]);
		}
		return min;
	}

	/** Render the plot to a BufferedImage */
	private BufferedImage renderPlot(){
		int width=800*scale;
		int height=600*scale;
		int margin=50*scale;

		BufferedImage img=new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g=img.createGraphics();

		// White background
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, width, height);

		// Draw points
		int plotWidth=width-2*margin;
		int plotHeight=height-2*margin;

		int numPoints=data[0].size();
		for(int i=0; i<numPoints; i++){
			float gcVal=data[0].array[i];    // X-axis
			float hhVal=data[1].array[i];    // Y-axis
			float cagaVal=data[2].array[i];  // Color

			// Map data coords to pixel coords
			float normX=(gcVal-xmin)/(xmax-xmin);
			float normY=(hhVal-ymin)/(ymax-ymin);

			int px=margin+(int)(normX*plotWidth);
			int py=height-margin-(int)(normY*plotHeight);  // Flip Y-axis

			Color c=cagaToColor(cagaVal);
			g.setColor(c);
			g.fillOval(px-pointsize, py-pointsize, pointsize*2, pointsize*2);
		}

		g.dispose();
		return img;
	}

	/** Map CAGA value to color: Red(0.0) → Blue(0.5) → Green(1.0) */
	private Color cagaToColor(float caga){
		if(caga<=0.5f){
			// Red (0.0) → Blue (0.5)
			float t=caga*2.0f;  // Map 0.0-0.5 to 0.0-1.0
			return interpolateColor(new Color(255, 0, 0), new Color(0, 0, 255), t);
		}else{
			// Blue (0.5) → Green (1.0)
			float t=(caga-0.5f)*2.0f;  // Map 0.5-1.0 to 0.0-1.0
			return interpolateColor(new Color(0, 0, 255), new Color(0, 255, 0), t);
		}
	}

	/** Interpolate between two colors */
	private Color interpolateColor(Color c1, Color c2, float t){
		int r=(int)(c1.getRed()+t*(c2.getRed()-c1.getRed()));
		int g=(int)(c1.getGreen()+t*(c2.getGreen()-c1.getGreen()));
		int b=(int)(c1.getBlue()+t*(c2.getBlue()-c1.getBlue()));
		return new Color(
			Math.max(0, Math.min(255, r)),
			Math.max(0, Math.min(255, g)),
			Math.max(0, Math.min(255, b))
		);
	}

	/** Write BufferedImage to PNG file */
	private void writeOutput(BufferedImage img) throws Exception{
		File outFile=new File(out1);
		ImageIO.write(img, "png", outFile);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Data storage: [0]=GC, [1]=HH, [2]=CAGA */
	private FloatList[] data=null;

	/*--------------------------------------------------------------*/

	/** Axis range parameters (negative = autoscale) */
	private float xmin=-1.0f;
	private float xmax=-1.0f;
	private float ymin=-1.0f;
	private float ymax=-1.0f;
	private float zmin=-1.0f;
	private float zmax=-1.0f;

	/** Rendering parameters */
	private int scale=1;
	private int pointsize=2;

	/*--------------------------------------------------------------*/

	private long pointsProcessed=0;
	private long bytesProcessed=0;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input File */
	private final FileFormat ffin1;
	/** Output File */
	private final FileFormat ffout1;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

}
