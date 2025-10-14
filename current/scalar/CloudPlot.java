package scalar;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;

import javax.imageio.ImageIO;

import clade.Clade;
import clade.SendClade;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SendSketch;
import structures.FloatList;
import tax.TaxTree;

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
			maxReads=parser.maxReads;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written

		ffout1=FileFormat.testOutput(out1, FileFormat.PNG, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
		if(useTree) {tree=TaxTree.sharedTree();}
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
			}else if(a.equals("xpct") || a.equals("xpercent")){
				xPercent=Float.parseFloat(b);
			}else if(a.equals("ypct") || a.equals("ypercent")){
				yPercent=Float.parseFloat(b);
			}else if(a.equals("zpct") || a.equals("zpercent")){
				zPercent=Float.parseFloat(b);
			}else if(a.equals("scale")){
				scale=Float.parseFloat(b);
			}else if(a.equals("pointsize")){
				pointsize=Float.parseFloat(b);
			}else if(a.equals("window")){
				window=Parse.parseIntKMG(b);
			}else if(a.equals("interval")){
				interval=Parse.parseIntKMG(b);
			}else if(a.equals("shred")){
				interval=window=Parse.parseIntKMG(b);
			}else if(a.equals("break")){
				breakOnContig=Parse.parseBoolean(b);
			}else if(a.equals("minlen")){
				minlen=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("gcHh") || a.equalsIgnoreCase("hhGc")){
				gcHhCorrelation=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("gcCaga") || a.equalsIgnoreCase("cagaGc")){
				gcCagaCorrelation=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("gcHhstrength") || a.equalsIgnoreCase("gcHhs")){
				gcHhStrength=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("gcCagastrength") || a.equalsIgnoreCase("gcCagas")){
				gcCagaStrength=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("cagaGcstrength") || a.equalsIgnoreCase("cagaGcs")){
				cagaGcStrength=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("hhGcstrength") || a.equalsIgnoreCase("hhGcs")){
				hhGcStrength=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("decorrelate")){
				decorrelate=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("autoscale")){
				autoscale=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("tree") || a.equalsIgnoreCase("usetree")){
				useTree=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("tax") || a.equalsIgnoreCase("colorbytax") || 
				a.equalsIgnoreCase("colorbytid")){
				colorByTax=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("parsetid")){
				ScalarData.parseTID=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("colorByName")){
				colorByName=ScalarIntervals.printName=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("level")){
				level=TaxTree.parseLevelExtended(b);
				useTree=level>1;
			}else if(a.equals("sketch") | a.equals("bbsketch")){
				ScalarData.makeSketch=Parse.parseBoolean(b);
			}else if(a.equals("clade") || a.equals("quickclade")){
				ScalarData.makeClade=Parse.parseBoolean(b);
			}else if(a.equals("mt")){
				ScalarIntervals.mt=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sendInThread")){
				ScalarIntervals.sendInThread=Parse.parseBoolean(b);
			}else if(a.equals("concurrency")){
				SendClade.maxConcurrency=SendSketch.maxConcurrency=Integer.parseInt(b);
			}else if(a.equals("order")){
				order=parseOrder(b);
			}else if(a.equals("concise")){
				Clade.CONCISE=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i>0 && Tools.looksLikeOutputStream(arg) && arg.endsWith(".png")){
				parser.out1=arg;
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
	
	public static int[] parseOrder(String b) {
		HashMap<String, Integer> map=new HashMap<String, Integer>();
		map.put("0", 0);
		map.put("1", 1);
		map.put("2", 2);
		map.put("gc", 0);
		map.put("hh", 1);
		map.put("caga", 2);
		String[] s=b.toLowerCase().split(",");
		int[] order=new int[b.length()];
		for(int i=0; i<order.length; i++) {order[i]=map.get(s[i]);}
		return order;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){

		// Read input data
		readData();

		if(decorrelate) {decorrelate();}

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
			// FASTA input - use ScalarIntervals
			data=ScalarIntervals.toIntervals(in1, window, interval, minlen, breakOnContig, maxReads);
		}else{
			// TSV input
			data=new ScalarData(true, -1).readTSV(ffin1);
		}
		bytesProcessed+=(data.bytesProcessed+data.basesProcessed);
	}
	
	private FloatList decorrelate(FloatList xList, FloatList yList, float correlation, float strength) {
		if(strength==0 || correlation==0) {return xList;}
		return modify(xList, yList, -correlation*strength);
	}
	
	private FloatList modify(FloatList xList, FloatList yList, float mult) {
		FloatList zList=new FloatList(xList.size());
		for(int i=0, lim=xList.size(); i<lim; i++) {
			float x=xList.get(i);
			float y=yList.get(i);
			float z=x+(y-0.5f)*mult;
			zList.add(z);
		}
		return zList;
	}
	
	void decorrelate() {
		//System.err.println("decorrelate");
		final FloatList gc0=data.gc, hh0=data.hh, caga0=data.caga;
//		printMinMax(gc0, hh0, caga0);
		
		FloatList hh=decorrelate(hh0, gc0, gcHhCorrelation, gcHhStrength);
//		printMinMax(gc0, hh, caga0);
		FloatList caga=decorrelate(caga0, gc0, gcCagaCorrelation, gcCagaStrength);
//		printMinMax(gc0, hh, caga);
		FloatList gc=decorrelate(gc0, caga0, gcCagaCorrelation, cagaGcStrength);
//		printMinMax(gc, hh, caga);
		gc=decorrelate(gc, hh0, gcHhCorrelation, hhGcStrength);
//		printMinMax(gc, hh, caga);
		
		data.gc=gc;
		data.hh=hh;
		data.caga=caga;
	}
	
	private static final void printMinMax(FloatList...data) {
		final FloatList gc0=data[0], hh0=data[1], caga0=data[2];
		System.err.println(gc0.min()+"-"+gc0.max()+", "+
			hh0.min()+"-"+hh0.max()+", "+caga0.min()+"-"+caga0.max());
	}
	
	/** Apply autoscaling to any axis with negative min/max values */
	private void autoscale(){
		if(!autoscale) {
			xmin=Tools.mid(xmin, 0, 1);
			xmax=Tools.mid(xmax, 0, 1);
			ymin=Tools.mid(ymin, 0, 1);
			ymax=Tools.mid(ymax, 0, 1);
			zmin=Tools.mid(zmin, 0, 1);
			zmax=Tools.mid(zmax, 0, 1);
		}else{
			FloatList[] lists=data.reorder(order);
			if(xmin<0){xmin=min(lists[0], xPercent);}
			if(xmax<0){xmax=max(lists[0], xPercent);}
			
			if(ymin<0){ymin=min(lists[1], yPercent);}
			if(ymax<0){ymax=max(lists[1], yPercent);}
			
			if(zmin<0){zmin=min(lists[2], zPercent);}
			if(zmax<0){zmax=max(lists[2], zPercent);}
			System.err.println(xmin+"-"+xmax+", "+ymin+"-"+ymax+", "+zmin+"-"+zmax);
		}
	}
	
	private float min(FloatList list, float percentile){
		if(percentile>=1 || percentile<=0) {return list.min();}
		list=list.copy();
		list.sort();
		return list.percentile(1-percentile);
	}
	
	private float max(FloatList list, float percentile){
		if(percentile>=1 || percentile<=0) {return list.max();}
		list=list.copy();
		list.sort();
		return list.percentile(percentile);
	}

	/** Render the plot to a BufferedImage */
	private BufferedImage renderPlot(){
		int width=(int)Math.round(1024*scale);
		int height=(int)Math.round(768*scale);
		int margin=(int)Math.round(50*scale);

		BufferedImage img=new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g=img.createGraphics();

		// White background
		g.setColor(Color.BLACK);
		g.fillRect(0, 0, width, height);

		// Draw points
		int plotWidth=width-2*margin;
		int plotHeight=height-2*margin;

		final FloatList[] lists=data.reorder(order);
		final FloatList xlist=lists[0], ylist=lists[1], zlist=lists[2];
		int numPoints=xlist.size();

//		System.err.println("zmin="+zmin+", zmax="+zmax);
		for(int i=0; i<numPoints; i++){
			int tid=data.tid(i);
			String name=data.name(i);
			float x=xlist.array[i];    // X-axis
			float y=ylist.array[i];    // Y-axis
			float z=zlist.array[i];    // Color

			// Map data coords to pixel coords
			float normX=(x-xmin)/(xmax-xmin);
			float normY=(y-ymin)/(ymax-ymin);

			// Normalize CAGA value to 0-1 range for color mapping
			float normZ;
			if(zmax-zmin<0.001f){
				normZ=0.5f;
			}else{
				normZ=(z-zmin)/(zmax-zmin);
			}
//			System.err.println(cagaVal+" -> "+normZ);

			int px=margin+(int)(normX*plotWidth);
			int py=height-margin-(int)(normY*plotHeight);  // Flip Y-axis

			// Convert normZ (0-1) to radians (0-2π)
			float angle = normZ * (float)(2 * Math.PI);

			// Create elongated ellipse
			float pointWidth=pointsize*0.8f;
			float pointLength=pointsize*(3.2f+0.5f*normY);
			Ellipse2D ellipse=new Ellipse2D.Float(px - pointWidth/2, py - pointLength/2, pointWidth, pointLength);

			// Rotate around center
			AffineTransform rotation = AffineTransform.getRotateInstance(angle, px, py);
			Shape rotatedEllipse = rotation.createTransformedShape(ellipse);

			final Color c;
			if(colorByTax) {
				if(tid<1) {c=new Color(200, 200, 200);}
				else {
					if(useTree && level>1) {tid=tree.getIdAtLevelExtended(tid, level);}
					long hash=Tools.hash64shift(tid);
					float[] rgb=new float[3];
					for(int color=0; color<3; color++) {
						rgb[color]=(((hash&1023L)/1024f)*0.95f)+0.05f;
						hash>>=10;
					}
					float max=Tools.max(rgb);
					float mult=1f/max;
					mult=1+0.6f*(mult-1);
					Tools.multiplyBy(rgb, mult);
					c=new Color(rgb[0], rgb[1], rgb[2]);
				}
			}else if(colorByName) {
					int hash=name.hashCode();
					float[] rgb=new float[3];
					for(int color=0; color<3; color++) {
						rgb[color]=(((hash&1023)/1024f)*0.95f)+0.05f;
						hash>>=10;
					}
					float max=Tools.max(rgb);
					float mult=1f/max;
					mult=1+0.6f*(mult-1);
					Tools.multiplyBy(rgb, mult);
					c=new Color(rgb[0], rgb[1], rgb[2]);
			}else {
				c=cagaToColor6(normZ);
			}
			g.setColor(c);
			g.fill(rotatedEllipse);
			pointsProcessed++;
		}

		g.dispose();
		return img;
	}

	/** Map CAGA value to color: Red(0.0) → Blue(0.5) → Green(1.0) */
	private Color cagaToColor6(float caga){
		
		//Handle out of range.
		if(caga<0) {return new Color(255, 64, 64);}
		if(caga>1) {return new Color(255, 255, 64);}
		
		if(caga<=0.2f){//Red → Purple
			float t=caga*5.0f;
			return interpolateColor(new Color(250, 0, 0), new Color(200, 0, 240), t);
		}else if(caga<=0.4f){//Purple → Blue
			float t=(caga-0.2f)*5.0f;
			return interpolateColor(new Color(200, 0, 240), new Color(32, 32, 255), t);
		}else if(caga<=0.6f){//Blue → Cyan
			float t=(caga-0.4f)*5.0f;
			return interpolateColor(new Color(32, 32, 255), new Color(0, 240, 240), t);
		}else if(caga<=0.8f){//Cyan → Green
			float t=(caga-0.6f)*5.0f;
			return interpolateColor(new Color(0, 240, 240), new Color(0, 200, 0), t);
		}else{//Green → Yellow
			float t=(caga-0.8f)*5.0f;
			return interpolateColor(new Color(0, 200, 0), new Color(250, 200, 0), t);
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
	
	private ScalarData data=null;
	
	private int[] order={2, 1, 0};

	/*--------------------------------------------------------------*/

	/** Axis range parameters (negative = autoscale) */
	private float xmin=-1.0f;
	private float xmax=-1.0f;
	private float ymin=-1.0f;
	private float ymax=-1.0f;
	private float zmin=-1.0f;
	private float zmax=-1.0f;
	private float xPercent=0.998f;
	private float yPercent=0.998f;
	private float zPercent=0.99f;
	
	/** Rendering parameters */
	private float scale=1.0f;
	private float pointsize=3.5f;

	/** FASTA processing parameters */
	private int window=50000;
	private int interval=10000;
	private int minlen=500;
	private boolean breakOnContig=true;
	private long maxReads=-1;
	
	boolean autoscale=true;
	boolean decorrelate=true;
	private float gcHhCorrelation=-0.5f;
	private float gcHhStrength=0.20f;
	private float hhGcStrength=1.40f;
	
	private float gcCagaCorrelation=0.1f;
	private float gcCagaStrength=0.5f;
	private float cagaGcStrength=0.0f;

	private boolean colorByTax=false;
	private boolean colorByName=false;
	private int level=1;
	private boolean useTree=false;
	private static TaxTree tree;

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
	/** Verbose output */
	private static boolean verbose=false;

}
