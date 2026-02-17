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
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * @contributor G11, Neptune
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
			}else if(a.equals("smin") || a.equals("minsize")){
				smin=Float.parseFloat(b);
			}else if(a.equals("smax") || a.equals("maxsize")){
				smax=Float.parseFloat(b);
			}else if(a.equals("spct") || a.equals("spercent") || a.equals("sizepercent")){
				sPercent=Float.parseFloat(b);
			}else if(a.equals("cpct") || a.equals("cpercent") || a.equals("colorpercent")){
				cPercent=Float.parseFloat(b);
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
			}else if(a.equals("cov") || a.equals("coverage") || a.equals("covfile")){
				covFile=b;
			}else if(a.equals("depth") || a.equals("depthfile") || a.equals("sam") || a.equals("bam")){
				depthFile=b;
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
			}else if(a.equals("colorby") || a.equals("color")){
				colorMetric=parseMetric(b);
			}else if(a.equals("concise")){
				Clade.CONCISE=Parse.parseBoolean(b);
			}else if(a.equals("logoffset")){
				logOffset=Float.parseFloat(b);
			}else if(a.equals("logshift")){
				logShift=Float.parseFloat(b);
			}else if(a.equals("logpower")){
				logPower=Float.parseFloat(b);
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

		// Enable name storage if we need to look up depth from coverage file
		if(covFile!=null || depthFile!=null){
			ScalarIntervals.printName=true;
		}

		// Enable TID parsing if we're coloring by taxonomy
		if(colorMetric==TAXONOMY || colorByTax || useTree){
			ScalarData.parseTID=true;
		}

		// Validate taxonomy placement (only color or rotation allowed)
		if(order[0]==TAXONOMY){  // X-axis
			throw new RuntimeException("Taxonomy cannot be assigned to X-axis. Use color or rotation (z) instead.");
		}
		if(order[1]==TAXONOMY){  // Y-axis
			throw new RuntimeException("Taxonomy cannot be assigned to Y-axis. Use color or rotation (z) instead.");
		}
		if(order[3]==TAXONOMY){  // Size
			throw new RuntimeException("Taxonomy cannot be assigned to size. Use color or rotation (z) instead.");
		}

		// Taxonomy on Z-axis (rotation) or color is OK
		// Warn if taxonomy not used for visualization
		boolean taxUsed=(order[2]==TAXONOMY || colorMetric==TAXONOMY);
		if(colorByTax && !taxUsed){
			outstream.println("Warning: colorByTax=true but taxonomy not assigned to any dimension");
		}

		return true;
	}
	
	/**
	 * Parse a single metric name to its index.
	 * @param b Metric name like "gc", "depth", "taxonomy"
	 * @return Metric index constant
	 */
	public static int parseMetric(String b) {
		if(b==null){return NONE;}
		b=b.toLowerCase().trim();
		switch(b){
			case "gc": case "0": return GC;
			case "hh": case "1": return HH;
			case "caga": case "2": return CAGA;
			case "depth": case "cov": case "coverage": case "3": return DEPTH;
			case "length": case "len": case "size": case "4": return LENGTH;
			case "tax": case "taxonomy": case "tid": case "5": return TAXONOMY;
			case "none": case "-1": return NONE;
			default: throw new RuntimeException("Unknown metric: "+b+". Valid: gc, hh, caga, depth, length, taxonomy, none");
		}
	}

	/**
	 * Parse order parameter to assign metrics to dimensions.
	 * Supports 3-element format (x,y,z) or 4-element format (x,y,z,size).
	 * @param b Order string like "gc,hh,caga" or "gc,hh,caga,depth"
	 * @return Array of metric indices [x, y, z, size]
	 */
	public static int[] parseOrder(String b) {
		String[] s=b.toLowerCase().split(",");
		int[] order=new int[4];  // Always return 4 elements: x, y, z, size

		// Parse provided dimensions
		for(int i=0; i<Math.min(s.length, 4); i++) {
			order[i]=parseMetric(s[i]);
		}

		// If only 3 elements provided (legacy format), default size to NONE
		if(s.length==3){
			order[3]=NONE;
		}

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
			// Set ScalarIntervals static fields so it knows about coverage files
			ScalarIntervals.covFile=covFile;
			ScalarIntervals.depthFile=depthFile;
			data=ScalarIntervals.toIntervals(in1, window, interval, minlen, breakOnContig, maxReads);
		}else{
			// TSV input - depth already loaded from TSV
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
			// Size autoscaling: default 0.8x to 3.0x pointsize
			if(smin<0){smin=pointsize*0.8f;}
			if(smax<0){smax=pointsize*3.0f;}
		}else{
			// Get metric lists for each dimension
			FloatList xList=getMetricList(order[0]);
			FloatList yList=getMetricList(order[1]);
			FloatList zList=getMetricList(order[2]);
			FloatList sList=(order[3]!=NONE ? getMetricList(order[3]) : null);

			if(xmin<0 && xList!=null){xmin=minLog(xList, xPercent, order[0]);}
			if(xmax<0 && xList!=null){xmax=maxLog(xList, xPercent, order[0]);}

			if(ymin<0 && yList!=null){ymin=minLog(yList, yPercent, order[1]);}
			if(ymax<0 && yList!=null){ymax=maxLog(yList, yPercent, order[1]);}

			if(zmin<0 && zList!=null){zmin=minLog(zList, zPercent, order[2]);}
			if(zmax<0 && zList!=null){zmax=maxLog(zList, zPercent, order[2]);}

			// Size: always pixel range, track data range separately
			if(smin<0){smin=pointsize*0.8f;}
			if(smax<0){smax=pointsize*3.0f;}

			// Get data range for size dimension (for log scaling)
			if(order[3]!=NONE && sList!=null){
				dataMin=min(sList, sPercent);
				dataMax=max(sList, sPercent);
				System.err.println(xmin+"-"+xmax+", "+ymin+"-"+ymax+", "+zmin+"-"+zmax+", size:"+smin+"-"+smax+" pixels, data:"+dataMin+"-"+dataMax);
			}else{
				System.err.println(xmin+"-"+xmax+", "+ymin+"-"+ymax+", "+zmin+"-"+zmax+", size:"+smin+"-"+smax+" pixels (fixed)");
			}
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

	/** Calculate min percentile - in log space for depth/length, linear for others */
	private float minLog(FloatList list, float percentile, int metric){
		if(metric==DEPTH || metric==LENGTH){
			// Log-transform, get percentile, transform back
			FloatList logList=new FloatList(list.size());
			for(int i=0; i<list.size(); i++){
				logList.add((float)(Math.log(list.get(i)+logOffset)+logShift));
			}
			float logMin=min(logList, percentile);
			return (float)(Math.exp(logMin-logShift)-logOffset);
		}else{
			return min(list, percentile);
		}
	}

	/** Calculate max percentile - in log space for depth/length, linear for others */
	private float maxLog(FloatList list, float percentile, int metric){
		if(metric==DEPTH || metric==LENGTH){
			// Log-transform, get percentile, transform back
			FloatList logList=new FloatList(list.size());
			for(int i=0; i<list.size(); i++){
				logList.add((float)(Math.log(list.get(i)+logOffset)+logShift));
			}
			float logMax=max(logList, percentile);
			return (float)(Math.exp(logMax-logShift)-logOffset);
		}else{
			return max(list, percentile);
		}
	}

	/**
	 * Get the FloatList for a given metric.
	 * @param metric Metric index (GC, HH, CAGA, DEPTH, LENGTH)
	 * @return FloatList containing values for that metric
	 */
	private FloatList getMetricList(int metric){
		switch(metric){
			case GC: return data.gc;
			case HH: return data.hh;
			case CAGA: return data.caga;
			case DEPTH: return data.depth;
			case LENGTH: return data.length;
			case TAXONOMY: return null;  // Categorical, not a FloatList
			case NONE: return null;
			default: throw new RuntimeException("Unknown metric: "+metric);
		}
	}

	/**
	 * Get metric value for a specific data point.
	 * @param index Data point index
	 * @param metric Metric constant
	 * @return Metric value
	 */
	private float getMetricValue(int index, int metric){
		switch(metric){
			case GC: return data.gc.get(index);
			case HH: return data.hh.get(index);
			case CAGA: return data.caga.get(index);
			case DEPTH: return data.depth.get(index);
			case LENGTH: return data.length.get(index);
			case TAXONOMY: return data.tid(index);  // Return TID as float
			case NONE: return 1.0f;  // Fixed value
			default: throw new RuntimeException("Unknown metric: "+metric);
		}
	}

	/** Render the plot to a BufferedImage */
	private BufferedImage renderPlot(){
		int width=(int)Math.round(1024*scale);
		int height=(int)Math.round(768*scale);
		int margin=(int)Math.round(50*scale);

		BufferedImage img=new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g=img.createGraphics();

		// Enable antialiasing
		g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(java.awt.RenderingHints.KEY_RENDERING, java.awt.RenderingHints.VALUE_RENDER_QUALITY);
		g.setRenderingHint(java.awt.RenderingHints.KEY_STROKE_CONTROL, java.awt.RenderingHints.VALUE_STROKE_PURE);

		// White background
		g.setColor(Color.BLACK);
		g.fillRect(0, 0, width, height);

		// Draw points
		int plotWidth=width-2*margin;
		int plotHeight=height-2*margin;

		int numPoints=data.gc.size();  // All FloatLists are same size
		boolean hasSizeDimension=(order[3]!=NONE);

//		System.err.println("zmin="+zmin+", zmax="+zmax);
		for(int i=0; i<numPoints; i++){
			int tid=data.tid(i);
			String name=data.name(i);

			// Get metric values for this point
			float x=getMetricValue(i, order[0]);    // X-axis
			float y=getMetricValue(i, order[1]);    // Y-axis
			float z=getMetricValue(i, order[2]);    // Z-axis (rotation)
			float s=(hasSizeDimension ? getMetricValue(i, order[3]) : 1.0f);  // Size

			// Map data coords to pixel coords - apply log scaling to depth/length
			float normX;
			if(order[0]==DEPTH || order[0]==LENGTH){
				float logX=(float)Math.log(x+logOffset)+logShift;
				float logXmin=(float)Math.log(xmin+logOffset)+logShift;
				float logXmax=(float)Math.log(xmax+logOffset)+logShift;
				normX=(logX-logXmin)/(logXmax-logXmin);
			}else{
				normX=(x-xmin)/(xmax-xmin);
			}

			float normY;
			if(order[1]==DEPTH || order[1]==LENGTH){
				float logY=(float)Math.log(y+logOffset)+logShift;
				float logYmin=(float)Math.log(ymin+logOffset)+logShift;
				float logYmax=(float)Math.log(ymax+logOffset)+logShift;
				normY=(logY-logYmin)/(logYmax-logYmin);
			}else{
				normY=(y-ymin)/(ymax-ymin);
			}

			// Normalize Z value to 0-1 range for rotation
			float normZ;
			if(order[2]==TAXONOMY){
				// Taxonomy rotation: hash TID to get consistent angle per taxon
				if(tid<1){
					normZ=0.5f;  // Default for unclassified
				}else{
					int rotTid=tid;
					if(useTree && level>1) {rotTid=tree.getIdAtLevelExtended(tid, level);}
					long hash=Tools.hash64shift(rotTid);
					normZ=((hash&1023L)/1024f);  // Extract 10 bits, normalize to 0-1
				}
			}else if(order[2]==DEPTH || order[2]==LENGTH){
				// Log scaling for depth/length on Z-axis
				float logZ=(float)Math.log(z+logOffset)+logShift;
				float logZmin=(float)Math.log(zmin+logOffset)+logShift;
				float logZmax=(float)Math.log(zmax+logOffset)+logShift;
				if(logZmax-logZmin<0.001f){
					normZ=0.5f;
				}else{
					normZ=(logZ-logZmin)/(logZmax-logZmin);
				}
			}else if(zmax-zmin<0.001f){
				normZ=0.5f;
			}else{
				normZ=(z-zmin)/(zmax-zmin);
			}

			// Normalize size value and scale to pixel range
			float pixelSize;
			if(hasSizeDimension){
				// Apply appropriate scaling based on size metric
				int sizeMetric=order[3];
				if(sizeMetric==DEPTH || sizeMetric==LENGTH){
					// Log scaling for unbounded metrics: (log(x+offset)+shift)^power
					float logS=(float)Math.pow(Math.log(s+logOffset)+logShift, logPower);
					float logMin=(float)Math.pow(Math.log(dataMin+logOffset)+logShift, logPower);
					float logMax=(float)Math.pow(Math.log(dataMax+logOffset)+logShift, logPower);
					float normS=(logS-logMin)/(logMax-logMin);
					pixelSize=smin+normS*(smax-smin);
				}else{
					// Scaling for 0-1 metrics (GC, HH, CAGA) - apply power for emphasis
					float normS=(s-dataMin)/(dataMax-dataMin);
					normS=(float)Math.pow(normS, logPower);  // Apply power to emphasize differences
					pixelSize=smin+normS*(smax-smin);
				}
				pixelSize=Tools.mid(pixelSize, smin, smax);  // Clamp to range
			}else{
				// Fixed size
				pixelSize=pointsize;
			}

			int px=margin+(int)(normX*plotWidth);
			int py=height-margin-(int)(normY*plotHeight);  // Flip Y-axis

			// Convert normZ (0-1) to radians (0-2π)
			float angle = normZ * (float)(2 * Math.PI);

			// Create elongated ellipse
			float pointWidth=pixelSize*0.8f;
			float pointLength;
			if(hasSizeDimension){
				// When size is a dimension, use proportional length (no Y-dependency)
				pointLength=pixelSize*3.2f;
			}else{
				// Legacy behavior: length varies with Y position
				pointLength=pixelSize*(3.2f+0.5f*normY);
			}
			Ellipse2D ellipse=new Ellipse2D.Float(px - pointWidth/2, py - pointLength/2, pointWidth, pointLength);

			// Rotate around center
			AffineTransform rotation = AffineTransform.getRotateInstance(angle, px, py);
			Shape rotatedEllipse = rotation.createTransformedShape(ellipse);

			// Determine color based on colorMetric
			final Color c;
			if(colorMetric==TAXONOMY) {
				// Color by taxonomy
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
			}else{
				// Color by metric gradient (GC, HH, CAGA, etc.)
				float colorValue=getMetricValue(i, colorMetric);
				float colorMin, colorMax;
				// Get appropriate min/max for color metric
				if(colorMetric==order[0]){colorMin=xmin; colorMax=xmax;}
				else if(colorMetric==order[1]){colorMin=ymin; colorMax=ymax;}
				else if(colorMetric==order[2]){colorMin=zmin; colorMax=zmax;}
				else{
					// Color metric not on any axis - use data range with percentile
					FloatList cList=getMetricList(colorMetric);
					if(cList!=null){
						colorMin=minLog(cList, cPercent, colorMetric);
						colorMax=maxLog(cList, cPercent, colorMetric);
					}else{
						colorMin=0;
						colorMax=1;
					}
				}

				// Apply log scaling to depth/length for color
				float normColor;
				if(colorMetric==DEPTH || colorMetric==LENGTH){
					float logValue=(float)Math.log(colorValue+logOffset)+logShift;
					float logMin=(float)Math.log(colorMin+logOffset)+logShift;
					float logMax=(float)Math.log(colorMax+logOffset)+logShift;
					normColor=(logMax-logMin<0.001f ? 0.5f : (logValue-logMin)/(logMax-logMin));
				}else{
					normColor=(colorMax-colorMin<0.001f ? 0.5f : (colorValue-colorMin)/(colorMax-colorMin));
				}
				c=cagaToColor6(normColor);
			}
			g.setColor(c);
			g.fill(rotatedEllipse);
			pointsProcessed++;
		}

		// Draw axis labels
		drawAxisLabels(g, width, height, margin);

		g.dispose();
		return img;
	}

	/** Draw axis labels showing min/max values */
	private void drawAxisLabels(Graphics2D g, int width, int height, int margin){
		g.setColor(Color.WHITE);
		java.awt.Font font=new java.awt.Font("SansSerif", java.awt.Font.PLAIN, (int)(12*scale));
		g.setFont(font);

		// X-axis labels
		String xLabel=getMetricName(order[0]);
		String xMinLabel=String.format("%s: %.2f", xLabel, xmin);
		String xMaxLabel=String.format("%.2f", xmax);
		g.drawString(xMinLabel, margin, height-margin/4);
		g.drawString(xMaxLabel, width-margin-50*(int)scale, height-margin/4);

		// Y-axis labels
		String yLabel=getMetricName(order[1]);
		String yMinLabel=String.format("%s: %.2f", yLabel, ymin);
		String yMaxLabel=String.format("%.2f", ymax);
		g.drawString(yMinLabel, margin/4, height-margin+15*(int)scale);
		g.drawString(yMaxLabel, margin/4, margin+15*(int)scale);
	}

	/** Get human-readable name for a metric */
	private String getMetricName(int metric){
		switch(metric){
			case GC: return "GC";
			case HH: return "HH";
			case CAGA: return "CAGA";
			case DEPTH: return "Depth";
			case LENGTH: return "Length";
			case TAXONOMY: return "Tax";
			default: return "?";
		}
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

	/** Dimension assignment: [x, y, z/rotation, size] */
	private int[] order={HH, CAGA, GC, NONE};  // Default: x=HH, y=CAGA, z=GC, size=none

	/** Metric to use for color */
	private int colorMetric=TAXONOMY;  // Default: color by taxonomy

	/*--------------------------------------------------------------*/
	/*----------------        Metric Constants      ----------------*/
	/*--------------------------------------------------------------*/

	static final int GC=0, HH=1, CAGA=2, DEPTH=3, LENGTH=4, TAXONOMY=5, NONE=-1;

	/*--------------------------------------------------------------*/

	/** Axis range parameters (negative = autoscale) */
	private float xmin=-1.0f;
	private float xmax=-1.0f;
	private float ymin=-1.0f;
	private float ymax=-1.0f;
	private float zmin=-1.0f;
	private float zmax=-1.0f;
	private float smin=-1.0f;  // Size min in pixels (autoscale if negative)
	private float smax=-1.0f;  // Size max in pixels (autoscale if negative)
	private float dataMin=-1.0f;  // Data range min for size dimension
	private float dataMax=-1.0f;  // Data range max for size dimension
	private float xPercent=0.998f;
	private float yPercent=0.998f;
	private float zPercent=0.99f;
	private float sPercent=0.99f;  // Size percentile for autoscaling
	private float cPercent=0.98f;  // Color percentile for autoscaling
	
	/** Rendering parameters */
	private float scale=1.0f;
	private float pointsize=3.5f;

	/** FASTA processing parameters */
	private int window=50000;
	private int interval=10000;
	private int minlen=500;
	private boolean breakOnContig=true;
	private long maxReads=-1;

	/** Depth/coverage support */
	private String covFile=null;  // Coverage file (pileup or covmaker format)
	private String depthFile=null;  // SAM/BAM file for depth calculation

	/** Log scaling parameters for depth/length */
	private float logOffset=0.25f;  // Add before log to handle zeros
	private float logShift=2.0f;    // Add after log to avoid negatives
	private float logPower=2.0f;    // Power to apply (0.5=sqrt, 2.0=squared for emphasis)
	
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
