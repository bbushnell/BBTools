package hiseq;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicLongArray;

import align2.QualityTools;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.LineParser2;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

public class TileDump {
	
	/*--------------------------------------------------------------*/
	/*----------------              CLI             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args) {
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		TileDump td=new TileDump(args);
		
		//Run the object
		td.process(t);
	}
	
	public void process(Timer t) {
		final FlowCell fc0=loadDump(in);
		FlowCell fc=fc0;

		targetX=Tools.max(targetX, Tile.xSize);
		targetY=Tools.max(targetY, Tile.ySize);
		if(targetX>Tile.xSize || targetY>Tile.ySize) {fc=fc.widen(targetX, targetY, true);}
		if(targetReads>0) {fc=fc.widenToTargetReads(targetReads);}
		if(blurTiles) {fc.blur();}
		
		//Temporarily widen to calculate a regression
		if(fc.avgAlignedReads<targetAlignedReads && fc.readsAligned>targetAlignedReads) {
			final int oldX=Tile.xSize, oldY=Tile.ySize;
			FlowCell temp=fc.widenToTargetAlignedReads(targetAlignedReads);
			temp.calcStats();
			fc.uniqueToReadErrorRateFormula=temp.uniqueToReadErrorRateFormula;
			fc.uniqueToBaseErrorRateFormula=temp.uniqueToBaseErrorRateFormula;
			Tile.xSize=oldX;
			Tile.ySize=oldY;
		}
		
		ArrayList<MicroTile> mtList=fc.calcStats();
		markTiles(fc, mtList, System.err);
//		long readsToDiscard=markTiles(mtList, fc.avgReads);//Non-static method...
		
		if(out!=null) {
			write(fc, out, overwrite);
		}
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public TileDump(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			in=parser.in1;
			out=parser.out1;
		}
		checkFileExistence(); //Ensure files can be read and written
	}
	
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

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("blur") || a.equals("blurtiles") || a.equals("smoothtiles")){
				blurTiles=Parse.parseBoolean(b);
			}else if(a.equals("x") || a.equals("xsize")){
				targetX=Integer.parseInt(b);
			}else if(a.equals("y") || a.equals("ysize")){
				targetY=Integer.parseInt(b);
			}else if(a.equals("target") || a.equals("targetreads") || a.equals("reads")){
				targetReads=Parse.parseIntKMG(b);
			}else if(a.equals("targetalignedreads") || a.equals("alignedreads")){
				targetAlignedReads=Parse.parseIntKMG(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parseStatic(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	public static boolean parseStatic(String arg, String a, String b) {
		if(a.equals("qdeviations") || a.equals("qd")){
			qDeviations=Float.parseFloat(b);
		}else if(a.equals("udeviations") || a.equals("ud")){
			uDeviations=Float.parseFloat(b);
		}else if(a.equals("edeviations") || a.equals("ed")){
			eDeviations=Float.parseFloat(b);
		}else if(a.equals("pgdeviations") || a.equals("pgd")){
			pgDeviations=Float.parseFloat(b);
		}
		else if(a.equals("qfraction") || a.equals("qf")){
			qualFraction=Float.parseFloat(b);
		}else if(a.equals("ufraction") || a.equals("uf")){
			uniqueFraction=Float.parseFloat(b);
		}else if(a.equals("efraction") || a.equals("ef")){
			errorFreeFraction=Float.parseFloat(b);
		}else if(a.equals("pgfraction") || a.equals("pgf")){
			polyGFraction=Float.parseFloat(b);
		}
		else if(a.equals("qabsolute") || a.equals("qa")){
			qualAbs=Float.parseFloat(b);
		}else if(a.equals("uabsolute") || a.equals("ua")){
			uniqueAbs=Float.parseFloat(b);
		}else if(a.equals("eabsolute") || a.equals("ea")){
			errorFreeAbs=Float.parseFloat(b);
		}else if(a.equals("pgabsolute") || a.equals("pga")){
			polyGAbs=Float.parseFloat(b);
		}
		else if(a.equals("mbf") || a.equals("mdf") || a.equals("maxbadfraction") || 
				a.equals("maxbadtilefraction") || a.equals("maxdiscardfraction")){
			maxDiscardFraction=Float.parseFloat(b);
		}else if(a.equals("impliederrorrate") || a.equals("inferrederrorrate") 
				|| a.equals("ier") || a.equals("maxier")) {
			maxImpliedErrorRate=Float.parseFloat(b);
		}else if(a.equals("inferredquality") || a.equals("impliedquality") || a.equals("miniq")) {
			float miniq=Float.parseFloat(b);
			maxImpliedErrorRate=(float)QualityTools.phredToProbError(miniq);
		}
		
		else {
			return false;
		}
		return true;
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			System.err.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+
					"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Writing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void write(FlowCell fc, String fname, boolean overwrite) {
		assert(fname!=null);
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();

		bsw.print("#Version\t").print(TileDump.VERSION_OUT).nl();
		
		bsw.print("#id\t").print(fc.name).nl();
		bsw.print("#lane\t");
		{
			String comma="";
			for(Lane lane : fc.lanes) {
				if(lane!=null && !lane.isEmpty()) {
					bsw.print(comma).print(lane.lane);
					comma=",";
				}
			}
			bsw.println();
		}
		assert(fc.k>0) : fc.k;
		bsw.print("#k\t").print(fc.k).nl();
		
		bsw.print("#xSize\t").print(Tile.xSize).nl();
		bsw.print("#ySize\t").print(Tile.ySize).nl();
		if(fc.xMax>=0 || fc.yMax>=0 || fc.xMin>=0 || fc.yMin>=0) {
			bsw.print("#xMin\t").print(fc.xMin).nl();
			bsw.print("#xMax\t").print(fc.xMax).nl();
			bsw.print("#yMin\t").print(fc.yMin).nl();
			bsw.print("#yMax\t").print(fc.yMax).nl();
			bsw.print("#tMin\t").print(fc.tMin).nl();
			bsw.print("#tMax\t").print(fc.tMax).nl();
		}
		bsw.print("#reads\t").print(fc.readsProcessed).nl();
		bsw.print("#avgReads\t").print(fc.avgReads, 1).nl();
		bsw.print("#avgLen\t").print(fc.basesProcessed/(1.0*fc.readsProcessed), 1).nl();
		bsw.print("#alignReads\t").print(fc.readsAligned, 1).nl();
		bsw.print("#alignBases\t").print(fc.basesAligned, 1).nl();
		bsw.print("#alignRate\t").print(fc.readsAligned/(double)fc.readsProcessed, 5).nl();
		bsw.print("#readErrRate\t").print(fc.readErrors/(double)fc.readsAligned, 5).nl();
		bsw.print("#baseErrRate\t").print(fc.baseErrors/(double)fc.basesAligned, 5).nl();
		
		bsw.print("#avgQuality\t").print(fc.avgQuality, 5).nl();
		bsw.print("#avgUnique\t").print(fc.avgUnique, 5).nl();
		bsw.print("#avgErrorFree\t").print(fc.avgErrorFree, 5).nl();
		bsw.print("#avgPolyG\t").print(fc.avgPolyG, 5).nl();
		if(fc.avgG>0) {bsw.print("#avgG\t").print(fc.avgG, 5).nl();}
		
		bsw.print("#stdQuality\t").print(fc.stdQuality, 6).nl();
		bsw.print("#stdUnique\t").print(fc.stdUnique, 6).nl();
		bsw.print("#stdErrorFree\t").print(fc.stdErrorFree, 6).nl();
		bsw.print("#stdPolyG\t").print(fc.stdPolyG, 6).nl();
		if(fc.avgG>0) {bsw.print("#stdG\t").print(fc.stdG, 6).nl();}
		
		if(fc.uniqueToReadErrorRateFormula!=null) {
			bsw.print("#uniqueToReadErrorRate\t").print(fc.uniqueToReadErrorRateFormula[0], 5).print('+').print(fc.uniqueToReadErrorRateFormula[1],5).print('X').nl();
			bsw.print("#uniqueToBaseErrorRate\t").print(fc.uniqueToBaseErrorRateFormula[0], 5).print('+').print(fc.uniqueToBaseErrorRateFormula[1],5).print('X').nl();
		}
		
		boolean printCycleStats=false;
		for(Lane lane : fc.lanes){
			if(lane!=null && printCycleStats){
				for(int pairnum=0; pairnum<2; pairnum++) {
					if(lane.depthCounts[pairnum].get(0)>0) {
						printTSV(bsw, lane.depthSums[pairnum], "#depthSums_"+lane.lane+"_"+pairnum);
						printTSV(bsw, lane.depthCounts[pairnum], "#depthCounts_"+lane.lane+"_"+pairnum);
					}
					if(lane.matchCounts[pairnum].get(0)>0) {
						printTSV(bsw, lane.matchCounts[pairnum], "#matchCounts_"+lane.lane+"_"+pairnum);
						printTSV(bsw, lane.subCounts[pairnum], "#subCounts_"+lane.lane+"_"+pairnum);
					}
				}
			}
		}
		
		bsw.print("#");
		bsw.println(MicroTile.header());
		
		for(Lane lane : fc.lanes){
			if(lane!=null){
				lane.print(bsw, fc.k, fc.uniqueToReadErrorRateFormula, fc.uniqueToBaseErrorRateFormula);
			}
		}
		bsw.poisonAndWait();
	}
	
	private static void printTSV(ByteStreamWriter bsw, AtomicLongArray list, String label) {
		bsw.print(label);
		int maxNonZero=-1;
		for(int i=0; i<list.length(); i++) {
			if(list.get(i)>0) {maxNonZero=i;}
		}
		for(int i=0; i<=maxNonZero; i++) {
			bsw.tab().print(list.get(i));
		}
		bsw.nl();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static FlowCell loadDump(String fname){
		return loadDump(fname, new FlowCell(31));
	}
	
	public static FlowCell loadDump(String fname, FlowCell fc){
		ByteFile bf=ByteFile.makeByteFile(fname, false);

		final int version=parseHeader(bf, fc);
		long tiles=0;
		if(version==1) {
			tiles=loadTiles1(bf, fc);
		}else if(version>1) {
			tiles=loadTiles2(bf, fc, version);
//		}else if(version==3) {
//			tiles=loadTiles3(bf, fc);
		}else {
			throw new RuntimeException("Unhandled dump file version "+version);
		}
		return fc;
	}

	public static int parseHeader(ByteFile bf, FlowCell fc) {
		LineParser1 lp=new LineParser1('\t');
		int version=1;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]=='#'){
				lp.set(line);
//				String[] split=new String(line).split("\t");
//				System.err.println(new String(line));
				if(Tools.startsWith(line, "#Version")){
					VERSION_IN=version=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#xSize")){
					Tile.xSize=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#ySize")){
					Tile.ySize=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#xMin")){
					fc.xMin=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#xMax")){
					fc.xMax=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#yMin")){
					fc.yMin=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#yMax")){
					fc.yMax=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#tMin")){
					fc.tMin=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#tMax")){
					fc.tMax=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#k")){
					fc.k=lp.parseInt(1);
				}else if(Tools.startsWith(line, "#id")){
					fc.name=lp.parseString(1);
				}else if(Tools.startsWith(line, "#lane")){
					//ignore
				}
				
				else if(Tools.startsWith(line, "#avgQuality")){
					fc.avgQuality=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#avgUnique")){
					fc.avgUnique=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#avgErrorFree")){
					fc.avgErrorFree=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#avgPolyG")){
					fc.avgPolyG=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#avgG")){
					fc.avgG=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#stdQuality")){
					fc.stdQuality=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#stdUnique")){
					fc.stdUnique=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#stdErrorFree")){
					fc.stdErrorFree=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#stdPolyG")){
					fc.stdPolyG=lp.parseDouble(1);
				}else if(Tools.startsWith(line, "#stdG")){
					fc.stdG=lp.parseDouble(1);
				}
				
				else if(Tools.startsWith(line, "#reads")){
					fc.readsProcessed=lp.parseLong(1);
				}else if(Tools.startsWith(line, "#avgReads")){
					fc.avgReads=lp.parseDouble(1);
				}else {}
			}else{
				bf.pushBack(line);
				break;
			}
		}
		return version;
	}
	
	public static long loadTiles1(ByteFile bf, FlowCell fc) {
		final LineParser2 lp=new LineParser2('\t');
		long lines=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			assert(line[0]!='#') : "Unexpected header line "+bf.lineNum()+":\n"+new String(line);
			lp.set(line);
			int lane=lp.parseInt();
			int tile=lp.parseInt();
			int x1=lp.parseInt();
			int x2=lp.parseInt();
			int y1=lp.parseInt();
			int y2=lp.parseInt();

			MicroTile mt=fc.getMicroTile(lane, tile, x1, y1);
			assert(mt.x1==x1 && mt.x2==x2) : 
				"Micro-tile size seems to be different:\n"+ 
				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
				+ "mt.x1="+mt.x1+", mt.x2="+mt.x2+", x1="+x1+", x2="+x2+"\n"
				+"line='"+new String(line)+"'";
			assert(mt.y1==y1 && mt.y2==y2) : 
				"Micro-tile size seems to be different:\n"+ 
				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
				+ "mt.y1="+mt.y1+", mt.y2="+mt.y2+", y1="+y1+",y2="+y2;

			final long reads=mt.readCount=lp.parseLong();
			mt.alignedReadCount=lp.parseLong();
			mt.alignedBaseCount=lp.parseLong();
			mt.readErrorCount=lp.parseLong();
			mt.baseErrorCount=lp.parseLong();
			mt.kmerReadErrorCount=lp.parseLong();
			mt.kmerBaseErrorCount=lp.parseLong();
			mt.readInsCount=lp.parseLong();
			mt.readDelCount=lp.parseLong();

			float uniquePercent=lp.parseFloat();
			float averageQualityProb=lp.parseFloat();
			float percentErrorFree=lp.parseFloat();
			float depth=lp.parseFloat();
			float alignmentRate=lp.parseFloat();
			float trueQuality=lp.parseFloat();
			float kErrRateR=lp.parseFloat();
			float kErrRateB=lp.parseFloat();
			float insRate=lp.parseFloat();
			float delRate=lp.parseFloat();

			mt.discard=lp.parseInt();

			long a=lp.parseLong();
			long c=lp.parseLong();
			long g=lp.parseLong();
			long t=lp.parseLong();
			long n=lp.parseLong();
			long hmpCount=lp.parseLong();
			long hmpSum=lp.parseLong();
			float hmpRate=lp.parseFloat();

			mt.misses=(long)(uniquePercent*reads*0.01);
			mt.hits=reads-mt.misses;
			mt.readQualityByProbSum=reads*averageQualityProb;
			mt.probErrorFreeSum=reads*percentErrorFree;
			mt.depthSum=(long)(reads*depth);

			mt.acgtn=new long[] {a,c,g,t,n};
			mt.homoPolyGCount=hmpCount;
			mt.homoPolyGSum=hmpSum;
		}
		return lines;
	}
	
	public static long loadTiles2(ByteFile bf, FlowCell fc, int dumpVersion) {
		final LineParser2 lp=new LineParser2('\t');
		long lines=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			assert(line[0]!='#') : "Unexpected header line "+bf.lineNum()+":\n"+new String(line);
			lp.set(line);
			int lane=lp.parseInt();
			int tile=lp.parseInt();
			int x1=lp.parseInt();
			int x2=lp.parseInt();
			int y1=lp.parseInt();
			int y2=lp.parseInt();

			MicroTile mt=fc.getMicroTile(lane, tile, x1, y1);
			assert(mt.x1==x1 && mt.x2==x2) : 
				"Micro-tile size seems to be different:\n"+ 
				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
				+ "mt.x1="+mt.x1+", mt.x2="+mt.x2+", x1="+x1+", x2="+x2+"\n"
				+"line='"+new String(line)+"'";
			assert(mt.y1==y1 && mt.y2==y2) : 
				"Micro-tile size seems to be different:\n"+ 
				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
				+ "mt.y1="+mt.y1+", mt.y2="+mt.y2+", y1="+y1+",y2="+y2;

			final long reads=mt.readCount=lp.parseLong();
			final long bases=mt.baseCount=lp.parseLong();
			mt.alignedReadCount=lp.parseLong();
			mt.alignedBaseCount=lp.parseLong();
			mt.readErrorCount=lp.parseLong();
			mt.baseErrorCount=lp.parseLong();
			mt.kmerReadErrorCount=lp.parseLong();
			mt.kmerBaseErrorCount=lp.parseLong();
			mt.readInsCount=lp.parseLong();
			mt.readDelCount=lp.parseLong();

			float rer=lp.parseFloat();
			float ber=lp.parseFloat();
			float uniquePercent=lp.parseFloat();
			float averageQualityProb=lp.parseFloat();
			float percentErrorFree=lp.parseFloat();
			float averageBaseErrorProb=lp.parseFloat();
			float depth=lp.parseFloat();
			float ierate1=lp.parseFloat();
			float ierate2=lp.parseFloat();
			float ierate3=lp.parseFloat();
			float iqscore=lp.parseFloat();
			float alignmentRate=lp.parseFloat();
			float trueQuality=lp.parseFloat();
			float kErrRateR=lp.parseFloat();
			float kErrRateB=lp.parseFloat();
			float insRate=lp.parseFloat();
			float delRate=lp.parseFloat();

			mt.discard=lp.parseInt();

			long a=lp.parseLong();
			long c=lp.parseLong();
			long g=lp.parseLong();
			long t=lp.parseLong();
			long n=lp.parseLong();
			long hmpCount=lp.parseLong();
			long hmpSum=lp.parseLong();
			
			mt.misses=(long)(uniquePercent*reads*0.01);
			mt.hits=reads-mt.misses;
			mt.readQualityByProbSum=reads*averageQualityProb;
			mt.probErrorFreeSum=reads*percentErrorFree;
			mt.baseErrorProbSum=bases*averageBaseErrorProb;
			mt.depthSum=(long)(reads*depth);
			
			mt.acgtn=new long[] {a,c,g,t,n};
			mt.homoPolyGCount=hmpCount;
			mt.homoPolyGSum=hmpSum;
			
			float polyRate=lp.parseFloat();
			mt.barcodes=lp.parseLong();
			mt.barcodeHDistSum=lp.parseLong();
			mt.barcodePolymers=lp.parseLong();
			
			float hdist=lp.parseFloat();
			float bcpoly=lp.parseFloat();
			
			if(dumpVersion>2) {
				mt.validKmerSum=lp.parseLong();
				mt.goodKmerSum=lp.parseDouble();
				float goodKmerFraction=lp.parseFloat();
			}
			
			if(lp.hasMore()) {
				mt.mergedReads=lp.parseLong();
				mt.insertSum=lp.parseLong();
				mt.overlapSum=lp.parseLong();
				mt.mergeErrorSum=lp.parseLong();
				float avgInsert=lp.parseFloat();
				float mergeRate=lp.parseFloat();
				float mergBaseErrorRate=lp.parseFloat();
			}
		}
		return lines;
	}
	
//	public static long loadTiles3(ByteFile bf, FlowCell fc) {
//		final LineParser2 lp=new LineParser2('\t');
//		long lines=0;
//		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
//			assert(line[0]!='#') : "Unexpected header line "+bf.lineNum()+":\n"+new String(line);
//			lp.set(line);
//			int lane=lp.parseInt();
//			int tile=lp.parseInt();
//			int x1=lp.parseInt();
//			int x2=lp.parseInt();
//			int y1=lp.parseInt();
//			int y2=lp.parseInt();
//
//			MicroTile mt=fc.getMicroTile(lane, tile, x1, y1);
//			assert(mt.x1==x1 && mt.x2==x2) : 
//				"Micro-tile size seems to be different:\n"+ 
//				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
//				+ "mt.x1="+mt.x1+", mt.x2="+mt.x2+", x1="+x1+", x2="+x2+"\n"
//				+"line='"+new String(line)+"'";
//			assert(mt.y1==y1 && mt.y2==y2) : 
//				"Micro-tile size seems to be different:\n"+ 
//				"xsize="+Tile.xSize+", ysize="+Tile.ySize+"\n"
//				+ "mt.y1="+mt.y1+", mt.y2="+mt.y2+", y1="+y1+",y2="+y2;
//
//			final long reads=mt.readCount=lp.parseLong();
//			final long bases=mt.baseCount=lp.parseLong();
//			mt.alignedReadCount=lp.parseLong();
//			mt.alignedBaseCount=lp.parseLong();
//			mt.readErrorCount=lp.parseLong();
//			mt.baseErrorCount=lp.parseLong();
//			mt.kmerReadErrorCount=lp.parseLong();
//			mt.kmerBaseErrorCount=lp.parseLong();
//			mt.readInsCount=lp.parseLong();
//			mt.readDelCount=lp.parseLong();
//
//			float rer=lp.parseFloat();
//			float ber=lp.parseFloat();
//			float uniquePercent=lp.parseFloat();
//			float averageQualityProb=lp.parseFloat();
//			float percentErrorFree=lp.parseFloat();
//			float averageBaseErrorProb=lp.parseFloat();
//			float depth=lp.parseFloat();
//			float ierate1=lp.parseFloat();
//			float ierate2=lp.parseFloat();
//			float ierate3=lp.parseFloat();
//			float iqscore=lp.parseFloat();
//			float alignmentRate=lp.parseFloat();
//			float trueQuality=lp.parseFloat();
//			float kErrRateR=lp.parseFloat();
//			float kErrRateB=lp.parseFloat();
//			float insRate=lp.parseFloat();
//			float delRate=lp.parseFloat();
//
//			mt.discard=lp.parseInt();
//
//			long a=lp.parseLong();
//			long c=lp.parseLong();
//			long g=lp.parseLong();
//			long t=lp.parseLong();
//			long n=lp.parseLong();
//			long hmpCount=lp.parseLong();
//			long hmpSum=lp.parseLong();
//
//			mt.misses=(long)(uniquePercent*reads*0.01);
//			mt.hits=reads-mt.misses;
//			mt.readQualityByProbSum=reads*averageQualityProb;
//			mt.probErrorFreeSum=reads*percentErrorFree;
//			mt.baseErrorProbSum=bases*averageBaseErrorProb;
//			mt.depthSum=(long)(reads*depth);
//
//			mt.acgtn=new long[] {a,c,g,t,n};
//			mt.homoPolyGCount=hmpCount;
//			mt.homoPolyGSum=hmpSum;
//			
//			float polyRate=lp.parseFloat();
//			mt.barcodes=lp.parseLong();
//			mt.barcodeHDistSum=lp.parseLong();
//			mt.barcodePolymers=lp.parseLong();
//			
//			float hdist=lp.parseFloat();
//			float bcpoly=lp.parseFloat();
//			
//			mt.validKmerSum=lp.parseLong();
//			mt.goodKmerSum=lp.parseDouble();
//			float goodKmerFraction=lp.parseFloat();
//			
//			if(lp.hasMore()) {
//				mt.mergedReads=lp.parseLong();
//				mt.insertSum=lp.parseLong();
//				mt.overlapSum=lp.parseLong();
//				mt.mergeErrorSum=lp.parseLong();
//				float avgInsert=lp.parseFloat();
//				float mergeRate=lp.parseFloat();
//				float mergBaseErrorRate=lp.parseFloat();
//			}
//		}
//		return lines;
//	}
	
	static final long markTiles(FlowCell flowcell, ArrayList<MicroTile> mtList, PrintStream outstream){
		for(MicroTile mt : mtList){
			mt.discard=0;
		}
		long readsToDiscard=0;
		
		long cDiscards=0, qDiscards=0, iqDiscards=0, eDiscards=0;
		long uDiscards=0, mtDiscards=0, gDiscards=0, mtRetained=0;
		final double alignmentRate=flowcell.alignmentRate();
		final double avgReads=flowcell.avgReads;
		for(MicroTile mt : mtList){
			double q=mt.averageReadQualityByProb();
			double e=mt.percentErrorFree();
			double u=mt.uniquePercent();
			double pg=mt.polyGPercent();
			double g=mt.maxG();
			double ier=mt.impliedErrorRate(flowcell.uniqueToBaseErrorRateFormula);
//			assert(false) : "\nier="+ier+", u="+u+", berf="+Arrays.toString(flowcell.uniqueToBaseErrorRateFormula)
//			+", maxImpliedErrorRate="+maxImpliedErrorRate+", alignmentRate="+alignmentRate;
			
			double dq=flowcell.avgQuality-q;
			double de=flowcell.avgErrorFree-e;
			double du=u-flowcell.avgUnique;
			double dpg=pg-flowcell.avgPolyG;
			double dg=g-flowcell.avgG;
			
			if(mt.readCount<10 && mt.readCount<0.02f*avgReads){
				mt.discard++;
				cDiscards++;
			}
			
			if(dq>qDeviations*flowcell.stdQuality && dq>flowcell.avgQuality*qualFraction && dq>qualAbs){
				mt.discard++;
				qDiscards++;
			}
			if(ier>maxImpliedErrorRate && alignmentRate>0.0001) {
				mt.discard++;
				iqDiscards++;
			}
			if(de>eDeviations*flowcell.stdErrorFree && 
					de>flowcell.avgErrorFree*errorFreeFraction && de>errorFreeAbs){
				mt.discard++;
				eDiscards++;
			}
//			assert(de<=0) : "e="+e+", de="+de+", eDeviations="+eDeviations+
//			", flowcell.stdErrorFree="+flowcell.stdErrorFree+", flowcell.avgErrorFree="+flowcell.avgErrorFree
//			+", errorFreeFraction="+errorFreeFraction+", errorFreeAbs="+errorFreeAbs;
			if(dpg>pgDeviations*flowcell.stdPolyG && 
					dpg>flowcell.avgPolyG*polyGFraction && dpg>polyGAbs){
				mt.discard++;
				gDiscards++;
//				assert(dpg<0) : "pg="+pg+", dpg="+dpg+", devs="+pgDeviations+", std="+flowcell.stdPolyG+"\n"
//				+ "avg="+flowcell.avgPolyG+", frac="+polyGFraction+", abs="+polyGAbs+
//				", t1="+ (pgDeviations*flowcell.stdPolyG)+", t2="+(flowcell.avgPolyG*polyGFraction);
			}
			if(flowcell.avgUnique>2 && flowcell.avgUnique<98){
				if(du>uDeviations*flowcell.stdUnique && du>flowcell.avgUnique*uniqueFraction && du>uniqueAbs){
					mt.discard++;
					uDiscards++;
				}
			}
			
//			if((discardG || gToN) && (dg>gDeviations*flowcell.stdG && dg>flowcell.avgG*gFraction && dg>gAbs)){
//				mt.discard++;
//				gDiscards++;
//			}
			if(mt.discard>0){
				mtDiscards++;
				readsToDiscard+=mt.readCount;
			}
			else{mtRetained++;}
		}
		
		long fullSize=mtList.size()-cDiscards;
		long fullSizeDiscards=mtDiscards-cDiscards;
		long maxDiscards=(long)(maxDiscardFraction*fullSize);
		if(fullSizeDiscards>maxDiscards) {
			Collections.sort(mtList);
			for(MicroTile mt : mtList) {
				if(mt.discard>0 && mt.readCount>=10 && mt.readCount>=0.02f*avgReads) {//full size discard
					mt.discard=0;
					mtDiscards--;
					fullSizeDiscards--;
					mtRetained++;
					readsToDiscard-=mt.readCount;
					if(fullSizeDiscards<=maxDiscards) {break;}
				}
			}
		}
		
		outstream.println();
		outstream.println("Flagged "+mtDiscards+" of "+(mtDiscards+mtRetained)+" micro-tiles, containing "+readsToDiscard+" reads:");
		outstream.println(uDiscards+" exceeded uniqueness thresholds.");
		outstream.println(qDiscards+" exceeded quality thresholds.");
		outstream.println(iqDiscards+" exceeded implied quality thresholds.");
		outstream.println(eDiscards+" exceeded error-free probability thresholds.");
		outstream.println(gDiscards+" contained G spikes.");
		outstream.println(cDiscards+" had too few reads to calculate statistics.");
		outstream.println();
		
		return readsToDiscard;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int VERSION_OUT=3;
	public static int VERSION_IN=1;
	public static boolean verbose;
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static float maxDiscardFraction=0.4f;
	static float maxImpliedErrorRate=0.012f;
	
	static float qDeviations=2.4f;
	static float uDeviations=1.5f;
	static float eDeviations=3.0f;
	static float pgDeviations=1.4f;
	static float gDeviations=3f;
	
	static float qualFraction=0.08f;
	static float uniqueFraction=0.01f;
	static float errorFreeFraction=0.2f;
	static float polyGFraction=0.2f;
	static float gFraction=.1f;
	
	static float qualAbs=2.0f;
	static float uniqueAbs=1f;
	static float errorFreeAbs=6f;
	static float polyGAbs=0.2f;
	static float gAbs=.1f;
	
	private String in=null;
	private String out=null;
	private int targetX=-1;
	private int targetY=-1;
	private int targetReads=-1;
	private int targetAlignedReads=250;
	private boolean blurTiles=false;
	private boolean overwrite=true;
	
}
