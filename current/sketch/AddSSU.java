package sketch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tax.TaxTree;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class AddSSU {
	
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
		AddSSU x=new AddSSU(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AddSSU(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
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
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.SKETCH, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.SKETCH, null, true, false);
		
		if(verbose){
			System.err.println("Set r16SFile="+r16SFile);
			System.err.println("Set r18SFile="+r18SFile);
		}
		
		tree=(treeFile!=null && (preferSSUMapEuks || preferSSUMapProks || clear16SEuks || clear18SEuks || 
				clear16SProks || clear18SProks || useSSUMapOnlyEuks || useSSUMapOnlyProks) ? TaxTree.loadTaxTree(treeFile, outstream, false, false) : null);
		
		if(preferSSUMapEuks || preferSSUMapProks || clear16SEuks || clear18SEuks || clear16SProks || clear18SProks || useSSUMapOnlyEuks || useSSUMapOnlyProks){
			assert(tree!=null) : "preferSSUMapForEuks, clear16SEuks, and clear18SEuks require a TaxTree.";
		}
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

			if(a.equalsIgnoreCase("16S") || a.equalsIgnoreCase("16Sfile")){
				r16SFile=b;
			}else if(a.equalsIgnoreCase("18S") || a.equalsIgnoreCase("18Sfile")){
				r18SFile=b;
			}else if(a.equalsIgnoreCase("tree") || a.equalsIgnoreCase("treefile")){
				treeFile=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				ByteFile1.verbose=verbose;
//				ByteFile2.verbose=verbose;
//				ReadWrite.verbose=verbose;
			}
			
			else if(a.equalsIgnoreCase("preferSSUMap")){
				preferSSUMap=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("preferSSUMapForEuks") || a.equalsIgnoreCase("preferSSUMapEuks")){
				preferSSUMapEuks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("useSSUMapOnly")){
				useSSUMapOnly=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("useSSUMapOnlyEuks") || a.equalsIgnoreCase("SSUMapOnlyEuks")){
				useSSUMapOnlyEuks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("useSSUMapOnlyProks") || a.equalsIgnoreCase("SSUMapOnlyProks")){
				useSSUMapOnlyProks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("preferSSUMapForProks") || a.equalsIgnoreCase("preferSSUMapProks")){
				preferSSUMapProks=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("clearAll")){
				clear16S=clear18S=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear16S")){
				clear16S=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear18S")){
				clear18S=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear16SEuks")){
				clear16SEuks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear18SEuks")){
				clear18SEuks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear16SProks")){
				clear16SProks=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("clear18SProks")){
				clear18SProks=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if("auto".equalsIgnoreCase(r16SFile)){r16SFile=TaxTree.default16SFile();}
		if("auto".equalsIgnoreCase(r18SFile)){r18SFile=TaxTree.default18SFile();}
		SSUMap.r16SFile=r16SFile;
		SSUMap.r18SFile=r18SFile;
		
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
		if(!Tools.testInputFiles(false, true, in1, r16SFile, r18SFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		assert(in1!=null) : "Input sketch file is required";
		assert(r16SFile!=null || r18SFile!=null) : "Input SSU file is required";
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, r16SFile, r18SFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		
		processInner(bf, bsw);
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println(Tools.linesBytesOut(linesProcessed, bytesProcessed, linesOut, bytesOut, 8, true));
		
		outstream.println();
		outstream.println(Tools.number("Sketches:", sketchCount, 8));
		outstream.println(Tools.number("16S In:", r16Sin, 8));
		outstream.println(Tools.number("18S In:", r18Sin, 8));
		outstream.println(Tools.number("16S Added:", r16SfromMap, 8));
		outstream.println(Tools.number("18S Added:", r18SfromMap, 8));
		outstream.println(Tools.numberPercent("16S Out:", r16Sout, r16Sout*100.0/sketchCount, 2, 8));
		outstream.println(Tools.numberPercent("18S Out:", r18Sout, r18Sout*100.0/sketchCount, 2, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
//	private void processInner_old(ByteFile bf, ByteStreamWriter bsw){
//		SSUMap.load(outstream);
//		
//		if(verbose){
//			System.err.println("Loaded SSUMap; |16S|="+SSUMap.r16SCount()+", |18S|="+SSUMap.r18SCount());
//		}
//		
//		byte[] line=bf.nextLine();
////		ByteBuilder bb=new ByteBuilder();
//
//		final byte[] ssuBytes="SSU:".getBytes();
//		final byte[] r16SBytes="16S:".getBytes();
//		final byte[] r18SBytes="18S:".getBytes();
//
//		while(line!=null){
//			if(line.length>0){
//				if(maxLines>0 && linesProcessed>=maxLines){break;}
//				linesProcessed++;
//				bytesProcessed+=(line.length+1);
//
//				final boolean header=(line[0]=='#');
//
//				linesOut++;
//				bytesOut+=(line.length+1);
//
//				if(header){
//					if(Tools.startsWith(line, "#SZ:")){
//						sketchCount++;
//
//						bsw.print(line);
//
//						final int tid=parseTaxID(line);
//						final boolean has16S=Tools.contains(line, ssuBytes, 0) || Tools.contains(line, r16SBytes, 0);
//						final boolean has18S=Tools.contains(line, r18SBytes, 0);
//
//						if(verbose){
//							System.err.println("For line "+new String(line)+":");
//							System.err.println("tid="+tid+", has16S="+has16S+", has18S="+has18S);
//						}
//
//						if(tid>0){
//							final byte[] r16S=has16S ? null : SSUMap.r16SMap.get(tid);
//							final byte[] r18S=has18S ? null : SSUMap.r18SMap.get(tid);
//							if(r16S!=null){bsw.print("\t16S:").print(r16S.length); ssuOut++;}
//							if(r18S!=null){bsw.print("\t18S:").print(r18S.length); ssuOut++;}
//							if(r16S!=null){bsw.print("\n#16S:").print(r16S);}
//							if(r18S!=null){bsw.print("\n#18S:").print(r18S);}
//
//							if(verbose){System.err.println("Found 16S: "+(r16S!=null)+"; found 18S: "+(r18S!=null));}
//						}
//						bsw.println();
//					}else if(Tools.startsWith(line, "#16S:") || Tools.startsWith(line, "#18S:") || Tools.startsWith(line, "#SSU:")){
//						bsw.println(line);
//						ssuIn++;
//						ssuOut++;
//					}else{
//						assert(Tools.startsWith(line, "##")) : new String(line);
//						bsw.println(line);
//					}
//				}else{
//					bsw.println(line);
//				}
//			}
//			line=bf.nextLine();
//		}
//	}
	
	private void processInner(ByteFile bf, ByteStreamWriter bsw){
		SSUMap.load(outstream);
		
		if(verbose){
			System.err.println("Loaded SSUMap; |16S|="+SSUMap.r16SCount()+", |18S|="+SSUMap.r18SCount());
		}
		
		byte[] line=bf.nextLine();
//		ByteBuilder bb=new ByteBuilder();

//		final byte[] ssuBytes="SSU:".getBytes();
//		final byte[] r16SBytes="16S:".getBytes();
//		final byte[] r18SBytes="18S:".getBytes();

		SketchHeader header=null;
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				final boolean isHeader=(line[0]=='#');

				if(isHeader){
					if(Tools.startsWith(line, "#SZ:")){
						assert(header==null) : "\nReplacing this:\n"+header.toBytes()+"\nWith this:\n"+new String(line)+"\n";
						header=new SketchHeader(line);
						sketchCount++;
					}else if(Tools.startsWith(line, "##")){
						bsw.println(line);

						linesOut++;
						bytesOut+=(line.length+1);
					}else{
						header.addLine(line);
					}
				}else{
					if(header!=null){
						try {
							processHeader(header);
						} catch (Throwable e) {
							e.printStackTrace();
							assert(false) : header.toBytes();
						}
						r16Sout+=(header.r16S==null ? 0 : 1);
						r18Sout+=(header.r18S==null ? 0 : 1);
						linesOut+=1+(header.r16S==null ? 0 : 1)+(header.r18S==null ? 0 : 1);
						ByteBuilder bb=header.toBytes();
						bytesOut+=(bb.length+1);
						bsw.println(bb);
						header=null;
					}
					bsw.println(line);

					linesOut++;
					bytesOut+=(line.length+1);
				}
			}
			line=bf.nextLine();
		}
	}
	
	void processHeader(SketchHeader header){
		
		if(verbose){System.err.println("Processing tid "+header.tid+":\n"+header.toBytes()+"\n");}

		final boolean euk=(tree!=null && header.tid>0 && header.tid<SketchObject.minFakeID) ? tree.isEukaryote(header.tid) : false;
		final boolean prok=(tree!=null && header.tid>0 && header.tid<SketchObject.minFakeID) ? tree.isProkaryote(header.tid) : false;
		if(useSSUMapOnly || (useSSUMapOnlyEuks && euk) || (useSSUMapOnlyProks && prok)){header.r16S=header.r18S=null;}
		if(header.tid>0){
			final boolean preferMap=(preferSSUMap || (preferSSUMapEuks && euk) || (preferSSUMapProks && prok));
			byte[] r16S=(SSUMap.r16SMap==null ? null : SSUMap.r16SMap.get(header.tid));
			byte[] r18S=(SSUMap.r18SMap==null ? null : SSUMap.r18SMap.get(header.tid));
			if(r16S!=null && (preferMap || header.r16S==null)){
				header.r16S=r16S;
				r16SfromMap++;
			}
			if(r18S!=null && (preferMap || header.r18S==null)){
				header.r18S=r18S;
				r18SfromMap++;
			}
		}
		if(clear16S || (clear16SEuks && euk) || (clear16SProks && prok)){header.r16S=null;}
		if(clear18S || (clear18SEuks && euk) || (clear18SProks && prok)){header.r18S=null;}
	}
	
	int parseTaxID(byte[] line){
		String[] split=Tools.tabPattern.split(new String(line));
		for(String s : split){
			if(s.startsWith("ID:") || s.startsWith("TAXID:")){
				final int colon=s.indexOf(':');
				final String sub=s.substring(colon+1);
				return Integer.parseInt(sub);
			}
		}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	
	//A very limited parser
	private class SketchHeader {
		
		SketchHeader(byte[] line){
			this(new String(line, 1, line.length-1));
		}
		
		SketchHeader(String line){
			if(line.charAt(0)=='#'){line=line.substring(1);}
			assert(line.startsWith("SZ:"));
			String[] split=Tools.tabPattern.split(line);
			fields=new ArrayList<String>(line.length()+2);
			int tid_=-1;
			for(String s : split){
				if(s.startsWith("16S:") || s.startsWith("18S:") || s.startsWith("SSU:")){
					//do nothing
				}else{
					if(s.startsWith("ID:") || s.startsWith("TAXID:")){
						final int colon=s.indexOf(':');
						final String sub=s.substring(colon+1);
						tid_=Integer.parseInt(sub);
					}
					fields.add(s);
				}
			}
			tid=tid_;
		}
		
		void addLine(byte[] line){
			assert(line[0]=='#');
			assert(line[1]=='1' || line[1]=='S') : new String(line);
			if(Tools.startsWith(line, "#16S:") || Tools.startsWith(line, "#SSU:")){
				assert(r16S==null);
				r16S=Arrays.copyOfRange(line, 5, line.length);
				r16Sin++;
			}else if(Tools.startsWith(line, "#18S:")){
				assert(r18S==null);
				r18S=Arrays.copyOfRange(line, 5, line.length);
				r18Sin++;
			}else{
				assert(false) : new String(line);
			}
		}
		
		ByteBuilder toBytes(){
			ByteBuilder bb=new ByteBuilder(1000);
			bb.append('#');
			for(int i=0; i<fields.size(); i++){
				if(i>0){bb.tab();}
				bb.append(fields.get(i));
			}
			if(r16S!=null){bb.tab().append("16S:").append(r16S.length);}
			if(r18S!=null){bb.tab().append("18S:").append(r18S.length);}
			
			if(r16S!=null){bb.nl().append("#16S:").append(r16S);}
			if(r18S!=null){bb.nl().append("#18S:").append(r18S);}
			return bb;
		}
		
		final int tid;
		ArrayList<String> fields;
		byte[] r16S;
		byte[] r18S;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String r16SFile="auto";
	private String r18SFile="auto";
	private String treeFile="auto";

	boolean preferSSUMap=false;
	boolean preferSSUMapEuks=false;
	boolean preferSSUMapProks=false;
	boolean useSSUMapOnly=false;
	boolean useSSUMapOnlyEuks=false;
	boolean useSSUMapOnlyProks=false;
	boolean clear16S=false;
	boolean clear18S=false;
	boolean clear16SEuks=false;
	boolean clear18SEuks=false;
	boolean clear16SProks=false;
	boolean clear18SProks=false;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long sketchCount=0;
	
	private long r16Sin=0;
	private long r16Sout=0;
	private long r16SfromMap=0;
	private long r18Sin=0;
	private long r18Sout=0;
	private long r18SfromMap=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	private final TaxTree tree;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
