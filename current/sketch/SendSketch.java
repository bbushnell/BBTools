package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import json.JsonObject;
import json.JsonParser;
import kmer.AbstractKmerTableSet;
import server.ServerTools;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.StringNum;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class SendSketch extends SketchObject {
	
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
		
		final boolean oldUnpigz=ReadWrite.USE_UNPIGZ;
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		SendSketch x=new SendSketch(args);
		
		//Run the object
		x.process(t);
		
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		Shared.setBufferLen(oldBufLen);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		assert(!x.errorState) : "This program ended in an error state.";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SendSketch(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
			silent=PreParser.silent;
			if(silent){AbstractKmerTableSet.DISPLAY_PROGRESS=false;}
		}
		
		//Set shared static variables
		ReadWrite.USE_UNPIGZ=true;
		KILL_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		defaultParams.inputVersion=Shared.BBTOOLS_VERSION_STRING;
		defaultParams.mode=PER_FILE;
		boolean setBlacklist=false;
		boolean setLocal=false;
		boolean setPrintDepth=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(a.equals("blacklist") || a.equalsIgnoreCase("bl")){
				setBlacklist=true;
				parseSketchFlags(arg, a, b);
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//Do nothing
			}else if(a.equals("local")){
				local=Parse.parseBoolean(b);
				setLocal=true;
			}else if(a.equals("refid") || a.equals("refids") || a.equals("refname") || a.equals("refnames")){
				refNames=b;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("address")){
				assert(b!=null) : "Bad parameter: "+arg;
				address=b;
			}
//			else if(a.equals("alternateaddress") || a.equals("altaddress") || a.equalsIgnoreCase("vm")){
//				boolean alt=Parse.parseBoolean(b);
//				switchDefaultAddresses(alt);
//			}
			

			else if(a.equalsIgnoreCase("nt") || a.equalsIgnoreCase("silva") || a.equalsIgnoreCase("ribo") || 
					a.equalsIgnoreCase("refseq") || a.equalsIgnoreCase("img") || a.equalsIgnoreCase("nr") || 
					a.equalsIgnoreCase("refseqprot") || a.equalsIgnoreCase("prokprot") || a.equalsIgnoreCase("protein") || 
					a.equalsIgnoreCase("protien") || a.equalsIgnoreCase("prot") || a.equalsIgnoreCase("mito") || a.equalsIgnoreCase("fungi")){
				address=a;
				blacklist=Blacklist.toBlacklist(a);
				System.err.println("Caught "+a);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxID=Integer.parseInt(b);
			}else if(a.equals("spid")){
				outSpid=Integer.parseInt(b);
			}else if(a.equals("imgid")){
				outImgID=Integer.parseInt(b);
			}else if((a.startsWith("meta_") || a.startsWith("mt_")) && b!=null){
				if(outMeta==null){outMeta=new ArrayList<String>();}
				int underscore=a.indexOf('_', 0);
				outMeta.add(a.substring(underscore+1)+":"+b);
			}
			
			else if(a.equals("outsketch") || a.equals("outs") || a.equals("sketchout") || a.equals("sketch")){
				outSketch=b;
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(b==null && new File(arg).isFile()){
				System.err.println("Adding "+arg);
				in.add(arg);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		address=toAddress(address);
		outMeta=SketchObject.fixMeta(outMeta);
		
		if(address!=null && !SET_AUTOSIZE_FACTOR){
			if(address.equals(refseqAddress())){
				AUTOSIZE_FACTOR=2.0f;
			}else if(address.equals(prokProtAddress())){
				AUTOSIZE_FACTOR=3.0f;
			}
		}
		
		while(address!=null && address.endsWith("/")){
			address=address.substring(0,  address.length()-1);
		}
		
		setFromAddress(address, setBlacklist);
		if(local){blacklist=null;}
		
		postParse();
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		//Ensure there is an input file
		if(in.isEmpty() && refNames==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		if(!ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, taxTreeFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(!local && refNames==null){
			if(!Tools.testInputFiles(true, false, in.toArray(new String[0]))){
				if(in.size()==1){
					String s=in.get(0);
					String s1=s.replaceFirst("#", "1"), s2=s.replaceFirst("#", "2");
					Tools.testInputFiles(true, false, s1, s2);
				}else{
					throw new RuntimeException("\nCan't read some input files.\n");  
				}
			}
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out, outSketch)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out+", "+outSketch+"\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		SSUMap.load(outstream);
		tool=new SketchTool(targetSketchSize, defaultParams);
		
//		assert(false) : defaultParams.toString()+"\n"+k+", "+amino+", "+HASH_VERSION;
		if(verbose){
			if(local){outstream.println("Running in local mode.");}
			if(useWhitelist){outstream.println("Using a whitelist.");}
			if(blacklist!=null){outstream.println("Using a blacklist.");}
		}
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile, silent ? null : outstream);}
		defaultParams.postParse(false, false);
		if(!defaultParams.printSSU){processSSU=false;}
		allowMultithreadedFastq=(in.size()==1 && Shared.threads()>2);
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
	}
	
	private void setFromAddress(String address, boolean setBlacklist){
		if(address.equals(nrAddress())){
			amino=true;
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
//			defaultParams.dbName="nr";
			assert(false) : "Need to set K.";
		}else if(address.equals(prokProtAddress())){
			translate=true;
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
//			defaultParams.dbName="nr";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.prokProtBlacklist();}
		}else if(address.equals(ntAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="nt";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.ntBlacklist();}
		}else if(address.equals(refseqAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="RefSeq";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.refseqBlacklist();}
		}else if(address.equals(silvaAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="Silva";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.silvaBlacklist();}
		}else if(address.equals(imgAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="IMG";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.imgBlacklist();}
		}else if(address.equals(mitoAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.mitoBlacklist();}
		}else if(address.equals(fungiAddress())){
			if(!setK){k=defaultK; k2=defaultK2;}
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.fungiBlacklist();}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		if(local){processLocal(t);}
		else if(refNames!=null){processRefMode(t);}
		else{processRemote(t);}
	}
	
	private void processRemote(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		inSketches=tool.loadSketches_MT(defaultParams, in);
		final int numLoaded=(inSketches.size());
		if(numLoaded>1 && defaultParams.mode==PER_FILE){
			inSketches.sort(SketchIdComparator.comparator);//Otherwise they come out disordered
		}
		
		t.stop();
		if(!silent){outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" in "+t);}
		assert(numLoaded<=MAX_ALLOWED_SKETCHES) : "\nSendSketch is configured to send at most "+MAX_ALLOWED_SKETCHES+" to prevent overwhelming the server.\n"
				+ "If you need to compare more than that, please use CompareSketch locally instead.\n"
				+ "References can be downloaded at http://portal.nersc.gov/dna/microbial/assembly/bushnell/\n";
		t.start();
//		outstream.println(inSketches.get(0));
		
		if(numLoaded>4000){
			SEND_BUFFER_MAX_BYTES*=4;
			SEND_BUFFER_MAX_SKETCHES*=4;
		}else if(numLoaded>1000){
			SEND_BUFFER_MAX_BYTES*=2;
			SEND_BUFFER_MAX_SKETCHES*=2;
		}
//		SEND_BUFFER_MAX_SKETCHES=1;//This is for testing json array format.
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		if(defaultParams.format==DisplayParams.FORMAT_QUERY_REF_ANI || defaultParams.format==DisplayParams.FORMAT_CONSTELLATION){tsw.println(defaultParams.header());}
		
		ByteBuilder bb=new ByteBuilder();
		
		int cntr=0;
		int chunks=0;
		boolean firstBlock=true; //Set to false after printing has started
		final Sketch lastSketch=Tools.getLast(inSketches); //Last input sketch to process
		for(Sketch sk : inSketches){
			
			if(sk.taxID<1 || sk.taxID>=minFakeID || outTaxID>0){sk.taxID=outTaxID;}
			if(defaultParams.printSSU()){sk.loadSSU();}
			
			if(outSpid>0){sk.spid=outSpid;}
			if(outImgID>0){sk.imgID=outImgID;}
			if(outTaxName!=null){sk.setTaxName(outTaxName);}
			if(outFname!=null){sk.setFname(outFname);}
			if(outName0!=null){sk.setName0(outName0);}
			sk.setMeta(outMeta);
			
			if(bb.length==0){
				bb.append(defaultParams.toString(chunks));
				chunks++;
			}
			sk.toBytes(bb);
			cntr++;
			if(cntr>=SEND_BUFFER_MAX_SKETCHES || bb.length>SEND_BUFFER_MAX_BYTES || sk==lastSketch){ //Don't allow too much data in a single transaction
				if(verbose){outstream.println("Sending:\n"+bb);}
//				outstream.println(cntr+", "+bb.length);
				
				byte[] message=bb.toBytes();
				bb.clear();
				try {
//					outstream.println("Sending to "+address+"\n"+message+"\n"); //123
					StringNum result=ServerTools.sendAndReceive(message, address);
					if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
						System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
						KillSwitch.kill();
					}
					errorState|=checkForError(result.s);
					if(!defaultParams.json()){
						tsw.print(result.s);
					}else{
						String s=result.s;
						
						{//Fixes rare case of multiple sketches broken into single sketches per block
							if(!firstBlock || sk!=lastSketch){
								if(s.charAt(0)!='['){s="["+s;}
								if(!s.endsWith("]")){;s=s+"]";}
							}
						}
						
						{//Case where there are more than one block of sketches, breaking JSON array format
							if(!firstBlock){
								s=","+s.substring(1);
							}
							if(sk!=lastSketch){
								s=s.substring(0, s.length()-1);
							}
						}
						
						tsw.print(s);
						firstBlock=false;
					}
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				cntr=0;
			}
		}
		
		if(bb.length>0){//This should no longer trigger unless there are no sketches since it's handled above.
			if(verbose){outstream.println("Sending:\n"+bb);}
			byte[] message=bb.toBytes();
			bb.clear();
			try {
				StringNum result=ServerTools.sendAndReceive(message, address);
				if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
					System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
					KillSwitch.kill();
				}
				errorState|=checkForError(result.s);
				tsw.print(result.s);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(!silent){tsw.println();}
		tsw.poison();
		
//		outstream.println("sending "+bb.toString());
		
		if(outSketch!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(outSketch, overwrite, append, true, FileFormat.SKETCH);
			bsw.start();
			for(Sketch sk : inSketches){
				sk.toBytes(bb);
				bsw.print(bb);
				bb.clear();
			}
			bsw.poisonAndWait();
			errorState|=bsw.errorState;
		}
		
		tsw.waitForFinish();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		if(!silent){outstream.println("Total Time: \t"+ttotal);}
	}
	

	public static ArrayList<JsonObject> sendSketches(ArrayList<Sketch> inSketches, String address, DisplayParams params){
		address=toAddress(address);
		final int numLoaded=(inSketches.size());
		assert(numLoaded<=MAX_ALLOWED_SKETCHES) : "\nSendSketch is configured to send at most "+MAX_ALLOWED_SKETCHES+" to prevent overwhelming the server.\n"
				+ "If you need to compare more than that, please use CompareSketch locally instead.\n"
				+ "References can be downloaded at http://portal.nersc.gov/dna/microbial/assembly/bushnell/\n";
		
		int SEND_BUFFER_MAX_BYTES=8000000;
		int SEND_BUFFER_MAX_SKETCHES=400;
		boolean errorState=false;
		
		ByteBuilder bb=new ByteBuilder();
		ArrayList<JsonObject> output=new ArrayList<JsonObject>(numLoaded);
		JsonParser jp=new JsonParser();
		
		int sketchesThisChunk=0;
		int chunkNum=0;
		boolean firstBlock=true; //Set to false after printing has started
		final Sketch lastSketch=Tools.getLast(inSketches); //Last input sketch to process
		for(Sketch sk : inSketches){
			if(bb.length==0){
				bb.append(params.toString(chunkNum));
				chunkNum++;
			}
			sk.toBytes(bb);//This does not handle null sketches
			sketchesThisChunk++;
			if(sketchesThisChunk>=SEND_BUFFER_MAX_SKETCHES || bb.length>SEND_BUFFER_MAX_BYTES || sk==lastSketch){ //Don't allow too much data in a single transaction
				if(verbose){System.err.println("Sending:\n"+bb);}
//				outstream.println(cntr+", "+bb.length);
				
				byte[] message=bb.toBytes();
				bb.clear();
				Object[] array=null;
				for(int i=0; i<10 && (array==null || array.length<sketchesThisChunk); i++) {
					String result=null;
					try {
						StringNum sn=ServerTools.sendAndReceive(message, address);
						if(sn.n<200 || sn.n>299){
							if(ServerTools.suppressErrors || i<9) {
								System.err.println("Unexpected server return code "+sn.n+"; retrying...");
								Tools.sleep(500*i);
								continue;
							}
							System.err.println("ERROR: Server returned code "+sn.n+" and this message:\n"+sn.s);
							KillSwitch.kill();
						}
						result=sn.s;
						errorState|=checkForError(result);
//						System.err.println("\nChunk "+chunkNum+":\n"+result);
						
						if(sketchesThisChunk>1) {array=jp.parseJsonArray(result);}
						else {array=new Object[] {jp.parseJsonObject(result)};}
						
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				assert(array.length==sketchesThisChunk);
				for(Object o : array) {
					JsonObject jo=(JsonObject)o;
					output.add(jo);
				}
				sketchesThisChunk=0;
			}
		}
		
		return output;
	}
	
	@Deprecated
	public static String sendSketch(Sketch sk, String address, int format, int chunkNum){
		DisplayParams params=defaultParams;
		if(format>=0){
			params=new DisplayParams();
			params.format=format;
		}
		return sendSketch(sk, address, params, chunkNum);
	}
	
	/** For programmatic use */
	public static String sendSketch(Sketch sk, String address, 
			DisplayParams params, int chunkNum){
		address=toAddress(address);
		ByteBuilder bb=new ByteBuilder();
		
		if(bb.length==0){bb.append(params.toString(chunkNum));}
		sk.toBytes(bb);

		byte[] message=bb.toBytes();
		try {
//			System.err.println("Sending to "+address+"\n"+new String(message)+"\n"); //123
			StringNum result=ServerTools.sendAndReceive(message, address);
			if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
				System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
				KillSwitch.kill();
			}
			return result.s;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	private static boolean checkForError(String s){
		if(s==null){return false;}
		return s.contains("java.io.IOException: Server returned HTTP response code:");
	}
	
	private void processLocal(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		
		final String message=defaultParams.toString(0);
		for(String fname : in){
			String address2=address+"/file/"+new File(fname).getAbsolutePath();
			
			if(verbose){outstream.println("Sending:\n"+message+"\nto "+address2);}
			try {
//				outstream.println("Sending to "+address2+"\n"+message+"\n"); //123
				StringNum result=ServerTools.sendAndReceive(message.getBytes(), address2);
				if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
					System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
					KillSwitch.kill();
				}
				tsw.print(result.s);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
		}
		tsw.println();
		
//		outstream.println("sending "+bb.toString());
		
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	private void processRefMode(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		
		final String message=defaultParams.toString(0);
		{
			String address2=address+"/ref/"+refNames;
			
			if(verbose){outstream.println("Sending:\n"+message+"\nto "+address2);}
			try {
//				outstream.println("Sending to "+address2+"\n"+message+"\n"); //123
				StringNum result=ServerTools.sendAndReceive(message.getBytes(), address2);
				if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
					System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
					KillSwitch.kill();
				}
				tsw.print(result.s);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
		}	
		tsw.println();
		
//		outstream.println("sending "+bb.toString());
		
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean addFiles(String a, Collection<String> list){
		int initial=list.size();
		if(a==null){return false;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){
				list.add(s);
			}
		}
		return list.size()>initial;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	private String outSketch=null;
	
	private String taxTreeFile=null;
	
	private final SketchTool tool;
	
	private ArrayList<Sketch> inSketches;

	private String address=null;
	private boolean local=false;
	/** List of reference names or TaxIDs to use as queries */
	private String refNames=null;
	
	/*Override metadata */
	private String outTaxName=null;
	private String outFname=null;
	private String outName0=null;
	private int outTaxID=-1;
	private long outSpid=-1;
	private long outImgID=-1;
	private ArrayList<String> outMeta=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
	private boolean silent=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final String toAddress(String b){
		String address=b;
		if(b==null){
			address=refseqAddress();//default
		}else if(b.equalsIgnoreCase("nt")){
			address=ntAddress();
		}else if(b.equalsIgnoreCase("refseq")){
			address=refseqAddress();
		}else if(b.equalsIgnoreCase("silva") || b.equalsIgnoreCase("ribo")){
			address=silvaAddress();
		}else if(b.equalsIgnoreCase("img")){
			address=imgAddress();
		}else if(b.equalsIgnoreCase("refseqprot") || b.equalsIgnoreCase("prokprot") 
				|| b.equalsIgnoreCase("protein") || b.equalsIgnoreCase("protien") || b.equalsIgnoreCase("prot")){
			address=prokProtAddress();
		}else if(b.equalsIgnoreCase("refseqmito") || b.equalsIgnoreCase("mito")){
			address=mitoAddress();
		}else if(b.equalsIgnoreCase("refseqfungi") || b.equalsIgnoreCase("fungi")){
			address=fungiAddress();
		}
		return address;
	}
	
	public int SEND_BUFFER_MAX_BYTES=8000000;
	public int SEND_BUFFER_MAX_SKETCHES=400;
	private static final int MAX_ALLOWED_SKETCHES=100000;
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
	private static String ntAddress(){return Shared.ntSketchServer()+"sketch";}
	private static String refseqAddress(){return Shared.refseqSketchServer()+"sketch";}
	private static String silvaAddress(){return Shared.riboSketchServer()+"sketch";}
	private static String imgAddress(){return null;}//Shared.SketchServer()+"sketch";
	private static String nrAddress(){return null;}//Shared.SketchServer()+"sketch";
	private static String prokProtAddress(){return Shared.proteinSketchServer()+"sketch";}
	private static String mitoAddress(){return null;}//Shared.SketchServer()+"sketch";
	private static String fungiAddress(){return null;}//Shared.SketchServer()+"sketch";
	
}
