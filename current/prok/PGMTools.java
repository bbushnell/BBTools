package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
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

/**
 * Static helpers for manipulating pgm files.
 * main() merges pgm files.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class PGMTools extends ProkObject {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Combines multiple pgm files into a single file */
	public static void main(String[] args){
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ProkObject.call18S=true;
		boolean overwrite=true;
		boolean allowDupes=false;
		String out=null;
		ArrayList<String> in=new ArrayList<String>(); 
		
		{
			Parser parser=new Parser();
			for(int i=0; i<args.length; i++){
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}
				
				if(a.equals("in")){
					assert(b!=null);
					Tools.addFiles(b, in);
				}else if(parseStatic(arg, a, b)){
					//do nothing
				}else if(a.equals("allowdupes") || a.equals("allowduplicates") || a.equals("dupes")){
					allowDupes=Parse.parseBoolean(b);
				}else if(a.equals("addcdsonly")){
					GeneModel.ADD_CDS_ONLY=Parse.parseBoolean(b);
				}else if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
					ReadWrite.verbose=verbose;
				}

				else if(parser.parse(arg, a, b)){
					//do nothing
				}else if(b==null && new File(arg.split("@")[0]).exists()){
					in.add(arg);
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
					//				throw new RuntimeException("Unknown parameter "+args[i]);
				}
			}
			overwrite=parser.overwrite;
			out=parser.out1;
		}
		
		checkFileExistence(in, out, overwrite, allowDupes);
		
		ArrayList<GeneModel> models=loadModels(in);
		GeneModel gm=mergeModels(models);
		boolean errorState=writeModel(gm, out, overwrite);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public static boolean parseStatic(String arg, String a, String b){
		if(a.equals("kinnercds")){
			int k=Integer.parseInt(b);
			GeneModel.setInnerK(k);
		}else if(a.equals("kstartcds")){
			int k=Integer.parseInt(b);
			GeneModel.setStartK(k);
		}else if(a.equals("kstopcds")){
			int k=Integer.parseInt(b);
			GeneModel.setStopK(k);
		}else if(a.equals("kinnerrna")){
			int k=Integer.parseInt(b);
			kInnerRNA=k;
		}else if(a.equals("kstartrna")){
			int k=Integer.parseInt(b);
			kStartRNA=k;
		}else if(a.equals("kstoprna")){
			int k=Integer.parseInt(b);
			kStopRNA=k;
		}else if(a.equals("startleftoffset")){
			int x=Integer.parseInt(b);
			GeneModel.setStartLeftOffset(x);
		}else if(a.equals("startrightoffset")){
			int x=Integer.parseInt(b);
			GeneModel.setStartRightOffset(x);
		}else if(a.equals("stopleftoffset")){
			int x=Integer.parseInt(b);
			GeneModel.setStopLeftOffset(x);
		}else if(a.equals("stoprightoffset")){
			int x=Integer.parseInt(b);
			GeneModel.setStopRightOffset(x);
		}else if(a.equalsIgnoreCase("callcdsonly") || a.equalsIgnoreCase("cdsonly")){
			callCDS=Parse.parseBoolean(b);
			calltRNA=call16S=call23S=call5S=call18S=!callCDS;
		}else if(a.equalsIgnoreCase("call16sonly") || a.equalsIgnoreCase("16sonly")){
			call16S=Parse.parseBoolean(b);
			calltRNA=call23S=call5S=call18S=callCDS=!call16S;
		}else if(a.equalsIgnoreCase("call23sonly") || a.equalsIgnoreCase("23sonly")){
			call23S=Parse.parseBoolean(b);
			calltRNA=call16S=call5S=call18S=callCDS=!call23S;
		}else if(a.equalsIgnoreCase("call5sonly") || a.equalsIgnoreCase("5sonly")){
			call5S=Parse.parseBoolean(b);
			calltRNA=call16S=call23S=call18S=callCDS=!call5S;
		}else if(a.equalsIgnoreCase("calltrnaonly") || a.equalsIgnoreCase("trnaonly")){
			calltRNA=Parse.parseBoolean(b);
			call16S=call23S=call5S=call18S=callCDS=!calltRNA;
		}else if(a.equalsIgnoreCase("call18sonly") || a.equalsIgnoreCase("18sonly")){
			call18S=Parse.parseBoolean(b);
			calltRNA=call16S=call23S=call5S=callCDS=!call18S;
		}
		
		else if(a.equalsIgnoreCase("callcds") || a.equalsIgnoreCase("cds")){
			callCDS=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("calltrna") || a.equalsIgnoreCase("trna")){
			calltRNA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("call16s") || a.equalsIgnoreCase("16s")){
			call16S=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("call18s") || a.equalsIgnoreCase("18s")){
			call18S=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("call23s") || a.equalsIgnoreCase("23s")){
			call23S=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("call5s") || a.equalsIgnoreCase("5s")){
			call5S=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("callrna") || a.equalsIgnoreCase("rna")){
			calltRNA=call16S=call18S=call5S=call23S=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("normalize")){
			normalize=Parse.parseBoolean(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public static ArrayList<GeneModel> loadModels(ArrayList<String> fnames){
		ArrayList<GeneModel> models=new ArrayList<GeneModel>(fnames.size());
		ArrayList<Double> mults=new ArrayList<Double>(fnames.size());
		for(String s : fnames){
			double mult=1;
			String fname=s;
			if(s.indexOf('@')>=0){
				String[] split=s.split("@");
				fname=split[0];
				mult=Double.parseDouble(split[1]);
			}
			GeneModel pgm=GeneModelParser.loadModel(fname);
			mults.add(mult);
			models.add(pgm);
		}
		if(normalize){
			long max=0;
			for(GeneModel gm : models){
				max=Tools.max(gm.basesProcessed, max);
			}
			for(GeneModel gm : models){
				if(max!=gm.basesProcessed){
					double mult=max/(double)(Tools.max(100, gm.basesProcessed));
					if(mult!=1) {
						gm.multiplyBy(mult);
						gm.stats5S.multiplyBy(0.01);
//						System.err.println("Multiplied by "+mult);
					}
				}
			}
		}
		for(int i=0; i<models.size(); i++){
			double mult=mults.get(i);
			GeneModel gm=models.get(i);
			if(mult!=1){
				gm.multiplyBy(mult);
//				System.err.println("Multiplied by "+mult);
			}
		}
		return models;
	}
	
	public static GeneModel mergeModels(ArrayList<GeneModel> models){
		if(models.size()==1){return models.get(0);}
		GeneModel pgmSum=new GeneModel(true);
		for(GeneModel pgm : models){
			pgmSum.add(pgm);
		}
		return pgmSum;
	}
	
	/** The first gene model is the dominant one */
	public static GeneModel mix(double cdsMult, double rRnaMult, double tRnaMult, boolean normalize, GeneModel...models) {
		assert(models.length>1);
		GeneModel pgm0=models[0];
		GeneModel pgm=new GeneModel(true);
		pgm.add(pgm0);
		long max=Tools.max(100, pgm0.basesProcessed);
		for(int i=1; i<models.length; i++) {
			GeneModel pgmI=models[i];
			GeneModel temp=new GeneModel(true);
			temp.add(pgmI);
			
			assert(temp.basesProcessed<max) : temp.basesProcessed+", "+max+"\nThis function is specifically for mixing a small model in with a large model.";
			double mult=max/(double)(Tools.max(100, temp.basesProcessed));
			temp.statsCDS.multiplyBy(mult*cdsMult);
			temp.stats5S.multiplyBy(mult*rRnaMult);
			temp.stats16S.multiplyBy(mult*rRnaMult);
			temp.stats18S.multiplyBy(mult*rRnaMult);
			temp.stats23S.multiplyBy(mult*rRnaMult);
			temp.statstRNA.multiplyBy(mult*tRnaMult);
			pgm.add(temp);
		}
		return pgm;
	}
	
	public static GeneModel loadAndMerge(ArrayList<String> in) {
		ArrayList<GeneModel> models=loadModels(in);
		return mergeModels(models);
	}
	
	public static boolean writeModel(GeneModel pgm, String out, boolean overwrite){
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, false, false);
		return writeModel(pgm, ffout);
	}
	
	public static boolean writeModel(GeneModel pgm, FileFormat ffout){
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);

		ByteBuilder bb=new ByteBuilder();
		pgm.appendTo(bb);

		boolean errorState=false;
		if(bsw!=null){
			bsw.addJob(bb);
			errorState|=bsw.poisonAndWait();
		}
		return errorState;
	}
	
	/** Ensure files can be read and written */
	private static void checkFileExistence(ArrayList<String> in, String out, boolean overwrite, boolean allowDupes){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		ArrayList<String> in2=new ArrayList<String>();
		for(String s : in){
			in2.add(s.split("@")[0]);
		}
		in=null;
		
		for(int i=0; i<in2.size(); i++){
			String s=in2.get(i);
			if(s.equalsIgnoreCase("auto") || s.equalsIgnoreCase("default")){
				in2.set(i, Data.findPath("?model.pgm"));
			}
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(in2);
		if(!Tools.testInputFiles(allowDupes, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!allowDupes){
			foo.add(out);
			if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
				throw new RuntimeException("\nSome file names were specified multiple times.\n");
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	/** Mix models equally */
	public static boolean normalize=false;
	
}
