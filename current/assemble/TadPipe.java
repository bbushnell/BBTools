package assemble;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import clump.Clumpify;
import dna.Data;
import fileIO.ReadWrite;
import jgi.BBDuk;
import jgi.BBMerge;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import stream.FASTQ;

/**
 * Merges and error-corrects reads, then assembles.
 * @author Brian Bushnell
 * @date August 16, 2016
 *
 */
public class TadPipe {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		TadPipe x=new TadPipe(args);
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public TadPipe(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set some shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(1+Shared.threads()/2);
		FASTQ.DETECT_QUALITY_OUT=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("phix") || a.equals("contam")){
				assert(false) : "Unsupported";
//				phix=b;
			}else if(a.equals("delete")){
				deleteTemp=Parse.parseBoolean(b);
			}else if(a.equals("gz")){
				gz=Parse.parseBoolean(b);
			}else if(a.equals("tmpdir") || a.equals("tempdir") || a.equals("temp")){
				tempdir=(b+"/").replaceAll("\\\\", "/").replaceAll("//", "/");
			}else if(a.startsWith("merge_")){
				mergeArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("ecco_")){
				eccoArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("ecc_") || a.startsWith("correct_")){
				eccArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("extend_") || a.startsWith("extend1_")){
				extendArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("extend2_")){
				extend2Args.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("filter_")){
				assert(false) : "Unsupported";
//				filterArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("clump_") || a.startsWith("clumpify_")){
				clumpifyArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("trim_")){
				trimArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("assemble_")){
				assembleArgs.add(arg.substring(arg.indexOf('_')+1));
			}else{
				assembleArgs.add(arg);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){
		
		File tmpfile=(tempdir==null ? null : new File(tempdir));
		if(tmpfile!=null && !tmpfile.exists()){tmpfile.mkdirs();}
		
//		String outFilter=null;
//
//		String outAdapter=null;

		String outTrimmed=null;
		String outEcco=null;
		String outClumped=null;

		String outMerged=null;
		String outUnmerged=null;

		String outMergedEcc=null;
		String outUnmergedEcc=null;

		String outMergedExtended=null;
		String outUnmergedExtended=null;

		String outMergedExtended2=null;
		String outUnmergedExtended2=null;
		
		String outMultiK=null;
		
		String ext=(gz ? ".fq.gz" : ".fq");
		
		try {
//			outFilter=File.createTempFile("filter_", ext, tmpfile).toString();
//			outAdapter=File.createTempFile("adapters_", ".fa", tmpfile).toString();
			outTrimmed=File.createTempFile("trimmed_", ext, tmpfile).toString();
			outEcco=File.createTempFile("ecco_", ext, tmpfile).toString();
			outClumped=File.createTempFile("clumped_", ext, tmpfile).toString();
			outMerged=File.createTempFile("merged_", ext, tmpfile).toString();
			outUnmerged=File.createTempFile("unmerged_", ext, tmpfile).toString();
			outMergedEcc=File.createTempFile("m_ecc_", ext, tmpfile).toString();
			outUnmergedEcc=File.createTempFile("u_ecc_", ext, tmpfile).toString();
			outMergedExtended=File.createTempFile("m_extended_", ext, tmpfile).toString();
			outUnmergedExtended=File.createTempFile("u_extended_", ext, tmpfile).toString();
			if(extend2){
				outMergedExtended2=File.createTempFile("m_extended2_", ext, tmpfile).toString();
				outUnmergedExtended2=File.createTempFile("u_extended2_", ext, tmpfile).toString();
			}
			outMultiK=File.createTempFile("multik_%_", ".fa", tmpfile).toString();
			delete(outMultiK);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		{
//			filterArgs.add("in="+in1);
//			if(in2!=null){filterArgs.add("in2="+in2);}
//
//			filterArgs.add("out="+outFilter);
//			filterArgs.add("ftm="+5);
//			filterArgs.add("ow");
//
//			if(phix!=null){
//				filterArgs.add("ref="+phix);
//				filterArgs.add("k="+31);
//			}
//
//			Collections.reverse(filterArgs);
//
//			BBDukF.main(filterArgs.toArray(new String[0]));
//		}
//
//		{
//			adapterArgs.add("in="+outFilter);
//			adapterArgs.add("outa="+outAdapter);
//
//			adapterArgs.add("ow");
//
//			Collections.reverse(adapterArgs);
//
//			BBMerge.main(adapterArgs.toArray(new String[0]));
//		}
		
		String adapterPath=Data.findPath("?adapters.fa", true);
		if(adapterPath!=null){
//			trimArgs.add("in="+outFilter);
//			trimArgs.add("ref="+outAdapter);
			trimArgs.add("in="+in1);
			trimArgs.add("ref="+adapterPath);
			trimArgs.add("out="+outTrimmed);
			
//			ArrayList<byte[]> alb=Tools.toAdapterList(outAdapter, 23);
//			boolean doTrim=false;
//			if(alb!=null){
//				for(byte[] array : alb){
//					if(array.length>=11){
//						doTrim=true;
//					}
//				}
//			}

//			if(doTrim){
			if(true){
				trimArgs.add("k=23");
				trimArgs.add("mink=11");
				trimArgs.add("hdist=1");
				trimArgs.add("ktrim=r");
			}

			trimArgs.add("qtrim=r");
			trimArgs.add("trimq=10");
			
			trimArgs.add("tbo");
			trimArgs.add("tpe");
			trimArgs.add("ow");
			trimArgs.add("minlen=62");
			
			Collections.reverse(trimArgs);
			
			BBDuk.main(trimArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outFilter);}
		}else{
			outTrimmed=in1;
		}
		
		{
			eccoArgs.add("in="+outTrimmed);
			
			eccoArgs.add("out="+outEcco);

			eccoArgs.add("strict");
			eccoArgs.add("ecco");
			eccoArgs.add("mix");
			eccoArgs.add("adapters=default");
			eccoArgs.add("ow");
			
			Collections.reverse(eccoArgs);
			
			BBMerge.main(eccoArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outTrimmed, outAdapter);}
			if(deleteTemp && outTrimmed!=in1){
				delete(outTrimmed);
//				assert(false) : outTrimmed+", "+in1;
			}
//			assert(false) : outTrimmed+", "+in1+", "+outTrimmed.equals(in1)+", "+(outTrimmed==in1);
		}
		
		{
			clumpifyArgs.add("in="+outEcco);
			
			clumpifyArgs.add("out="+outClumped);
			
			clumpifyArgs.add("ow");
			clumpifyArgs.add("shortname");
			clumpifyArgs.add("ecc");
			clumpifyArgs.add("passes=8");
			clumpifyArgs.add("unpair");
			clumpifyArgs.add("repair");
			
			Collections.reverse(clumpifyArgs);
			
			Clumpify.main(clumpifyArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outTrimmed, outAdapter);}
			if(deleteTemp){delete(outEcco);}
		}
		
		{
			mergeArgs.add("in="+outClumped);
			
			mergeArgs.add("out="+"j.fq");
			mergeArgs.add("outu="+"q.fq");

			mergeArgs.add("k=75");
			mergeArgs.add("extend2=120");
			mergeArgs.add("rem");
			mergeArgs.add("ecct");
//			mergeArgs.add("loose");
//			mergeArgs.add("adapters="+outAdapter);
			mergeArgs.add("adapters=default");
			mergeArgs.add("ow");
			mergeArgs.add("ordered");
			
			Collections.reverse(mergeArgs);
			
			BBMerge.main(mergeArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outTrimmed, outAdapter);}
			if(deleteTemp){delete(outClumped);}
		}
		
		{
			eccArgs.add("in="+outMerged+","+outUnmerged);
			eccArgs.add("out="+outMergedEcc+","+outUnmergedEcc);

			eccArgs.add("k=50");
			eccArgs.add("ecc");
			eccArgs.add("tossjunk");
			eccArgs.add("deadzone=2");
			eccArgs.add("ow");
			eccArgs.add("ordered");
			
			Collections.reverse(eccArgs);
			
			Tadpole.main(eccArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMerged, outUnmerged);}
		}
		
		{
			extendArgs.add("in="+outMergedEcc+","+outUnmergedEcc);
			extendArgs.add("out="+outMergedExtended+","+outUnmergedExtended);

			extendArgs.add("k=81");
			extendArgs.add("mode=extend");
			extendArgs.add("el=100");
			extendArgs.add("er=100");
			extendArgs.add("ow");
			extendArgs.add("tossjunk");
			extendArgs.add("deadzone=0");
			extendArgs.add("ecc");
			extendArgs.add("ordered");
			
			Collections.reverse(extendArgs);
			
			Tadpole.main(extendArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMergedEcc, outUnmergedEcc);}
		}
		
		if(extend2){
			extend2Args.add("in="+outMergedExtended+","+outMergedExtended);
			extend2Args.add("out="+outMergedExtended2+","+outUnmergedExtended2);

			extend2Args.add("k=124");
			extend2Args.add("mode=extend");
			extend2Args.add("el=60");
			extend2Args.add("er=60");
			extend2Args.add("ow");
			extend2Args.add("tossjunk");
			extend2Args.add("deadzone=0");
			extend2Args.add("ecc");
			extend2Args.add("ordered");
			
			Collections.reverse(extend2Args);
			
			Tadpole.main(extend2Args.toArray(new String[0]));
			if(deleteTemp){delete(outMergedExtended, outMergedExtended);}
		}
		
		{
			if(extend2){
				assembleArgs.add("in="+outMergedExtended2+","+outUnmergedExtended2);
			}else{
				assembleArgs.add("in="+outMergedExtended+","+outUnmergedExtended);
			}
			assembleArgs.add("out="+outMultiK);
			assembleArgs.add("outfinal="+out);

			assembleArgs.add("k=210,250,290");
//			assembleArgs.add("quitearly");
			assembleArgs.add("expand");
			assembleArgs.add("bisect");
			assembleArgs.add("shave");
			assembleArgs.add("rinse");
			assembleArgs.add("pop");
			if(deleteTemp){assembleArgs.add("delete");}
			
			Collections.reverse(assembleArgs);
			
			int best=TadpoleWrapper.process(assembleArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMergedExtended, outUnmergedExtended);}
		}
		
//		outstream.println();
//		AssemblyStats2.main(new String[] {"in="+out});
		
//		{
//			assembleArgs.add("in="+outMergedExtended+","+outUnmergedExtended);
//			assembleArgs.add("out="+out);
//
//			assembleArgs.add("k=250");
//
//			Collections.reverse(assembleArgs);
//
//			Tadpole.main(assembleArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outMergedExtended, outUnmergedExtended);}
//		}
		
	}
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String...names){
		for(String s : names){
			if(s!=null){
				assert(!s.equals(in1)) : s;
				if(verbose){outstream.println("Trying to delete "+s);}
				File f=new File(s);
				if(f.exists()){
					f.delete();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
//	private ArrayList<String> filterArgs=new ArrayList<String>();
//
//	private ArrayList<String> adapterArgs=new ArrayList<String>();
	
	private ArrayList<String> trimArgs=new ArrayList<String>();
	
	private ArrayList<String> eccoArgs=new ArrayList<String>();
	
	private ArrayList<String> clumpifyArgs=new ArrayList<String>();
	
	private ArrayList<String> assembleArgs=new ArrayList<String>();
	
	private ArrayList<String> mergeArgs=new ArrayList<String>();
	
	private ArrayList<String> eccArgs=new ArrayList<String>();
	
	private ArrayList<String> extendArgs=new ArrayList<String>();
	
	private ArrayList<String> extend2Args=new ArrayList<String>();
	
	private String in1, in2;
	
	private String out="contigs.fa";
	
	private String tempdir=Shared.tmpdir();

	private boolean deleteTemp=true;
	private boolean gz=false;
	
	private boolean extend2=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
