package driver;

import java.io.PrintStream;
import java.util.Collection;
import java.util.LinkedHashSet;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;

/**
 * Filters text lines by exact match or substring.
 * @author Brian Bushnell
 * @date Jul 6, 2015
 *
 */
public class CountSharedLines {

	public static void main(String[] args){
		Timer t=new Timer();
		CountSharedLines x=new CountSharedLines(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public CountSharedLines(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						in1.add(s);
					}
				}
			}else if(a.equals("names") || a.equals("in2")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						in2.add(s);
					}
				}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("lines") || a.equals("maxlines")){
				maxLines=Parse.parseKMG(b);
			}else if(a.equals("substrings") || a.equals("substring")){
				if(b==null){b="t";}
				if(b.equals("header")){
					lineSubstringOfName=true;
				}else if(b.equals("name")){
					nameSubstringOfLine=true;
				}else{
					nameSubstringOfLine=lineSubstringOfName=Parse.parseBoolean(b);
				}
			}else if(a.equals("prefix") || a.equals("prefixmode")){
				prefixMode=Parse.parseBoolean(b);
			}else if(a.equals("replace")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				assert(split2.length==2);
				replace1=split2[0];
				replace2=split2[1];
			}else if(a.equals("casesensitive") || a.equals("case")){
				ignoreCase=!Parse.parseBoolean(b);
			}else if(a.equals("include") || a.equals("retain")){
				exclude=!Parse.parseBoolean(b);
			}else if(a.equals("exclude") || a.equals("remove")){
				exclude=Parse.parseBoolean(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

//		{
//			String[] x=in1.toArray(new String[in1.size()]);
//			in1.clear();
//			for(String s : x){
//				Tools.addNames(s, in1);
//			}
//			x=in2.toArray(new String[in2.size()]);
//			in2.clear();
//			for(String s : x){
//				Tools.addNames(s, in2);
//			}
//		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		if(in1==null || in2==null){throw new RuntimeException("Error - at least one input file is required from each set.");}
	}
	
	final static String getOutputName(String fname){
		fname=fname.replaceAll("\\\\", "/");
		if(!fname.contains("/")){fname="./"+fname;}
		int idx=fname.lastIndexOf('/');
		final String out=fname.substring(0, idx+1)+"out_"+fname.substring(idx+1);
		return out;
	}
	
	void process(Timer t){

		for(String fname : in1){
			processInner(fname, getOutputName(fname), in2);
		}
		for(String fname : in2){
			processInner(fname, getOutputName(fname), in1);
		}
		
		t.stop();
		
		outstream.println("\nTime:               "+t);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
			
	}
	
	LinkedHashSet<String> getContents(String fname){
		final FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, true);
		final LinkedHashSet<String> set=new LinkedHashSet<String>();
		final TextFile tf=new TextFile(ff);
		
		for(String line0=tf.readLine(true); line0!=null; line0=tf.readLine(true)){
			String line=(ignoreCase ? line0.toLowerCase() : line0);
			if(replace1!=null){line=line.replace(replace1, replace2);}
			if(prefixMode){
				for(int x=1; x<line.length(); x++){
					char c=line.charAt(x-1);
					char next=line.charAt(x);
					if(Character.isWhitespace(c)){
						line=line.substring(0, x).trim();
						break;
					}
				}
			}
			set.add(line);
		}
		errorState|=tf.close();
		return set;
	}
	
	void processInner(String fnameIn, String fnameOut, Collection<String> list){
		
		final LinkedHashSet<String> set1=getContents(fnameIn);
		final FileFormat ffout=FileFormat.testOutput(fnameOut, FileFormat.TEXT, null, true, overwrite, append, false);
		
		final TextStreamWriter tsw;
		if(ffout!=null){
			tsw=new TextStreamWriter(ffout);
			tsw.start();
		}else{tsw=null;}
		
		for(String fname2 : list){
			long shared=0;
			final LinkedHashSet<String> set2=getContents(fname2);
			for(String s : set1){
				if(set2.contains(s)){
					shared++;
				}
			}
			if(tsw!=null){
				tsw.print(ReadWrite.stripToCore(fname2)+"\t"+shared+"\n");
			}
		}

		if(tsw!=null){
			errorState|=tsw.poisonAndWait();
		}
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/

	private LinkedHashSet<String> in1=new LinkedHashSet<String>();
	private LinkedHashSet<String> in2=new LinkedHashSet<String>();
	
	/*--------------------------------------------------------------*/
	
	private boolean exclude=true;
	private boolean nameSubstringOfLine=false;
	private boolean lineSubstringOfName=false;
	private boolean ignoreCase=true;
	private boolean prefixMode=false;
	private long maxLines=-1;

	private String replace1=null;
	private String replace2=null;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean useSharedHeader=false;
	
}
