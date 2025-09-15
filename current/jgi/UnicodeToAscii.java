package jgi;

import java.io.File;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Apr 21, 2015
 *
 */
public class UnicodeToAscii {
	
	public static void main(String[] args){

		
		Timer t=new Timer();
		UnicodeToAscii x=new UnicodeToAscii(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public UnicodeToAscii(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			out2=parser.out2;
			
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2) || !ReadStats.testFiles(false)){
			throw new RuntimeException("Duplicate filenames are not allowed.");
		}
	}
	
	private void process(Timer t){

		if(in1!=null && out1!=null){process(in1, out1);}
		if(in2!=null && out2!=null){process(in2, out2);}
		
	}
		
	private void process(String infile, String outfile){
		TextFile tf=new TextFile(infile, true);
		TextStreamWriter tsw=new TextStreamWriter(outfile, overwrite, append, true);
		tsw.start();
		for(String line=tf.readLine(false); line!=null; line=tf.readLine(false)){
			String line2=line;
			try {
				line2=new String(line.getBytes(), "UTF-8");
			} catch (UnsupportedEncodingException e) {
				try {
					line2=new String(line.getBytes(), "UTF-16");
				} catch (UnsupportedEncodingException e1) {}
			}
			tsw.println(Tools.fixHeader(line2, false, true));
//			tsw.println(Normalizer.normalize(line, Normalizer.Form.NFD));
		}
		tf.close();
		tsw.poisonAndWait();
	}
	
	private PrintStream outstream=System.err;

	private String in1, in2;
	private String out1, out2;
	@SuppressWarnings("unused")
	private boolean verbose=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
