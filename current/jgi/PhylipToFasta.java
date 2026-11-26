package jgi;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Converts Phylip format sequence files to FASTA format.
 * Reads the legacy interleaved layout handled here and outputs FASTA.
 * General PHYLIP name/layout support is incomplete; sequential input is unsupported.
 * Parses sequence names and concatenates sequence data across multiple lines.
 * Output honors overwrite and append settings; append takes precedence.
 *
 * @author Brian Bushnell
 * @date Oct 3, 2014
 */
public class PhylipToFasta {
	


	public static void main(String[] args){
		Timer t=new Timer();
		PhylipToFasta x=new PhylipToFasta(args);
		try{
			x.process(t);
		}finally{
			//JGI-162: close adopted diagnostics on every exit, preserving System streams.
			Shared.closeStream(x.outstream);
		}
	}
	
	public PhylipToFasta(String[] args){
		boolean success=false;
		try{
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
				//JGI-159: preserve complete values, but retain the old empty-value null behavior.
				final int equals=arg.indexOf('=');
				String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
				String b=(equals<0 || equals==arg.length()-1 ? null : arg.substring(equals+1));

				if(parser.parse(arg, a, b)){
					//do nothing
				}else if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
					ReadWrite.verbose=verbose;
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
				Parser.processQuality();

				in1=parser.in1;

				out1=parser.out1;
				//JGI-158: parsing these flags is not enough; the output descriptor needs them.
				overwrite=parser.overwrite;
				append=parser.append;
			}

			if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

			if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}

			if(!Tools.testOutputFiles(overwrite, append, false, out1)){
				outstream.println((out1==null)+", "+out1);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
			}

			ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);

			ffin1=FileFormat.testInput(in1, FileFormat.PHYLIP, ".phylip", true, true);
			success=true;
		}finally{
			//The caller cannot close our adopted diagnostics if construction does not return.
			if(!success){Shared.closeStream(outstream);}
		}
	}

	void process(Timer t){
		
		ArrayList<StringBuilder> data=new ArrayList<StringBuilder>();
		long bases=0;
		
		//JGI-160: close owned input and decoder resources without closing borrowed System.in.
		try(final InputStream input=ReadWrite.getManagedInputStream(ffin1.name(), true, ffin1.allowSubprocess());
			final BufferedReader reader=new BufferedReader(new InputStreamReader(input), 32768)){
			String s=nextLogicalLine(reader); //first line is some numbers

			//TODO: Probable bug JGI-161 - standard ten-character names and unindented
			//continuation blocks are not handled; count validation and gap/symbol policy remain open.
			for(s=nextLogicalLine(reader); s!=null; s=nextLogicalLine(reader)){
				if(s.startsWith("        ")){break;}
				StringBuilder sb=new StringBuilder();
				data.add(sb);
				sb.append('>');
				int pos=0;
				for(pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isWhitespace(c)){break;}
					sb.append(c);
				}
				sb.append('\n');
				while(pos<s.length() && Character.isWhitespace(s.charAt(pos))){pos++;}
				while(pos<s.length()){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
					pos++;
				}
			}

			final int mod=data.size();
			for(int i=0; s!=null; i++){
				StringBuilder sb=data.get(i%mod);
				for(int pos=0; pos<s.length(); pos++){
					char c=s.charAt(pos);
					if(Character.isLetter(c)){
						sb.append(c);
						bases++;
					}
				}
				s=nextLogicalLine(reader);
			}
		}catch(final IOException e){
			throw new RuntimeException("Error reading or closing "+ffin1.name(), e);
		}
		final long reads=data.size();
		
		if(ffout1!=null){
			TextStreamWriter tsw=new TextStreamWriter(ffout1);
			tsw.start();
			for(int i=0; i<data.size(); i++){
				StringBuilder sb=data.set(i, null);
				sb.append('\n');
				tsw.print(sb);
			}
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Returns the next line using TextFile's blank-line policy, without trimming its contents. */
	private static String nextLogicalLine(final BufferedReader reader) throws IOException{
		String line;
		while((line=reader.readLine())!=null){
			if(line.length()==0){continue;}
			if(Character.isWhitespace(line.charAt(0)) && Character.isWhitespace(line.charAt(line.length()-1))
					&& line.trim().length()==0){continue;}
			return line;
		}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;

	private String out1=null;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
