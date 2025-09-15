package tax;

import java.io.File;
import java.io.PrintStream;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date April 4, 2017
 *
 */
public class ShrinkAccession {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ShrinkAccession x=new ShrinkAccession(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ShrinkAccession(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		if(Data.PIGZ()){
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("gi")){
				KEEP_GI_NUMBERS=Parse.parseBoolean(b);
			}else if(a.equals("outgi") || a.equals("giout") || a.equals("gi")){
				giOut=b;
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;

			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutGi=FileFormat.testOutput(giOut, FileFormat.TXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		
	}
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		long linesProcessed=0;
		long charsProcessed=0;
		long badLines=0;
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder(10000);
		int columns=4;
		while(line!=null){
			if(Tools.startsWith(line, "accession\t")){
				bb.append(line);
				bb.nl();
			}else if(Tools.startsWith(line, "accession.version\ttaxid")){
				columns=2;
				bb.append("accession\t\ttaxid\t");//dummy header
				bb.nl();
			}else{
				charsProcessed+=line.length+1;
				linesProcessed++;
				
				final int tid=(columns==4 ? AccessionToTaxid.parseLineToTaxid(line, (byte)'\t') : 
					AccessionToTaxid.parseLineToTaxid_2col(line, (byte)'\t'));
				if(tid<1){
					badLines++;
				}else{
					int i=0;
					
					while(i<line.length){//Accession
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(columns==4){
						while(i<line.length){//Accession with decimal
							byte b=line[i];
							//						bb.append(b);
							i++;
							if(b=='\t'){break;}
						}
					}
					bb.append('\t');
					
					while(i<line.length){//Taxid
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(KEEP_GI_NUMBERS){
						if(line.length>i && Tools.isDigit(line[i])){//GI number or "na"
							while(i<line.length){
								byte b=line[i];
								bb.append(b);
								i++;
//								if(b=='\t'){break;}
							}
						}
					}
					bb.nl();
				}
				
//				String[] split=new String(line).split("\t");
//				bb.append(split[0]);
//				bb.tab();
//				bb.tab();
//				bb.append(split[2]);
//				bb.tab();
//				bb.nl();
			}
			if(bb.length()>8000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		outstream.println("Discarded "+badLines+" lines.\n");
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	private String in=null;
	private String out=null;
	private String giOut=null;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin;
	private final FileFormat ffout;
	private final FileFormat ffoutGi;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean KEEP_GI_NUMBERS=true;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
