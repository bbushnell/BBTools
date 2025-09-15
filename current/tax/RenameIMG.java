package tax;

import java.io.File;
import java.io.PrintStream;

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
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.IntHashSet;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class RenameIMG {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameIMG x=new RenameIMG(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RenameIMG(String[] args){
		
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

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("img")){
				imgFile=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}//TODO: why are these set to the same default?
		if("auto".equalsIgnoreCase(in1)){in1=TaxTree.defaultImgFile();}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FA, null, true, overwrite, append, false);
	}
	
	void copyFiles(ImgRecord[] array){
		if(useSet){set=new IntHashSet(10000);}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		for(ImgRecord ir : array){
			if(ir.taxID>0){set.add(ir.taxID);}
			else{unknownTaxid++;}
			FileFormat ffin=FileFormat.testInput(ir.path(), FileFormat.FA, null, true, true);
			process_inner(ffin, bsw, ir.imgID);
		}
		knownTaxid=set.size();
		set=null;
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	void process(Timer t){
		ImgRecord[] array=ImgRecord.toArray(in1, TaxTree.IMG_HQ);
		if(imgFile==null){
			TaxTree.loadIMG(array);
		}else{
			ImgRecord[] array2=ImgRecord.toArray(imgFile, TaxTree.IMG_HQ);
			TaxTree.loadIMG(array2);
		}
		
		copyFiles(array);
		
		t.stop();

		final int spaces=8;
		String fpstring=""+filesProcessed;
		String cpstring=Tools.padKMB(sequencesProcessed, spaces);
		String bapstring=Tools.padKMB(basesProcessed, spaces);
		String tpstring=""+knownTaxid;
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Files Processed:    "+fpstring);
		outstream.println("Contigs Processed:  "+cpstring);
		outstream.println("Bases Processed:    "+bapstring);
		if(useSet){outstream.println("TaxIDs Processed:   "+tpstring+" \t"+"("+unknownTaxid+" unknown)");}
		outstream.println(Tools.linesBytesProcessed(t.elapsed, linesProcessed, bytesProcessed, spaces));
		
		outstream.println();
		outstream.println("Valid Files:       \t"+filesValid);
		outstream.println("Invalid Files:     \t"+(filesProcessed-filesValid));
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	void process_inner(final FileFormat ffin, final ByteStreamWriter bsw, final long img){
		
		filesProcessed++;
		{
			File f=new File(ffin.name());
			if(!f.exists() || !f.canRead()){
				System.err.println("Can't find "+f);
				errorState=true;
				return;
			}
		}
		final int tid=TaxTree.imgToTaxid(img);
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;

				linesValid++;
				if(line[0]=='>'){
					sequencesProcessed++;
					bb.append('>');
					if(tid>=0){
						bb.append("tid|");
						bb.append(tid);
						bb.append('|');
					}
					bb.append("img|");
					bb.append(img);
					bb.append(' ');
					for(int i=1; i<line.length; i++){
						bb.append(line[i]);
					}
				}else{
					basesProcessed+=line.length;
					bb.append(line);
				}
				bb.nl();
				bsw.print(bb.toBytes());
				bb.clear();
			}
			line=bf.nextLine();
		}
		
		filesValid++;
		errorState|=bf.close();
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String imgFile=null;
	
	/*--------------------------------------------------------------*/
	
	private IntHashSet set=null;
	private int knownTaxid=0;
	private int unknownTaxid=0;
	private boolean useSet=true;
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;

	private long basesProcessed=0;
	private long sequencesProcessed=0;
	private long filesProcessed=0;
	private long filesValid=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
