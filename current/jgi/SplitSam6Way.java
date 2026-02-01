package jgi;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SamLine;

/**
 * Splits SAM files into six separate output files based on read characteristics.
 * Separates reads by pair number (R1/R2), strand orientation (plus/minus),
 * and mapping status (mapped/unmapped). Processes header lines and alignment
 * records to organize data for downstream analysis.
 *
 * @author Brian Bushnell
 * @date Jun 15, 2017
 */
public class SplitSam6Way {
	
	/** Program entry point: parses args, constructs SplitSam6Way, and closes redirected streams.
	 * Args: input SAM, then six outputs (r1plus, r1minus, r1unmapped, r2plus, r2minus, r2unmapped); use \"null\" to skip. */
	public static void main(String[] args){
		SplitSam6Way x=new SplitSam6Way(args);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Prints usage syntax and expected arguments to stderr. */
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx128m -cp <path> jgi.SplitSam6Way <input> <r1plus> <r1minus> <r1unmapped> <r2plus> <r2minus> <r2unmapped>");
		outstream.println("If you do not want one of the output files, use the word 'null'.\n");
	}
	
	/**
	 * Parses command-line arguments, sets up output writers, and streams SAM records into six buckets based on pair/strand/mapping.
	 * Supports optional maxReads argument.
	 * @param args [input, r1plus, r1minus, r1unmapped, r2plus, r2minus, r2unmapped, (optional) maxReads]
	 */
	public SplitSam6Way(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		if(args==null || args.length<7){
			printOptions();
			System.exit(1);
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		Timer t=new Timer();
		long reads=0, bases=0;
		long r1preads=0, r1mreads=0, r1ureads=0;
		long r2preads=0, r2mreads=0, r2ureads=0;
		
		String fin=args[0];
		String fr1plus=args[1];
		String fr1minus=args[2];
		String fr1unmapped=args[3];
		String fr2plus=args[4];
		String fr2minus=args[5];
		String fr2unmapped=args[6];
		
		long maxReads=Long.MAX_VALUE;
		if(args.length>7){
			maxReads=Parse.parseKMG(args[7]);
		}
		
//		ByteFile.FORCE_MODE_BF1=true;
		ByteFile tf=ByteFile.makeByteFile(fin, true);
		ByteStreamWriter r1plus=("null".equalsIgnoreCase(fr1plus) ? null : new ByteStreamWriter(fr1plus, true, false, true, FileFormat.SAM));
		ByteStreamWriter r1minus=("null".equalsIgnoreCase(fr1minus) ? null : new ByteStreamWriter(fr1minus, true, false, true, FileFormat.SAM));
		ByteStreamWriter r1unmapped=("null".equalsIgnoreCase(fr1unmapped) ? null : new ByteStreamWriter(fr1unmapped, true, false, true, FileFormat.SAM));
		ByteStreamWriter r2plus=("null".equalsIgnoreCase(fr2plus) ? null : new ByteStreamWriter(fr2plus, true, false, true, FileFormat.SAM));
		ByteStreamWriter r2minus=("null".equalsIgnoreCase(fr2minus) ? null : new ByteStreamWriter(fr2minus, true, false, true, FileFormat.SAM));
		ByteStreamWriter r2unmapped=("null".equalsIgnoreCase(fr2unmapped) ? null : new ByteStreamWriter(fr2unmapped, true, false, true, FileFormat.SAM));

		if(r1plus!=null){r1plus.start();}
		if(r1minus!=null){r1minus.start();}
		if(r1unmapped!=null){r1unmapped.start();}
		if(r2plus!=null){r2plus.start();}
		if(r2minus!=null){r2minus.start();}
		if(r2unmapped!=null){r2unmapped.start();}
		
		final LineParser1 lp=new LineParser1('\t');
		for(byte[] line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line[0]=='@'){
				if(r1plus!=null){r1plus.println(line);}
				if(r1minus!=null){r1minus.println(line);}
				if(r1unmapped!=null){r1unmapped.println(line);}
				if(r2plus!=null){r2plus.println(line);}
				if(r2minus!=null){r2minus.println(line);}
				if(r2minus!=null){r2minus.println(line);}
			}else{
				if(reads>=maxReads){break;}
				
				SamLine sl=new SamLine(lp.set(line));
				reads++;
				bases+=sl.seq.length;
				
				if(sl.pairnum()==0){
					if(sl.mapped()){
						if(sl.strand()==Shared.PLUS){
							if(r1plus!=null){r1plus.println(line);}
							r1preads++;
						}else{
							if(r1minus!=null){r1minus.println(line);}
							r1mreads++;
						}
					}else{
						if(r1unmapped!=null){r1unmapped.println(line);}
						r1ureads++;
					}
				}else{
					if(sl.mapped()){
						if(sl.strand()==Shared.PLUS){
							if(r2plus!=null){r2plus.println(line);}
							r2preads++;
						}else{
							if(r2minus!=null){r2minus.println(line);}
							r2mreads++;
						}
					}else{
						if(r2unmapped!=null){r2unmapped.println(line);}
						r2ureads++;
					}
				}
			}
		}
		
		tf.close();
		
		if(r1plus!=null){r1plus.poisonAndWait();}
		if(r1minus!=null){r1minus.poisonAndWait();}
		if(r1unmapped!=null){r1unmapped.poisonAndWait();}
		if(r2plus!=null){r2plus.poisonAndWait();}
		if(r2minus!=null){r2minus.poisonAndWait();}
		if(r2unmapped!=null){r2unmapped.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
		outstream.println("R1 Plus Reads:      "+r1preads);
		outstream.println("R1 Minus Reads:     "+r1mreads);
		outstream.println("R1 Unmapped Reads:  "+r1ureads);
		outstream.println("R1 Plus Reads:      "+r2preads);
		outstream.println("R1 Minus Reads:     "+r2mreads);
		outstream.println("R1 Unmapped Reads:  "+r2ureads);
	}
	
	private PrintStream outstream=System.err;
	
}
