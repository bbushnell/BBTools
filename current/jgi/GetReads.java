package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;

import dna.Data;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaQualWriterST;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Extracts reads with requested zero-based numeric IDs.  IDs may be listed
 * individually, as inclusive ranges, or both, comma-delimited with {@code id=}.
 * If a requested root has a mate, the mate is emitted to {@code out2} when set,
 * otherwise immediately after the root in the same output stream.
 *
 * @author Brian Bushnell
 * @date Jul 10, 2013
 */
public class GetReads{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the extraction and closes the redirected status stream, if any. */
	public static void main(String[] args){
		GetReads x=new GetReads(args);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses command-line arguments, streams the input for the requested passes,
	 * and writes matching roots to the primary output, plus selected mates to the
	 * same stream by default or to {@code out2} when a mate output is requested.
	 */
	public GetReads(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Timer t=new Timer();

		Parser parser=new Parser();
		String in1=null;
		String in2=null;

		String qfin1=null;
		String qfin2=null;

		String out1=null;
		String out2=null;

		String qfout1=null;
		String qfout2=null;

		boolean errorState=false;
		long maxReads=-1;
		int passes=1;
		boolean testsize=false;
		boolean overwrite=false, append=false;
		float samplerate=1f;
		long sampleseed=1;



		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;

		HashSet<Long> table=new HashSet<Long>();

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("id") || a.equals("number")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] b2=b.split(",");
				for(String c : b2){
					final long x, y;
					if(c.indexOf('-')>=0){
						String[] c2=c.split("-");
						assert(c2.length==2) : c;
						x=Long.parseLong(c2[0]);
						y=Long.parseLong(c2[1]);
					}else{
						x=y=Long.parseLong(c);
					}
					for(long z=x; z<=y; z++){
						table.add(z);
					}
				}
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				assert(b!=null) : "Bad parameter: "+arg;
				in1=b;
				if(b.indexOf('#')>-1 && !new File(b).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				assert(b!=null) : "Bad parameter: "+arg;
				out1=b;
				if(b.indexOf('#')>-1){
					out1=b.replace("#", "1");
					out2=b.replace("#", "2");
				}
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("testsize")){
				testsize=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=arg.replace("#", "1");
					in2=arg.replace("#", "2");
				}
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();
		}

		assert(FastaReadInputStream.settingsOK());

		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		if(out1==null){
			if(out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			out1="stdout";
		}

		if(!parser.setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		final boolean useQualityInput=(qfin1!=null || qfin2!=null);
		final boolean useQualityOutput=(qfout1!=null || qfout2!=null);

		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(useQualityOutput && !Tools.testOutputFiles(overwrite, append, false, qfout1, qfout2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to quality output files "
				+qfout1+", "+qfout2+"\n");
		}
		if(useQualityInput && !Tools.testInputFiles(false, true, qfin1, qfin2)){
			throw new RuntimeException("\nCan't read quality input files "+qfin1+", "+qfin2+"\n");
		}
		if((useQualityInput || useQualityOutput) &&
				!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, qfin1, qfin2, qfout1, qfout2)){
			throw new RuntimeException("\nInput, output, and quality sidecar names must be unique.\n");
		}


		FileFormat ffin=FileFormat.testInput(in1, 0, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out1, 0, null, true, overwrite, append, false);
		FileFormat ffout2=(useQualityOutput ? FileFormat.testOutput(out2, 0, null, true, overwrite, append, false) : null);
		validateQualityOutput(ffout, ffout2, qfout1, qfout2);


		final boolean useSharedHeader=(ffin!=null && ffout!=null && ffin.samOrBam() && ffout.samOrBam());

		if(ffin!=null && ffout!=null && ffin.samOrBam() && (ffout.samOrBam() || ffout.bread())){
			throw new RuntimeException("\nDirect conversion of sam to sam or bread are not currently supported.\nAll other conversions are possible.");
		}


		ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			validateQualityInput(ff1, ff2, qfin1, qfin2);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ff1, ff2, qfin1, qfin2);
		}

		cris.setSampleRate(samplerate, sampleseed);
		outstream.println("Input is "+(cris.paired() ? "paired" : "unpaired"));
		cris.start(); //4567

		TextStreamWriter tsw=new TextStreamWriter(out1, overwrite, append, false);
		TextStreamWriter tsw2=(out2==null ? null : new TextStreamWriter(out2, overwrite, append, false));
		TextStreamWriter qtsw=(qfout1==null ? null : new TextStreamWriter(qfout1, overwrite, append, false, FileFormat.QUAL));
		TextStreamWriter qtsw2=(qfout2==null ? null : new TextStreamWriter(qfout2, overwrite, append, false, FileFormat.QUAL));
		tsw.start();
		if(tsw2!=null){tsw2.start();}
		if(qtsw!=null){qtsw.start();}
		if(qtsw2!=null){qtsw2.start();}


		long readsProcessed=0;
		long basesProcessed=0;

		for(int pass=1; pass<=passes; pass++){
//			outstream.println("pass="+pass);
			if(pass>1){
				FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ff1, ff2, qfin1, qfin2);
				cris.setSampleRate(samplerate, sampleseed);
				cris.start();
			}

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin==null || ffin.samOrBam() || (r.mate!=null)==cris.paired());//ffin cannot be null
			}

			while(reads!=null && reads.size()>0 && !table.isEmpty() && ln!=null){//ln!=null is implied
				ByteBuilder qbb1=null, qbb2=null;

				for(Read r1 : reads){
					{
						readsProcessed++;
						basesProcessed+=r1.length();
					}
					Read r2=r1.mate;
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=r2.length();
					}

					if(table.remove(r1.numericID)){
						tsw.println(r1);
						if(qtsw!=null){
							if(qbb1==null){qbb1=new ByteBuilder();}
							FastaQualWriterST.appendQual(r1, qbb1);
						}
						if(r2!=null){
							if(tsw2==null){
								tsw.println(r2);
								if(qtsw!=null){
									if(qbb1==null){qbb1=new ByteBuilder();}
									FastaQualWriterST.appendQual(r2, qbb1);
								}
							}else{
								tsw2.println(r2);
								if(qtsw2!=null){
									if(qbb2==null){qbb2=new ByteBuilder();}
									FastaQualWriterST.appendQual(r2, qbb2);
								}
							}
						}
						if(table.isEmpty()){break;}
					}
				}
				if(qbb1!=null && qbb1.length()>0){qtsw.print(qbb1);}
				if(qbb2!=null && qbb2.length()>0){qtsw2.print(qbb2);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			errorState|=ReadWrite.closeStream(cris);
		}

		if(tsw!=null){
			errorState|=tsw.poisonAndWait();
		}
		if(tsw2!=null){
			errorState|=tsw2.poisonAndWait();
		}
		if(qtsw!=null){
			errorState|=qtsw.poisonAndWait();
		}
		if(qtsw2!=null){
			errorState|=qtsw2.poisonAndWait();
		}

		errorState|=(cris.errorState());

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length()))*passes;
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+Tools.format("%.2fm bytes/sec", xpnano*1000));
		}

		if(errorState){
			throw new RuntimeException("GetReads terminated in an error state; the output may be corrupt.");
		}

	}

	/**
	 * Validate optional FASTA+QUAL sidecars before CRIS opens readers.
	 * FastaQualReadInputStream bypasses the ordinary FastaReadInputStream path for
	 * forced interleaving, amino reads, min-length filtering, and read splitting.
	 */
	private static void validateQualityInput(FileFormat ff1, FileFormat ff2, String qfin1, String qfin2){
		if(qfin1==null && qfin2==null){return;}

		if(qfin2!=null && ff2==null){
			throw new RuntimeException("qfin2 requires in2; forced interleaved FASTA+QUAL input is not supported.");
		}
		if(FASTQ.FORCE_INTERLEAVED){
			throw new RuntimeException("qfin/qfin2 do not support forced interleaved input.");
		}
		if(FastaReadInputStream.MIN_READ_LEN!=1 || FastaReadInputStream.SPLIT_READS || FastaReadInputStream.FORCE_SECTION_NAME){
			throw new RuntimeException("qfin/qfin2 require FASTA input with minscaf=1, no FASTA read splitting, and forceSectionName=f.");
		}
		if(Shared.AMINO_IN){
			throw new RuntimeException("qfin/qfin2 require nucleotide FASTA input.");
		}
		if(ff1==null || !ff1.fasta() || ff1.preferShreds() || ff1.amino()){
			throw new RuntimeException("qfin/qfin2 require the primary input to be non-shredding nucleotide FASTA.");
		}
		if(ff2!=null && (!ff2.fasta() || ff2.preferShreds() || ff2.amino())){
			throw new RuntimeException("qfin/qfin2 require the secondary input to be non-shredding nucleotide FASTA.");
		}
	}

	/** Validate optional QUAL sidecars before any output writer is opened. */
	private static void validateQualityOutput(FileFormat ffout1, FileFormat ffout2, String qfout1, String qfout2){
		if(qfout1==null && qfout2==null){return;}

		if(isConsoleOutputName(qfout1) || isConsoleOutputName(qfout2)){
			throw new RuntimeException("qfout/qfout2 require file output, not stdout or stderr.");
		}
		if(qfout1!=null && (ffout1==null || !ffout1.fasta())){
			throw new RuntimeException("qfout requires primary FASTA output.");
		}
		if(qfout2!=null && ffout2==null){
			throw new RuntimeException("qfout2 requires out2 FASTA output.");
		}
		if(qfout2!=null && !ffout2.fasta()){
			throw new RuntimeException("qfout2 requires secondary FASTA output.");
		}
	}

	/** Returns true for borrowed console-output aliases, including suffixed names. */
	private static boolean isConsoleOutputName(String name){
		if(name==null){return false;}
		String lower=name.trim().toLowerCase(Locale.ROOT);
		return lower.equals("stdout") || lower.equals("stderr") ||
				lower.equals("standardout") || lower.equals("standarderr") ||
				lower.startsWith("stdout.") || lower.startsWith("stderr.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for status and timing messages; PreParser may redirect it. */
	private PrintStream outstream=System.err;

	/** Enables verbose input diagnostics. */
	public static boolean verbose=false;

}
