package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import align2.QualityTools;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Mar 16, 2014
 *
 */
public class AddAdapters {
	
	public static void main(String[] args){
		Timer t=new Timer();
		AddAdapters x=new AddAdapters(args);
		if(x.writeMode){
			x.write(t);
		}else{
			x.read(t);
		}
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public AddAdapters(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.PARSE_CUSTOM=true;
		
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
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
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("adapter") || a.equals("adapters") || a.equals("ref")){
				adapterFile=b;
			}else if(a.equals("literal") || a.equals("literals")){
				literals=(b==null ? null : b.split(","));
			}else if(a.equals("rate") || a.equals("prob")){
				adapterProb=Float.parseFloat(b);
			}else if(a.equals("minlength") || a.equals("minlen") || a.equals("ml")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("3'") || a.equalsIgnoreCase("3prime") || a.equalsIgnoreCase("3-prime") || a.equalsIgnoreCase("right") || a.equalsIgnoreCase("r")){
				right=Parse.parseBoolean(b);
			}else if(a.equals("5'") || a.equalsIgnoreCase("5prime") || a.equalsIgnoreCase("5-prime") || a.equalsIgnoreCase("left") || a.equalsIgnoreCase("l")){
				right=!Parse.parseBoolean(b);
			}else if(a.equals("end")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.equals("3'") || b.equalsIgnoreCase("3prime") || b.equalsIgnoreCase("3-prime") || b.equalsIgnoreCase("right") || a.equalsIgnoreCase("r")){
					right=true;
				}else if(b.equals("5'") || b.equalsIgnoreCase("5prime") || b.equalsIgnoreCase("5-prime") || b.equalsIgnoreCase("left") || a.equalsIgnoreCase("l")){
					right=true;
				}
			}else if(a.equals("addslash")){
				addslash=Parse.parseBoolean(b);
			}else if(a.equals("adderrors")){
				adderrors=Parse.parseBoolean(b);
			}else if(a.equals("addreversecomplement") || a.equals("arc")){
				addRC=Parse.parseBoolean(b);
			}else if(a.equals("addpaired")){
				addPaired=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("write")){
				writeMode=Parse.parseBoolean(b);
			}else if(a.equals("grade")){
				writeMode=!Parse.parseBoolean(b);
			}else if(a.equals("mode")){
				if("grade".equalsIgnoreCase(b) || "read".equalsIgnoreCase(b)){
					writeMode=false;
				}else if("generate".equalsIgnoreCase(b) || "write".equalsIgnoreCase(b) || "add".equalsIgnoreCase(b)){
					writeMode=true;
				}else{
					throw new RuntimeException("Unknown mode "+b);
				}
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(writeMode && out1==null){
			if(out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			System.err.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
		}
		
		if(!parser.setInterleaved){
			assert(in1!=null && (!writeMode || out1!=null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else if(writeMode){ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		ffa=FileFormat.testInput(adapterFile, FileFormat.FASTA, null, true, true);
		
		adapters=makeAdapterList();
		
		if(writeMode){
			if(adapters==null || adapters.isEmpty()){
				throw new RuntimeException("\n\nPlease specify adapters with 'adapters=file.fa' or 'literal=AGCTACGT'\n");
			}
			randy=Shared.threadLocalRandom();
		}
	}
	
	private final ArrayList<byte[]> makeAdapterList(){
		boolean oldTI=FASTQ.TEST_INTERLEAVED;
		boolean oldFI=FASTQ.FORCE_INTERLEAVED;
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		ArrayList<byte[]> x=makeAdapterList2();
		FASTQ.TEST_INTERLEAVED=oldTI;
		FASTQ.FORCE_INTERLEAVED=oldFI;
		return x;
	}
	
	private final ArrayList<byte[]> makeAdapterList2(){
		if(ffa==null && literals==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		if(ffa!=null){
			FastaReadInputStream fris=new FastaReadInputStream(ffa, false, false, -1);
			for(Read r=fris.next(); r!=null; r=fris.next()){
				if(r.bases!=null){
					list.add(r.bases);
				}
			}
			fris.close();
		}
		if(literals!=null){
			for(String s : literals){
				if(s!=null && !"null".equalsIgnoreCase(s)){
					list.add(s.getBytes());
				}
			}
		}
		
		if(addRC){
			int x=list.size();
			for(int i=0; i<x; i++){
				list.add(AminoAcid.reverseComplementBases(list.get(i)));
			}
		}
		
		return list.size()>0 ? list : null;
	}
	
	void write(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
			ros.start();
		}

		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					addAdapter(r1, addPaired);
					if(r2!=null && !addPaired){
						addAdapter(r2, addPaired);
					}
					
					if(r2==null){
						r1.id=r1.numericID+"_"+r1.id;
					}else{
						String base=r1.numericID+"_"+r1.id+"_"+r2.id;
						if(addslash){
							r1.id=base+" /1";
							r2.id=base+" /2";
						}else{
							r1.id=base;
							r2.id=base;
						}
					}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
//		System.err.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			System.err.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			System.err.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
//		assert(false);
		
		t.stop();

		outstream.println("Adapters Added:         \t"+adaptersAdded+" reads ("+Tools.format("%.2f",adaptersAdded*100.0/readsProcessed)+"%) \t"+
				adapterBasesAdded+" bases ("+Tools.format("%.2f",adapterBasesAdded*100.0/basesProcessed)+"%)");

		outstream.println("Valid Output:           \t"+validReads+" reads ("+Tools.format("%.2f",validReads*100.0/readsProcessed)+"%) \t"+
				validBases+" bases ("+Tools.format("%.2f",validBases*100.0/basesProcessed)+"%)");
		
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void addAdapter(Read r, final int loc){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final int remaining, initial=(bases==null ? 0 : bases.length);
		final byte[] adapter;
		int ab=0, rb=0;
		
		readsProcessed++;
		basesProcessed+=initial;
		
		if(bases==null){assert(false); return;}
		if(initial>0 && loc>=0 && loc<initial){
			adapter=adapters.get(randy.nextInt(adapters.size()));
			adaptersAdded++;

			if(right){
				final int lim=Tools.min(initial, adapter.length+loc);
				for(int i=loc, j=0; i<lim; i++, j++){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=adapter[j];
						if(adderrors){
							byte q=(quals==null ? 30 : quals[i]);
							if(randy.nextFloat()<QualityTools.PROB_ERROR[q]){
								int old=AminoAcid.baseToNumber[bases[i]];
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt(3)+1)&3];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i<initial; i++){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];
					}
					rb++;
				}
				remaining=loc;
			}else{
				final int lim=Tools.max(-1, loc-adapter.length);
				for(int i=loc, j=adapter.length-1; i>lim; i--, j--){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=adapter[j];
						if(adderrors){
							byte q=(quals==null ? 30 : quals[i]);
							if(randy.nextFloat()<QualityTools.PROB_ERROR[q]){
								int old=AminoAcid.baseToNumber[bases[i]];
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt(3)+1)&3];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i>-1; i--){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];
					}
					rb++;
				}
				remaining=initial-loc-1;
			}
			assert(remaining<initial) : "\nremaining="+remaining+", initial="+initial+", rb="+rb+", ab="+ab+
				", loc="+loc+", adapter.length="+(adapter==null ? 0 : adapter.length)+"\n";
		}else{
			adapter=null;
			remaining=initial;
		}
		
		assert(remaining==initial-(rb+ab));
		assert(remaining>=0);

		adapterBasesAdded+=ab;
		randomBasesAdded+=rb;
		r.id=initial+"_"+remaining;
		if(remaining>=minlen){
			validReads++;
			validBases+=remaining;
		}
	}
	
	private void addAdapter(Read r, boolean addPaired){
		final byte[] bases=r.bases;
		final int initial=(bases==null ? 0 : bases.length);
		final int loc;
		
		if(initial>0 && randy.nextFloat()<adapterProb){
			loc=randy.nextInt(initial);
		}else{
			loc=-1;
		}
		
		addAdapter(r, loc);
		if(addPaired && r.mate!=null){addAdapter(r.mate, loc);}
	}
	
	/*--------------------------------------------------------------*/
	
	void read(Timer t){

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					grade(r1, r2);
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
		
		t.stop();
		
		long validBasesRemoved=validBasesExpected-validBasesCounted;
		long incorrect=readsProcessed-correct;
		long incorrectBases=basesProcessed-correctBases;
		
		outstream.println("Total output:                        \t"+readsProcessed+" reads                  \t"+basesProcessed+" bases          ");
		outstream.println("Perfectly Correct (% of output):     \t"+correct+" reads ("+Tools.format("%.3f",correct*100.0/readsProcessed)+
				"%)        \t"+correctBases+" bases ("+Tools.format("%.3f",correctBases*100.0/basesProcessed)+"%)");
		outstream.println("Incorrect (% of output):             \t"+incorrect+" reads ("+Tools.format("%.3f",incorrect*100.0/readsProcessed)+
				"%)        \t"+incorrectBases+" bases ("+Tools.format("%.3f",incorrectBases*100.0/basesProcessed)+"%)");
		outstream.println();
//		outstream.println("Too Short:              \t"+tooShort+" reads ("+Tools.format("%.3f",tooShort*100.0/readsProcessed)+"%) \t"+
//				tooShortBases+" bases ("+Tools.format("%.3f",tooShortBases*100.0/basesProcessed)+"%)");
//		outstream.println("Too Long:               \t"+tooLong+" reads ("+Tools.format("%.3f",tooLong*100.0/readsProcessed)+"%) \t"+
//				tooLongBases+" bases ("+Tools.format("%.3f",tooLongBases*100.0/basesProcessed)+"%)");
		
		outstream.println("Adapters Remaining (% of adapters):  \t"+(adapterReadsRemaining)+" reads ("+Tools.format("%.3f",adapterReadsRemaining*100.0/adapterReadsTotal)+
				"%)        \t"+adapterBasesRemaining+" bases ("+Tools.format("%.3f",adapterBasesRemaining*100.0/basesProcessed)+"%)");
		outstream.println("Non-Adapter Removed (% of valid):    \t"+tooShort+" reads ("+Tools.format("%.4f",tooShort*100.0/readsProcessed)+
				"%)        \t"+validBasesRemoved+" bases ("+Tools.format("%.4f",validBasesRemoved*100.0/validBasesExpected)+"%)");
		
		if(broken>0 || mispaired>0){
			outstream.println("Broken:                              \t"+broken+" reads ("+Tools.format("%.2f",broken*100.0/readsProcessed)+"%)");
			outstream.println("Mispaired:                           \t"+mispaired+" reads ("+Tools.format("%.2f",mispaired*100.0/readsProcessed)+"%)");
		}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
//	private void grade_old(Read r1, Read r2){
//		
//		final String a=r1.id.split(" ")[0];
//		final String b=(r2==null ? a : r2.id.split(" ")[0]);
//		final int len=a.split("_").length;
//		
//		if(r2!=null){
//			if(r1.id.endsWith(" /2") || r2.id.endsWith(" /1") || !a.equals(b)){
//				mispaired+=2;
//			}
//			if(len==3){
//				r2.setPairnum(0);
//			}else if(len==5){
//				if(r1.id.endsWith(" /2")){r1.setPairnum(1);}
//				if(r2.id.endsWith(" /1")){r2.setPairnum(0);}
//			}else{
//				throw new RuntimeException("Headers are corrupt. They must be generated by AddAdapters or RenameReads.");
//			}
//		}else{
//			if(len!=3){
//				throw new RuntimeException("Headers are corrupt, or paired reads are being processed as unpaired.  Try running with 'int=t' or with 'in1=' and 'in2='");
//			}
//		}
//		grade(r1);
//		grade(r2);
//	}
	
	private void grade(Read r1, Read r2){
		grade(r1);
		grade(r2);
	}
	
	private void grade(Read r){
		if(r==null){return;}
		
		int insert=r.insert();
		int originalLength=r.stop-r.start+1;
		int length=r.length();
		
		final int offset=(2*r.pairnum());
		
//		String[] sa=r.id.split(" ")[0].split("_");
//		final long id=Long.parseLong(sa[0]);
		final int initial=originalLength;
		final int remaining=Tools.min(initial, insert);
		final int actual=length;
		
		readsProcessed++;
		basesProcessed+=actual;
		
		assert(initial>=remaining);
		
		if(actual>initial){broken++;}
		
		validBasesExpected+=remaining;
		

//		System.err.println("initial="+initial+", remaining="+remaining+", actual="+actual);
		
		if(initial==remaining){//Should not have trimmed
			if(actual==remaining || (actual<2 && (remaining<1 || remaining<minlen))){
				correct++;
				correctBases+=remaining;
				validBasesCounted+=remaining;
				trueNeg++;
			}else if(actual<remaining){
				tooShort++;
				tooShortReadBases+=actual;
				tooShortBases+=(remaining-actual);
				validBasesCounted+=actual;
				falsePos++;
			}else if(actual>remaining){
				tooLong++;
				tooLongReadBases+=remaining;
				tooLongBases+=(actual-remaining);
				validBasesCounted+=remaining;
				falseNeg++;
			}
		}else{//Should have trimmed
			
			adapterBasesTotal+=(initial-remaining);
			adapterReadsTotal++;
			
			if(actual==remaining || (actual<2 && (remaining<1 || remaining<minlen))){
				correct++;
				correctBases+=remaining;
				validBasesCounted+=remaining;
				truePos++;
			}else if(actual<remaining){
				tooShort++;
				tooShortReadBases+=actual;
				tooShortBases+=(remaining-actual);
				validBasesCounted+=actual;
				truePos++;
			}else if(actual>remaining){
				tooLong++;
				tooLongReadBases+=actual;
				tooLongBases+=(actual-remaining);
				adapterBasesRemaining+=(actual-remaining);
				validBasesCounted+=remaining;
				falseNeg++;
				adapterReadsRemaining++;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public boolean errorState=false;
	
	private String in1=null;
	private String in2=null;

	private String out1=null;
	private String out2=null;
	
	private String extin=null;
	private String extout=null;

	private String adapterFile=null;
	private String[] literals=null;
	
	private boolean overwrite=true;
	private boolean append=false;

	/** Add /1 and /2 to paired reads */
	private boolean addslash=true;
	/** Encode correct answer in read ID field */
	private boolean changename=true;
	/** Add errors from quality value */
	private boolean adderrors=true;

	/** Add adapters to the same location for read 1 and read 2 */
	private boolean addPaired=true;
	/** Add reverse-complemented adapters also */
	private boolean addRC=false;
	/** aka 3' */
	private boolean right=true;
	
	private long maxReads=-1;
	private int minlen=1;
	
	private boolean writeMode=true;
	private float adapterProb=0.5f;
	
	private long readsProcessed=0;
	private long basesProcessed=0;
	private long adaptersAdded=0;
	private long adapterBasesAdded=0;
	private long randomBasesAdded=0;
	private long validReads=0;
	private long validBases=0;

	private long truePos=0;
	private long trueNeg=0;
	private long falsePos=0;
	private long falseNeg=0;
	private long broken=0;
	private long mispaired=0;
	
	private long tooShort=0;
	private long tooLong=0;
	private long correct=0;
	private long fullyRemoved=0;

	private long tooShortBases=0;
	private long tooLongBases=0;
	private long tooShortReadBases=0;
	private long tooLongReadBases=0;
	private long correctBases=0;

	private long validBasesCounted=0;
	private long validBasesExpected=0;
	
//	private long invalidBasesCounted=0;
	private long adapterBasesTotal=0;
	private long adapterReadsTotal=0;
	private long adapterReadsRemaining=0;
	private long adapterBasesRemaining=0;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	private final FileFormat ffa;
	
	private final ArrayList<byte[]> adapters;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	
	private java.util.Random randy;
	
}
