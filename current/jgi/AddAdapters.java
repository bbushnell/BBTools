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
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * Generates insert-unaware adapter-contaminated read data, with a legacy length-grading mode.
 * Generation replaces bases without changing read lengths or quality arrays,
 * then generates {@code id_original_remaining} or paired
 * {@code id_orig1_rem1_orig2_rem2} headers. Paired headers may append
 * {@code " /1"} or {@code " /2"} tags.
 * Mates may share a replacement position; they still choose adapters independently.
 * The minimum-length setting classifies statistics; it does not filter output.
 * Grading compares surviving read lengths with annotations in their IDs, not sequence identity.
 * It cannot account for reads discarded upstream and does not run trimming tools.
 * Its classifications are historical length-based rules, not adapter-sequence detection.
 * This is the legacy utility behind the deprecated {@code addadapters.sh} launcher.
 * Construction changes process-wide parser and I/O settings; instances are not
 * isolated for concurrent library use.
 * @author Brian Bushnell
 * @date Mar 16, 2014
 */
public class AddAdapters{

	/**
	 * Runs generation or grading, then closes redirected diagnostics on normal return.
	 * Diagnostic closure follows successful processing; standard output/error
	 * remain open through Shared.closeStream.
	 * @param args Command-line arguments
	 */
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

	/**
	 * Parses options, configures shared settings, prepares formats and loads adapters.
	 * Defaults to ordinary header handling; explicit parsecustom options retain their
	 * normal precedence. Generation asserts case-insensitive output/read-input name
	 * differences before writer acquisition; this is not filesystem alias resolution.
	 * Main streams are acquired by write/read; adapter-reference input is consumed
	 * here. Shared parser/I/O changes are not rolled back as a whole.
	 * @param args Command-line arguments array
	 */
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
		FASTQ.PARSE_CUSTOM=false;

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
				if(b.equals("3'") || b.equalsIgnoreCase("3prime") || b.equalsIgnoreCase("3-prime") || b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){
					right=true;
				}else if(b.equals("5'") || b.equalsIgnoreCase("5prime") || b.equalsIgnoreCase("5-prime") || b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){
					right=false;
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
		if(writeMode && out1!=null){
			// Reject aliases before starting producers or opening a writer that could truncate an input.
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in2)) :
				"Input file and output file have same name; refusing to overwrite input: "+out1;
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) :
				"Input file and output file have same name; refusing to overwrite input: "+out2;
		}

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

	/**
	 * Loads adapters with interleaving detection/forcing temporarily disabled.
	 * Restores those two inherited flags on return or failure; this does not restore
	 * every process-global setting changed by construction or support concurrent use.
	 * @return Mutable adapter list, or null when no sequences remain
	 */
	private final ArrayList<byte[]> makeAdapterList(){
		final boolean oldTI=FASTQ.TEST_INTERLEAVED;
		final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		try{return makeAdapterList2();}
		finally{
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
		}
	}

	/**
	 * Builds adapter sequence list from both file-based and literal sources.
	 * Retains FASTA base arrays without defensive copying, converts literals with the
	 * default charset and optionally appends reverse-complement copies. Entries are
	 * not deduplicated; input quality arrays are not used.
	 * Closes an acquired reader even after failure; read/close errors reject the adapter
	 * file before literals are appended, preserving any earlier unchecked failure.
	 * @return Mutable adapter list, or null when no sequences remain
	 */
	private final ArrayList<byte[]> makeAdapterList2(){
		if(ffa==null && literals==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		if(ffa!=null){
			final FastaReadInputStream fris=new FastaReadInputStream(ffa, false, false, -1);
			Throwable primary=null;
			try{
				for(ArrayList<Read> reads=fris.nextList(); reads!=null; reads=fris.nextList()){
					for(Read r : reads){
						if(r.bases!=null){list.add(r.bases);}
					}
				}
			}catch(RuntimeException | Error e){
				primary=e;
				throw e;
			}finally{
				try{
					if(fris.close()){throw new RuntimeException("Error reading adapter file: "+ffa.name());}
				}catch(RuntimeException | Error e){
					if(primary==null){throw e;}
					if(e!=primary){primary.addSuppressed(e);}
				}
			}
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

	/**
	 * Generates fixed-length contaminated reads with annotated length IDs.
	 * Edits arrays/IDs and accumulates mate-inclusive counters on this caller while
	 * stream wrappers perform background I/O. minlen classifies expected remainders;
	 * every normally processed read is submitted when output is configured.
	 * Once input acquisition returns, attempts a pending-batch return and independent
	 * shutdown in finally. Processing/return failures abort queued output and remain
	 * primary if cleanup also fails. Factory acquisition, partial startup, failed
	 * close requests and blocked I/O retain their own limits. Reported close errors
	 * are rejected after normal reporting; counters do not confirm physical writes.
	 *
	 * @param t Running timer stopped after normal processing and cleanup
	 */
	void write(Timer t){

		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
		ConcurrentReadOutputStream ros=null;
		boolean outputStarted=false;
		ListNum<Read> ln=null;
		Throwable failure=null;
		try{
			if(verbose){System.err.println("Started cris");}
			cris.start();
			final boolean paired=cris.paired();
			if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}
			if(out1!=null){
				final int buff=4;
				if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
					outstream.println("Writing interleaved.");
				}
				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
				ros.start();
				outputStarted=true;
			}
			ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

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

				final ListNum<Read> spent=ln;
				ln=null;
				cris.returnList(spent);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			finishGenerationStreams(cris, ros, outputStarted, ln, failure);
		}

		t.stop();

		outstream.println("Adapters Added:         \t"+adaptersAdded+" reads ("+Tools.format("%.2f",adaptersAdded*100.0/readsProcessed)+"%) \t"+
				adapterBasesAdded+" bases ("+Tools.format("%.2f",adapterBasesAdded*100.0/basesProcessed)+"%)");

		outstream.println("Valid Output:           \t"+validReads+" reads ("+Tools.format("%.2f",validReads*100.0/readsProcessed)+"%) \t"+
				validBases+" bases ("+Tools.format("%.2f",validBases*100.0/basesProcessed)+"%)");

		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		if(errorState){
			throw new RuntimeException("AddAdapters terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Replaces the suffix from loc through the end, or the prefix from zero through loc.
	 * Right replacement starts with the adapter's first base; left replacement uses
	 * its last base at loc and walks backward. One adapter entry is sampled, followed
	 * by random filler where it does not cover the selected region. Undefined input
	 * bases remain unchanged, but their selected positions still enter counters.
	 * Optional substitutions affect adapter-covered, defined input positions, using
	 * the existing quality or Q30 when quality is null. Qualities/length stay unchanged.
	 * Rewrites the ID as {@code initial_remaining}, even when no region is selected.
	 *
	 * @param r Nonnull read with a nonnull base array; modified in place
	 * @param loc Inclusive replacement boundary; an out-of-read position means no adapter
	 */
	private void addAdapter(Read r, final int loc){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final int remaining, initial=(bases==null ? 0 : bases.length);
		final byte[] adapter;
		int ab=0, rb=0;

		readsProcessed++;
		basesProcessed+=initial;

		if(bases==null){assert(false) : "Adapter generation requires read bases; read id="+r.id; return;}
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
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt3()+1)&3];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i<initial; i++){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt()&3];
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
								bases[i]=AminoAcid.numberToBase[(old+randy.nextInt3()+1)&3];
							}
						}
					}
					ab++;
				}
				for(int i=lim; i>-1; i--){
					if(AminoAcid.isFullyDefined(bases[i])){
						bases[i]=AminoAcid.numberToBase[randy.nextInt()&3];
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

	/**
	 * Samples whether and where to replace bases, using this read's length for the position.
	 * Selection uses this instance's adapterProb and random source. When requested,
	 * the mate receives the same position but independently chooses its adapter.
	 * Insert size is ignored. An out-of-range position on a shorter mate leaves that
	 * mate's bases unchanged while still updating its counters and length ID.
	 *
	 * @param r Nonnull read supplying the placement length; must have sequence bases
	 * @param addPaired Whether to apply the sampled position to a linked mate as well
	 */
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

	/**
	 * Compares surviving read lengths with expected lengths recorded in their IDs.
	 * Attempts one outstanding-batch return and input closure after acquisition,
	 * preserving an earlier processing exception. Does not compare sequences or infer
	 * discarded reads. A zero/one-base result may count as correct under the explicit
	 * expected-remainder/minlen tolerance. Totals describe surviving annotations;
	 * their credit and denominator rules are not direct adapter detection.
	 * Prints normal summaries before rejecting any reported input-error flag.
	 *
	 * @param t Running timer stopped after normal processing and cleanup
	 */
	void read(Timer t){

		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
		ListNum<Read> ln=null;
		Throwable failure=null;
		try{
			if(verbose){System.err.println("Started cris");}
			cris.start();
			final boolean paired=cris.paired();
			if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}
			ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

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

				final ListNum<Read> spent=ln;
				ln=null;
				cris.returnList(spent);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}catch(RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{
			finishGradingInput(cris, ln, failure);
		}

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

		outstream.println("Adapters Remaining (% of adapters):  \t"+(adapterReadsRemaining)+" reads ("+Tools.format("%.3f",adapterReadsRemaining*100.0/adapterReadsTotal)+
				"%)        \t"+adapterBasesRemaining+" bases ("+Tools.format("%.3f",adapterBasesRemaining*100.0/basesProcessed)+"%)");
		outstream.println("Non-Adapter Removed (% of valid):    \t"+tooShort+" reads ("+Tools.format("%.4f",tooShort*100.0/readsProcessed)+
				"%)        \t"+validBasesRemoved+" bases ("+Tools.format("%.4f",validBasesRemoved*100.0/validBasesExpected)+"%)");

		if(broken>0 || mispaired>0){
			outstream.println("Broken:                              \t"+broken+" reads ("+Tools.format("%.2f",broken*100.0/readsProcessed)+"%)");
			outstream.println("Mispaired:                           \t"+mispaired+" reads ("+Tools.format("%.2f",mispaired*100.0/readsProcessed)+"%)");
		}

		if(errorState){
			throw new RuntimeException("AddAdapters terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Applies the single-read length grader to both reads without changing mate links.
	 * @param r1 First surviving read, or null
	 * @param r2 Surviving mate, or null
	 */
	private void grade(Read r1, Read r2){
		grade(r1);
		grade(r2);
	}

	/**
	 * Tallies one surviving read using its annotated original and retained lengths.
	 * Decodes before changing counters and borrows all Read/mate state without modification.
	 * A linked mate supplies fallback mate context; an unpaired read's default pairnum does not.
	 * Accepted short-read tolerance credits the expected remainder, not necessarily
	 * the observed base count. Overtrimmed annotated positives still count as true
	 * positives; non-tolerated overlong results enter the historical false-negative bucket.
	 *
	 * @param r Surviving annotated read, or null to skip
	 * @throws IllegalArgumentException If length/mate annotations cannot be decoded
	 */
	private void grade(Read r){
		if(r==null){return;}

		final long lengths=decodeExpectedLengths(r.id, r.mate==null ? -1 : r.pairnum());
		final int initial=(int)(lengths>>>32);
		final int remaining=(int)(lengths&0xFFFFFFFFL);
		final int actual=r.length();

		readsProcessed++;
		basesProcessed+=actual;

		assert(initial>=remaining);

		if(actual>initial){broken++;}

		validBasesExpected+=remaining;

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

	/** Attempts a pending generation return, input close and independent output shutdown.
	 * Aborts output before input closure when processing or batch return failed.
	 * Joins only after output start returned successfully: a partially started backend
	 * can remain NEW, which the output join contract does not support.
	 * Underlying blocking I/O and partial startup retain their own lifecycle limits.
	 * A failed output close request is not followed by a separate abort fallback here.
	 * @param cris Captured input wrapper
	 * @param ros Captured output wrapper, or null
	 * @param outputStarted Whether the output's outer start returned successfully
	 * @param pending Borrowed batch still needing one return attempt, or null
	 * @param primary Existing unchecked failure, or null on normal processing
	 */
	private void finishGenerationStreams(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros,
			final boolean outputStarted, final ListNum<Read> pending, final Throwable primary){
		Throwable failure=primary;
		try{cris.returnList(pending);}
		catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		if(ros!=null && failure!=null){
			try{ros.abort();}
			catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		}
		try{errorState|=ReadWrite.closeStream(cris);}
		catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		if(ros!=null){
			try{
				if(failure==null && outputStarted){ros.close();}
				else{ros.abort();}
			}catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
			if(outputStarted){
				try{
					ros.join();
					errorState|=ros.errorState() || !ros.finishedSuccessfully();
				}catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
			}else{errorState=true;}
		}
		if(primary==null && failure!=null){
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw (Error)failure;
		}
	}

	/** Attempts a pending return once and input closure even if either step throws.
	 * Keeps an existing primary failure and suppresses a distinct cleanup failure.
	 * Factory-internal acquisition and blocking I/O are outside this helper's control.
	 * @param cris Captured input wrapper, or null to do nothing
	 * @param pending Borrowed batch needing one return attempt, or null
	 * @param primary Existing unchecked failure, or null on normal processing */
	private void finishGradingInput(final ConcurrentReadInputStream cris, final ListNum<Read> pending, final Throwable primary){
		if(cris==null){return;}
		Throwable failure=primary;
		try{cris.returnList(pending);}
		catch(RuntimeException|Error e){failure=retainFailure(failure, e);}
		try{errorState|=ReadWrite.closeStream(cris);}
		catch(RuntimeException|Error e){failure=retainFailure(failure, e);}
		if(primary==null && failure!=null){
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw (Error)failure;
		}
	}

	/** Keeps the first failure and attaches a different secondary object.
	 * @param primary Earlier failure, or null
	 * @param secondary Nonnull cleanup failure
	 * @return Primary when present, otherwise secondary; never self-suppresses */
	private static Throwable retainFailure(final Throwable primary, final Throwable secondary){
		if(primary==null){return secondary;}
		if(primary!=secondary){primary.addSuppressed(secondary);}
		return primary;
	}

	/** Decodes original length (upper32 bits) and expected remainder (lower32 bits).
	 * Accepts three-field per-read or five-field paired length tags, and SYN_ bodies.
	 * Explicit {@code /1}/{@code /2} tags, attached or after spaces, or space-delimited
	 * {@code 1:}/{@code 2:} markers identify the original mate before paired-context
	 * fallback. Unknown mate context is -1; ambiguous paired tags fail.
	 * Validates length fields and body structure, not the full unused SYN_ metadata schema.
	 * Does not modify parser globals or allocate on the successful path.
	 * @param id Original header annotation, retained by upstream trimming
	 * @param mate Known paired position 0 or 1, or -1 for an unpaired read
	 * @return Packed nonnegative int lengths, with remaining no greater than original
	 * @throws IllegalArgumentException If required fields, lengths or mate context are invalid
	 */
	private static long decodeExpectedLengths(final String id, int mate){
		if(id==null || id.isEmpty()){throw badHeader(id, "missing length annotation");}
		if(mate<-1 || mate>1){throw badHeader(id, "mate context must be -1, 0 or 1");}
		final int space=id.indexOf(' ');
		int end=space<0 ? id.length() : space;
		int suffix=-1;
		if(end>=2 && id.charAt(end-2)=='/'){
			suffix=headerMate(id.charAt(end-1), id);
			end-=2;
		}
		if(space>=0){
			int p=space;
			while(p<id.length() && id.charAt(p)==' '){p++;}
			int other=-1;
			if(p+1<id.length() && id.charAt(p)=='/' && id.charAt(p+1)>='0' && id.charAt(p+1)<='9'){
				other=headerMate(id.charAt(p+1), id);
				if(p+2<id.length() && id.charAt(p+2)!=' '){throw badHeader(id, "invalid slash mate suffix");}
			}else if(p+1<id.length() && id.charAt(p+1)==':' && id.charAt(p)>='0' && id.charAt(p)<='9'){
				other=headerMate(id.charAt(p), id);
			}
			if(other>=0){
				if(suffix>=0 && suffix!=other){throw badHeader(id, "conflicting mate suffixes");}
				suffix=other;
			}
		}
		if(suffix>=0){mate=suffix;}
		if(id.startsWith("SYN_")){
			final int middle=headerDelimiter(id, '&', 4, end);
			if(middle<0){return syntheticLengths(id, 4, end);}
			if(headerDelimiter(id, '&', middle+1, end)>=0){throw badHeader(id, "more than two synthetic bodies");}
			if(mate<0){throw badHeader(id, "two synthetic bodies require a mate suffix or paired input");}
			return mate==0 ? syntheticLengths(id, 4, middle) : syntheticLengths(id, middle+1, end);
		}
		final int a=headerDelimiter(id, '_', 0, end);
		if(a<0){throw badHeader(id, "expected id_original_remaining or a paired length tag");}
		headerNumber(id, 0, a, true);//Validate the numeric identity without storing it.
		final int b=headerDelimiter(id, '_', a+1, end);
		if(b<0){throw badHeader(id, "missing remaining length");}
		final int c=headerDelimiter(id, '_', b+1, end);
		final int initial=headerInt(id, a+1, b, false);
		if(c<0){return packLengths(id, initial, headerInt(id, b+1, end, false));}
		final int d=headerDelimiter(id, '_', c+1, end);
		if(d<0 || headerDelimiter(id, '_', d+1, end)>=0){throw badHeader(id, "length tags require three or five fields");}
		final long first=packLengths(id, initial, headerInt(id, b+1, c, false));
		final long second=packLengths(id, headerInt(id, c+1, d, false), headerInt(id, d+1, end, false));
		if(mate<0){throw badHeader(id, "paired length tags require a mate suffix or paired input");}
		return mate==0 ? first : second;
	}

	/** Selects the original mate; other markers must not select a length tuple.
	 * @param marker Original mate digit, '1' or '2'
	 * @param id Header retained for diagnostics
	 * @return Zero for read1, one for read2
	 * @throws IllegalArgumentException If the marker is not '1' or '2' */
	private static int headerMate(final char marker, final String id){
		if(marker=='1'){return 0;}
		if(marker=='2'){return 1;}
		throw badHeader(id, "mate suffix must identify read1 or read2");
	}

	/** Finds a separator inside one selected body, excluding suffixes/descriptions.
	 * @param id Whole header
	 * @param delimiter Separator to find
	 * @param from Inclusive slice start
	 * @param end Exclusive slice end
	 * @return Separator index within the slice, or -1 */
	private static int headerDelimiter(final String id, final char delimiter, final int from, final int end){
		assert(from>=0 && from<=end && end<=id.length()) :
			"Header field scans must stay inside their selected body: from="+from+", end="+end+", length="+id.length();
		final int found=id.indexOf(delimiter, from);
		return found>=0 && found<end ? found : -1;
	}

	/** Reads the required numeric fields of one nine-field SYN_ body without allocating a CustomHeader.
	 * Original length comes from the inclusive span. Requires known positive insert
	 * metadata, clamped to that span; nonpositive insert is rejected as unknown,
	 * unlike a valid explicit zero remainder in a legacy length tag.
	 * @param id Whole header
	 * @param from Inclusive selected-body start, after SYN_ or a body separator
	 * @param end Exclusive selected-body end
	 * @return Packed original/expected lengths
	 * @throws IllegalArgumentException If the required numeric fields or span/insert are invalid */
	private static long syntheticLengths(final String id, int from, final int end){
		int start=0, stop=0, insert=0;
		for(int field=0; field<9; field++){
			final int separator=headerDelimiter(id, '_', from, end);
			final int to=separator<0 ? end : separator;
			if(from>=to || (field<8 && separator<0) || (field==8 && separator>=0)){
				throw badHeader(id, "synthetic bodies require nine nonempty fields");
			}
			if(field==0){headerNumber(id, from, to, true);}
			else if(field==1){start=headerInt(id, from, to, true);}
			else if(field==2){stop=headerInt(id, from, to, true);}
			else if(field==3){insert=headerInt(id, from, to, true);}
			if(field<8){from=to+1;}
		}
		final long span=(long)stop-start+1;
		if(span<0 || span>Integer.MAX_VALUE){throw badHeader(id, "synthetic span is outside the supported length range");}
		if(insert<=0){throw badHeader(id, "synthetic insert must be known and positive");}
		final int initial=(int)span;
		return packLengths(id, initial, Tools.min(initial, insert));
	}

	/** Packs validated lengths without allowing the two int components to overlap.
	 * @param id Header retained for diagnostics
	 * @param initial Nonnegative original length
	 * @param remaining Expected length between zero and initial, inclusive
	 * @return Original length in the upper 32 bits, expected length in the lower 32
	 * @throws IllegalArgumentException If the lengths violate those bounds */
	private static long packLengths(final String id, final int initial, final int remaining){
		if(initial<0 || remaining<0 || remaining>initial){throw badHeader(id, "require 0 <= remaining <= original length");}
		return ((long)initial<<32)|(long)remaining;
	}

	/** Narrows only after checked decimal conversion and int-range validation.
	 * @param id Whole header
	 * @param from Inclusive numeric-field start
	 * @param end Exclusive numeric-field end
	 * @param signed Whether an initial plus/minus sign is allowed
	 * @return Parsed int
	 * @throws IllegalArgumentException If the slice is malformed or outside int range */
	private static int headerInt(final String id, final int from, final int end, final boolean signed){
		final long value=headerNumber(id, from, end, signed);
		if(value<Integer.MIN_VALUE || value>Integer.MAX_VALUE){throw badHeader(id, "integer field overflow");}
		return (int)value;
	}

	/** Parses a decimal slice without allocation, using negative accumulation to represent MIN_VALUE.
	 * Checks bounds before each multiply/subtract so long fields cannot wrap around.
	 * @param id Whole header
	 * @param from Inclusive numeric-field start
	 * @param end Exclusive numeric-field end
	 * @param signed Whether an initial plus/minus sign is allowed
	 * @return Parsed long
	 * @throws IllegalArgumentException If the field is empty, malformed or outside long range */
	private static long headerNumber(final String id, int from, final int end, final boolean signed){
		assert(from>=0 && from<=end && end<=id.length()) :
			"Numeric header slices must stay inside the selected field: from="+from+", end="+end+", length="+id.length();
		if(from==end){throw badHeader(id, "empty numeric field");}
		boolean negative=false;
		final char first=id.charAt(from);
		if(first=='-' || first=='+'){
			if(!signed){throw badHeader(id, "length fields must contain unsigned decimal digits");}
			negative=first=='-'; from++;
			if(from==end){throw badHeader(id, "sign without digits");}
		}
		final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
		final long multiplyLimit=limit/10;
		long value=0;
		for(; from<end; from++){
			final int digit=id.charAt(from)-'0';
			if(digit<0 || digit>9){throw badHeader(id, "nondecimal character in numeric field");}
			if(value<multiplyLimit){throw badHeader(id, "numeric field overflow");}
			value*=10;
			if(value<limit+digit){throw badHeader(id, "numeric field overflow");}
			value-=digit;
		}
		return negative ? value : -value;
	}

	/** Builds a diagnostic only on rejected input; valid records allocate no parser objects.
	 * @param id Original header, possibly null
	 * @param reason Rejection explanation
	 * @return Diagnostic exception containing the header and reason */
	private static IllegalArgumentException badHeader(final String id, final String reason){
		return new IllegalArgumentException("Cannot grade read header '"+id+"': "+reason+".");
	}

	/** Aggregated stream-error flags; thrown processing failures are a separate channel. */
	public boolean errorState=false;

	/** Primary input file path. */
	private String in1=null;
	/** Secondary input file path for paired reads. */
	private String in2=null;

	/** Primary output file path. */
	private String out1=null;
	/** Secondary output file path for paired reads. */
	private String out2=null;

	/** Optional input format extension override passed to FileFormat. */
	private String extin=null;
	/** Optional output format extension override passed to FileFormat. */
	private String extout=null;

	/** Path to the FASTA adapter reference, consumed during construction. */
	private String adapterFile=null;
	/** Literal adapter strings, converted with the default charset during construction. */
	private String[] literals=null;

	/** Whether existing output files may be overwritten. */
	private boolean overwrite=true;
	/** Whether output is appended to existing files. */
	private boolean append=false;

	/** Append space-/1 and space-/2 tags to the shared generated paired ID. */
	private boolean addslash=true;
	/** Dormant setting; current generation always rewrites read IDs. */
	private boolean changename=true;
	/** Simulate substitutions at adapter-covered defined positions, using Q30 for null quality. */
	private boolean adderrors=true;

	/** Share the sampled replacement position between mates, not the chosen adapter. */
	private boolean addPaired=true;
	/** Whether reverse-complement copies are appended to the adapter entries. */
	private boolean addRC=false;
	/** Replace the right (3-prime) suffix when true, or the left (5-prime) prefix otherwise. */
	private boolean right=true;

	/** Limit forwarded to CRIS; normally fragments with reader-specific units; -1 is unlimited. */
	private long maxReads=-1;
	/** Generation validity threshold; also enters grading's zero/one-base tolerance.
	 * Never filters generation output. */
	private int minlen=1;

	/** True for generation, false for length grading of surviving annotated reads. */
	private boolean writeMode=true;
	/** Probability of replacement per root read, or per mate when addPaired is false. */
	private float adapterProb=0.5f;

	/** Generation read attempts including mates, or successfully decoded surviving graded reads. */
	private long readsProcessed=0;
	/** Original read lengths in generation; observed surviving lengths in grading. */
	private long basesProcessed=0;
	/** Reads selecting an in-range replacement region, even when undefined bases are retained. */
	private long adaptersAdded=0;
	/** Positions covered by adapters, including undefined bases left unchanged. */
	private long adapterBasesAdded=0;
	/** Positions assigned random filler beyond adapters, including undefined bases left unchanged. */
	private long randomBasesAdded=0;
	/** Generated reads whose expected unmodified remainder meets minlen. */
	private long validReads=0;
	/** Sum of expected unmodified remainder lengths for validReads, not emitted read lengths. */
	private long validBases=0;

	/** Annotated positives accepted or overtrimmed by the length grader, without sequence checking. */
	private long truePos=0;
	/** Annotated negative controls accepted by the exact-length/short-read tolerance rule. */
	private long trueNeg=0;
	/** Annotated negative controls shortened beyond the accepted tolerance. */
	private long falsePos=0;
	/** Historical false-negative bucket: non-tolerated overlong results from either annotation category. */
	private long falseNeg=0;
	/** Graded reads longer than their recorded original length. */
	private long broken=0;
	/** Historical counter retained for reporting; dormant, and the active grader never increments it. */
	private long mispaired=0;

	/** Non-tolerated graded reads shorter than their expected remainder. */
	private long tooShort=0;
	/** Non-tolerated overlong graded reads, including controls expected to remain untrimmed. */
	private long tooLong=0;
	/** Graded reads at expected length, including the accepted empty/one-base short-read tolerance. */
	private long correct=0;
	/** Dormant counter; the active grader does not track missing reads. */
	private long fullyRemoved=0;

	/** Expected minus observed lengths for non-tolerated short graded reads. */
	private long tooShortBases=0;
	/** Excess observed bases above the expected remainder in non-tolerated overlong records. */
	private long tooLongBases=0;
	/** Observed bases in non-tolerated short graded reads. */
	private long tooShortReadBases=0;
	/** Historical aggregate: expected bases for overlong untrimmed controls, actual bases for other overlong reads. */
	private long tooLongReadBases=0;
	/** Expected valid bases credited to accepted reads, including short-read tolerance cases. */
	private long correctBases=0;

	/** Expected-length credit for accepted/overlong reads; observed-length credit for other short reads. */
	private long validBasesCounted=0;
	/** Expected remainder lengths summed over surviving graded records only. */
	private long validBasesExpected=0;

	/** Annotated original-minus-remaining bases among surviving positive records. */
	private long adapterBasesTotal=0;
	/** Surviving graded reads annotated as requiring trimming. */
	private long adapterReadsTotal=0;
	/** Non-tolerated overlong annotated positives; not adapter-sequence detection. */
	private long adapterReadsRemaining=0;
	/** Observed-minus-expected lengths for those overlong annotated positives. */
	private long adapterBasesRemaining=0;

	/** Primary read-input descriptor. */
	private final FileFormat ffin1;
	/** Optional mate-input descriptor. */
	private final FileFormat ffin2;

	/** Optional primary generation-output descriptor; unused by grading. */
	private final FileFormat ffout1;
	/** Optional mate-output descriptor; unused by grading. */
	private final FileFormat ffout2;

	/** Optional adapter-reference descriptor. */
	private final FileFormat ffa;

	/** Mutable adapter entries; may be null in grading mode. Generation samples entries uniformly. */
	private final ArrayList<byte[]> adapters;

	/*--------------------------------------------------------------*/

	/** Diagnostic/report stream; main closes redirects on normal return only. */
	private PrintStream outstream=System.err;
	/** Process-wide verbose switch; the constructor also sets several shared I/O switches. */
	public static boolean verbose=false;

	/** Generation random source, initialized only in write mode. */
	private shared.Random randy;

}
