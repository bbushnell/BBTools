package stream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.bam.BgzfInputStreamMT2;
import stream.bam.BgzfSettings;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/**
 * Counts sequence records and field lengths with reusable byte-scan buffers.
 * Dedicated scanners handle FASTQ, FASTA, SAM, SCARF, GFA and the supported
 * single-line FASTG layout; BAM and other sequence formats use Streamer fallbacks.
 * This is a lightweight counter with selected structural checks, not a validator
 * of every field. CLI and counting API behavior differ: only the single-parser CLI
 * calls corruption(), while FASTQ requests above one parsing thread delegate to
 * FastqScanMT and its distinct checking/counting contract.
 *
 * Instances accumulate counters and are intended for one scan, not concurrent use.
 * Static entry points also modify shared parser/decompression settings; they do
 * not isolate those settings from concurrent calls. Text bytes count input bytes
 * after decompression, excluding the synthetic newline used to finish a final line.
 * BAM and generic fallback byte/header counts follow their own rules below.
 * @author Brian Bushnell
 * @contributor Collei
 * @date November 22, 2025
 */
public final class FastqScan{

	/** Scans a path and prints time, Records, Bases, Quals and Bytes on standard output.
	 * SAM/BAM also prints Headers. The single-parser path prints detected corruption
	 * and exits with status one; FASTQ with more than one parsing thread delegates
	 * the complete command to FastqScanMT instead. Parser options affect global state.
	 * @param args Path (optionally in=), then a thread count or flag=value options retaining embedded '=' */
	public static void main(String[] args){
		Timer t=new Timer(System.out);
		if(args.length<1){throw new RuntimeException("Usage: fastqscan.sh filename");}
		String fname=args[0];
		while(fname.startsWith("-")){fname=fname.substring(1);}
		if(fname.startsWith("in=")){fname=fname.substring(3);}
		int threads=1;
		BgzfSettings.READ_THREADS=Tools.mid(1, 18, Shared.threads());
		for(int i=1; i<args.length; i++){
			String arg=args[i];
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			String b=(equals<0 || equals==arg.length()-1 ? null : arg.substring(equals+1));
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equalsIgnoreCase("simd")){Shared.SIMD&=Parse.parseBoolean(b);}//&= can DISABLE simd but never force-enable it on Java-8/non-AVX2 (Shared.SIMD is already false there) — the safe gate; no NoClassDefFoundError.
			else if(Tools.isNumeric(arg)){threads=Integer.parseInt(arg);}
			else if(Parser.parseCommonStatic(arg, a, b)){}
			else if(Parser.parseZip(arg, a, b)){}
			else{assert(false) : "Unknown parameter "+arg;}
		}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		if(threads>1 && ff.fastq()){
			FastqScanMT.main(args);
			return;
		}
		if(ff.stdin()){
			//Do nothing
		}else{
			File f=new File(fname);
			if(!f.isFile() || !f.canRead()){
				throw new RuntimeException("Can't read "+fname);
			}
		}
		final int rt=BgzfSettings.READ_THREADS=Tools.mid(1, BgzfSettings.READ_THREADS, Shared.threads());
		FastqScan fqs=new FastqScan(ff);
		try{fqs.read();}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.out.println("Records:\t"+fqs.totalRecords);
		System.out.println("Bases:  \t"+fqs.totalBases);
		System.out.println("Quals:  \t"+fqs.totalQuals);
		System.out.println("Bytes:  \t"+fqs.totalBytes);
		if(ff.samOrBam()){System.out.println("Headers:\t"+fqs.totalHeaders);}
		ByteBuilder bb=fqs.corruption();
		if(fqs.slashrLines>0){
			System.out.println("Contained Windows-style \r\n");
		}
		if(bb!=null){
			System.out.print(bb);
			System.exit(1);
		}
	}

	/** Resolves a FASTQ-default descriptor and delegates to the counting API.
	 * @param fname Input path
	 * @param halveInterleaved Probe FASTQ pairing and halve the first result when interleaved
	 * @param readThreads Above one selects FastqScanMT for FASTQ only
	 * @param zipThreads BGZF thread request; values at most one select an automatic count
	 * @return Four counters with the descriptor overload's format-specific meaning,
	 * or null if its read attempt throws IOException */
	public static long[] countReadsAndBases(String fname, boolean halveInterleaved, int readThreads, int zipThreads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countReadsAndBases(ff, halveInterleaved, readThreads, zipThreads);
	}

	/** Returns counts without calling corruption() or rejecting its structural flags.
	 * The first value is the record count divided by two only when a requested FASTQ
	 * interleaving probe detects pairs; integer division truncates an odd count.
	 * The next values are records, bases and headers. Here headers means SAM header
	 * lines or non-record GFA/FASTG lines, and is zero for FASTQ/FASTA/SCARF and fallbacks.
	 * Delegated FastqScanMT instead returns FASTQ record count as its header value.
	 * Generic fallback records count returned Read entries, not expanded mate counts.
	 *
	 * Compressed input temporarily changes BgzfSettings.READ_THREADS around read();
	 * the earlier pairing probe and construction are outside that restoration block.
	 * Other backend-global changes, including readBam's SamLine flags, are not restored.
	 * @param ff Nonnull resolved input descriptor
	 * @param halveInterleaved Probe FASTQ pairing to determine the first counter's divisor
	 * @param readThreads Above one delegates FASTQ to FastqScanMT; otherwise one parser
	 * @param zipThreads BGZF request for compressed input; at most one uses a shared-thread cap of 18
	 * @return {record count/divisor, record count, bases, format-specific headers},
	 * or null after printing an IOException from the read attempt */
	public static long[] countReadsAndBases(FileFormat ff, boolean halveInterleaved, int readThreads, int zipThreads){
		if(readThreads>1 && ff.fastq()){return FastqScanMT.countReadsAndBases(ff, halveInterleaved, readThreads, zipThreads);}
		final int oldZT=BgzfSettings.READ_THREADS;
		if(ff.compressed()){
			BgzfSettings.READ_THREADS=(zipThreads>1 ? zipThreads : Tools.mid(1, Shared.threads(), 18));
		}
		int recordsPerRead=1;
		if(ff.fastq() && halveInterleaved){
			int[] iq=FileFormat.testInterleavedAndQuality(ff.name(), false);
			recordsPerRead=(iq[1]==FileFormat.INTERLEAVED ? 2 : 1);
		}
		FastqScan fqs=new FastqScan(ff);
		try{fqs.read();}
		catch(IOException e){
			e.printStackTrace();
			//throw new RuntimeException(e);
			return null;
		}finally{BgzfSettings.READ_THREADS=oldZT;}
		long[] ret=new long[]{fqs.totalRecords/recordsPerRead, fqs.totalRecords,
			fqs.totalBases, fqs.totalHeaders};
		return ret;
	}

	/** Retains the descriptor; counters and reusable storage initialize without opening input.
	 * @param ff_ Nonnull input descriptor for the later read dispatch */
	FastqScan(FileFormat ff_){ff=ff_;}

	/** Formats detected structural problems without changing counters or flags.
	 * Reports lower bounds, not an exact count of distinct damaged records.
	 * Missing terminal newline is intentionally a CLI error, even when counts exist.
	 * @return Newly allocated diagnostic text, or null when no tracked flag is set */
	public ByteBuilder corruption(){
		//[stream/FastqScan#001 RESOLVED - Brian] missingTerminalNewline → exit 1 is INTENTIONAL and correct: sequence files are ALWAYS supposed to end in a newline. Its absence flags truncation and eases parsers, and a file lacking it is a mistake — usually someone making a fastq/fasta/sam by hand rather than programmatically. Crash-loud per BBTools philosophy. CLI-only (the countReadsAndBases API never calls this).
		if(partialRecords<1 && missingSequences<1 && !qualMismatch && !missingTerminalNewline && !missingPlus && !missingAt){
			return null;
		}
		ByteBuilder bb=new ByteBuilder();
		if(partialRecords>0 || missingAt || missingPlus || qualMismatch){
			bb.appendln("At least "+Math.max(partialRecords, 1)+" corrupt records.");
		}
		if(partialRecords>0){bb.appendln("At least "+partialRecords+" incomplete records.");}
		if(qualMismatch){bb.appendln("At least "+1+" base/quality mismatches.");}
		if(missingAt){bb.appendln("At least "+1+" missing @ symbols.");}
		if(missingPlus){bb.appendln("At least "+1+" missing + symbols.");}
		if(missingSequences>0){
			bb.append("GFA segments without stored sequence: ").append(missingSequences);
			bb.append(" (first segment record ").append(firstMissingSequence).appendln("). Review before use as sequence data.");
		}
		if(missingTerminalNewline){bb.appendln("Missing terminal newline.");}
		assert(bb.length()>0);
		return bb;
	}

	/** Adds counts from the selected format; unsupported non-sequence formats use FASTQ.
	 * Text scanners reuse a growable buffer and slide incomplete residue; BAM uses a
	 * streamer. Other sequence formats reach readOther, whose current factory rejects
	 * them before opening a reader. Does not reset state or call corruption().
	 * @throws IOException If the selected reading path propagates an I/O error */
	//TODO: Probable bug - the six text scanners ignore finishReading's boolean error
	//result (ReadWrite.finishReading), so a close-reported error does not reach the
	//caller. They also lack finally-based input closure; review separately from docs.
	void read() throws IOException{
		if(ff.fastq()){readFastq();}
		else if(ff.fasta()){readFasta();}
		else if(ff.sam()){readSam();}
		else if(ff.scarf()){readScarf();}
		else if(ff.gfa()){readGfa();}
		else if(ff.fastg()){readFastg();}
		else if(ff.bam()){readBam();}
		else if(ff.isSequence()){readOther();}
		else{readFastq();}
	}

	/** Counts complete four-line records, sequence/quality lengths and selected markers.
	 * CR before LF is removed from sequence/quality lengths. EOF residue is reported
	 * as one partial record; missing final LF is flagged after a synthetic LF is added.
	 * No alphabet, quality-value or header-identity validation is performed here.
	 * @throws IOException If reading input fails */
	void readFastq() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			//4 newlines per FASTQ record; integer-divide leaves any trailing partial record in residue, recovered next read or counted as partial at EOF (L195).
			final int records=newlines.size/4;
			totalRecords+=records;
			int recordStart=0;
			for(int i=0, j=0; i<records; i++, j+=4){
				final int headerEnd=newlines.get(j);
				final int basesEnd=newlines.get(j+1);
				final int plusEnd=newlines.get(j+2);
				final int recordEnd=newlines.get(j+3);
				int slashr1=(buffer[basesEnd-1]=='\r') ? 1 : 0;
				int slashr2=(buffer[recordEnd-1]=='\r') ? 1 : 0;
				slashrLines+=slashr1+slashr2;
				final int bases=basesEnd-headerEnd-1-slashr1;
				final int quals=recordEnd-plusEnd-1-slashr2;
				totalBases+=bases;
				totalQuals+=quals;
				bstart=recordEnd+1;
				qualMismatch|=(quals!=bases);
				missingAt|=(buffer[recordStart]!='@');
				missingPlus|=(buffer[basesEnd+1]!='+');
				recordStart=recordEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){
				//FASTQ-only deviation vs the line-based readers below: a 4-line record can be left incomplete in residue at EOF. The line formats self-complete via the missing-newline injection above, so they need no partial count here.
				if(residue>0){partialRecords++;}
				break;
			}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	/** Counts lines beginning with greater-than as records and other line lengths as bases.
	 * Excludes one CR before LF, but otherwise counts sequence-line bytes verbatim.
	 * Flags a non-header first buffered byte before any record has been counted.
	 * @throws IOException If reading input fails */
	void readFasta() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			if(totalRecords==0 && buffer[0]!='>'){partialRecords++;}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++){
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='>');
				int slashr=(lineEnd>0 && buffer[lineEnd-1]=='\r') ? 1 : 0;
				slashrLines+=slashr;
				if(header){
					totalRecords++;
				}else{
					int bases=lineEnd-bstart-slashr;
					totalBases+=bases;
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	/** Counts at-sign lines as headers and remaining lines as SAM records.
	 * Reads SEQ/QUAL lengths from tab positions, treating a leading asterisk as absent;
	 * fewer than ten tabs marks a partial record. Only nonzero QUAL lengths are compared
	 * with SEQ lengths. No alignment or optional-field parsing occurs here.
	 * @throws IOException If reading input fails */
	void readSam() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++){
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='@');
				if(header){
					totalHeaders++;
				}else{
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)'\t', symbols.clear());
					if(symbols.size>=10){
						int basesStartTab=symbols.get(8);
						int basesStopTab=symbols.get(9);
						int qualsStopSymbol=(symbols.size>10 ? symbols.get(10) : lineEnd);
						int slashr=(symbols.size==10 && buffer[qualsStopSymbol-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=(buffer[basesStartTab+1]=='*' ? 0 : basesStopTab-basesStartTab-1);
						int quals=(buffer[basesStopTab+1]=='*' ? 0 : qualsStopSymbol-basesStopTab-1-slashr);
						qualMismatch|=(quals>0 && quals!=bases);
						totalBases+=bases;
						totalQuals+=quals;
					}else{
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	/** Counts each line as a record and uses the last two colons for sequence/quality lengths.
	 * Fewer than two colons marks a partial record. Longer quality text beginning with
	 * a digit is tolerated as decimal qualities; totalQuals still counts text bytes.
	 * @throws IOException If reading input fails */
	void readScarf() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++){
				int lineEnd=newlines.array[i];
				totalRecords++;
				Vector.findSymbols(buffer, bstart, lineEnd, (byte)':', symbols.clear());
				final int colonCount=symbols.size;
				if(colonCount>=2){
					int basesStartSym=symbols.get(colonCount-2);
					int basesStopSym=symbols.get(colonCount-1);
					int qualsStopSym=lineEnd;
					int slashr=(buffer[qualsStopSym-1]=='\r') ? 1 : 0;
					slashrLines+=slashr;
					int bases=basesStopSym-basesStartSym-1;
					int quals=qualsStopSym-basesStopSym-1-slashr;
					qualMismatch|=(quals<bases ||
						(quals>bases && !Tools.isDigit(buffer[basesStopSym+1])));// Quals can be decimal
					totalBases+=bases;
					totalQuals+=quals;
				}else{
					partialRecords++;
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	/** Counts lines beginning with S as records and all other lines as headers.
	 * Counts stored third-field sequence bytes, ignoring optional fields including LN.
	 * An exact asterisk placeholder or empty sequence contributes zero bases and is
	 * flagged for review by corruption(); declared lengths never supply missing bases.
	 * Fewer than two tabs marks a partial record. No graph or alphabet validation occurs.
	 * @throws IOException If reading input fails */
	void readGfa() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++){
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]!='S');
				if(header){
					totalHeaders++;
				}else{
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)'\t', symbols.clear());
					final int size=symbols.size;
					if(size>=2){
						int basesStartSym=symbols.get(1);
						int basesStopSym=(size>2 ? symbols.get(2) : lineEnd);
						int slashr=(size==2 && buffer[basesStopSym-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=basesStopSym-basesStartSym-1-slashr;
						//STR384: count actual sequence, never trust LN to supply absent bases.
						if(bases==0 || (bases==1 && buffer[basesStartSym+1]=='*')){
							bases=0;
							if(missingSequences==0){firstMissingSequence=totalRecords;}
							missingSequences++;
						}
						totalBases+=bases;
					}else{
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

//	>NODE_1:NODE_2\'; ATGCGTACGTTAG
//	>NODE_1\'; CTAACGTACGCAT
	/** Counts greater-than lines in the supported single-line header/sequence layout.
	 * Sequence length is measured after the final semicolon and one separator byte;
	 * other lines count as headers, not continued sequence. Missing semicolons mark
	 * partial records. This path does not implement a general multiline FASTG parser.
	 * @throws IOException If reading input fails */
	void readFastg() throws IOException{
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			totalBytes+=r;
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++){
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]!='>');
				if(header){
					totalHeaders++;
				}else{
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)';', symbols.clear());
					final int size=symbols.size;
					if(size>=1){
						//Assumes the FASTG "...'; SEQ" separator includes the space: basesStartSym lands ON the space, and the -1 in the bases calc below skips it. A separator lacking the space undercounts by 1 (format-trusting; FASTG writers emit "; ").
						int basesStartSym=symbols.get(size-1)+1;
						int basesStopSym=lineEnd;
						int slashr=(buffer[basesStopSym-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=basesStopSym-basesStartSym-1-slashr;
						totalBases+=bases;
					}else{
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	/** Counts returned SamLine entries using a factory-selected BAM input streamer.
	 * Disables selected global SamLine parsing fields and FLIP_ON_LOAD without restoring
	 * them. Header count stays unchanged; bytes come from exact-class BamStreamer only.
	 * The streamer is closed after normal consumption, without checking its error flag.
	 * @throws IOException Declared for the common dispatcher contract */
	void readBam() throws IOException{
		SamLine.PARSE_0=SamLine.PARSE_2=SamLine.PARSE_5=SamLine.PARSE_6=false;
		SamLine.PARSE_7=SamLine.PARSE_8=SamLine.PARSE_OPTIONAL=false;
		SamLine.FLIP_ON_LOAD=false;
		Streamer st=StreamerFactory.makeStreamer(ff, 0, false, -1, false, false, -1);
		st.start();
		for(ListNum<SamLine> ln=st.nextLines(); ln!=null && !ln.poison(); ln=st.nextLines()){
			for(SamLine sl : ln){
				int bases=sl.seq==null ? 0 : sl.seq.length;
				int quals=sl.qual==null ? 0 : sl.qual.length;
				totalRecords++;
				totalBases+=bases;
				totalQuals+=quals;
				qualMismatch|=(quals>0 && quals!=bases);
			}
		}
		//TODO: Probable bug - neither this fallback nor readOther checks st.errorState()
		//after terminal consumption; a backend-reported error can leave partial counts.
		st.close();
		if(st.getClass()==BamStreamer.class){
			totalBytes=((BamStreamer)st).bytesProcessed();
		}
	}

	/** Attempts to count Read entries through the no-QUAL streamer factory.
	 * The ordinary read() dispatcher handles all currently supported factory formats
	 * before this fallback; other formats are rejected before a reader is created.
	 * Direct package-level calls with supported formats can still use this helper.
	 * Counts each root's bases/qualities and countFastqBytes estimate without expanding
	 * mates; header count stays unchanged. Unlike the BAM path, this does not close
	 * the streamer after consumption or check its error flag.
	 * @throws IOException Declared for the common dispatcher contract */
	void readOther() throws IOException{
		//TODO: STR275 - supported direct calls or future factory extensions need explicit
		//cleanup review: this helper has no close call. Current read() dispatch cannot open a reader here.
		Streamer st=StreamerFactory.makeStreamer(ff, 0, false, -1, false, true, -1);
		st.start();
		for(ListNum<Read> ln=st.nextList(); ln!=null && !ln.poison(); ln=st.nextList()){
			for(Read r : ln){
				int bases=r.length();
				int quals=r.quality==null ? 0 : r.quality.length;
				totalRecords++;
				totalBases+=bases;
				totalQuals+=quals;
				totalBytes+=r.countFastqBytes();
				qualMismatch|=(quals>0 && quals!=bases);
			}
		}
	}

	/** Grows the reusable text buffer up to MAX_ARRAY_LEN, preserving buffered bytes.
	 * Called for a full unconsumed residue or to append a synthetic EOF newline.
	 * With assertions enabled, rejects a buffer that can no longer grow. */
	//One reusable buffer plus residue sliding avoids per-record allocation in text scans.
	private void expand(){
		long newlen=Math.min(buffer.length*2L, Shared.MAX_ARRAY_LEN);
		assert(newlen>buffer.length) : "Record "+totalRecords+" is too long.";
		buffer=Arrays.copyOf(buffer, (int)newlen);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained input descriptor; read() selects the scan path. */
	private final FileFormat ff;

	/** Initial buffer allocation size; not updated when the buffer expands. */
	private int bufferLen=262144;
	/** Reusable text scan buffer, expanded when a record/line does not fit. */
	private byte[] buffer=new byte[bufferLen];
	/** SAM header lines or non-record GFA/FASTG lines; other paths do not increment it. */
	long totalHeaders;
	/** Format-specific records or returned fallback entries, without mate expansion. */
	long totalRecords;
	/** Sequence-field lengths or fallback sequence lengths. */
	long totalBases;
	/** Quality-field text lengths or fallback quality-array lengths. */
	long totalQuals;
	/** Decompressed text bytes, BamStreamer byte count or generic FASTQ-size estimates. */
	long totalBytes=0;

	/** Lower-bound count of incomplete/malformed record structures. */
	long partialRecords;
	/** GFA segment records with an exact asterisk placeholder or empty sequence field. */
	long missingSequences;
	/** One-based position among segment records of the first missing GFA sequence; zero if none. */
	long firstMissingSequence;
	/** CR-before-LF occurrences at positions checked by each format, not all file lines. */
	long slashrLines;
	/** A format-specific sequence/quality-length mismatch was observed. */
	boolean qualMismatch;
	/** A text scan inserted a final LF that was absent from the input. */
	boolean missingTerminalNewline;
	/** At least one FASTQ separator line lacked its initial plus sign. */
	boolean missingPlus;
	/** At least one FASTQ record lacked its initial at sign. */
	boolean missingAt;
}
