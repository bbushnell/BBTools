package stream;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Counts reads and bases in a sequence file, with low overhead.
 * @author Brian Bushnell
 * @contributor Collei
 * @date November 22, 2025
 */
public class FastqScan{

	public static void main(String[] args) {
		Timer t=new Timer();
		if(args.length!=1) {
			System.err.println("Usage: fastqscan.sh filename");
		}
		String fname=args[0];
		while(fname.startsWith("-")) {fname=fname.substring(1);}
		if(fname.startsWith("in=")) {fname=fname.substring(3);}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		FastqScan fqs=new FastqScan(ff);
		try{fqs.read();}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.err.println("Records:\t"+fqs.totalRecords);
		System.err.println("Bases:  \t"+fqs.totalBases);
		if(ff.samOrBam()) {System.err.println("Headers:\t"+fqs.totalHeaders);}
		ByteBuilder bb=fqs.corruption();
		if(fqs.slashrLines>0) {
			System.err.println("Contained Windows-style \r\n");
		}
		if(bb!=null) {
			System.err.print(bb);
			System.exit(1);
		}
	}

	public static long[] countReadsAndBases(String fname, boolean halveInterleaved) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countReadsAndBases(ff, halveInterleaved);
	}

	/** Returns molecules, reads, bases, file headers */
	public static long[] countReadsAndBases(FileFormat ff, boolean halveInterleaved) {
		int recordsPerRead=1;
		if(ff.fastq() && halveInterleaved) {
			int[] iq=FileFormat.testInterleavedAndQuality(ff.name(), false);
			recordsPerRead=(iq[1]==FileFormat.INTERLEAVED ? 2 : 1);
		}
		FastqScan fqs=new FastqScan(ff);
		try{fqs.read();}
		catch(IOException e){
			e.printStackTrace();
			//throw new RuntimeException(e);
			return null;
		}
		long[] ret=new long[] {fqs.totalRecords/recordsPerRead, fqs.totalRecords, 
			fqs.totalBases, fqs.totalRecords};
		return ret;
	}

	FastqScan(FileFormat ff_) {ff=ff_;}
	
	public ByteBuilder corruption() {
		if(partialRecords<1 && !qualMismatch && !missingTerminalNewline && !missingPlus && !missingAt) {
			return null;
		}
		ByteBuilder bb=new ByteBuilder();
		if(partialRecords>0 || missingAt || missingPlus || qualMismatch) {
			bb.appendln("At least "+Math.max(partialRecords, 1)+" corrupt records.");
		}
		if(partialRecords>0) {bb.appendln("At least "+partialRecords+" incomplete records.");}
		if(qualMismatch) {bb.appendln("At least "+1+" base/quality mismatches.");}
		if(missingAt) {bb.appendln("At least "+1+" missing @ symbols.");}
		if(missingPlus) {bb.appendln("At least "+1+" missing + symbols.");}
		if(missingTerminalNewline) {bb.appendln("Missing terminal newline.");}
		assert(bb.length()>0);
		return bb;
	}

	void read() throws IOException {
		if(ff.fastq()) {readFastq();}
		else if(ff.fasta()) {readFasta();}
		else if(ff.sam()) {readSam();}
		else if(ff.scarf()) {readScarf();}
		else if(ff.gfa()) {readGfa();}
		else {readFastq();}
	}

	void readFastq() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			final int records=newlines.size/4;
			totalRecords+=records;
			int recordStart=0;
			for(int i=0, j=0; i<records; i++, j+=4) {
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
				bstart=recordEnd+1;
				qualMismatch|=(quals!=bases);
				missingAt|=(buffer[recordStart]!='@');
				missingPlus|=(buffer[basesEnd+1]!='+');
				recordStart=recordEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {
				if(residue>0) {partialRecords++;}
				break;
			}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	void readFasta() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			if(totalRecords==0 && buffer[0]!='>') {partialRecords++;}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='>');
				int slashr=(lineEnd>0 && buffer[lineEnd-1]=='\r') ? 1 : 0;
				slashrLines+=slashr;
				if(header) {
					totalRecords++;
				}else {
					int bases=lineEnd-bstart-slashr;
					totalBases+=bases;
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	void readSam() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='@');
				if(header) {
					totalHeaders++;
				}else {
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)'\t', symbols.clear());
					if(symbols.size>=10) {
						int basesStartTab=symbols.get(8);
						int basesStopTab=symbols.get(9);
						int qualsStopSymbol=(symbols.size>10 ? symbols.get(10) : lineEnd);
						int slashr=(symbols.size==10 && buffer[qualsStopSymbol-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=(buffer[basesStartTab+1]=='*' ? 0 : basesStopTab-basesStartTab-1);
						int quals=(buffer[basesStopTab+1]=='*' ? 0 : qualsStopSymbol-basesStopTab-1-slashr);
						qualMismatch|=(quals>0 && quals!=bases);
						totalBases+=bases;
					}else {
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	void readScarf() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				totalRecords++;
				Vector.findSymbols(buffer, bstart, lineEnd, (byte)':', symbols.clear());
				final int colonCount=symbols.size;
				if(colonCount>=2) {
					int basesStartSym=symbols.get(colonCount-2);
					int basesStopSym=symbols.get(colonCount-1);
					int qualsStopSym=lineEnd;
					int slashr=(buffer[qualsStopSym-1]=='\r') ? 1 : 0;
					slashrLines+=slashr;
					int bases=basesStopSym-basesStartSym-1;
					int quals=qualsStopSym-basesStopSym-1-slashr;
					qualMismatch|=(quals<bases || 
						(quals>bases && !Tools.isDigit(buffer[basesStopSym+1]))); //Quals can be decimal
					totalBases+=bases;
				}else {
					partialRecords++;
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}
	
	void readGfa() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]!='S');
				if(header) {
					totalHeaders++;
				}else {
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)'\t', symbols.clear());
					final int size=symbols.size;
					if(size>=2) {
						int basesStartSym=symbols.get(1);
						int basesStopSym=(size>2 ? symbols.get(2) : lineEnd);
						int slashr=(size==2 && buffer[basesStopSym-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=basesStopSym-basesStartSym-1-slashr;
						totalBases+=bases;
					}else {
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}
	
//	>NODE_1:NODE_2\'; ATGCGTACGTTAG
//	>NODE_1\'; CTAACGTACGCAT
	void readFastg() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList symbols=new IntList(128);
		int bstop=0, bstart=0;
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]!='>');
				if(header) {
					totalHeaders++;
				}else {
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)';', symbols.clear());
					final int size=symbols.size;
					if(size>=1) {
						int basesStartSym=symbols.get(size-1)+1;
						int basesStopSym=lineEnd;
						int slashr=(buffer[basesStopSym-1]=='\r') ? 1 : 0;
						slashrLines+=slashr;
						int bases=basesStopSym-basesStartSym-1-slashr;
						totalBases+=bases;
					}else {
						partialRecords++;
					}
				}
				bstart=lineEnd+1;
			}

			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {break;}
		}
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}

	private void expand() {
		long newlen=Math.min(buffer.length*2L, Shared.MAX_ARRAY_LEN);
		assert(newlen>buffer.length) : "Record "+totalRecords+" is too long.";
		buffer=Arrays.copyOf(buffer, (int)newlen);
	}

	private final FileFormat ff;
	
	private int bufferLen=262144;
	private byte[] buffer=new byte[bufferLen];
	long totalHeaders;
	long totalRecords;
	long totalBases;
	
	long partialRecords;
	long slashrLines;
	boolean qualMismatch;
	boolean missingTerminalNewline;
	boolean missingPlus;
	boolean missingAt;
}
