package stream;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Vector;
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
		String fname=args[0];
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		FastqScan fqs=new FastqScan(ff);
		try{fqs.read();}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:");
		System.err.println("Records:\t"+fqs.totalRecords);
		System.err.println("Bases:  \t"+fqs.totalBases);
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

	void read() throws IOException {
		if(ff.fasta()) {readFasta();}
		else if(ff.sam()) {readSam();}
		else {readFastq();}
	}

	void readFastq() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, residue=0, bstart=0;
		for(int r=is.read(buffer); r>0; r=is.read(buffer, residue, buffer.length-residue)) {
			assert(bstop==residue);
			assert(bstart==0);
			bstop+=r;
			//			if(r<1 && buffer[bstop-1]!='\n') {buffer[bstop-1]='\n'; bstop++;}//Files without ending newline
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			int records=newlines.size/4;
			totalRecords+=records;
			for(int i=0, j=0; i<records; i++, j+=4) {
				int headerEnd=newlines.get(j);
				int basesEnd=newlines.get(j+1);
				int recordEnd=newlines.get(j+3);
				int bases=basesEnd-headerEnd-1-(buffer[basesEnd-1]=='\r' ? 1 : 0);
				totalBases+=bases;
				bstart=recordEnd+1;
			}

			residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			newlines.clear();
		}

		if(residue>0) {//For files missing terminal newlines
			newlines.clear();
			// Scan the residue for internal newlines
			Vector.findSymbols(buffer, 0, residue, (byte)'\n', newlines);

			// A valid FASTQ record without a trailing newline will have exactly 3 newlines
			// (Header\nSeq\n+\nQual[EOF])
			if(newlines.size==3) {
				int headerEnd=newlines.get(0);
				int basesEnd=newlines.get(1);

				// Calculate bases length
				int bases=basesEnd-headerEnd-1;
				if(bases>0 && buffer[basesEnd-1]=='\r') {bases--;}

				totalBases += bases;
				totalRecords++;
			}
		}
		is.close();
	}

	void readFasta() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, residue=0, bstart=0;
		for(int r=is.read(buffer); r>0; r=is.read(buffer, residue, buffer.length-residue)) {
			assert(bstop==residue);
			assert(bstart==0);
			bstop+=r;
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='>');
				if(header) {
					totalRecords++;
				}else {
					int bases=lineEnd-bstart-(buffer[lineEnd-1]=='\r' ? 1 : 0);
					totalBases+=bases;
				}
				bstart=lineEnd+1;
			}

			residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			newlines.clear();
		}
		if(residue>0) {//For files missing terminal newlines
			// Verify it's not just a trailing newline that became residue
			if(buffer[0]=='>') {
				totalRecords++;
			} else {
				int bases=residue;
				// Handle Windows \r at EOF
				if(bases>0 && buffer[bases-1]=='\r') {bases--;}
				totalBases += bases;
			}
		}
		is.close();
	}

	void readSam() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		IntList tabs=new IntList(128);
		int bstop=0, residue=0, bstart=0;
		for(int r=is.read(buffer); r>0; r=is.read(buffer, residue, buffer.length-residue)) {
			assert(bstop==residue);
			assert(bstart==0);
			bstop+=r;
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(buffer[bstart]=='@');
				if(header) {
					totalHeaders++;
				}else {
					totalRecords++;
					Vector.findSymbols(buffer, bstart, lineEnd, (byte)'\t', tabs.clear());
					int basesStart=tabs.get(8);
					int basesStop=tabs.get(9);
					int bases=basesStop-basesStart-1-(buffer[lineEnd-1]=='\r' ? 1 : 0);
					totalBases+=bases;
				}
				bstart=lineEnd+1;
			}

			residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			newlines.clear();
		}

		if(residue>0) {//For files missing terminal newlines
			if(buffer[0]=='@') {
				totalHeaders++;
			} else {
				totalRecords++;
				// We must find the tabs in the residue to get the sequence length
				tabs.clear();
				Vector.findSymbols(buffer, 0, residue, (byte)'\t', tabs);

				// SAM format: Sequence is the 10th column (index 9)
				// It is located between the 9th tab (index 8) and 10th tab (index 9)
				if(tabs.size>=9) {
					int basesStart=tabs.get(8)+1;
					// If there is no 10th tab (EOF in sequence), use residue length
					int basesStop=(tabs.size>9) ? tabs.get(9) : residue;

					int bases=basesStop-basesStart;
					// Handle Windows \r at EOF
					if(basesStop==residue && bases>0 && buffer[bases-1]=='\r') {bases--;}

					totalBases += Math.max(0, bases);
				}
			}
		}
		is.close();
	}

	private void expand() {
		long newlen=Math.min(buffer.length*2L, Shared.MAX_ARRAY_LEN);
		assert(newlen>buffer.length) : "Record "+totalRecords+" is too long.";
		buffer=Arrays.copyOf(buffer, (int)newlen);
	}

	private byte[] buffer=new byte[262144];
	private final FileFormat ff;
	long totalHeaders;
	long totalRecords;
	long totalBases;
}
