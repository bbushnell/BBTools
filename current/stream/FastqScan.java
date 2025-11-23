package stream;

import java.io.IOException;
import java.io.InputStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Vector;
import structures.IntList;

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
	
	FastqScan(FileFormat ff_) {ff=ff_;}
	
	void read() throws IOException {
		if(ff.fasta()) {readFasta();}
		else {readFastq();}
	}
	
	void readFastq() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, residue=0, bstart=0;
		for(int r=is.read(buffer); r>0 || residue>0; r=is.read(buffer, residue, buffer.length-residue)) {
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
				int bases=basesEnd-headerEnd-1;
				totalBases+=bases;
				bstart=recordEnd+1;
			}
			
			residue=bstop-bstart;
			if(residue>0) {
				System.arraycopy(buffer, bstart, buffer, 0, residue);
			}
			bstart=0;
			bstop=residue;
			newlines.clear();
		}
		is.close();
	}
	
	void readFasta() throws IOException {
		InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		int bstop=0, residue=0, bstart=0;
		for(int r=is.read(buffer); r>0 || residue>0; r=is.read(buffer, residue, buffer.length-residue)) {
			assert(bstop==residue);
			assert(bstart==0);
			bstop+=r;
//			if(r<1 && buffer[bstop-1]!='\n') {buffer[bstop-1]='\n'; bstop++;}//Files without ending newline
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			int lines=newlines.size;
			for(int i=0; i<lines; i++) {
				int lineEnd=newlines.array[i];
				boolean header=(bstart>0 && buffer[bstart-1]=='>');
				if(header) {
					totalRecords++;
				}else {
					int bases=lineEnd-bstart+1;
					totalBases+=bases;
				}
				bstart=lineEnd+1;
			}
			
			residue=bstop-bstart;
			if(residue>0) {
				System.arraycopy(buffer, bstart, buffer, 0, residue);
			}
			bstart=0;
			bstop=residue;
			newlines.clear();
		}
		is.close();
	}
	
	//TODO - buffer expansion logic 
	
	private final int buflen=262144;//65536
	private final byte[] buffer=new byte[262144];
	private final FileFormat ff;
	long totalRecords;
	long totalBases;
}
