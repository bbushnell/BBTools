package jgi;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;
import structures.IntList;

public class Assembly {
	
	public Assembly(String fname_) {
		fname=fname_;
		load();
	}
	
	void load() {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		assert(ff.fasta());
		ByteFile bf=ByteFile.makeByteFile(ff);

		clear();
		acgtnio=new long[7];
		int contigLen=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line[0]=='>') {
				headerLength+=(line.length-1);
				if(firstHeader==null && line.length>1) {
					firstHeader=new String(line, 1, line.length-1);
				}
				if(contigLen>0) {
					contigs.add(contigLen);
				}
				contigLen=0;
			}else {
				addToACGTNIO(line);
				contigLen+=line.length;
			}
		}
		bf.close();
		contigs.sort();
		contigs.reverse();
		length=contigs.sumLong();
	}
	
	void addToACGTNIO(byte[] line) {
		for(byte b : line) {
			byte x=baseToACGTNIO[b];
			acgtnio[x]++;
		}
	}
	
	void clear() {
		contigs.clear();
		length=0;
		headerLength=0;
		firstHeader=null;
		acgtnio=null;
	}
	
	float gc() {
		float AT=acgtnio[A]+acgtnio[T]+acgtnio[U];
		float GC=acgtnio[G]+acgtnio[C];
		return GC/Tools.max(1, AT+GC);
	}
	
	long lengthAtLeast(int minimum) {
		long sum=0;
		for(int i=0; i<contigs.size; i++) {
			int len=contigs.get(i);
			if(len<minimum) {break;}
			sum+=len;
		}
		return sum;
	}
	
	final String fname;
	IntList contigs=new IntList();
	long length=0;
	long headerLength=0;
	String firstHeader=null;
	long[] acgtnio;
	
	public static final byte[] baseToACGTNIO=makeBaseToACGTUNIO();
	private static final byte A=0, C=1, G=2, T=3, U=4, N=5, IUPAC=6, OTHER=7;
	
	private static final byte[] makeBaseToACGTUNIO() {
		final byte[] array=new byte[128];
		Arrays.fill(array, OTHER);
		array['a']=array['A']=A;
		array['c']=array['C']=C;
		array['g']=array['G']=G;
		array['t']=array['T']=T;
		array['u']=array['U']=U;
		array['n']=array['N']=N;
		for(int i=0; i<array.length; i++) {
			if(AminoAcid.baseToNumberExtended[i]>=0 && array[i]==OTHER) {
				array[i]=IUPAC;
			}
		}
		return array;
	}
	
}
