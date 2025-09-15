package prok;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Stores frame-relative kmer counts for a type of genomic feature, such as a coding start site.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 */
public class FrameStats {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public FrameStats(String name_, int k_, int frames_, int leftOffset_){
		name=name_;
		k=k_;
		mask=~((-1)<<(2*k));
		frames=frames_;
		kMax=1<<(2*k);
		invFrames=1.0f/frames;
		leftOffset=leftOffset_;
		
		probs=new float[frames][kMax];
		countsTrue=new long[frames][kMax];
		countsFalse=new long[frames][kMax];
		counts=new long[][][] {countsFalse, countsTrue};
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(int kmer, int frame, int valid){
		counts[valid][frame][kmer]++;
		validSums[valid]++;
	}
	
	public boolean compatibleWith(FrameStats fs) {
		return fs.name.equals(name) && fs.leftOffset==leftOffset && fs.k==k && fs.frames==frames;
	}
	
	public void clear() {
		Tools.fill(counts, 0);
		Arrays.fill(validSums, 0);
	}
	
	public void setFrom(FrameStats fs) {
		assert(compatibleWith(fs)) : name+", "+frames+", "+fs.name+", "+fs.frames;
		clear();
		add(fs);
	}
	
	public void add(FrameStats fs){
		assert(fs.name.equals(name));
		assert(fs.leftOffset==leftOffset);
		assert(fs.k==k);
		assert(fs.frames==frames) : name+", "+frames+", "+fs.name+", "+fs.frames;
//		for(int x=0; x<counts.length; x++) {
//			for(int y=0; y<counts[x].length; y++) {
//				for(int z=0; z<counts[x][y].length; z++) {
//					assert(fs.counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z]+", "+fs.name;
//					assert(counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//				}
//			}
//		}

		Tools.add(counts, fs.counts);
		Tools.add(validSums, fs.validSums);
//		for(int x=0; x<counts.length; x++) {
//			for(int y=0; y<counts[x].length; y++) {
//				for(int z=0; z<counts[x][y].length; z++) {
//					assert(fs.counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//					assert(counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//				}
//			}
//		}
//		calculate();
	}
	
	public void multiplyBy(double mult) {
		Tools.multiplyBy(counts, mult);
		Tools.multiplyBy(validSums, mult);
	}
	
	void calculate(){
		average=(float)((validSums[1]+1.0)/(validSums[0]+validSums[1]+1.0));
		invAvg=1.0f/average;
		
		for(int a=0; a<frames; a++){
			for(int b=0; b<kMax; b++){
				long t=countsTrue[a][b];
				long f=countsFalse[a][b];
				probs[a][b]=(float)(t/(t+f+1.0))*invAvg;
			}
		}
	}
	
	public float scorePoint(int point, byte[] bases){
		final int mask=~((-1)<<(2*k));
		
		int kmer=0;
		int len=0;
		float score=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-leftOffset;
		for(int i=start, frame=0-k+1; i<bases.length && frame<frames; i++, frame++){
			byte b=(i>=0 ? bases[i] : (byte)'A');
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=probs[frame][kmer];
					float dif=prob-0.99f;
					score+=dif;
					
//					if(name.equals("startStats")){
//						System.err.println("frame="+frame+" kmer="+AminoAcid.kmerToString(kmer, k)+
//								Tools.format(" prob=%.4f\tdif=%.4f\tscore=%.4f", prob, dif, score)+
//								"\tvalid="+counts[1][frame][kmer]+"\tinvalid="+counts[0][frame][kmer]);
//					}
				}
			}else{len=0;}
		}
		
		return score*invFrames;
	}
	
	void processCDSFrames(byte[] bases, byte[] validFrames){
		if(!ProkObject.callCDS){return;}
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int vf=validFrames[i];
					for(int frame=0; frame<frames; frame++){
						int valid=vf&1;
						add(kmer, frame, valid);
						//For CDS start (0-based) of 189, i=192, j=189, vf=1, frame=0 - all as expected.
//						assert(valid==0) : "vf="+vf+", frame="+frame+", len="+len+", kmer="+AminoAcid.kmerToString(kmer, k)+", i="+i+", j="+j;
						vf=(vf>>1);
					}
				}
			}else{len=0;}
		}
	}
	
	void processPoint(byte[] bases, int point, int valid){
		
		//Degenerate cases where the point is at the end, possibly from a truncated gene.
		if(point<3){return;}
		if(point>=bases.length-3){return;}
		
		int kmer=0;
		int len=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-leftOffset;
		
		int i=start, frame=0-k+1;
		while(i<0){i++; frame++;}
		for(; i<bases.length && frame<frames; i++, frame++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					add(kmer, frame, valid);
				}
			}else{len=0;}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Text Methods         ----------------*/
	/*--------------------------------------------------------------*/

	public void parseData(byte[] line) {
		int a=0, b=0;
		final int valid, frame;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		valid=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		frame=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		assert(valid==0 || valid==1);
		assert(frame>=0 && frame<frames);
		try {
			final long[] row=counts[valid][frame];
			long sum=0;
			for(int kmer=0; kmer<row.length; kmer++){
				while(b<line.length && line[b]!='\t'){b++;}
				assert(b>a) : "Missing field 1: "+new String(line);
				long count=Parse.parseLong(line, a, b);
				b++;
				a=b;
				row[kmer]=count;
				sum+=count;
			}
			validSums[valid]+=sum;
		} catch (Exception e) {
			System.err.println(new String(line)+"\n"+name);
			assert(false) : e;
		}
	}
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append("#name\t").append(name).nl();
		bb.append("#k\t").append(k).nl();
		bb.append("#frames\t").append(frames).nl();
		bb.append("#offset\t").append(leftOffset).nl();
		bb.append("#valid\tframe");
		for(int i=0; i<kMax; i++){bb.tab().append(AminoAcid.kmerToString(i, k));}
		bb.nl();
		for(int a=0; a<2; a++){//valid
			for(int b=0; b<frames; b++){
				bb.append(a);
				bb.tab().append(b);
				for(int c=0; c<kMax; c++){
					bb.tab().append(counts[a][b][c]);
				}
				bb.nl();
			}
		}
		return bb;
	}
	
	public ByteBuilder append0(ByteBuilder bb){
		bb.append("#name\t").append(name).nl();
		bb.append("#k\t").append(k).nl();
		bb.append("#frames\t").append(frames).nl();
		bb.append("#offset\t").append(leftOffset).nl();
		bb.append("#valid\tframe");
		for(int i=0; i<kMax; i++){bb.tab().append(AminoAcid.kmerToString(i, k));}
		bb.nl();
		for(int a=0; a<2; a++){//valid
			for(int b=0; b<frames; b++){
				bb.append(a);
				bb.tab().append(b);
				for(int c=0; c<kMax; c++){
					bb.tab().append(0);
				}
				bb.nl();
			}
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final String name;
	public final int k;
	public final int mask;
	public final int frames;
	public final int kMax;
	public final float invFrames;
	public final int leftOffset;
	public int rightOffset() {return frames-leftOffset-1;}
	
	public final float[][] probs;
	public final long[][] countsTrue;
	public final long[][] countsFalse;
	public final long[][][] counts;

	public final long[] validSums=KillSwitch.allocLong1D(2);
	private float average=-1;
	private float invAvg=-1;
	
}
