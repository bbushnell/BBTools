package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import hiseq.IlluminaHeaderParser2;
import shared.Shared;
import structures.ByteBuilder;

public class FastqReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		FastqReadInputStream fris=new FastqReadInputStream(args[0], true);
		
		Read r=fris.next();
		System.out.println(r.toText(false));
		
	}
	
	public FastqReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}
	
	public FastqReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		shrinkHeaders=FASTQ.SHRINK_HEADERS;
		if(!ff.fastq()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		
		if(FASTQ.PARSE_CUSTOM){
			try {
				String s[]=ff.name().split("_");
//				maxSnps=toNumber(s[3]);
//				maxInss=toNumber(s[4]);
//				maxDels=toNumber(s[5]);
//				maxSubs=toNumber(s[6]);
				
//				s=s[8].split("\\.");
//
//				s=s[0].split("-");
				
//				if(s.length!=8 && s.length!=9){
//					if(Shared.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
//				}
				
//				minChrom=Gene.toChromosome(s[0]);
//				maxChrom=Gene.toChromosome(s[1]);

			} catch (Exception e) {
				// TODO Auto-generated catch block
				//			e.printStackTrace();
//				if(Shared.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
			}
		}
//		interleaved=false;
		interleaved=ff.interleaved();//(ff.stdio()) ? FASTQ.FORCE_INTERLEAVED : FASTQ.isInterleaved(ff.name(), false);
		tf=ByteFile.makeByteFile(ff);
//		assert(false) : interleaved;
	}

	@Override
	public void start() {
//		if(cris!=null){cris.start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	@Override
	public Read next() {
		if(!hasMore()){return null;}
		Read r=buffer.set(next, null);
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=FASTQ.toReadList(tf, BUF_LEN, nextReadID, interleaved, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}else if(shrinkHeaders) {//Kind of slow to regenerate strings but oh well
			for(Read r : buffer) {shrinkHeader(r, r.mate);}
		}
	}
	
	private void shrinkHeader(Read r1, Read r2) {
		ihp.parse(r1.id);
		if(!ihp.canShrink()) {return;}
		bbh.clear().colon().colon().colon();
		ihp.appendCoordinates(bbh).space().append(1).colon();
		r1.id=bbh.toString();
		if(r2!=null) {
			bbh.set(bbh.length-2, (byte)'2');
			r2.id=bbh.toString();
		}
	}
	
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	@Override
	public boolean paired() {return interleaved;}
	
	@Override
	public String fname(){return tf.name();}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}

	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final ByteFile tf;
	private final boolean interleaved;
	private final int flag;
	private boolean shrinkHeaders;
	
	private final ByteBuilder bbh=new ByteBuilder(128);
	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

	private final int BUF_LEN=Shared.bufferLen();;
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean stdin;
	public static boolean verbose=false;

}
