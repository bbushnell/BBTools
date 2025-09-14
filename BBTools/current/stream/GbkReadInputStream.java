package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class GbkReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		GbkReadInputStream fris=new GbkReadInputStream(args[0], true);
		
		Read r=fris.next();
		System.out.println(r.toText(false));
		
	}
	
	public GbkReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.GBK, null, allowSubprocess_, false));
	}
	
	public GbkReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.gbk()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		bf=ByteFile.makeByteFile(ff);
//		assert(false) : interleaved;
	}

	@Override
	public void start() {
//		if(cris!=null){cris.start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(bf.isOpen()){
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
		buffer=toReadList(bf, BUF_LEN, nextReadID, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){bf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}
	}
	

	
	public static ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn, long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		
		int added=0;
		
		String idLine=null;
		ByteBuilder bb=new ByteBuilder();
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
//			if(Tools.startsWith(s, "ID")){
//				idLine=new String(s, 2, s.length-2).trim();
//			}else 
			if(Tools.startsWith(s, "ORIGIN")){
//				System.err.println(new String(s));
				byte[] line=null;
				for(line=bf.nextLine(); line!=null && line[0]!='/'; line=bf.nextLine()){
					for(byte b : line){
						if(Tools.isLetter(b)){
							bb.append(Tools.toUpperCase(b));
						}
					}
				}
				assert(line==null || Tools.startsWith(line, "//")) : new String(line);
				
				Read r=new Read(bb.toBytes(), null, idLine==null ? ""+numericID : idLine, numericID, flag);
				list.add(r);
				added++;
				numericID++;
				
				bb.clear();
				idLine=null;
				
				if(added>=maxReadsToReturn){break;}
			}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		errorState|=bf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		bf.reset();
	}

	@Override
	public boolean paired() {return false;}
	
	@Override
	public String fname(){return bf.name();}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}

	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final ByteFile bf;
	private final int flag;

	private final int BUF_LEN=Shared.bufferLen();;
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean stdin;
	public static boolean verbose=false;

}
