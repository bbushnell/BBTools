package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class FastaQualReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		FastaQualReadInputStream fris=new FastaQualReadInputStream(args[0], args[1], true);
		
		Read r=fris.next();
		int i=0;
		while(r!=null){
			System.out.println(r.toText(false));
			r=fris.next();
			if(i++>3){break;}
		}
		
	}
	
	public FastaQualReadInputStream(String fname, String qfname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false), qfname);
	}
	
	public FastaQualReadInputStream(FileFormat ff, String qfname){
		
		if(!ff.fasta() && !ff.stdio()){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+ff.name());
		}
		
		btf=ByteFile.makeByteFile(ff);
		qtf=ByteFile.makeByteFile(FileFormat.testInput(qfname, FileFormat.QUAL, null, ff.allowSubprocess(), false));
		interleaved=false;
		
	}

	@Override
	public void start() {
//		if(cris!=null){cris.start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(btf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	@Override
	public Read next() {
		if(!hasMore()){
			if(verbose){System.err.println("hasMore() returned false;  buffer="+(buffer==null ? null : buffer.size())+", next="+next+", consumed="+consumed);}
			return null;
		}
		Read r=buffer.set(next, null);
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return r;
	}
	
	private synchronized void fillBuffer(){
		if(builder==null){builder=new ByteBuilder(2000);}
		if(verbose){System.err.println("Filling buffer.  buffer="+(buffer==null ? null : buffer.size()));}
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReads(BUF_LEN, nextReadID, interleaved);
		final int count=(buffer==null ? 0 : buffer.size());

		if(verbose){System.err.println("Filled buffer.  size="+count);}
		
		nextReadID+=count;
		if(count<BUF_LEN){
			if(verbose){System.err.println("Closing tf");}
			errorState|=close();
		}
		
		generated+=count;
		if(verbose){System.err.println("generated="+generated);}
	}
	
	private ArrayList<Read> toReads(int maxReadsToReturn, long numericID, boolean interleaved){
		ArrayList<Read> list=toReadList(maxReadsToReturn, numericID, interleaved);
		if(list==null){assert(finished);}
		else{assert(list.size()<=maxReadsToReturn);}
		return list;
	}
	
	private ArrayList<Read> toReadList(int maxReadsToReturn, long numericID, boolean interleaved){
		if(finished){return null;}
		if(verbose){System.err.println("FastaQualRIS fetching a list.");}
		
		if(currentHeader==null && numericID==0){
//			assert(numericID==0) : numericID;
			nextBases(btf, builder);
			nextQualities(qtf, builder);
			if(nextHeaderB==null){
				finish();
				return null;
			}
			assert(Tools.equals(nextHeaderB, nextHeaderQ)) : "Quality and Base headers differ for read "+numericID;
			currentHeader=nextHeaderB;
			nextHeaderB=nextHeaderQ=null;
			if(currentHeader==null){
				finish();
				return null;
			}
		}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(1000, maxReadsToReturn));
		
		int added=0;
		
		Read prev=null;
		
		while(added<maxReadsToReturn){
			Read r=makeRead(numericID);
			if(verbose){System.err.println("Made "+r);}
			if(r==null){
				finish();
				if(verbose){System.err.println("makeRead returned null.");}
				break;
			}
			if(interleaved){
				if(prev==null){prev=r;}
				else{
					prev.mate=r;
					r.mate=prev;
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}
		}
		
		assert(list.size()<=maxReadsToReturn);
		if(verbose){System.err.println("FastaQualRIS returning a list.  Size="+list.size());}
		return list;
	}
	
	private final byte[] nextBases(ByteFile btf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=btf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			bb.append(line);
			line=btf.nextLine();
		}
		
		
		if(line==null){//end of file
			//do nothing
		}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderB=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	private final byte[] nextQualities(ByteFile qtf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=qtf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			if(NUMERIC_QUAL && line.length>0){
				int x=0;
				for(int i=0; i<line.length; i++){
					byte b=line[i];
					if(b==space){
						assert(i>0);
						bb.append((byte)x);
						x=0;
					}else{
						x=10*x+(b-zero);
					}
				}
				bb.append((byte)x);
			}else{
				for(byte b : line){bb.append((byte)(b-FASTQ.ASCII_OFFSET));}
			}
			line=qtf.nextLine();
		}
//		assert(bb.length()<1) : "'"+Arrays.toString(bb.toBytes())+"'";
		
		if(line==null){//end of file
			//do nothing
		}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderQ=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	private Read makeRead(long numericID){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		if(currentHeader==null){return null;}
		assert(nextHeaderB==null);
		assert(nextHeaderQ==null);
		
		final byte[] bases=nextBases(btf, builder);
		final byte[] quals=nextQualities(qtf, builder);
		final byte[] header=currentHeader;
		
		currentHeader=nextHeaderB;
		nextHeaderB=nextHeaderQ=null;
		
		if(bases==null){
			if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
			return null;
		}
		
		assert(bases.length==quals.length) :
			"\nFor sequence "+numericID+", name "+new String(header)+":\n" +
					"The bases and quality scores are different lengths, "+bases.length+" and "+quals.length;
		
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Tools.toUpperCase(bases[i]);
		}
//		for(int i=0; i<quals.length; i++){
//			quals[i]=(byte)(quals[i]-FASTQ.ASCII_OFFSET);
//		}
		
		assert(bases[0]!=carrot) : new String(bases)+"\n"+numericID+"\n"+header[0];
		String hd=new String(header, 1, header.length-1);
		Read r=new Read(bases, quals, hd, numericID);
		return r;
	}
	
	@Override
	public synchronized boolean close(){
		if(closed){return errorState;}
		if(verbose){System.err.println("FastaQualRIS closing.");}
//		if(verbose){new Exception().printStackTrace(System.err);}
		builder=null;
		finish();
		boolean a=btf.close();
		boolean b=qtf.close();
		closed=true;
		return a|b;
	}

	@Override
	public synchronized void restart() {
		if(verbose){System.err.println("FastaQualRIS restarting.");}
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		
		finished=false;
		closed=false;

		buffer=null;
		nextHeaderB=null;
		nextHeaderQ=null;
		currentHeader=null;
		builder=null;

		btf.reset();
		qtf.reset();
	}

	@Override
	public boolean paired() {
		return interleaved;
	}
	
	private synchronized void finish(){
		if(verbose){System.err.println("FastaQualRIS setting finished "+finished+" -> "+true);}
		finished=true;
	}
	
	@Override
	public String fname(){return btf.name()+","+qtf.name();}

	private ArrayList<Read> buffer=null;
	private int next=0;

	private final ByteFile btf;
	private final ByteFile qtf;
	private final boolean interleaved;

	private final int BUF_LEN=Shared.bufferLen();;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public static boolean NUMERIC_QUAL=true;
	
	public static boolean verbose=false;

	private byte[] nextHeaderB=null;
	private byte[] nextHeaderQ=null;
	
	private byte[] currentHeader=null;
	
	private ByteBuilder builder=null;
	
	private boolean finished=false, closed=false;
	private final byte carrot='>', space=' ', zero='0';
	
}
