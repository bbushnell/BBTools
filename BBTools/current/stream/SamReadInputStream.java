package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class SamReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		SamReadInputStream sris=new SamReadInputStream(args[0], false, false, true);
		
		Read r=sris.next();
		System.out.println(r.toText(false));
		System.out.println();
		System.out.println(r.obj.toString());
		System.out.println();
	}
	
	public SamReadInputStream(String fname, boolean loadHeader_, boolean interleaved_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.SAM, null, allowSubprocess_, false), loadHeader_, interleaved_);
	}
		
	public SamReadInputStream(FileFormat ff, boolean loadHeader_, boolean interleaved_){
		loadHeader=loadHeader_;
//		assert(loadHeader);
//		interleaved=((tf.is==System.in || stdin) ? FASTQ.FORCE_INTERLEAVED : true);
		interleaved=interleaved_;
		
		stdin=ff.stdio();
		if(!ff.samOrBam()){
			System.err.println("Warning: Did not find expected sam file extension for filename "+ff.name());
		}
		
		tf=ByteFile.makeByteFile(ff);
		header=new ArrayList<byte[]>();
		
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
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return list;
	}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReadList(tf, BUF_LEN, nextReadID, FASTQ.PARSE_CUSTOM);
		nextReadID+=buffer.size();
		generated+=buffer.size();
		
		if(buffer.size()<BUF_LEN){tf.close();}
	}
	
	/**
	 * @param tf2
	 * @param bUF_LEN2
	 * @param nextReadID2
	 * @param interleaved2
	 * @return
	 */
	private final ArrayList<Read> toReadList(ByteFile tf2, int buflen, long nextReadID2, boolean parseCustom) {
		ArrayList<Read> list=new ArrayList<Read>(buflen);
		while(list.size()<buflen){
			byte[] line=tf2.nextLine();
//			System.out.println("A: Read line "+new String(line));
			while(line!=null && line[0]=='@'){
//				System.out.println(">"+new String(line));
				if(Shared.TRIM_RNAME){line=trimHeaderSQ(line);}
				if(loadHeader){header.add(line);}
				line=tf2.nextLine();
//				assert(false) : new String(line)+"\n"+header.size()+", "+SHARED_HEADER;
//				System.out.println("B: Read line "+new String(line));
			}
			if(loadHeader && nextReadID2==0){setSharedHeader(header);}
			if(line==null){return list;}
			SamLine sl1=new SamLine(line);
			
			Read r1=sl1.toRead(parseCustom);
			r1.samline=sl1;
			r1.numericID=nextReadID2;
			list.add(r1);
			if(interleaved && (sl1.flag&0x1)!=0){
				assert((sl1.flag&0x40)!=0) : r1+"\n\n"+sl1;
				byte[] line2=tf2.nextLine();
				SamLine sl2=null;
				Read r2=null;
				if(line2!=null){
					sl2=new SamLine(line2);
					r2=sl2.toRead(parseCustom);
					r2.numericID=nextReadID2;
				}else{
					assert(false) : r1+"\n\n"+sl1;
				}
				if(sl2!=null){
					assert((sl2.flag&0x1)!=0);
					assert((sl2.flag&0x80)!=0) : r2+"\n\n"+sl2+"\nflag="+Integer.toBinaryString(sl2.flag)+"\n";
					r1.mate=r2;
					r2.mate=r1;
					
					int lim=Tools.min(sl1.qname.length(), sl2.qname.length());
					for(int i=0; i<lim; i++){
						char a=sl1.qname.charAt(i);
						char b=sl2.qname.charAt(i);
						if(a=='/' || b=='/' || Character.isWhitespace(a) || Character.isWhitespace(b)){break;}
						assert(a==b) : "Name mismatch for paired reads: '"+sl1.qname+"' != '"+sl2.qname+"'\n\n"+sl1+"\n\n"+sl2;
					}
					
				}
			}
			nextReadID2++;
		}
		return list;
	}

	@Override
	public boolean close(){
		return tf.close();
	}
	
	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		header=new ArrayList<byte[]>();
		tf.reset();
	}
	
	public static synchronized ArrayList<byte[]> getSharedHeader(boolean wait){
		if(!wait || SHARED_HEADER!=null){return SHARED_HEADER;}
		System.err.println("Waiting on header to be read from a sam file.");
		while(SHARED_HEADER==null){
			try {
				SamReadInputStream.class.wait(100);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return SHARED_HEADER;
	}
	
	public static synchronized void setSharedHeader(ArrayList<byte[]> list){
//		assert(false) : list.size();
//		if(list!=null && Shared.TRIM_RNAME){
//			for(int i=0; i<list.size(); i++){
//				byte[] line=list.get(i);
//				list.set(i, trimHeaderSQ(line));
//			}
//		}
		SHARED_HEADER=list;
		SamReadInputStream.class.notifyAll();
	}
	
	static byte[] trimHeaderSQ(byte[] line){
		if(line==null || !Tools.startsWith(line, "@SQ")){return line;}
		
		final int idx=Tools.indexOfDelimited(line, "SN:", 2, (byte)'\t');
		if(idx<0){
			assert(false) : "Bad header: "+new String(line);
			return line;
		}
		
		int trimStart=-1;
		for(int i=idx; i<line.length; i++){
			final byte b=line[i];
			if(b=='\t'){return line;}
			if(Character.isWhitespace(b)){
				trimStart=i;
				break;
			}
		}
		if(trimStart<0){return line;}
		
		final int trimStop=Tools.indexOf(line, (byte)'\t', trimStart+1);
		final int bbLen=trimStart+(trimStop<0 ? 0 : line.length-trimStop);
		final ByteBuilder bb=new ByteBuilder(bbLen);
		for(int i=0; i<trimStart; i++){bb.append(line[i]);}
		if(trimStop>=0){
			for(int i=trimStop; i<line.length; i++){bb.append(line[i]);}
		}
		assert(bb.length==bbLen) : bbLen+", "+bb.length+", idx="+idx+", trimStart="+trimStart+", trimStop="+trimStop+"\n\n"+new String(line)+"\n\n"+bb+"\n\n";
		
		return bb.array;
	}
	
	@Override
	public String fname(){return tf.name();}
	
	private static volatile ArrayList<byte[]> SHARED_HEADER;
//	private static boolean SET_SHARED_HEADER;
	
	@Override
	public boolean paired() {return interleaved;}

	private ArrayList<Read> buffer=null;
	private ArrayList<byte[]> header=null;
	private int next=0;
	
	private final ByteFile tf;
	private final boolean interleaved;
	private final boolean loadHeader;

	private final int BUF_LEN=Shared.bufferLen();;
	
	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean stdin;

}
