package stream;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * This class allows multiple files as input.
 * These files are synchronized, so a read will be created by merging the sitescores from the same line of each file.
 * @author Brian Bushnell
 * @date Jul 16, 2013
 *
 */
public class RTextInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		RTextInputStream rtis=new RTextInputStream(args, 0);
		ArrayList<Read> list=rtis.nextList();
		while(list!=null){
			for(Read r : list){
				System.out.println(r.toText(true));
			}
			list=rtis.nextList();
		}
	}

	public RTextInputStream(FileFormat ff1, FileFormat ff2, long crisReadLimit){
		this(ff1.name(), (ff2==null ? null : ff2.name()), crisReadLimit);
	}

	public RTextInputStream(String fname1, String fname2, long crisReadLimit){
		this(new String[] {fname1}, (fname2==null || "null".equalsIgnoreCase(fname2)) ? null : new String[] {fname2}, crisReadLimit);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
	}
	public RTextInputStream(String[] fnames_, long crisReadLimit){this(fnames_, null, crisReadLimit);}
	
	public RTextInputStream(String[] fnames_, String[] mate_fnames_, long crisReadLimit){
		fnames=fnames_;
		textfiles=new TextFile[fnames.length];
		for(int i=0; i<textfiles.length; i++){
			textfiles[i]=new TextFile(fnames[i], true);
		}
		
		readLimit=(crisReadLimit<0 ? Long.MAX_VALUE : crisReadLimit);
		if(readLimit==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
		interleaved=(mate_fnames_!=null ? false :
			(!FASTQ.TEST_INTERLEAVED || textfiles[0].is==System.in) ? FASTQ.FORCE_INTERLEAVED : isInterleaved(fnames[0]));
		
//		assert(false) : (mate_fnames_!=null)+", "+(textfiles[0].is==System.in)+", "+interleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+isInterleaved(fnames[0]);
		
		mateStream=(mate_fnames_==null ? null : new RTextInputStream(mate_fnames_, null, crisReadLimit));
		cris=((!USE_CRIS || mateStream==null) ? null : new ConcurrentLegacyReadInputStream(mateStream, crisReadLimit));
		if(cris!=null){cris.start();}
	}
	
	public static boolean isInterleaved(String fname){
		File f=new File(fname);
		assert(f.exists() && f.isFile());
		TextFile tf=new TextFile(fname, false);
		String s=tf.nextLine();
		tf.close();
		return "#INTERLEAVED".equals(s);
	}

	@Override
	public void start() {
		assert(false); //Not fully implemented everywhere...
		if(cris!=null){cris.start();}
	}
	
//	@Override
//	public synchronized Read[] nextBlock(){
//		ArrayList<Read> list=readList();
//		if(list==null || list.size()==0){return null;}
//		return list.toArray(new Read[list.size()]);
//	}
	
	@Override
	public synchronized ArrayList<Read> nextList(){
//		System.out.println((mateStream==null ? "F5: " : "F3: ")+"Grabbing a list: finished="+finished);
		if(finished){return null;}
		return readList();
	}
	
	private synchronized ArrayList<Read> readList(){
		assert(buffer==null);
//		System.out.println((mateStream==null ? "F5: " : "F3: ")+" Entering readList");
		if(finished){return null;}
		
		ArrayList<Read> merged=getListFromFile(textfiles[0]);
		
		if(textfiles.length>1){
			ArrayList<Read>[] temp=new ArrayList[textfiles.length];
			temp[0]=merged;
			for(int i=0; i<temp.length; i++){
				temp[i]=getListFromFile(textfiles[i]);
			}
			
			for(int i=0; i<merged.size(); i++){
				Read r=merged.get(i);
				for(int j=1; j<temp.length; j++){
					Read r2=temp[j].get(i);
					assert(r2.numericID==r.numericID);
					assert(r2.id.equals(r.id));
					if(r.sites==null){r.sites=r2.sites;}
					else if(r2.sites!=null){r.sites.addAll(r2.sites);}
				}
			}
		}
		
//		System.out.println((mateStream==null ? "F5: " : "F3: ")+"Merged: "+merged==null ? "null" : ""+merged.size());
		
		if(cris!=null){
			//				System.out.println((mateStream==null ? "F5: " : "F3: ")+"Grabbing a mate list: finished="+mateStream.finished);
			ListNum<Read> mates0=cris.nextList();
			ArrayList<Read> mates=mates0.list;
			assert((mates==null || mates.size()==0) == (merged==null || merged.size()==0)) : (merged==null)+", "+(mates==null);
			if(merged!=null && mates!=null){
				
				assert(mates.size()==merged.size()) : "\n"+mates.size()+", "+merged.size()+", "+paired()+"\n"+
						merged.get(0).toText(false)+"\n"+mates.get(0).toText(false)+"\n\n"+
						merged.get(merged.size()-1).toText(false)+"\n"+mates.get(mates.size()-1).toText(false)+"\n\n"+
						merged.get(Tools.min(merged.size(), mates.size())-1).toText(false)+"\n"+
						mates.get(Tools.min(merged.size(), mates.size())-1).toText(false)+"\n\n";

				for(int i=0; i<merged.size(); i++){
					Read r1=merged.get(i);
					Read r2=mates.get(i);
					r1.mate=r2;
					assert(r1.pairnum()==0);
					
					if(r2!=null){
						r2.mate=r1;
						r2.setPairnum(1);
						assert(r2.numericID==r1.numericID) : "\n\n"+r1.toText(false)+"\n\n"+r2.toText(false)+"\n";
//						assert(r2.id.equals(r1.id)) : "\n\n"+r1.toText(false)+"\n\n"+r2.toText(false)+"\n";
					}
					
				}
			}
			cris.returnList(mates0.id, mates0.list.isEmpty());
		}else if(mateStream!=null){
			//			System.out.println((mateStream==null ? "F5: " : "F3: ")+"Grabbing a mate list: finished="+mateStream.finished);
			ArrayList<Read> mates=mateStream.readList();
			assert((mates==null || mates.size()==0) == (merged==null || merged.size()==0)) : (merged==null)+", "+(mates==null);
			if(merged!=null && mates!=null){
				assert(mates.size()==merged.size()) : mates.size()+", "+merged.size();

				for(int i=0; i<merged.size(); i++){
					Read r1=merged.get(i);
					Read r2=mates.get(i);
					r1.mate=r2;
					r2.mate=r1;
					
					assert(r1.pairnum()==0);
					r2.setPairnum(1);

					assert(r2.numericID==r1.numericID) : "\n\n"+r1.toText(false)+"\n\n"+r2.toText(false)+"\n";
					assert(r2.id.equals(r1.id)) : "\n\n"+r1.toText(false)+"\n\n"+r2.toText(false)+"\n";
				}
			}
		}
		
		if(merged.size()<READS_PER_LIST){
			if(merged.size()==0){merged=null;}
			shutdown();
		}
		
		return merged;
	}
	
	private ArrayList<Read> getListFromFile(TextFile tf){
		
		int len=READS_PER_LIST;
		if(readLimit-readCount<len){len=(int)(readLimit-readCount);}
		
		ArrayList<Read> list=new ArrayList<Read>(len);
		
		for(int i=0; i<len; i++){
			String s=tf.nextLine();
			while(s!=null && s.charAt(0)=='#'){s=tf.nextLine();}
			if(s==null){break;}
			Read r=Read.fromText(s);
//			assert(r.toString().equals(s)) : "\n\n"+s+"\n!=\n"+r.toString()+"\n\n";
//			assert(r.chrom>0 == r.mapScore>0) : r.toText(false);
			if(interleaved){
				s=tf.nextLine();
				assert(s!=null) : "Odd number of reads in interleaved file "+tf.name;
				if(s!=null){
					Read r2=Read.fromText(s);
					assert(r2.numericID==r.numericID) : "Different numeric IDs for paired reads in interleaved file "+tf.name;
					r2.numericID=r.numericID;
					r2.mate=r;
					r.mate=r2;
				}
			}
			list.add(r);
		}
		readCount+=list.size();
		
		if(list.size()<len){
			assert(tf.nextLine()==null);
			tf.close();
		}
		return list;
	}

	@Override
	public boolean paired() {
		return mateStream!=null || interleaved;
	}
	
	public final void shutdown(){
		finished=true;
		if(mateStream!=null){mateStream.shutdown();}
		if(cris!=null){cris.shutdown();}
	}
	
	@Override
	public String fname(){return Arrays.toString(fnames);}
	
	public boolean finished=false;
	public String[] fnames;
	public TextFile[] textfiles;
	
	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private long readCount;
	private final long readLimit;
	private final boolean interleaved;
	
	public static final int READS_PER_LIST=Shared.bufferLen();;

	private final RTextInputStream mateStream;
	private final ConcurrentLegacyReadInputStream cris;
	public static boolean USE_CRIS=true; //Doubles read speed for zipped paired files

	@Override
	/** This is optimistic and may return "true" incorrectly. */
	public boolean hasMore() {
		if(buffer!=null && next<buffer.size()){return true;}
		return !finished;
	}
	
	
	@Override
	/** ONLY CALL FROM A SINGLE THREAD! */
	public Read next() {
		if(buffer==null || next>=buffer.size()){
			buffer=null;
			next=0;
			if(!finished){
				buffer=nextList();
			}
		}
		
		if(buffer==null || next>=buffer.size()){
			assert(finished);
			return null;
		}
		Read r=buffer.get(next);
		buffer.set(next, null);
		next++;
		return r;
	}
	
	
	@Override
	public synchronized void restart() {
		finished=false;
		next=0;
		buffer=null;
		for(TextFile tf : textfiles){tf.reset();}
		if(cris!=null){
			cris.restart();
			cris.start();
		}else if(mateStream!=null){mateStream.restart();}
	}

	@Override
	public synchronized boolean close() {
		boolean error=false;
		for(TextFile tf : textfiles){error|=tf.close();}
		if(cris!=null){
			error|=ReadWrite.closeStream(cris);;
		}else if(mateStream!=null){
			mateStream.close();
			error|=mateStream.errorState();
		}
		return error;
	}
	
}
