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
 * List reader for BBTools serialized Read text, decoded by Read.fromText.
 * Multiple primary files contribute site lists to corresponding first-file
 * records. Corresponding entries must have matching names and numeric IDs;
 * list sizes must align. Bases and qualities come from the first file. Files can contain interleaved pairs or
 * have a separate mate reader, optionally wrapped by the legacy concurrent adapter.
 * Filename arrays and several state fields are exposed without defensive copies;
 * method-level synchronization is not a general concurrent-use guarantee.
 *
 * @author Brian Bushnell
 * @date Jul 16, 2013
 */
public class RTextInputStream extends ReadInputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Diagnostic multi-file merge printing serialized reads to stdout.
	 * Uses the constructor's negative-limit convention to scan all input entries.
	 * @param args Nonempty synchronized serialized-read filename array
	 */
	public static void main(String[] args){
		//Resolved #002: scan all entries using the constructor's negative/unlimited convention.
		//Zero is a rejected limit when assertions are enabled, not an unlimited sentinel.
		RTextInputStream rtis=new RTextInputStream(args, -1);
		ArrayList<Read> list=rtis.nextList();
		while(list!=null){
			for(Read r : list){
				System.out.println(r.toText(true));
			}
			list=rtis.nextList();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens inputs using only format names; other format options are not forwarded.
	 * @param ff1 Nonnull primary format
	 * @param ff2 Optional mate format
	 * @param crisReadLimit Positive entry limit, or negative for unlimited; zero is rejected
	 * when assertions are enabled
	 */
	public RTextInputStream(FileFormat ff1, FileFormat ff2, long crisReadLimit){
		this(ff1.name(), (ff2==null ? null : ff2.name()), crisReadLimit);
	}

	/**
	 * Opens one primary file and an optional separate mate file.
	 * @param fname1 Primary serialized-read filename
	 * @param fname2 Mate filename, null or the case-insensitive literal "null"
	 * @param crisReadLimit Positive entry limit, or negative for unlimited
	 */
	public RTextInputStream(String fname1, String fname2, long crisReadLimit){
		this(new String[]{fname1}, (fname2==null || "null".equalsIgnoreCase(fname2)) ? null : new String[]{fname2}, crisReadLimit);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
	}
	/** Opens synchronized primary files without a separate mate-file array.
	 * @param fnames_ Nonempty filename array, retained by reference
	 * @param crisReadLimit Positive entry limit, or negative for unlimited */
	public RTextInputStream(String[] fnames_, long crisReadLimit){this(fnames_, null, crisReadLimit);}
	
	/**
	 * Opens primary TextFiles with subprocess support and selects pairing once.
	 * Separate mate files suppress primary interleaving. Otherwise FASTQ's testing/
	 * forcing flags and stdin status select forced pairing or marker testing.
	 * A separate mate reader is constructed recursively; USE_CRIS starts its optional
	 * concurrent adapter during construction. The entry quota counts each merged
	 * primary entry once, independent of the number of contributing files.
	 * @param fnames_ Nonempty synchronized primary filename array, retained by reference
	 * @param mate_fnames_ Optional synchronized mate filename array
	 * @param crisReadLimit Positive list-entry limit, or negative for unlimited
	 */
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
	
		mateStream=(mate_fnames_==null ? null : new RTextInputStream(mate_fnames_, null, crisReadLimit));
		cris=((!USE_CRIS || mateStream==null) ? null : new ConcurrentLegacyReadInputStream(mateStream, crisReadLimit));
		if(cris!=null){cris.start();}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests the first nonblank line of a regular file for the exact #INTERLEAVED marker.
	 * Uses a separate TextFile without subprocess support, closed after the test.
	 * @param fname Existing regular-file path
	 * @return true only when the first returned line equals the marker
	 */
	public static boolean isInterleaved(String fname){
		File f=new File(fname);
		assert(f.exists() && f.isFile());
		TextFile tf=new TextFile(fname, false);
		String s=tf.nextLine();
		tf.close();
		return "#INTERLEAVED".equals(s);
	}
	
	/** Returns the next merged batch, or null after termination; entries may have mates. */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(finished){return null;}
		return readList();
	}
	
	/**
	 * Merges aligned site lists, associates mates and terminates on a short batch.
	 * Primary-file entries supply all fields except added sites and mate links.
	 * Separate mates require equal batch sizes and numeric IDs; direct mate mode
	 * also checks names. Interleaved records retain serialized pair-number flags.
	 * @return Merged primary entries, or null when no entries remain
	 */
	private synchronized ArrayList<Read> readList(){
		assert(buffer==null);
		if(finished){return null;}
		
		ArrayList<Read> merged=getListFromFile(textfiles[0]);
		
		if(textfiles.length>1){
			ArrayList<Read>[] temp=new ArrayList[textfiles.length];
			temp[0]=merged;
			//[stream/RTextInputStream#001] start at 1: temp[0] is already 'merged' (textfiles[0] read above).
			//Starting at 0 re-read textfiles[0] a 2nd time, advancing its pointer 2 batches/call (silently
			//dropping every other batch of file 0 + desyncing the merge). Repository caller search
			//found named-file forms externally; array forms occur in local main and mate recursion.
			for(int i=1; i<temp.length; i++){
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
		// Fixed #004: charge the quota once per merged entry, after all files have
		//read the same remaining quota. Charging each file truncated finite merges.
		readCount+=merged.size();
	
		if(cris!=null){
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
					}
					
				}
			}
			cris.returnList(mates0.id, mates0.list.isEmpty());
		}else if(mateStream!=null){
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
	
	/**
	 * Reads up to the remaining shared entry quota from one TextFile.
	 * TextFile skips blank lines; # comments are skipped before each primary record.
	 * An interleaved mate is the next nonblank line without comment skipping.
	 * The caller counts merged entries once, including each pair anchor. Short reads
	 * close this TextFile; quota exhaustion itself does not guarantee closure.
	 * @param tf Open input containing normal serialized Read records
	 * @return Nonnull list, possibly empty
	 */
	private ArrayList<Read> getListFromFile(TextFile tf){
		
		int len=READS_PER_LIST;
		if(readLimit-readCount<len){len=(int)(readLimit-readCount);}
		
		ArrayList<Read> list=new ArrayList<Read>(len);
		
		for(int i=0; i<len; i++){
			String s=tf.nextLine();
			while(s!=null && s.charAt(0)=='#'){s=tf.nextLine();}
			if(s==null){break;}
			Read r=Read.fromText(s);
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
		
		if(list.size()<len){
			assert(tf.nextLine()==null);
			tf.close();
		}
		return list;
	}

	/** Returns whether separate mates or interleaved input were selected. */
	@Override
	public boolean paired(){return mateStream!=null || interleaved;}
	
	/** Marks this reader and mate components finished; does not close primary TextFiles. */
	public final void shutdown(){
		finished=true;
		if(mateStream!=null){mateStream.shutdown();}
		if(cris!=null){cris.shutdown();}
	}
	
	/** Returns Arrays.toString of the borrowed primary filename array. */
	@Override
	public String fname(){return Arrays.toString(fnames);}
	
	/** Returns an optimistic state check without looking ahead or verifying open files. */
	@Override
	public boolean hasMore(){
		if(buffer!=null && next<buffer.size()){return true;}
		return !finished;
	}
	
	/**
	 * Reopens primary files and restarts mate components with existing pairing choices.
	 * Clears local iteration state and renews the original entry quota.
	 * TextFile.reset closes and reopens each file.
	 */
	@Override
	public synchronized void restart(){
		// Fixed #003: rewinding files must also renew a finite quota; otherwise a fully
		//consumed reader returns no records after restart despite the file rewind.
		finished=false;
		next=0;
		buffer=null;
		readCount=0;
		for(TextFile tf : textfiles){tf.reset();}
		if(cris!=null){
			cris.restart();
			cris.start();
		}else if(mateStream!=null){mateStream.restart();}
	}

	/**
	 * Closes primary TextFiles and the selected mate reader or adapter.
	 * Does not mark this reader finished or update its inherited error flag. TextFile
	 * close returns false while storing its own error state, not collected here;
	 * the direct mate close result is also discarded (#005).
	 * @return OR of collected close return values and selected mate status
	 */
	@Override
	public synchronized boolean close(){
		//TODO: Probable bug #005 - TextFile.close returns false and keeps errors in its own flag.
		//Those flags are not collected here; direct mate close also discards its return value.
		boolean error=false;
		for(TextFile tf : textfiles){error|=tf.close();}
		if(cris!=null){
			error|=ReadWrite.closeStream(cris);
		}else if(mateStream!=null){
			mateStream.close();
			error|=mateStream.errorState();
		}
		return error;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Termination flag checked by list iteration and optimistic hasMore. */
	public boolean finished=false;
	/** Borrowed primary filename array, exposed for legacy callers. */
	public String[] fnames;
	/** Primary file readers opened during construction. */
	public TextFile[] textfiles;
	
	/** Legacy buffer checked by hasMore; the list path requires it to be null. */
	private ArrayList<Read> buffer=null;
	/** Legacy buffer index, cleared by restart. */
	private int next=0;
	
	/** Consumed merged primary entries, counted once per entry regardless of file count. */
	private long readCount;
	/** Entry quota; negative constructor arguments become Long.MAX_VALUE. */
	private final long readLimit;
	/** Pairing choice captured at construction, not recomputed by restart. */
	private final boolean interleaved;
	
	/** Maximum batch entries, captured at class initialization. */
	public static final int READS_PER_LIST=Shared.bufferLen();

	/** Optional recursively constructed mate reader. */
	private final RTextInputStream mateStream;
	/** Optional legacy concurrent adapter around the mate reader. */
	private final ConcurrentLegacyReadInputStream cris;
	/** Whether construction wraps separate mates in the legacy concurrent adapter. */
	public static boolean USE_CRIS=true;//Historical note reports faster zipped paired input; not remeasured.

}
