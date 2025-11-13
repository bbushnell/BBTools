package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.bam.BgzfSettings;
import structures.IntList;
import structures.ListNum;

/**
 * ByteFile variant that returns FASTA records with newlines already stripped.
 * Maintains an IntList of newline positions for the consumer.
 * 
 * @author Brian Bushnell & Isla
 * @date November 12, 2025
 */
public final class ByteFile1Fb extends ByteFile {

	public ByteFile1Fb(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false));
	}

	public ByteFile1Fb(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1Fb("+ff+")");}
		is=open();
	}

	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
		firstRecord=true;
		listPos=0;
		recordPositions.clear();
		allNewlines.clear();
	}

	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		errorState|=ReadWrite.finishReading(is, name(), (allowSubprocess() || FileFormat.isBamFile(name())));

		is=null;
		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}

	/** Compatibility wrapper - creates temporary IntList */
	@Override
	public final byte[] nextLine(){
		IntList temp=new IntList();
		return nextLine(temp);
	}

	/** 
	 * Return next FASTA record with internal newlines stripped.
	 * @param newlines Will be populated with newline positions in the returned array
	 */
	public final byte[] nextLine(IntList newlines){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		// For the first record, find the first '>'
		if(firstRecord){
			if(condensedStop==0){fillAndCondense();}
			while(condensedStart<condensedStop && condensed[condensedStart]!=carrot){condensedStart++;}
			if(condensedStart>=condensedStop){
				close();
				return null;
			}
			firstRecord=false;
		}

		// Need more data?
		if(listPos>=recordPositions.size()){
			fillAndCondense();
			if(listPos>=recordPositions.size() && condensedStart>=condensedStop){
				// EOF
				close();
				return null;
			}
		}

		lineNum++;

		// Determine record boundaries
		final int limit;
		if(listPos<recordPositions.size()){
			// Next boundary is at recordPositions[listPos] (points to \n after record)
			limit=recordPositions.get(listPos++);
		}else{
			// No more boundaries, use rest of buffer
			limit=condensedStop;
		}

		if(condensedStart>=limit){
			// Empty record
			condensedStart=limit+1;
			newlines.clear();
			return blankLine;
		}

		// Extract record
		byte[] record=KillSwitch.copyOfRange(condensed, condensedStart, limit+1); // Include the \n

		// Populate newlines IntList with positions relative to this record
		newlines.clear();
		// Find newlines within this record's range
		for(int i=0; i<allNewlines.size(); i++){
			int pos=allNewlines.get(i);
			if(pos>=condensedStart && pos<=limit){
				newlines.add(pos-condensedStart); // Relative position
			}else if(pos>limit){
				break; // Past this record
			}
		}

		// Move past this record
		condensedStart=limit+1;

		return record;
	}
	
	/** Compatibility wrapper */
	@Override
	public final ListNum<byte[]> nextList(){
		IntList temp=new IntList();
		return nextList(temp);
	}

	/**
	 * Return list of FASTA records.
	 * @param newlines Will be populated with newline positions across all records
	 */
	public final ListNum<byte[]> nextList(IntList newlines){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		
		newlines.clear();
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);
		for(int i=0; i<listSize; i++){
			IntList recordNewlines=new IntList();
			byte[] record=nextLine(recordNewlines);
			if(record==null){
				break;
			}
			list.add(record);
			// Accumulate newline positions
			int offset=(newlines.size()>0 ? newlines.get(newlines.size()-1)+1 : 0);
			for(int j=0; j<recordNewlines.size(); j++){
				newlines.add(offset+recordNewlines.get(j));
			}
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Fill buffer, find newlines, strip internal ones, build condensed array */
	private void fillAndCondense(){
		// Shift remaining data to start of buffer
		if(bstart>0 && bstart<bstop){
			int extra=bstop-bstart;
			System.arraycopy(buffer, bstart, buffer, 0, extra);
			bstop=extra;
			bstart=0;
		}else if(bstart>=bstop){
			bstart=0;
			bstop=0;
		}

		// Clear state
		recordPositions.clear();
		allNewlines.clear();

		// Read data until we have at least one complete record
		while(recordPositions.isEmpty()){
			if(bstop==buffer.length){
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			
			int r=-1;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}catch(NullPointerException e){
				if(!Shared.anomaly){e.printStackTrace();}
			}
			
			if(r<=0){break;} // EOF
			
			final int from=Math.max(0, bstop-1);
			bstop+=r;

			// Find all \n> boundaries
			Vector.findFastaHeaders(buffer, from, bstop, recordPositions);
		}

		// Find all newlines in the buffer
		IntList tempNewlines=new IntList();
		Vector.findSymbols(buffer, 0, bstop, slashn, tempNewlines);

		// Determine where complete records end
		int completeEnd=bstop;
		if(!recordPositions.isEmpty()){
			completeEnd=recordPositions.get(recordPositions.size()-1)+2; // After last \n>
		}

		// Condense: strip internal newlines, keep structural ones
		condensedStart=0;
		condensedStop=0;
		if(condensed.length<bstop){
			condensed=new byte[bstop*2];
		}

		int writePos=0;
		int readPos=0;
		int nlIdx=0;

		while(readPos<completeEnd){
			if(nlIdx<tempNewlines.size() && readPos==tempNewlines.get(nlIdx)){
				// At a newline
				int nlPos=tempNewlines.get(nlIdx);
				boolean isStructural=false;

				// Check if this is a structural newline (before >, or last in record)
				if(nlPos+1<completeEnd && buffer[nlPos+1]==carrot){
					isStructural=true; // \n> boundary
				}
				// Could also check if it's the last newline in a record, but \n> should catch those

				if(isStructural){
					// Keep this newline
					// But strip \r if present
					if(writePos>0 && condensed[writePos-1]=='\r'){
						writePos--; // Remove the \r
					}
					condensed[writePos++]=slashn;
					allNewlines.add(writePos-1); // Record position in condensed array
				}
				// Skip this newline in input (already copied or skipped)
				readPos++;
				nlIdx++;
			}else{
				// Regular character - copy it (unless it's \r before \n)
				if(buffer[readPos]!='\r' || (readPos+1<completeEnd && buffer[readPos+1]!=slashn)){
					condensed[writePos++]=buffer[readPos];
				}
				readPos++;
			}
		}

		// Add terminal newline if last record doesn't have one
		if(writePos>0 && condensed[writePos-1]!=slashn){
			condensed[writePos++]=slashn;
			allNewlines.add(writePos-1);
		}

		condensedStop=writePos;

		// Shift incomplete data back to buffer start
		if(completeEnd<bstop){
			int remaining=bstop-completeEnd;
			System.arraycopy(buffer, completeEnd, buffer, 0, remaining);
			bstop=remaining;
			bstart=0;
		}else{
			bstart=0;
			bstop=0;
		}

		// Reset list position
		listPos=0;
	}

	@Override
	public void pushBack(byte[] record){
		throw new UnsupportedOperationException("pushBack not supported for ByteFile1Fb");
	}

	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile1Fb "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=0;
		bstop=0;
		condensedStart=0;
		condensedStop=0;
		firstRecord=true;
		listPos=0;
		recordPositions.clear();
		allNewlines.clear();
		return is;
	}

	@Override
	public boolean isOpen(){return open;}

	@Override
	public final InputStream is(){return is;}

	@Override
	public final long lineNum(){return lineNum;}

	private boolean open=false;
	private byte[] buffer=new byte[bufferlen];
	private byte[] condensed=new byte[bufferlen];
	private static final byte[] blankLine=new byte[0];
	private int bstart=0, bstop=0;
	private int condensedStart=0, condensedStop=0;
	public InputStream is;
	public long lineNum=-1;
	private IntList recordPositions=new IntList(); // Positions of \n before >
	private IntList allNewlines=new IntList(); // All structural newlines in condensed
	private int listPos=0;
	private boolean firstRecord=true;

	public static boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=65536;

	private boolean errorState=false;
}