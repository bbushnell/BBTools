package synth;

import java.util.ArrayList;

import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import template.BBTool_ST;

/**
 * Fuses sequences together, with N-padding in between.
 * @author Brian Bushnell
 * @date Jan 20, 2015
 *
 */
public final class FuseSequence extends BBTool_ST {
	
	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		FuseSequence fs=new FuseSequence(args);
		fs.process(t);
	}
	
	/**
	 * Constructs FuseSequence with command-line arguments.
	 * Sets default padding symbol and enables amino acid mode if specified.
	 * @param args Command-line arguments for configuration
	 */
	public FuseSequence(String[] args){
		super(args);
		reparse(args);
		if(PAD_SYMBOL==0){PAD_SYMBOL='N';}
		if(amino){
			Shared.AMINO_IN=true;
			if(PAD_SYMBOL=='N'){PAD_SYMBOL='X';}
		}
	}
	
	@Override
	protected void setDefaults(){
		npad=300;
		defaultQuality=30;
		fusePairs=false;
		PAD_SYMBOL=0;
		amino=false;
	}
	
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		if(a.equals("pad") || a.equals("npad") || a.equals("ns")){
			npad=Integer.parseInt(b);
			return true;
		}else if(a.equals("q") || a.equals("quality")){
			defaultQuality=Byte.parseByte(b);
			return true;
		}else if(a.equals("fp") || a.equals("fusepairs")){
			fusePairs=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("maxlen")){
			maxlen=(int)Tools.min(Parse.parseKMG(b), Shared.MAX_ARRAY_LEN);
			return true;
		}else if(a.equals("rename") || a.equals("name") || a.equals("prefix")){
			name=(b==null ? "" : b);
			return true;
		}else if(a.equals("addnumber")){
			addNumber=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("padsymbol")){
			PAD_SYMBOL=b.charAt(0);
			return true;
		}else if(a.equals("amino")){
			amino=Parse.parseBoolean(b);
			return true;
		}
		return false;
	}
	
	@Override
	protected void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;

		long outNum=0;
		long lastListID=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				lastListID=ln.id;
				
				ArrayList<Read> readsOut=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
					if(!fusePairs && maxlen>0 && bases.length>0 && bases.length+initialLength1+initialLength2+npad>maxlen){
						Read r=bufferToRead(outNum);
						outNum++;
						readsOut.add(r);
					}
					processReadPair(r1, r2);
					if(fusePairs){
						readsOut.add(r1);
					}else if(maxlen>=0 && bases.length>=maxlen){
						Read r=bufferToRead(outNum);
						outNum++;
						readsOut.add(r);
					}
				}
				
				if(ros!=null && (fusePairs || maxlen>0)){ros.add(readsOut, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(ros!=null && bases.length()>0){
			ArrayList<Read> readsOut=new ArrayList<Read>(1);
			Read r=bufferToRead(outNum);
			outNum++;
			readsOut.add(r);
			ros.add(readsOut, lastListID);
		}
	}
	
	/**
	 * Converts accumulated sequences in buffer to a Read object.
	 * Clears the buffer after conversion and sets the read ID.
	 * @param id Numeric identifier for the output read
	 * @return Read object containing the concatenated sequence
	 */
	Read bufferToRead(long id){
		Read r=new Read(bases.toBytes(), (quals==null ? null : quals.toBytes()), 0);
		bases.clear();
		if(quals!=null){quals.clear();}
		r.id=(name==null || name.length()==0 ? ""+id : (addNumber ? name+" "+id : name));
		addNumber=true;
		return r;
	}
	
	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#processReadPair(stream.Read, stream.Read)
	 */
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
		if(fusePairs){
			fusePair(r1, r2);
			return true;
		}
		if(r1!=null && r1.length()>0){processRead(r1);}
		if(r2!=null && r2.length()>0){processRead(r2);}
		return false;
	}
	
	/**
	 * Fuses paired-end reads into a single sequence.
	 * Reverse complements R2 and joins it to R1 with N-padding between them.
	 * Quality scores are preserved if available from both reads.
	 *
	 * @param r1 First read (kept as-is)
	 * @param r2 Second read (reverse complemented before joining)
	 */
	private void fusePair(Read r1, Read r2) {
		if(r2==null){return;}
		
		r2.reverseComplementFast();
		final int len=r1.length()+r2.length()+npad;
		byte[] bases=new byte[len];
		byte[] quals=(r1.quality==null || r2.quality==null ? null : new byte[len]);
		
		for(int i=0, max=r1.length(); i<max; i++){
			bases[i]=r1.bases[i];
			if(quals!=null){quals[i]=r1.quality[i];}
		}
		for(int i=0, j=r1.length(); i<npad; i++, j++){
			bases[j]=(byte)PAD_SYMBOL;
		}
		for(int i=0, j=r1.length()+npad, max=r2.length(); i<max; i++, j++){
			bases[j]=r2.bases[i];
			if(quals!=null){quals[j]=r2.quality[i];}
		}
		
		r1.mate=r2.mate=null;
		r1.bases=bases;
		r1.quality=quals;
	}
	
	/**
	 * Adds a single read to the accumulation buffer.
	 * Appends padding symbols between reads and handles quality score assignment.
	 * Uses the first read's ID as the base name if none is specified.
	 * @param r Read to add to the buffer
	 */
	private void processRead(Read r) {
		if(name==null){
			name=r.id;
		}
		if(bases.length>0){
			for(int i=0; i<npad; i++){
				bases.append(PAD_SYMBOL);
				quals.append((byte)0);
			}
		}
		bases.append(r.bases);
		if(r.quality!=null){
			quals.append(r.quality);
		}else{
			for(int i=0, max=r.length(); i<max; i++){
				quals.append(defaultQuality);
			}
		}
	}
	
	@Override
	protected void startupSubclass() {}
	
	@Override
	protected void shutdownSubclass() {}
	
	@Override
	protected final boolean useSharedHeader(){return false;}
	
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn) {}
	
	/** Maximum length for output sequences before flushing buffer */
	int maxlen=Shared.MAX_ARRAY_LEN;
	/** Number of padding characters to insert between sequences */
	int npad;
	/** Default quality score for bases when input lacks quality information */
	byte defaultQuality;
	/** Whether to fuse paired-end reads into single sequences */
	boolean fusePairs;
	/** Buffer for accumulating sequence bases */
	ByteBuilder bases=new ByteBuilder();
	/** Buffer for accumulating quality scores */
	ByteBuilder quals=new ByteBuilder();
	/** Name prefix for output sequences */
	String name;
//	boolean prefix=true;
	/** Whether to append numeric suffixes to output sequence names */
	boolean addNumber=false;
	/**
	 * Character used for padding between sequences ('N' for DNA, 'X' for amino acids)
	 */
	char PAD_SYMBOL='N';
	/** Whether input sequences are amino acids rather than nucleotides */
	boolean amino=false;
	
}
