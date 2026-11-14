package synth;

import java.util.ArrayList;

import parse.Parse;
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
 * Fuses sequences together with configurable padding between them.
 * Can operate in two modes: concatenating reads into long sequences with N-padding,
 * or fusing paired-end reads by reverse complementing R2 and joining with R1.
 * Supports amino acid sequences and customizable padding characters.
 *
 * @author Brian Bushnell
 * @date Jan 20, 2015
 */
public final class FuseSequence extends BBTool_ST {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FuseSequence fs=new FuseSequence(args);
		fs.process(t);
	}
	
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
	public boolean parseArgument(String arg, String a, String b){
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
	
	/**
	 * Consumes input reads and emits either fused pairs or buffered concatenations.
	 * Pair-fusing mode keeps one output read per retained pair, while concatenation
	 * mode appends reads to a shared buffer and flushes it at positive length limits
	 * or at end of input.
	 *
	 * @param cris Input stream for reading sequences
	 * @param ros Output stream for writing processed sequences
	 */
	@Override
	protected void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		if(discardOutputQualities()){
			quals=null;
		}

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
					
					// Two nonempty mates add two gaps when joining a nonempty buffer.
					final long padding=(long)npad*(initialLength1>0 && initialLength2>0 ? 2 : 1);
					if(!fusePairs && maxlen>0 && bases.length>0 && (long)bases.length+initialLength1+initialLength2+padding>maxlen){
						Read r=bufferToRead(outNum);
						outNum++;
						readsOut.add(r);
					}
					processReadPair(r1, r2);
					if(fusePairs){
						readsOut.add(r1);
					}else if(maxlen>0 && bases.length>=maxlen){
						// Nonpositive limits retain all reads for the final unlimited-mode flush.
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
	 * Converts the accumulated bases and qualities into an output read.
	 * @param id numeric suffix to append when output numbering is enabled
	 * @return a read containing the previous buffer contents
	 */
	Read bufferToRead(long id){
		Read r=new Read(bases.toBytes(), (quals==null ? null : quals.toBytes()), 0);
		bases.clear();
		if(quals!=null){quals.clear();}
		r.id=(name==null || name.length()==0 ? ""+id : (addNumber ? name+" "+id : name));
		addNumber=true;
		return r;
	}
	
	/**
	 * Handles one primary read and its optional mate according to the active mode.
	 * @param r1 first read in the pair
	 * @param r2 second read in the pair, or null for unpaired input
	 * @return true when the primary read should stay in the current output list
	 */
	@Override
	protected boolean processReadPair(Read r1, Read r2){
		if(fusePairs){
			fusePair(r1, r2);
			return true;
		}
		if(r1!=null && r1.length()>0){processRead(r1);}
		if(r2!=null && r2.length()>0){processRead(r2);}
		return false;
	}
	
	/**
	 * Joins a primary read, padding and the reverse-complemented mate in place.
	 * @param r1 primary read to receive the fused bases and qualities
	 * @param r2 mate to reverse-complement and append, or null to leave r1 unchanged
	 */
	private void fusePair(Read r1, Read r2){
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
	 * Appends one read to the concatenation buffer, inserting padding first when
	 * previous reads are already buffered.
	 * @param r read to append
	 */
	private void processRead(Read r){
		if(name==null){
			name=r.id;
		}
		if(bases.length>0){
			for(int i=0; i<npad; i++){
				bases.append(PAD_SYMBOL);
				if(quals!=null){quals.append((byte)0);}
			}
		}
		bases.append(r.bases);
		if(quals!=null){
			if(r.quality!=null){
				quals.append(r.quality);
			}else{
				for(int i=0, max=r.length(); i<max; i++){
					quals.append(defaultQuality);
				}
			}
		}
	}
	
	/**
	 * Indicates whether concatenated reads need an inline quality buffer.
	 * FASTA primary output discards read qualities unless separate qfout output is active.
	 * @return true if FuseSequence can skip synthetic quality construction
	 */
	private boolean discardOutputQualities(){
		return ffout1!=null && ffout1.fasta() && qfout1==null;
	}

	/** No startup work is needed after the shared streams are allocated. */
	@Override
	protected void startupSubclass(){}
	
	/** No cleanup is needed before the shared streams are closed. */
	@Override
	protected void shutdownSubclass(){}
	
	/** Indicates whether to pass a shared SAM/BAM header to output streams.
	 * @return false; FuseSequence writes ordinary sequence output */
	@Override
	protected final boolean useSharedHeader(){return false;}
	
	/** FuseSequence adds no tool-specific summary statistics. */
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn){}
	
	int maxlen=Shared.MAX_ARRAY_LEN;
	int npad;
	byte defaultQuality;
	boolean fusePairs;
	ByteBuilder bases=new ByteBuilder();
	ByteBuilder quals=new ByteBuilder();
	String name;
//	boolean prefix=true;
	boolean addNumber=false;
	char PAD_SYMBOL='N';
	boolean amino=false;
	
}
