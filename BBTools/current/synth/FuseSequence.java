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
	
	private void fusePair(Read r1, Read r2) {
		if(r2==null){return;}
		
		r2.reverseComplement();
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
