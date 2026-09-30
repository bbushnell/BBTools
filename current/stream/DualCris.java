package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;
import structures.ListNum;

/**
 * Combines batches from two optional ConcurrentReadInputStream backends.
 * Pairs overlapping entries by position, sets second-source pair numbers to one,
 * and retains source Read references. Coordinates backend calls without running
 * a producer thread of its own; callers must coordinate lifecycle and consumption.
 * SplitPairsAndSingles uses this wrapper for its two-input repair path, then clears
 * temporary mate links and rematches by name. Positional links are not name validation.
 *
 * CAVEAT: pairing is POSITIONAL (i-th R1 &lt;-&gt; i-th R2) and assumes the two streams stay
 * buffer-aligned; that can break for long/variable-length reads - see historical #001.
 * Return batches with the three-argument returnList and per-source presence flags.
 * The inherited ListNum return overload reaches the unsupported two-argument method.
 *
 * @author Brian Bushnell
 */
public class DualCris extends ConcurrentReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostic driver that prints batch IDs while tracking source presence.
	 * Retains its conditional null-wrapper accesses (#002); not a general validation tool.
	 * @param args Primary input path and optional secondary input path
	 */
	public static void main(String[] args){
		String a=args[0];
		String b=args.length>1 ? args[1] : null;
		FileFormat ff1=FileFormat.testInput(a, null, false);
		FileFormat ff2=(b==null ? null : FileFormat.testInput(b, null, false));
		DualCris cris=getReadInputStream(-1, false, ff1, ff2, null, null);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		boolean foundR1=false, foundR2=false;
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r1 : reads){
				Read r2=r1.mate;
				if(r1.pairnum()==0){foundR1=true;}else{foundR2=true;}
				if(r2!=null){
					if(r2.pairnum()==0){foundR1=true;}else{foundR2=true;}
				}
			}
			
			System.err.print(ln.id);
			
			cris.returnList(ln.id, foundR1, foundR2);
			foundR1=foundR2=false;
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			System.err.print(",");
		}
		System.err.print("Finished.");
		//TODO: Possible bug [stream/DualCris#002] - LOW (dev-only main): the while loop exits when ln becomes null (~L47), so this ln.id can NPE. Test driver only; no shell-script reach.
		cris.returnList(ln.id, foundR1, foundR2);
		ReadWrite.closeStreams(cris);
	}

	/** Creates independent single-input backends and retains them in a wrapper.
	 * @param maxReads Limit passed separately to each backend factory
	 * @param keepSamHeader Whether backend factories should retain SAM headers
	 * @param ff1 Optional first input descriptor
	 * @param ff2 Optional second input descriptor
	 * @param qf1 Optional separate quality input for the first backend
	 * @param qf2 Optional separate quality input for the second backend
	 * @return Unstarted wrapper; null descriptors leave the corresponding backend absent
	 */
	public static DualCris getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		ConcurrentReadInputStream cris1=(ff1==null ? null : ConcurrentReadInputStream.getReadInputStream(maxReads, keepSamHeader, ff1, null, qf1, null));
		ConcurrentReadInputStream cris2=(ff2==null ? null : ConcurrentReadInputStream.getReadInputStream(maxReads, keepSamHeader, ff2, null, qf2, null));
		return new DualCris(cris1, cris2);
	}
	
	/** Retains backends without starting them or checking their pairing compatibility.
	 * Builds the inherited filename from available backend names; active flags start false.
	 * @param cris1_ Optional first input stream
	 * @param cris2_ Optional second input stream
	 */
	public DualCris(ConcurrentReadInputStream cris1_, ConcurrentReadInputStream cris2_){
		super((cris1_==null ? "null" : cris1_.fname())+(cris2_==null ? "null" : ","+cris2_.fname()));
		cris1=cris1_;
		cris2=cris2_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Obtains one batch from each active backend and combines their Read references.
	 * A null backend result marks that source inactive; an empty wrapper alone does not.
	 * Sets all second-source reads to pair number one and links overlapping positions;
	 * first-source pair numbers are not changed.
	 * Extra second-source entries are appended to the first list; unmatched entries retain
	 * existing mate fields. No names, numeric IDs or cross-batch offsets are checked.
	 * @return First wrapper when present, otherwise the second; null when neither is present
	 */
	@Override
	public ListNum<Read> nextList(){
		
		ListNum<Read> ln1=null, ln2=null;
		if(cris1Active && cris1!=null){
			ln1=cris1.nextList();
			if(ln1==null){
				synchronized(this){
					cris1Active=false;
					System.err.println("\nSet cris1Active="+cris1Active);
				}
			}
		}
		if(cris2Active && cris2!=null){
			ln2=cris2.nextList();
			if(ln2!=null){
				for(Read r : ln2.list){r.setPairnum(1);}
			}else{
				synchronized(this){
					cris2Active=false;
					System.err.println("\nSet cris2Active="+cris2Active);
				}
			}
		}
		
		if(ln1!=null && ln2!=null){
			final int size1=ln1.size(), size2=ln2.size();
			final int min=Tools.min(size1, size2);
			//Historical #001 below concerns positional links, not a proved final-output failure.
			//SplitPairsAndSingles.process3_repair saves both references, then repair clears each
			//temporary mate link and rematches by name. Its behavior qualifies the older claim.
			//TODO: Possible bug [stream/DualCris#001] - LOW (latent; re-graded from MEDIUM after adversarial verify): positional pairing (ln1[i]<->ln2[i]) assumes cris1 and cris2 deliver BUFFER-ALIGNED lists. Holds whenever the 200-read cap (Shared.READ_BUFFER_LENGTH) binds = all short-read usage. Can break ONLY when the 400000-base cap (READ_BUFFER_MAX_DATA, applied in each sub-stream's crisG.readLists repack) binds first - avg read >2000bp - AND R1/R2 length distributions differ enough to fire it at DIFFERENT read counts; SYMMETRIC long reads (the norm) still stay aligned. When it does fire: buffers desync, reads mispair across boundaries, the size-mismatch handling below orphans the surplus mates - silent pair corruption. Trigger (asymmetric long paired reads through two-file repair.sh/bbsplitpairs.sh) essentially never occurs => LOW, but silent-when-fired => documented + FLAGGED FOR BRIAN (design assumption; see package summary). NOT auto-fixed (structural - needs coordinated/realigned buffering).
			for(int i=0; i<min; i++){
				Read r1=ln1.get(i);
				Read r2=ln2.get(i);
				r1.mate=r2;
				r2.mate=r1;
			}
			if(size2>size1){
				for(int i=size1; i<size2; i++){
					ln1.add(ln2.get(i));
				}
			}
		}else if(ln2!=null){
			ln1=ln2;
		}
		
		return ln1;
	}
	
	/** Rejects a single shared terminal flag; this wrapper requires per-source flags.
	 * @param listNum Unused batch ID
	 * @param poison Unused shared terminal flag
	 * @throws RuntimeException Always; use the three-argument overload
	 */
	@Override
	public void returnList(long listNum, boolean poison){//By design: DualCris needs the 3-arg returnList(listNum, foundR1, foundR2) to independently poison each sub-stream; the standard 2-arg form can't express that, so it throws. SplitPairsAndSingles.process3_repair calls the 3-arg version.
		throw new RuntimeException("Unsupported.");
	}
	
	/** Returns the supplied ID to each active backend, with independently derived terminal flags.
	 * Marks a source inactive when its presence flag is false. Determine flags from the
	 * consumed reads' pair numbers before changing them, as the repair caller does.
	 * @param listNum Same batch ID forwarded to each active backend
	 * @param foundR1 Whether the consumed batch contained first-source reads
	 * @param foundR2 Whether the consumed batch contained second-source reads
	 */
	public void returnList(long listNum, boolean foundR1, boolean foundR2){
		if(cris1!=null && cris1Active){
			cris1.returnList(listNum, !foundR1);
			if(!foundR1){cris1Active=false;}
		}
		if(cris2!=null && cris2Active){
			cris2.returnList(listNum, !foundR2);
			if(!foundR2){cris2Active=false;}
		}
	}
	
	/** Sets inherited started, starts each present backend and marks it active.
	 * Creates no producer thread for this wrapper; coordinate against other lifecycle calls.
	 */
	@Override
	public void start(){
		started=true;
		if(cris1!=null){
			cris1.start();
			cris1Active=true;
		}
		if(cris2!=null){
			cris2.start();
			cris2Active=true;
		}
	}
	
	/** Rejects direct execution with assertions enabled; otherwise performs no work.
	 * Use start to launch the backends, not a Thread wrapping this object.
	 */
	@Override
	public void run(){assert(false);}//DualCris has no producer thread of its own; start delegates to cris1/cris2 and never invokes this guard.
	
	/** Delegates shutdown to present backends, then marks both sources inactive. */
	@Override
	public void shutdown(){
		if(cris1!=null){cris1.shutdown();}
		if(cris2!=null){cris2.shutdown();}
		cris1Active=cris2Active=false;
	}
	
	/** Delegates restart and marks present backends active; does not itself call start.
	 * Retains this wrapper's accumulated error flag and inherited started state.
	 * Backend-specific restart contracts still apply.
	 */
	@Override
	public void restart(){
		if(cris1!=null){
			cris1.restart();
			cris1Active=true;
		}
		if(cris2!=null){
			cris2.restart();
			cris2Active=true;
		}
	}
	
	/** Closes present backends and marks both inactive; does not reset cached errors. */
	@Override
	public void close(){
		if(cris1!=null){cris1.close();}
		if(cris2!=null){cris2.close();}
		cris1Active=cris2Active=false;
	}
	
	/** Infers paired mode from the second backend's presence, otherwise asks the first.
	 * Requires at least one backend with assertions enabled; does not inspect active flags.
	 * @return true when a second backend exists, otherwise the first backend's paired mode
	 */
	@Override
	public boolean paired(){
		assert(cris1!=null || cris2!=null);
		if(cris2!=null){return true;}
		if(cris1!=null){return cris1.paired();}
		return false;
	}
	
	/** Flattens present backends' producer arrays, first source before second.
	 * @return Newly allocated array retaining each backend's producer references
	 */
	@Override
	public Object[] producers(){
		ArrayList<Object> list=new ArrayList<Object>();
		if(cris1!=null){
			for(Object o : cris1.producers()){list.add(o);}
		}
		if(cris2!=null){
			for(Object o : cris2.producers()){list.add(o);}
		}
		return list.toArray();
	}
	
	/** Accumulates each present backend's current error flag into this wrapper's flag.
	 * A prior true result remains true across close and restart; no completion barrier.
	 * @return Whether an error has been observed by this method
	 */
	@Override
	public boolean errorState(){
		if(cris1!=null){errorState|=cris1.errorState();}
		if(cris2!=null){errorState|=cris2.errorState();}
		return errorState;
	}
	
	/** Rejects wrapper-level sampling without changing either backend.
	 * @param rate Sampling rate (unused)
	 * @param seed Random seed (unused)
	 * @throws RuntimeException Always thrown as this method is invalid
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		throw new RuntimeException("Invalid.");
	}
	
	/** Returns the sum of present backends' observed base counters.
	 * @return Sum of current snapshots; no completion barrier
	 */
	@Override
	public long basesIn(){
		return (cris1==null ? 0 : cris1.basesIn())+(cris2==null ? 0 : cris2.basesIn());
	}
	
	/** Returns the sum of present backends' observed read counters.
	 * @return Sum in the backends' own units; no completion barrier
	 */
	@Override
	public long readsIn(){
		return (cris1==null ? 0 : cris1.readsIn())+(cris2==null ? 0 : cris2.readsIn());
	}
	
	/** Returns the local verbose flag, initialized false and not changed here. */
	@Override
	public boolean verbose(){return verbose;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained optional first backend. */
	private final ConcurrentReadInputStream cris1;
	/** Retained optional second backend. */
	private final ConcurrentReadInputStream cris2;
	/** Whether each corresponding source remains eligible for nextList/returnList calls. */
	private boolean cris1Active, cris2Active;
	/** Accumulated backend error observations; retained across lifecycle calls. */
	private boolean errorState=false;
	/** Local diagnostic flag with no setter in this implementation. */
	private boolean verbose=false;

}
