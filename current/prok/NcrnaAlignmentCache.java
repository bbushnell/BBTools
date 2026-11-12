package prok;

import java.util.Arrays;

import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.ObjectMap;
import shared.Shared;

/** Bounded, worker-owned reuse of the ordinary ncRNA identity alignments.
 * One cache binds one immutable model library. Exact window bytes, model index,
 * aligner/direction and Quantum's incoming offset identify a result. Neither
 * accepted calls nor mutable AlignmentStats objects are cached. Missing entries
 * execute the same production kernel, including its original operand order.
 *
 * A retained window owns primitive result arrays for the entire library; there
 * is no object or sequence copy per model attempt. A reusable lookup key avoids
 * allocations on hits. FIFO replacement bounds retained windows and each window
 * is limited to maxWindowLength bases. Larger requests execute without retention.
 * Clear between input contigs to release old sequence references. This class is
 * not thread-safe; all use, including clear(), belongs to one worker.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaAlignmentCache{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows an immutable library; zero windows provides the uncached control.
	 * @param library_ Model arrays, kept unchanged for the cache's entire lifetime
	 * @param maxWindows Maximum retained distinct sequence windows
	 * @param maxWindowLength_ Largest sequence eligible for retention, in bases
	 */
	public NcrnaAlignmentCache(final byte[][] library_, final int maxWindows,
			final int maxWindowLength_){
		if(library_==null || library_.length<1 || library_.length>Integer.MAX_VALUE/7
				|| maxWindows<0 || maxWindowLength_<1){
			throw new IllegalArgumentException("Alignment reuse needs a model library and nonnegative window budget");
		}
		for(byte[] model:library_){
			if(model==null || model.length<1){throw new IllegalArgumentException("Cached alignment models must be nonempty");}
		}
		library=library_;
		maxWindowLength=maxWindowLength_;
		windows=new Window[maxWindows];
		simd=Shared.SIMD;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Alignment           ----------------*/
	/*--------------------------------------------------------------*/

	/** Same model-to-window Quantum call as NcrnaScavenger's primary verifier.
	 * positions is both input and output: its incoming first element affects the
	 * bandwidth through IDAlignerStatics.decideBandwidth. The four output integers
	 * are copied, never lent to the caller.
	 */
	public float quantum(final int model, final byte[] sequence, final int[] positions){
		if(positions==null || positions.length!=4){
			throw new IllegalArgumentException("The primary Quantum cache requires its production four-position buffer");
		}
		return align(model, sequence, positions, QUANTUM);
	}

	/** Null positions preserves Scrabble's sequence-to-model call and shorter-query
	 * swapping. A two-element buffer preserves the explicit model-to-window call.
	 */
	public float scrabble(final int model, final byte[] sequence, final int[] positions){
		if(positions!=null && positions.length!=2){
			throw new IllegalArgumentException("Observed Scrabble coordinates require the production two-position buffer");
		}
		return align(model, sequence, positions, positions==null ? SCRABBLE : SCRABBLE_POSITIONED);
	}

	private float align(final int model, final byte[] sequence, final int[] positions, final int mode){
		if(model<0 || model>=library.length || sequence==null || sequence.length<1){
			throw new IllegalArgumentException("Alignment reuse requires a known model and a nonempty sequence");
		}
		if(simd!=Shared.SIMD){
			clear();
			simd=Shared.SIMD;
		}
		requests++;
		final boolean visualizing=mode==QUANTUM ? QuantumAligner.output!=null : ScrabbleAligner.output!=null;
		final Window window=visualizing ? null : window(sequence);
		final int incomingOffset=mode==QUANTUM ? positions[0] : 0;
		final int bit=1<<mode;
		final int valueIndex=mode*library.length+model;
		final int positionIndex=model*7;
		if(window!=null && (window.valid[model]&bit)!=0
				&& (mode!=QUANTUM || window.positions[positionIndex]==incomingOffset)){
			hits++;
			final NcrnaAlignmentWork work=NcrnaAlignmentWork.current();
			if(work!=null){work.reused(library[model].length, sequence.length);}
			if(mode==QUANTUM){System.arraycopy(window.positions, positionIndex+1, positions, 0, 4);}
			else if(mode==SCRABBLE_POSITIONED){System.arraycopy(window.positions, positionIndex+5, positions, 0, 2);}
			return window.identities[valueIndex];
		}
		final long start=NcrnaCpuProfile.cpuNow();
		final float identity;
		if(mode==QUANTUM){identity=NcrnaAlignmentWork.quantum(library[model], sequence, positions, library[model].length);}
		else if(mode==SCRABBLE_POSITIONED){identity=NcrnaAlignmentWork.scrabble(library[model], sequence, positions, library[model].length);}
		else{identity=NcrnaAlignmentWork.scrabble(sequence, library[model], null, library[model].length);}
		if(start>=0){kernelCpuNanos+=NcrnaCpuProfile.cpuNow()-start;}
		fills++;
		if(window!=null){
			window.identities[valueIndex]=identity;
			if(mode==QUANTUM){
				window.positions[positionIndex]=incomingOffset;
				System.arraycopy(positions, 0, window.positions, positionIndex+1, 4);
			}else if(mode==SCRABBLE_POSITIONED){System.arraycopy(positions, 0, window.positions, positionIndex+5, 2);}
			window.valid[model]|=bit;
		}
		assert(requests==hits+fills) : "Every completed cache request either reuses one result or executes one production kernel";
		return identity;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Bounded Storage        ----------------*/
	/*--------------------------------------------------------------*/

	private Window window(final byte[] sequence){
		if(windows.length==0 || sequence.length>maxWindowLength){return null;}
		probe.set(sequence);
		final Window existing=index.get(probe);
		if(existing!=null){return existing;}
		Window next=windows[cursor];
		if(next==null){
			next=new Window(library.length);
			windows[cursor]=next;
		}else{
			if(next.bases!=null){
				final Window removed=index.remove(next);
				assert(removed==next) : "FIFO eviction must remove the old immutable key before replacing its sequence";
				evictions++;
			}
			Arrays.fill(next.valid, (byte)0);
		}
		next.set(sequence.clone());
		index.put(next, next);
		cursor++;
		if(cursor==windows.length){cursor=0;}
		return next;
	}

	/** Releases sequences and invalidates results; retains primitive scratch and counters. */
	public void clear(){
		index.clear();
		for(Window window:windows){
			if(window!=null){
				window.bases=null;
				Arrays.fill(window.valid, (byte)0);
			}
		}
		probe.bases=null;
		cursor=0;
	}

	/** Identity binding avoids applying model indices from another loaded library. */
	public boolean usesLibrary(final byte[][] models){return library==models;}
	public long requests(){return requests;}
	public long hits(){return hits;}
	public long fills(){return fills;}
	public long evictions(){return evictions;}
	/** Physically executed kernel CPU only; zero when CPU profiling is disabled. */
	public long kernelCpuNanos(){return kernelCpuNanos;}

	/*--------------------------------------------------------------*/
	/*----------------         Window State         ----------------*/
	/*--------------------------------------------------------------*/

	/** Equality uses every byte; a hash collision never supplies alignment evidence. */
	private static final class Window{
		Window(final int models){
			identities=models==0 ? null : new float[models*3];
			positions=models==0 ? null : new int[models*7];
			valid=models==0 ? null : new byte[models];
		}
		void set(final byte[] sequence){
			assert(sequence!=null) : "Stored and probing keys require the exact candidate sequence";
			bases=sequence;
			hash=Arrays.hashCode(sequence);
		}
		@Override
		public int hashCode(){return hash;}
		@Override
		public boolean equals(final Object other){
			return other==this || other instanceof Window && Arrays.equals(bases, ((Window)other).bases);
		}
		byte[] bases;
		int hash;
		final float[] identities;
		final int[] positions;
		final byte[] valid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final byte[][] library;
	private final int maxWindowLength;
	private final Window[] windows;
	private final ObjectMap<Window,Window> index=new ObjectMap<Window,Window>(Window.class, Window.class);
	private final Window probe=new Window(0);
	private int cursor;
	private boolean simd;
	private long requests, hits, fills, evictions, kernelCpuNanos;
	private static final int QUANTUM=0, SCRABBLE=1, SCRABBLE_POSITIONED=2;
}
