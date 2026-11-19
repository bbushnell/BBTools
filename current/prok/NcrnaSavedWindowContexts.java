package prok;

import java.util.ArrayList;
import java.util.Comparator;

import dna.AminoAcid;
import map.LongObjectMap;
import map.ObjectSet;
import stream.Read;

/** One exact original input file's worker-owned saved contexts. Candidate bytes
 * are always present; zero marks unavailable flank bytes, distinct from real Ns.
 * Overlapping saved candidates supply either strand before work starts. Only an
 * actual missing refinement read opens the original file, once, streaming each
 * contig and checking every retained candidate/context byte before filling holes.
 * No original contig or whole genome is retained after its sequence callback.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSavedWindowContexts implements NcrnaScavenger.SavedContextSource{

	public NcrnaSavedWindowContexts(final String input_, final String fasta_){
		require(input_!=null && !input_.isEmpty() && fasta_!=null && !fasta_.isEmpty(),
			"A context owner must retain the exact original input identity and FASTA path");
		input=input_;
		fasta=fasta_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Metadata Owner        ----------------*/
	/*--------------------------------------------------------------*/

	/** The caller has already joined this event to its exact input/frame/path.
	 * Candidate bytes are borrowed during this call, then released by the owner. */
	public Window add(final String event, final long ordinal, final String name,
			final int contigLength, final int strand, final int start, final int stop,
			final int flank, final byte[] candidate){
		require(!prepared && event!=null && events.add(event) && ordinal>=0 && name!=null
			&& name.equals(NcrnaDdlOracle.seqid(name)) && contigLength>0 && (strand==0 || strand==1)
			&& start>=0 && stop>=start && stop<contigLength && flank>=0 && candidate!=null
			&& candidate.length==stop-start+1,
			"Each retained event needs unique input-bound identity, original contig geometry and complete candidate bytes");
		Contig contig=contigs.get(ordinal);
		if(contig==null){
			contig=new Contig(ordinal, name, contigLength);
			contigs.put(ordinal, contig);
			ordered.add(contig);
		}else{require(contig.name.equals(name) && contig.length==contigLength,
			"One original input/contig ordinal cannot change its native name or length");}
		final Window window=new Window(this, event, contig, strand, start, stop, flank);
		for(int i=0; i<candidate.length; i++){
			final byte base=upper(candidate[i]);
			require(base!=0, "A captured candidate cannot contain an unavailable-byte marker");
			window.bases[start-window.origin+i]=base;
		}
		contig.windows.add(window);
		windows.add(window);
		return window;
	}

	/** Optional new-capture sidecar, bound to the same exact candidate/context.
	 * Old captures omit it and recover actual missing reads lazily instead. */
	public void supply(final Window window, final int origin, final byte[] context){
		require(!prepared && window!=null && window.owner==this && origin==window.origin
			&& context!=null && context.length==window.bases.length,
			"A saved context sidecar must bind the original event and complete declared refinement interval");
		for(int i=0; i<context.length; i++){put(window, i, upper(context[i]));}
	}

	/** Freeze metadata, then fill context only from saved candidates in this exact
	 * input/contig. Prefix maximum ends make overlap lookup bounded by candidates
	 * that can actually reach the requested interval, including nested windows. */
	public void prepare(){
		require(!prepared, "Prepare each original-file context collection once before its worker starts");
		for(final Contig contig:ordered){
			contig.windows.sort(FORWARD_ORDER);
			final int[] ends=new int[contig.windows.size()];
			for(int i=0; i<ends.length; i++){
				ends[i]=Math.max(i==0 ? -1 : ends[i-1], contig.windows.get(i).forwardStop);
			}
			for(final Window target:contig.windows){
				final int from=target.strand==0 ? target.origin : contig.length-target.origin-target.bases.length;
				final int to=from+target.bases.length-1;
				int low=0, high=ends.length;
				while(low<high){
					final int middle=(low+high)>>>1;
					if(ends[middle]<from){low=middle+1;}else{high=middle;}
				}
				for(int i=low; i<ends.length; i++){
					final Window source=contig.windows.get(i);
					if(source.forwardStart>to){break;}
					final int left=Math.max(from, source.forwardStart), right=Math.min(to, source.forwardStop);
					for(int pos=left; pos<=right; pos++){
						final int targetPosition=(target.strand==0 ? pos : contig.length-1-pos)-target.origin;
						final int sourcePosition=(source.strand==0 ? pos : contig.length-1-pos)-source.origin;
						byte base=source.bases[sourcePosition];
						assert(base!=0) : "Only the source candidate interval, never its unavailable flank, supplies overlap context";
						if(source.strand!=target.strand){base=complement(base);}
						put(target, targetPosition, base);
					}
				}
			}
		}
		prepared=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input Worker         ----------------*/
	/*--------------------------------------------------------------*/

	/** Activate synchronously immediately before this window's verifier call.
	 * The same source object is reused; no per-window task or callback allocation. */
	public NcrnaScavenger.SavedContextSource sourceFor(final Window window){
		final long thread=Thread.currentThread().getId();
		if(owner<0){owner=thread;}
		require(prepared && owner==thread && window!=null && window.owner==this,
			"One input-file worker owns its exact saved windows and all later source reads");
		active=window;
		return this;
	}

	@Override
	public void fill(final byte[] context, final int origin, final int from, final int to){
		require(owner==Thread.currentThread().getId() && active!=null && active.bases==context
			&& active.origin==origin && from>=origin && to>=from && (long)to<origin+(long)context.length,
			"A refinement request must belong to the active worker, original input and bounded context array");
		boolean missing=false;
		for(int i=from-origin; i<=to-origin; i++){if(context[i]==0){missing=true;break;}}
		if(!missing){return;}
		if(!sourceAttempted){
			sourceAttempted=true;
			sourceLoads++;
			NcrnaDdlSavedData.sequences(fasta, this::original);
			for(final Contig contig:ordered){require(contig.seen,
				"The exact original FASTA lacks captured contig ordinal "+contig.ordinal+" for input "+input);}
		}
		for(int i=from-origin; i<=to-origin; i++){
			require(context[i]!=0, "The exact original source did not supply every requested refinement base");
		}
	}

	/** Streamer owns the Read and releases it after this callback. Even contigs
	 * without saved windows contribute to the explicit input-I/O accounting. */
	private void original(final Read read){
		assert(sourceAttempted && owner==Thread.currentThread().getId() && read!=null && read.bases!=null)
			: "Only the active original-file worker may consume a borrowed native source contig";
		final Contig contig=contigs.get(sourceContigs++);
		sourceBases+=read.length();
		if(contig==null){return;}
		require(!contig.seen && contig.length==read.length() && contig.name.equals(NcrnaDdlOracle.seqid(read.id)),
			"Original source ordinal/name/length disagree with captured input "+input+", contig ordinal "+contig.ordinal);
		for(final Window window:contig.windows){
			for(int i=0; i<window.bases.length; i++){
				final int oriented=window.origin+i;
				byte base=upper(read.bases[window.strand==0 ? oriented : contig.length-1-oriented]);
				if(window.strand==1){base=complement(base);}
				put(window, i, base);
			}
		}
		contig.seen=true;
	}

	private static void put(final Window window, final int position, final byte base){
		assert(position>=0 && position<window.bases.length)
			: "Orientation conversion must preserve the retained bounded context interval";
		if(base==0 || window.bases[position]!=0 && window.bases[position]!=base){
			throw new IllegalArgumentException("Saved/original context bytes disagree for input "+window.owner.input
				+", event "+window.event+", oriented position "+(window.origin+position));
		}
		window.bases[position]=base;
	}

	private static byte upper(final byte base){return base>='a' && base<='z' ? (byte)(base-32) : base;}
	private static byte complement(final byte base){
		require(base>0 && base<AminoAcid.baseToComplementExtended.length,
			"Original and captured FASTA symbols must use the native extended DNA complement table");
		return upper(AminoAcid.baseToComplementExtended[base]);
	}

	public long sourceLoads(){return sourceLoads;}
	public long sourceBases(){return sourceBases;}
	public int windows(){return windows.size();}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Metadata      ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Window{
		private Window(final NcrnaSavedWindowContexts owner_, final String event_, final Contig contig,
				final int strand_, final int start, final int stop, final int flank){
			owner=owner_;
			event=event_;
			strand=strand_;
			origin=Math.max(0, start-flank);
			final int end=(int)Math.min((long)contig.length-1, (long)stop+flank);
			bases=new byte[end-origin+1];
			forwardStart=strand==0 ? start : contig.length-1-stop;
			forwardStop=strand==0 ? stop : contig.length-1-start;
		}
		public final String event;
		public final int strand, origin;
		/** Mutable only by this input's context owner; borrowed by its one verifier. */
		public final byte[] bases;
		private final int forwardStart, forwardStop;
		private final NcrnaSavedWindowContexts owner;
	}

	private static final class Contig{
		Contig(final long ordinal_, final String name_, final int length_){ordinal=ordinal_;name=name_;length=length_;}
		final long ordinal;
		final String name;
		final int length;
		final ArrayList<Window> windows=new ArrayList<Window>();
		boolean seen;
	}
	private static final Comparator<Window> FORWARD_ORDER=(a, b)->{
		final int order=Integer.compare(a.forwardStart, b.forwardStart);
		return order!=0 ? order : a.event.compareTo(b.event);
	};
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	private final String input, fasta;
	private final LongObjectMap<Contig> contigs=new LongObjectMap<Contig>(Contig.class);
	private final ObjectSet<String> events=new ObjectSet<String>(String.class);
	private final ArrayList<Contig> ordered=new ArrayList<Contig>();
	private final ArrayList<Window> windows=new ArrayList<Window>();
	private boolean prepared, sourceAttempted;
	private long owner=-1, sourceLoads, sourceContigs, sourceBases;
	private Window active;
}
