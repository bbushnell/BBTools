package stream;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Set;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Lazily creates and retains one output stream per name until explicit teardown.
 * Names substitute into primary/optional-mate patterns; there is no destination cap
 * or rotation. The maxSize option configures each child stream's buffering.
 * Creation and cached lookup share the map's monitor; collection traversal does not.
 * Child submissions occur outside that monitor. Lifecycle/status traversal
 * requires callers to stop creating streams and coordinate use of the exposed collections.
 * Ordered children require complete per-destination list IDs, including empty batches.
 *
 * @author Brian Bushnell
 * @date Apr 12, 2015
 */
public class MultiCros{
	
	/** Copies every input root to each supplied name using an unordered bucket set.
	 * @param args Positional input file, output pattern and zero or more destination names
	 * @throws RuntimeException If input or output teardown reports an error */
	public static void main(String[] args){
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){
			names.add(args[i]);
		}
		final int buff=Tools.max(16, 2*Shared.threads());
		MultiCros mcros=new MultiCros(pattern, null, false, false, false, false, false, FileFormat.FASTQ, buff);
		
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		ArrayListSet als=new ArrayListSet(false);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				als.add(r1, names);
			}
			cris.returnList(ln);
			if(mcros!=null){mcros.add(als, ln.id);}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		//STR250: all buckets were drained in the loop; no terminal submission needs ln.id.
		//STR251: evaluate both teardown results before deciding whether to fail.
		final boolean inputError=ReadWrite.closeStreams(cris);
		final boolean outputError=ReadWrite.closeStreams(mcros);
		if(inputError || outputError){
			throw new RuntimeException("MultiCros encountered an input or output error; output may be incomplete.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captures output policy without creating or starting child streams.
	 * A hash expands into primary/mate patterns when pattern2 is absent.
	 * @param pattern1_ Primary pattern containing a percent placeholder
	 * @param pattern2_ Optional mate pattern, also requiring a percent placeholder
	 * @param ordered_ Request child ordering; caller supplies complete IDs per destination
	 * @param overwrite_ Forwarded to each child output descriptor
	 * @param append_ Forwarded to each child output descriptor
	 * @param allowSubprocess_ Output descriptor subprocess permission
	 * @param useSharedHeader_ Request shared headers when creating each child
	 * @param defaultFormat_ Fallback output format
	 * @param maxSize_ Child stream buffering parameter, not a destination-count limit */
	public MultiCros(String pattern1_, String pattern2_,
			boolean ordered_, boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, int maxSize_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern2_.indexOf('%')>=0);//#003 FIXED: was pattern1_.indexOf (redundant, never validated pattern2_) - a non-null pattern2_ without '%' collapses all R2 output to one file. Twin of BufferedMultiCros#001.
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}
		
		ordered=ordered_;
		overwrite=overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		useSharedHeader=useSharedHeader_;
		
		defaultFormat=defaultFormat_;
		maxSize=maxSize_;

		streamList=new ArrayList<ConcurrentReadOutputStream>();
		streamMap=new LinkedHashMap<String, ConcurrentReadOutputStream>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Detaches and submits each nonnull named bucket, retaining the set's name metadata.
	 * Does not synthesize empty submissions for absent buckets or renumber list IDs.
	 * @param set Caller-owned set; each producer must coordinate its own bucket mutation
	 * @param listnum ID forwarded unchanged to every submitted destination */
	public void add(ArrayListSet set, long listnum){
		for(String s : set.getNames()){
			ArrayList<Read> list=set.getAndClear(s);
			//TODO: Compatibility issue - STR249: sparse buckets omit list IDs, while ordered
			//generic children require contiguous IDs starting at zero for every destination.
			if(list!=null){
				add(list, listnum, s);
			}
		}
	}
		
	/** Submits a list to the lazily opened named child; Read payloads remain borrowed.
	 * @param list Nonnull read list passed to the child
	 * @param listnum Child ordering ID; ordered output needs a complete sequence per name
	 * @param name Replacement string used for the first percent in each pattern */
	public void add(ArrayList<Read> list, long listnum, String name){
		ConcurrentReadOutputStream ros=getStream(name);
		ros.add(list, listnum);
	}
	
	/** Requests closure of every created child without joining; submissions/creation must stop first. */
	public void close(){
		for(ConcurrentReadOutputStream cros : streamList){cros.close();}
	}
	
	/** Joins every created child; does not request closure or prevent new stream creation. */
	public void join(){
		for(ConcurrentReadOutputStream cros : streamList){cros.join();}
	}
	
	/** Delegates ordering reset to current children, which may wait for pending ordered data.
	 * Requires ordered children, a quiescent destination collection and coordinated submissions. */
	public void resetNextListID(){
		for(ConcurrentReadOutputStream cros : streamList){cros.resetNextListID();}
	}
	
	/** @return Primary output pattern, rather than one resolved destination filename */
	public String fname(){return pattern1;}
	
	/** Aggregates current child error flags and this object's retained flag without waiting.
	 * @return true if any observed flag indicates an error */
	public boolean errorState(){
		//#001 FIX [stream/MultiCros#001]: was `b=b&&cros.errorState()` seeded from errorState (a field never set true in this class), so it ALWAYS returned false - silently masking any sub-stream error. An error aggregator must OR, not AND (cf. ConcurrentGenericReadOutputStream.errorState; finishedSuccessfully() below correctly uses &&). The live demux/seal callers detect errors via ReadWrite.closeStreams(mc), which ORs the sub-streams directly, so this broken method was NOT biting them - but it is public and must be correct for any direct caller. Latent LOW; fixed anyway (one char, zero perf, called only at teardown).
		boolean b=errorState;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b||cros.errorState();
		}
		return b;
	}

	/** Combines current child completion flags without closing or joining them.
	 * Does not include this object's separate errorState field; an empty collection returns true.
	 * @return true if every observed child reports successful completion */
	public boolean finishedSuccessfully(){
		boolean b=true;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b&&cros.finishedSuccessfully();
		}
		return b;
	}
	
	/** @return Live name view; treat as read-only and traverse after creation has stopped */
	public Set<String> getKeys(){return streamMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates descriptors and an unstarted child for one replacement name.
	 * @param name Replacement string used by String.replaceFirst, not a quoted literal
	 * @return New child stream with forwarded output policy */
	private ConcurrentReadOutputStream makeStream(String name){
		//Substitute the name into the '%' placeholder of the output pattern(s) to form this name's actual file(s), then open a cros for it.
		String s1=pattern1.replaceFirst("%", name);
		String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
		final FileFormat ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, overwrite, append, ordered);
		final FileFormat ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, overwrite, append, ordered);
		ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, maxSize, null, useSharedHeader);
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	//Lazy per-name output stream: created once on first sight of a name, cached in streamMap (and streamList for teardown).
	//STR248 (legacy stream/MultiCros#002): Seal shares this object across workers.
	//Lookup and publication must use the same monitor for the plain LinkedHashMap;
	//locking only creation left cache hits racing with insertion. Child add stays outside.
	/** Returns a cached child or creates, starts and registers a new one.
	 * Lookup and creation hold the map monitor; callers still coordinate lifecycle operations
	 * and must not mutate the exposed map/list independently.
	 * @param name Nonnull replacement name suitable for both output patterns
	 * @return Started output stream retained for later lifecycle operations */
	public ConcurrentReadOutputStream getStream(String name){
		synchronized(streamMap){
			ConcurrentReadOutputStream ros=streamMap.get(name);
			if(ros==null){
				ros=makeStream(name);
				ros.start();
				streamList.add(ros);
				streamMap.put(name, ros);
			}
			return ros;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resolved primary/optional-mate patterns after any hash expansion. */
	public final String pattern1, pattern2;
	/** Live creation-order child list; callers must not mutate it independently of the map. */
	public final ArrayList<ConcurrentReadOutputStream> streamList;
	/** Live name-to-child map; exposed for compatibility, with no external mutation protection. */
	public final LinkedHashMap<String, ConcurrentReadOutputStream> streamMap;
	/** Child ordering policy; does not add missing per-destination IDs. */
	public final boolean ordered;
	
	/** Retained local error flag; this class does not currently set it true. */
	boolean errorState=false;
	/** Unused retained lifecycle field; children start individually in getStream. */
	boolean started=false;
	/** Child descriptor overwrite policy. */
	final boolean overwrite;
	/** Child descriptor append policy. */
	final boolean append;
	/** Child descriptor subprocess permission. */
	final boolean allowSubprocess;
	/** Fallback format for each destination. */
	final int defaultFormat;
	/** Buffering parameter forwarded to each child factory invocation. */
	final int maxSize;
	/** Shared-header option forwarded when each child is constructed. */
	final boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Unused diagnostic option retained for compatibility. */
	public static boolean verbose=false;

}
