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
 * Allows output of reads to multiple different output streams.
 * @author Brian Bushnell
 * @date Apr 12, 2015
 *
 */
public class MultiCros {
	
	/**
	 * Command-line entry point for MultiCros functionality.
	 * Processes reads from input file and distributes them to multiple output
	 * streams based on the provided pattern and names.
	 * @param args Command line arguments: input_file pattern name1 name2 ...
	 */
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
		if(mcros!=null){mcros.add(als, ln.id);}
		ReadWrite.closeStreams(cris);
		ReadWrite.closeStreams(mcros);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a MultiCros with specified output patterns and configuration.
	 * Handles pattern processing including automatic paired-end pattern generation
	 * when '#' placeholder is used instead of separate patterns.
	 * @param pattern1_ Primary output file pattern with '%' placeholder for names
	 * @param pattern2_ Secondary output file pattern (may be null for single-end)
	 * @param ordered_ Whether to maintain read order in output streams
	 * @param overwrite_ Whether to overwrite existing output files
	 * @param append_ Whether to append to existing output files
	 * @param allowSubprocess_ Whether to allow subprocess creation for compression
	 * @param useSharedHeader_ Whether output streams share common headers
	 * @param defaultFormat_ Default file format for output streams
	 * @param maxSize_ Maximum size for output stream buffers
	 */
	public MultiCros(String pattern1_, String pattern2_,
			boolean ordered_, boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, int maxSize_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0);
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}
		
		ordered=ordered_;
		overwrite= overwrite_;
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
	
	
	/**
	 * Adds read sets from ArrayListSet to appropriate output streams.
	 * Iterates through all names in the set and distributes reads to
	 * corresponding streams based on name-pattern matching.
	 * @param set Container holding read lists organized by names
	 * @param listnum List identifier for maintaining read order
	 */
	public void add(ArrayListSet set, long listnum){
		for(String s : set.getNames()){
			ArrayList<Read> list=set.getAndClear(s);
			if(list!=null){
				add(list, listnum, s);
			}
		}
	}
		
	/**
	 * Adds a read list to the output stream corresponding to the given name.
	 * Creates the stream if it doesn't exist using the configured patterns.
	 * @param list Read list to add to output stream
	 * @param listnum List identifier for maintaining read order
	 * @param name Name used to identify target output stream
	 */
	public void add(ArrayList<Read> list, long listnum, String name){
		ConcurrentReadOutputStream ros=getStream(name);
		ros.add(list, listnum);
	}
	
	/** Closes all managed output streams.
	 * Should be called to ensure proper cleanup and file finalization. */
	public void close(){
		for(ConcurrentReadOutputStream cros : streamList){cros.close();}
	}
	
	/** Waits for all output streams to complete their operations.
	 * Blocks until all background writing threads have finished. */
	public void join(){
		for(ConcurrentReadOutputStream cros : streamList){cros.join();}
	}
	
	/** Resets the next list ID counter for all managed output streams.
	 * Used to restart list numbering sequence. */
	public void resetNextListID(){
		for(ConcurrentReadOutputStream cros : streamList){cros.resetNextListID();}
	}
	
	/** Returns the primary output file pattern */
	public String fname(){return pattern1;}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){
		boolean b=errorState;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b&&cros.errorState();
		}
		return b;
	}

	/** Checks if all managed output streams finished successfully.
	 * @return true if all streams completed without errors */
	public boolean finishedSuccessfully(){
		boolean b=true;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b&&cros.finishedSuccessfully();
		}
		return b;
	}
	
	/** Returns the set of stream names currently managed */
	public Set<String> getKeys(){return streamMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a new ConcurrentReadOutputStream for the specified name.
	 * Substitutes the name into the configured patterns and creates
	 * appropriate FileFormat objects for the output streams.
	 * @param name Name to substitute into output patterns
	 * @return New configured ConcurrentReadOutputStream
	 */
	private ConcurrentReadOutputStream makeStream(String name){
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
	
	/**
	 * Gets or creates an output stream for the specified name.
	 * Uses double-checked locking to ensure thread-safe stream creation.
	 * Newly created streams are automatically started and added to management lists.
	 * @param name Name identifying the desired output stream
	 * @return ConcurrentReadOutputStream for the specified name
	 */
	public ConcurrentReadOutputStream getStream(String name){
		ConcurrentReadOutputStream ros=streamMap.get(name);
		if(ros==null){
			synchronized(streamMap){
				ros=streamMap.get(name);
				if(ros==null){
					ros=makeStream(name);
					ros.start();
					streamList.add(ros);
					streamMap.put(name, ros);
				}
			}
		}
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String pattern1, pattern2;
	/** List of all managed output streams for iteration and bulk operations */
	public final ArrayList<ConcurrentReadOutputStream> streamList;
	/** Map from stream names to output streams for fast lookup */
	public final LinkedHashMap<String, ConcurrentReadOutputStream> streamMap;
	/** Whether to maintain read order in output streams */
	public final boolean ordered;
	
	/** Whether this MultiCros instance has encountered an error */
	boolean errorState=false;
	/** Whether this MultiCros instance has been started */
	boolean started=false;
	/** Whether to overwrite existing output files */
	final boolean overwrite;
	/** Whether to append to existing output files */
	final boolean append;
	/** Whether to allow subprocess creation for compression */
	final boolean allowSubprocess;
	/** Default file format for created output streams */
	final int defaultFormat;
	/** Maximum buffer size for output streams */
	final int maxSize;
	/** Whether output streams should share common headers */
	final boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global verbosity flag for debug output */
	public static boolean verbose=false;

}
