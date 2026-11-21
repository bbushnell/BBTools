package jgi;

import shared.Timer;
import stream.Read;
import stream.SamLine;
import template.BBTool_ST;

/**
 * Remaps numeric Phred quality arrays while retaining every read.
 *
 * <p>The default map leaves all scores unchanged except 2..41, which become
 * {@code 43-score}.  A nonempty {@code map=} value skips that reversal and overlays
 * semicolon-separated {@code from,to} assignments on an otherwise identity
 * table, leaving unlisted scores unchanged. Destinations must fit a signed byte
 * ({@code -128..127}); representable values are retained without further
 * quality-range clamping.</p>
 *
 * @author Brian Bushnell
 * @date Apr 27, 2015
 */
public class RemapQuality extends BBTool_ST{
	
	/**
	 * Runs quality remapping on the configured input and output streams.
	 *
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RemapQuality bbt=new RemapQuality(args);
		bbt.process(t);
	}
	
	/** Leaves inherited defaults unchanged before argument parsing. */
	@Override
	protected void setDefaults(){}

	/**
	 * Parses arguments, enables SamLine reuse without a loaded reference, and builds
	 * the numeric Phred remapping table.
	 *
	 * <p>The default table is identity except for scores 2 through 41; score
	 * {@code i} in that inclusive range maps to {@code 43-i}.  With a nonempty mapping value,
	 * that reversal is skipped and the supplied semicolon-separated pairs overlay
	 * direct {@code from,to} assignments on the identity table. Destinations outside
	 * {@code -128..127} are rejected before read processing begins.</p>
	 *
	 * <p>{@link SamLine#SET_FROM_OK} is set globally for this JVM so SAM records
	 * can be reused from retained SamLine data when no genome is loaded.  This
	 * tool does not restore the previous global value.</p>
	 *
	 * @param args Command-line arguments, including optional {@code map=from,to;...}
	 * @throws NumberFormatException If a parsed token is not a decimal integer,
	 * or a destination cannot be represented as a signed byte
	 */
	public RemapQuality(String[] args){
		super(args);
		SamLine.SET_FROM_OK=true;
		map=new byte[256];
		for(int i=0; i<map.length; i++){
			map[i]=(byte)i;
		}

		if(mapString==null){//reverse quality
			for(int i=2; i<=41; i++){
				map[i]=(byte)(43-i);
			}
		}else{
			String[] pairs=mapString.split(";");
			for(String pair : pairs){
				String[] split=pair.split(",");
				int a=Integer.parseInt(split[0]);
				byte b=Byte.parseByte(split[1]);
				map[a]=b;
			}
		}
	}

	/**
	 * Parses RemapQuality-specific arguments.
	 *
	 * @param arg Full argument string
	 * @param a Argument key, lower-case left side of {@code =}
	 * @param b Argument value, or null when {@code =} is absent or terminal
	 * @return true if this method recognized {@code arg}
	 */
	@Override
	public boolean parseArgument(String arg, String a, String b){
		if(a.equals("map")){
			mapString=b;
			return true;
		}else if(false){
			return true;
		}
		return false;
	}

	/** Performs no additional startup before the common streams are opened. */
	@Override
	protected void startupSubclass(){}

	/** Performs no additional shutdown before the common streams are closed. */
	@Override
	protected void shutdownSubclass(){}

	/**
	 * Adds no subclass-specific statistics beyond the common read/base counts.
	 *
	 * @param t Completed processing timer
	 * @param readsIn Total reads processed
	 * @param basesIn Total bases processed
	 */
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn){}
	
	/** Requests a shared input SAM header for output stream setup. */
	@Override
	protected final boolean useSharedHeader(){return true;}

	/**
	 * Mutates each borrowed numeric quality array in place by direct table lookup.
	 * Each present quality must be nonnegative to index the remapping table.
	 *
	 * @param r1 First read in the pair, or null
	 * @param r2 Second read in the pair, or null
	 * @return Always true; this tool never filters reads
	 */
	@Override
	protected boolean processReadPair(Read r1, Read r2){
		if(r1!=null && r1.quality!=null){
			final byte[] qual=r1.quality;
			for(int i=0; i<qual.length; i++){qual[i]=map[qual[i]];}
		}
		if(r2!=null && r2.quality!=null){
			final byte[] qual=r2.quality;
			for(int i=0; i<qual.length; i++){qual[i]=map[qual[i]];}
		}
		return true;
	}
	
	/** Optional semicolon-separated list of {@code from,to} map assignments. */
	public String mapString;
	/** Numeric Phred score remapping table built after argument parsing. */
	public final byte[] map;

}
