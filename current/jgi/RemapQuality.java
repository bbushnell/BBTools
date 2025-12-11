package jgi;

import shared.Timer;
import stream.Read;
import stream.SamLine;
import template.BBTool_ST;

/**
 * Remaps quality scores from one encoding to another using configurable mapping.
 * Supports reversing quality scores or custom mappings via semicolon-separated pairs.
 * Default behavior reverses quality scores in the range 2-41 to create inverted quality.
 *
 * @author Brian Bushnell
 * @date Apr 27, 2015
 */
public class RemapQuality extends BBTool_ST {
	
	/** Program entry point for quality score remapping.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//Example:
		Timer t=new Timer();
		RemapQuality bbt=new RemapQuality(args);
		bbt.process(t);
	}
	
	@Override
	protected void setDefaults(){}

	/**
	 * Constructs RemapQuality with command-line arguments and initializes mapping table.
	 * Creates identity mapping by default, then applies reverse mapping (43-i for i=2-41)
	 * unless custom mapping string is provided.
	 * @param args Command-line arguments including optional map parameter
	 */
	public RemapQuality(String[] args) {
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
				int b=Integer.parseInt(split[1]);
				map[a]=(byte)b;
			}
		}
	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#parseArgument(java.lang.String, java.lang.String, java.lang.String)
	 */
	/**
	 * Parses tool-specific command-line arguments.
	 * Recognizes 'map' parameter for custom quality score mapping.
	 *
	 * @param arg Full argument string
	 * @param a Argument key (left side of =)
	 * @param b Argument value (right side of =)
	 * @return true if argument was recognized and parsed
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

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#startupSubclass()
	 */
	@Override
	protected void startupSubclass() {
		// TODO Auto-generated method stub

	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#shutdownSubclass()
	 */
	@Override
	protected void shutdownSubclass() {
		// TODO Auto-generated method stub

	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#showStatsSubclass(dna.Timer, long, long)
	 */
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn) {
		// TODO Auto-generated method stub

	}
	
	@Override
	protected final boolean useSharedHeader(){return true;}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#processReadPair(stream.Read, stream.Read)
	 */
	/**
	 * Processes a read pair by remapping quality scores using the configured mapping table.
	 * Applies the byte mapping to every quality score in both reads if quality arrays exist.
	 *
	 * @param r1 First read in the pair (may be null)
	 * @param r2 Second read in the pair (may be null)
	 * @return Always returns true to retain all reads
	 */
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
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
	
	public String mapString;
	public final byte[] map;

}
