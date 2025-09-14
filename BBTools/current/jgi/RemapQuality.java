package jgi;

import shared.Timer;
import stream.Read;
import stream.SamLine;
import template.BBTool_ST;

/**
 * Changes quality scores to other quality scores.
 * @author Brian Bushnell
 * @date Apr 27, 2015
 *
 */
public class RemapQuality extends BBTool_ST {
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Example:
		Timer t=new Timer();
		RemapQuality bbt=new RemapQuality(args);
		bbt.process(t);
	}
	
	@Override
	protected void setDefaults(){}

	/**
	 * @param args
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
