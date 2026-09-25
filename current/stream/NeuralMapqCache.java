package stream;

/** Runtime-transient loose and strict neural MAPQs packed in one Read short. */
public final class NeuralMapqCache {
	private NeuralMapqCache(){}
	public static void set(final Read read,final int mapq){setLoose(read,mapq);}
	public static void setLoose(final Read read,final int mapq){
		validate(read,mapq,"Loose");
		read.neuralMapqs=(short)((read.neuralMapqs&~LOOSE_MASK)|(mapq+1));
	}
	public static void setStrict(final Read read,final int mapq){
		validate(read,mapq,"Strict");
		read.neuralMapqs=(short)((read.neuralMapqs&~STRICT_MASK)|((mapq+1)<<STRICT_SHIFT));
	}
	public static int get(final Read read){return getLoose(read);}
	public static int getLoose(final Read read){return decode(read,0);}
	public static int getStrict(final Read read){return decode(read,STRICT_SHIFT);}
	public static void clear(final Read read){if(read!=null){read.neuralMapqs=0;}}
	private static int decode(final Read read,final int shift){
		if(read==null){return -1;}final int encoded=(read.neuralMapqs>>>shift)&COMPONENT_MASK;return encoded-1;
	}
	private static void validate(final Read read,final int mapq,final String label){
		if(read==null||mapq<0||mapq>63){throw new IllegalArgumentException(label+" neural MAPQ outside [0,63]: "+mapq);}
	}
	private static final int COMPONENT_MASK=127,LOOSE_MASK=COMPONENT_MASK,STRICT_SHIFT=7,STRICT_MASK=COMPONENT_MASK<<STRICT_SHIFT;
}
