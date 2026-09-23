package stream;

/** Runtime-transient primary neural MAPQ stored in currently unused Read flag bits. */
public final class NeuralMapqCache {
	private NeuralMapqCache(){}
	public static void set(final Read read,final int mapq){
		if(read==null||mapq<0||mapq>63){throw new IllegalArgumentException("Neural MAPQ outside [0,63]: "+mapq);}
		read.flags=(read.flags&~CACHE_MASK)|VALID_MASK|(mapq<<SHIFT);
	}
	public static int get(final Read read){return read==null||(read.flags&VALID_MASK)==0 ? -1 : (read.flags&VALUE_MASK)>>>SHIFT;}
	public static void clear(final Read read){if(read!=null){read.flags&=~CACHE_MASK;}}
	public static final int VALID_MASK=1<<22,SHIFT=23,VALUE_MASK=63<<SHIFT,CACHE_MASK=VALID_MASK|VALUE_MASK;
}
