package ddl;

/** Position tags and bucket placement are part of protein sketch package identity.
 * Relative coordinates use the sequence's residue length, not its count of word starts.
 * Multiscale retrieval is a query policy over relative sketches, not another stored layout.
 * @author Brian, Keqing
 */
public enum ProteinPositionMode{

	PLAIN("plain", "plain"),
	FIXED("fixed", "positional-v1"),
	RELATIVE("relative", "positional-relative-v1"),
	FIXED_QUADRANT("fixed-quadrant", "positional-fixed-quadrant-v1");

	/** Binds the CLI spelling and serialized sketch-layout token. */
	ProteinPositionMode(String flag_, String sketch_){
		flag=flag_;
		sketch=sketch_;
	}

	/** Parses a case-insensitive CLI spelling; null/unknown values throw rather than default. */
	public static ProteinPositionMode parse(String value){
		for(ProteinPositionMode mode : values()){
			if(mode.flag.equalsIgnoreCase(value)){return mode;}
		}
		throw new IllegalArgumentException("Unknown protein position mode: "+value);
	}

	/** Parses an exact case-sensitive wire token; null/unknown values throw. */
	public static ProteinPositionMode fromSketch(String value){
		for(ProteinPositionMode mode : values()){
			if(mode.sketch.equals(value)){return mode;}
		}
		throw new IllegalArgumentException("Unknown protein sketch position layout: "+value);
	}

	/**
	 * Resolves CLI aliases in order: explicit mode, sketch token, legacy boolean, default.
	 * An explicit mode must agree with any other supplied option. A missing boolean is
	 * distinct from false; true accepts any positional explicit mode, or selects FIXED
	 * when used alone. The sketch spelling is lowercased here before wire-token parsing.
	 * @throws IllegalArgumentException if the default is null, aliases conflict, a token
	 * is unknown, or both sketch and positional are supplied
	 */
	public static ProteinPositionMode resolve(String explicit, String sketch, Boolean positional, ProteinPositionMode defaultMode){
		if(defaultMode==null){throw new IllegalArgumentException("Position mode requires a declared default");}
		if(sketch!=null && positional!=null){throw new IllegalArgumentException("Use either sketch or positional, not both");}
		final ProteinPositionMode stored=sketch==null ? null : fromSketch(sketch.toLowerCase(java.util.Locale.ROOT));
		final ProteinPositionMode mode=explicit!=null ? parse(explicit) : stored!=null ? stored :
			positional==null ? defaultMode : positional ? FIXED : PLAIN;
		if(stored!=null && mode!=stored || positional!=null && positional.booleanValue()!=mode.positional()){
			throw new IllegalArgumentException("Contradictory protein position mode, sketch or positional options");
		}
		return mode;
	}

	/** Returns whether this mode uses position tags. */
	public boolean positional(){return this!=PLAIN;}
	/** Both fixed layouts retain the same tagged words and stop-relative boundaries. */
	public boolean fixedCoordinates(){return this==FIXED || this==FIXED_QUADRANT;}
	/** CLI spelling and serialized layout token, respectively. */
	public final String flag, sketch;
}
