package icecream;

import parse.Parse;

/**
 * Mutable parser for PacBio coordinate headers of the form
 * {@code movie/zmw/start_stop}, such as
 * {@code m64021_190821_100154/102/32038_35649}. Stores the supplied movie
 * coordinates without checking their order or agreement with a read length.
 * The full parser requires the underscore coordinate suffix; use
 * {@link #parseZMW(String)} when only the numeric well ID is needed from a
 * header with another suffix, such as {@code movie/zmw/ccs}.
 * Public fields and parsing are unsynchronized. A failed parse can leave
 * partially updated fields; do not share an instance across concurrent updates.
 * @author Brian Bushnell, Sucrose
 */
public class PBHeader{

	/** Parses a coordinate header immediately.
	 * @param s Non-null header in movie/zmw/start_stop form
	 */
	public PBHeader(String s){
		setFrom(s);
	}
	
	/**
	 * Replaces fields from the first two slashes and the following underscore.
	 * Stop-coordinate parsing ends at the first nondigit after the underscore,
	 * allowing a trailing description without treating it as part of the coordinate.
	 * The original header is stored before validation. Delimiter checks use Java
	 * assertions; malformed numeric fields are left to the integer parser.
	 * @param header Non-null coordinate header; numeric fields are decimal integers
	 * @throws AssertionError If assertions are enabled and a required delimiter is missing
	 */
	public void setFrom(final String header){
		original=header;
		int slash1=header.indexOf('/');
		assert(slash1>=0) : "Misformatted PBHeader: "+header;
		int slash2=header.indexOf('/', slash1+1);
		assert(slash2>=0) : "Misformatted PBHeader: "+header;
		int under3=header.indexOf('_', slash2+1);
		assert(under3>=0) : "Misformatted PBHeader: "+header;
		int terminator=under3+1;
		while(terminator<header.length() && Character.isDigit(header.charAt(terminator))){terminator++;}
		runID=header.substring(0, slash1);
		zmwID=Parse.parseInt(header, slash1+1, slash2);
		start=Parse.parseInt(header, slash2+1, under3);
		stop=Parse.parseInt(header, under3+1, terminator);
	}
	
	/**
	 * Extracts the numeric field between the first two slashes without constructing
	 * a PBHeader. The remaining suffix is not parsed or checked, so coordinate and
	 * CCS suffixes are both accepted here. Delimiter checks use Java assertions.
	 * @param header Non-null header containing movie/numeric-well/suffix
	 * @return Parsed well ID, without range or movie-identity validation
	 * @throws AssertionError If assertions are enabled and either slash is missing
	 */
	public static int parseZMW(final String header){
		int slash1=header.indexOf('/');
		assert(slash1>=0) : "Misformatted PBHeader: "+header;
		int slash2=header.indexOf('/', slash1+1);
		assert(slash2>=0) : "Misformatted PBHeader: "+header;
		int zmwID=Parse.parseInt(header, slash1+1, slash2);
		return zmwID;
	}
	
	/** Most recently supplied header, including any description; assigned before parsing. */
	public String original;
	/** Movie/run prefix before the first slash. */
	public String runID;
	/** Numeric well ID; the movie prefix is needed to distinguish wells across movies. */
	public int zmwID;
	/** Encoded movie start coordinate, with no ordering or read-length validation. */
	public int start;
	/** Encoded movie stop coordinate, with no ordering or read-length validation. */
	public int stop;
	
}
