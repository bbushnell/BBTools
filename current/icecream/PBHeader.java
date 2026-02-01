package icecream;

import parse.Parse;

/**
 * Parses and extracts structured metadata from Pacific Biosciences
 * sequencing run headers.
 * Handles complex header strings with format: runID/zmwID/start_stop
 * (e.g., "m64021_190821_100154/102/32038_35649").
 *
 * @author Brian Bushnell
 */
public class PBHeader {

	/** Constructs a PBHeader by parsing the provided header string.
	 * @param s PacBio header string to parse */
	public PBHeader(String s){
		setFrom(s);
	}
	
	//m64021_190821_100154/102/32038_35649
	//m64021_190821_100154 is run id
	//102 is zmw id
	//32038_35649 are movie coordinates
	//This function allows extra stuff like whitespace after the end of the formal header.
	/**
	 * Parses PacBio header string into component fields.
	 * Extracts run ID, ZMW ID, start and stop coordinates from header format:
	 * runID/zmwID/start_stop. Handles trailing whitespace after coordinates.
	 * @param header Header string with format like "m64021_190821_100154/102/32038_35649"
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
//		System.err.println("\n"+header+"\n"+slash1+","+slash2);
		zmwID=Parse.parseInt(header, slash1+1, slash2);
		start=Parse.parseInt(header, slash2+1, under3);
		stop=Parse.parseInt(header, under3+1, terminator);
	}
	
	/**
	 * Extracts ZMW ID from a PacBio header string without creating an object.
	 * @param header PacBio header string to parse
	 * @return ZMW ID as integer
	 */
	public static int parseZMW(final String header){
		int slash1=header.indexOf('/');
		assert(slash1>=0) : "Misformatted PBHeader: "+header;
		int slash2=header.indexOf('/', slash1+1);
		assert(slash2>=0) : "Misformatted PBHeader: "+header;
		int zmwID=Parse.parseInt(header, slash1+1, slash2);
		return zmwID;
	}
	
	/** Original unparsed header string */
	public String original;
	/** Run identifier extracted from header (e.g., "m64021_190821_100154") */
	public String runID;
	/** Zero-Mode Waveguide identifier extracted from header */
	public int zmwID;
	/** Movie start coordinate extracted from header */
	public int start;
	/** Movie stop coordinate extracted from header */
	public int stop;
	
}
