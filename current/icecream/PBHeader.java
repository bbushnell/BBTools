package icecream;

import shared.Parse;

public class PBHeader {

	public PBHeader(String s){
		setFrom(s);
	}
	
	//m64021_190821_100154/102/32038_35649
	//m64021_190821_100154 is run id
	//102 is zmw id
	//32038_35649 are movie coordinates
	//This function allows extra stuff like whitespace after the end of the formal header.
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
	
	public static int parseZMW(final String header){
		int slash1=header.indexOf('/');
		assert(slash1>=0) : "Misformatted PBHeader: "+header;
		int slash2=header.indexOf('/', slash1+1);
		assert(slash2>=0) : "Misformatted PBHeader: "+header;
		int zmwID=Parse.parseInt(header, slash1+1, slash2);
		return zmwID;
	}
	
	public String original;
	public String runID;
	public int zmwID;
	public int start;
	public int stop;
	
}
