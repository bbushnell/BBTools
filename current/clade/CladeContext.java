package clade;

import java.util.Arrays;

import shared.LineParser1;
import shared.Parse;
import shared.Tools;

/**
 * Context object for connection-specific parameters in CladeServer.
 * Ensures thread-safe handling of concurrent requests by giving
 * each connection its own isolated parameter set.
 *
 * @author Chloe
 * @date September 16, 2025
 */
public class CladeContext {

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor with default values */
	public CladeContext() {
		// Default values
	}

	/** Copy constructor for cloning contexts */
	public CladeContext(CladeContext other) {
		this.format = other.format;
		this.hits = other.hits;
		this.heap = other.heap;
		this.printQTID = other.printQTID;
		this.banSelf = other.banSelf;
	}
	


	/** Parse parameters from request body (SendClade format) */
	void parseRequestBody(byte[] requestBody){
		if(requestBody==null || requestBody.length<1){return;}

		// SendClade sends parameters in first line: format=oneline/hits=10/heap=1/
		int idx=Tools.indexOf(requestBody, '\n');
		byte[] line=(idx<0 ? requestBody : Arrays.copyOf(requestBody, idx+1));

		LineParser1 lp=new LineParser1('/');
		lp.set(line);
		
		for(int i=0; i<lp.terms(); i++){
			//TODO: lp.termContains('=') would be nice
			String part=lp.parseString(i);
			if(part.isEmpty()){continue;}
			String[] kv=part.split("=");
			if(kv.length!=2){continue;}

			String key=kv[0].toLowerCase();
			String value=kv[1];
			
			if(key.equals("format")){
				if(value.equals("oneline") || value.equals("machine")){
					format=CladeSearcher.MACHINE;
				}else{
					format=CladeSearcher.HUMAN;
				}
			}else if(key.equals("hits")){
				hits=Integer.parseInt(value);
			}else if(key.equals("heap")){
				heap=Integer.parseInt(value);
			}else if(key.equals("printqtid")){
				printQTID=Parse.parseBoolean(value);
			}else if(key.equals("banself")){
				banSelf=Parse.parseBoolean(value);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output format (HUMAN or ONELINE) */
	public int format = CladeSearcher.HUMAN;

	/** Number of hits to return per query */
	public int hits = 1;

	/** Heap size for comparisons */
	public int heap = 1;

	/** Print query TaxID */
	public boolean printQTID = false;

	/** Ban self-matches */
	public boolean banSelf = false;
}