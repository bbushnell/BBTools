package jgi;

import dna.AminoAcid;
import shared.Timer;
import stream.Read;
import template.BBTool_ST;

/**
 * @author Brian Bushnell
 * @date Mar 16, 2015
 *
 */
public class RemoveBadBarcodes extends BBTool_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RemoveBadBarcodes bbt=new RemoveBadBarcodes(args);
		bbt.process(t);
	}

	/**
	 * @param args
	 */
	public RemoveBadBarcodes(String[] args) {
		super(args);
	}
	
	@Override
	protected void setDefaults(){}
	
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		return false;
	}
	
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
		String id=r1.id;
		int loc=(id==null ? -1 : id.lastIndexOf(':'));
		if(loc<0 || loc>=id.length()-1){
			noBarcode++;
			return false;
		}
		for(int i=loc+1; i<id.length(); i++){
			char c=id.charAt(i);
			boolean ok=(c=='+' || AminoAcid.isFullyDefined(c));
			if(!ok){
				bad++;
				return false;
			}
		}
		good++;
		return true;
	}
	
	@Override
	protected void startupSubclass() {}
	
	protected @Override
	void shutdownSubclass() {}
	
	@Override
	protected final boolean useSharedHeader(){return true;}
	
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn) {
		
		outstream.println();
		outstream.println("Good:               "+good);
		outstream.println("Bad:                "+bad);
		outstream.println("No Barcode:         "+noBarcode);
	}
	
	long good=0;
	long bad=0;
	long noBarcode=0;
	
}
