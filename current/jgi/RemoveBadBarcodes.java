package jgi;

import dna.AminoAcid;
import shared.Timer;
import stream.Read;
import template.BBTool_ST;

/**
 * Filters sequence reads based on barcode quality in read identifiers.
 * Removes reads with missing, malformed, or invalid barcodes from sequencing data.
 * Barcodes are expected to be located after the last colon in the read ID.
 *
 * @author Brian Bushnell
 * @date Mar 16, 2015
 */
public class RemoveBadBarcodes extends BBTool_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		RemoveBadBarcodes bbt=new RemoveBadBarcodes(args);
		bbt.process(t);
	}

	/** Constructs a RemoveBadBarcodes instance with command line arguments.
	 * @param args Command line arguments for configuration */
	public RemoveBadBarcodes(String[] args) {
		super(args);
	}
	
	/** Sets default parameter values (no parameters for this tool) */
	@Override
	protected void setDefaults(){}
	
	/**
	 * Parses command line arguments (this tool accepts no specific arguments).
	 *
	 * @param arg Full argument string
	 * @param a Argument key
	 * @param b Argument value
	 * @return Always false as this tool has no custom arguments
	 */
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		return false;
	}
	
	/**
	 * Processes a read pair to determine if barcodes are valid.
	 * Extracts barcode from read ID (after last colon) and validates characters.
	 * Valid barcode characters are fully-defined DNA bases (A,C,G,T) or '+'.
	 *
	 * @param r1 First read in pair
	 * @param r2 Second read in pair (may be null)
	 * @return true if barcode is valid, false if invalid or missing
	 */
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
	
	/** Performs startup initialization (no action required for this tool) */
	@Override
	protected void startupSubclass() {}
	
	/** Performs cleanup during shutdown (no action required for this tool) */
	protected @Override
	void shutdownSubclass() {}
	
	/** Determines if shared header should be used for output files.
	 * @return Always true to use shared header */
	@Override
	protected final boolean useSharedHeader(){return true;}
	
	/**
	 * Displays filtering statistics after processing.
	 * Shows counts of good reads, bad barcodes, and reads without barcodes.
	 *
	 * @param t Timer used for execution timing
	 * @param readsIn Total number of input reads processed
	 * @param basesIn Total number of input bases processed
	 */
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
