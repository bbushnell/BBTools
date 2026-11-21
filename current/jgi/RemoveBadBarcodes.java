package jgi;

import dna.AminoAcid;
import shared.Timer;
import stream.Read;
import template.BBTool_ST;

/**
 * Filters reads by the final colon-delimited barcode in R1 read names.
 * Missing or empty suffixes are rejected.  A/C/G/T/U in either case and '+'
 * are retained; there is no separator-count or barcode-length validation.
 * For paired input, only r1.id is inspected and the whole fragment is kept
 * or discarded from that R1 suffix.
 *
 * @author Brian Bushnell
 * @date Mar 16, 2015
 */
public class RemoveBadBarcodes extends BBTool_ST{
	
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
	public RemoveBadBarcodes(String[] args){
		super(args);
	}
	
	/** Sets default parameter values. */
	@Override
	protected void setDefaults(){}
	
	/**
	 * Parses tool-specific command line arguments.
	 * This tool has no custom arguments beyond BBTool_ST's common parser.
	 *
	 * @param arg Full argument string
	 * @param a Argument key
	 * @param b Argument value
	 * @return always false
	 */
	@Override
	public boolean parseArgument(String arg, String a, String b){
		return false;
	}
	
	/**
	 * Tests the final colon-delimited suffix of r1.id.
	 * The suffix must be nonempty and all characters must be ASCII table
	 * entries: '+' or a fully defined DNA base.  AminoAcid.isFullyDefined
	 * accepts A/C/G/T/U in either case.  r2.id is not checked.
	 *
	 * @param r1 first read in a fragment
	 * @param r2 second read in a fragment, or null
	 * @return true if the fragment should be retained
	 */
	@Override
	protected boolean processReadPair(Read r1, Read r2){
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
	
	/** Performs startup initialization. */
	@Override
	protected void startupSubclass(){}
	
	/** Performs cleanup during shutdown. */
	@Override
	protected void shutdownSubclass(){}
	
	/** Determines if shared headers should be used for output files.
	 * @return true */
	@Override
	protected final boolean useSharedHeader(){return true;}
	
	/**
	 * Displays filtering statistics after processing.
	 * Good, Bad, and No Barcode count fragments: single reads, or pairs once.
	 *
	 * @param t Timer used for execution timing
	 * @param readsIn Total number of input reads processed
	 * @param basesIn Total number of input bases processed
	 */
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn){
		
		outstream.println();
		outstream.println("Good:               "+good);
		outstream.println("Bad:                "+bad);
		outstream.println("No Barcode:         "+noBarcode);
	}
	
	/** Fragments whose R1 barcode suffix passed. */
	long good=0;
	/** Fragments with a nonempty R1 suffix containing an invalid character. */
	long bad=0;
	/** Fragments with no final colon-delimited R1 suffix. */
	long noBarcode=0;
	
}
