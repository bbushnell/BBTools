package barcode.stub;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import barcode.Barcode;
import barcode.PCRMatrix;
import structures.ByteBuilder;
/**
 * Abstract superclass for PCRMatrixProb variants.
 * 
 * @author Brian Bushnell
 * @date May 15, 2024
 *
 */
public abstract class PCRMatrixProbAbstract extends PCRMatrix {
	
	public PCRMatrixProbAbstract() {super(0, 0, 0, false);}
	
	public static final boolean parseStatic(String arg, String a, String b){return false;}
	
	public final static void postParseStatic(){}
	
	@Override
	public final boolean parse(String arg, String a, String b) {return false;}
	
	@Override
	public final void refine(Collection<Barcode> cb, long c) {}
	
	@Override
	public final HashMap<String, String> makeAssignmentMap(Collection<Barcode> cb, long x) {return null;}

	@Override
	public final void populateCounts(ArrayList<Barcode> list, long minCount) {}
	
	@Override
	public final void makeProbs() {}

	@Override
	public final void initializeData() {}
	
	@Override
	public final void populateUnexpected() {}

	@Override
	protected final boolean valid() {return false;}
	
	@Override
	public final Barcode findClosest(String s) {return null;}
	
	@Override
	public final ByteBuilder toBytesProb(ByteBuilder bb) {return null;}
	
	public static final boolean clientside() {return true;}
	
}
