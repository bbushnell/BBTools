package test.dna;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.ChromosomeArray;
import shared.Shared;

/** Finite array-boundary oracles; each case runs in a fresh JVM.
 * Contract-edge failures are evidence for review, not authorization to change callers.
 * @author Sucrose
 */
public final class ChromosomeArrayContractTest{

	/** Runs one explicit case with checks active under both assertion modes. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Require case=offset|padding|empty|resize");
		}
		final String mode=args[0].substring(5);
		ChromosomeArray.CHANGE_U_TO_T=ChromosomeArray.CHANGE_DEGENERATE_TO_N=true;
		if(mode.equals("offset")){offset();}
		else if(mode.equals("padding")){padding();}
		else if(mode.equals("empty")){empty();}
		else if(mode.equals("resize")){resize();}
		else{throw new IllegalArgumentException(mode);}
		System.out.println("PASS "+mode);
	}

	/** Both bulk setters must retain the first and last written coordinates. */
	private static void offset(){
		final ChromosomeArray text=array(), bytes=array();
		text.set(5, "aCg");
		bytes.set(5, "aCg".getBytes(StandardCharsets.US_ASCII));
		for(ChromosomeArray ca : new ChromosomeArray[]{text, bytes}){
			System.out.println("offset min="+ca.minIndex+" max="+ca.maxIndex+
					" bases="+ca.getString(0, 9));
			require(ca.minIndex==5 && ca.maxIndex==7, "Bulk setter lost its written bounds");
			require(ca.getString(0, 9).equals("NNNNNACGNN"), "Unexpected base normalization or padding");
			require(Arrays.equals(ca.getBytes(5, 7), "ACG".getBytes(StandardCharsets.US_ASCII)),
					"In-range byte extraction differs from literal bases");
		}
	}

	/** Tests the documented outside-range N padding against getString's contract. */
	private static void padding(){
		final ChromosomeArray ca=array();
		ca.set(3, "ACG");
		final String expected="NNNNACGNN";
		require(ca.getString(-1, 7).equals(expected), "Reference per-base padding control failed");
		try{
			final String actual=new String(ca.getBytes(-1, 7), StandardCharsets.US_ASCII);
			System.out.println("padding actual="+actual+" expected="+expected);
			require(actual.equals(expected), "getBytes outside-range padding differs from getString");
		}catch(RuntimeException failure){
			System.out.println("padding threw="+failure.getClass().getName());
			throw failure;
		}
		final int[][] bounds={{-4, -2}, {0, 2}, {2, 4}, {4, 7}, {1001, 1003}, {3, 2}};
		final String[] padded={"NNN", "NNN", "NAC", "CGNN", "NNN", ""};
		for(int i=0; i<bounds.length; i++){
			require(Arrays.equals(ca.getBytes(bounds[i][0], bounds[i][1]),
					padded[i].getBytes(StandardCharsets.US_ASCII)), "Padding boundary control failed at "+i);
		}
		ca.array[0]='T'; //Storage outside logical bounds must not escape through getBytes.
		ca.array[6]='T';
		require(Arrays.equals(ca.getBytes(0, 6), "NNNACGN".getBytes(StandardCharsets.US_ASCII)),
				"Extraction exposed storage outside the logical interval");
		final ChromosomeArray gap=array();
		gap.set(3, 'A');
		gap.set(5, 'G');
		require(Arrays.equals(gap.getBytes(2, 6), new byte[]{'N', 'A', 0, 'G', 'N'}),
				"Padding changed a stored byte inside the logical interval");
	}

	/** An empty bulk write should not manufacture a newly populated interval. */
	private static void empty(){
		boolean unchanged=true;
		for(boolean bytes : new boolean[]{false, true}){
			final ChromosomeArray ca=array();
			ca.set(3, "ACG");
			final byte[] before=ca.array.clone();
			if(bytes){ca.set(20, new byte[0]);}else{ca.set(20, "");}
			System.out.println("empty bytes="+bytes+" min="+ca.minIndex+" max="+ca.maxIndex);
			unchanged&=ca.minIndex==3 && ca.maxIndex==5 && Arrays.equals(before, ca.array);
		}
		require(unchanged, "Empty bulk write changed logical bounds or stored bytes");
	}

	/** Shrinking to an inclusive maximum would discard that last valid base. */
	private static void resize(){
		final ChromosomeArray ca=array();
		ca.set(0, "ACG");
		final byte[] before=ca.array.clone();
		Throwable rejection=null;
		try{ca.resize(2);}
		catch(RuntimeException | AssertionError failure){rejection=failure;}
		System.out.println("resize rejection="+(rejection==null ? "none" : rejection.getClass().getName())+
				" capacity="+ca.array.length+" max="+ca.maxIndex);
		require(rejection!=null, "Resize accepted capacity equal to the inclusive maximum index");
		require(ca.minIndex==0 && ca.maxIndex==2 && Arrays.equals(before, ca.array),
				"Rejected resize changed the array or logical bounds");
		ca.resize(3);
		require(ca.array.length==3 && ca.getString(0, 2).equals("ACG"), "Exact-fit resize control failed");
	}

	/** Makes an empty array with the usual sentinel bounds. */
	private static ChromosomeArray array(){return new ChromosomeArray(7, Shared.PLUS);}
	private static void require(boolean ok, String message){
		if(!ok){throw new AssertionError(message);}
	}
}
