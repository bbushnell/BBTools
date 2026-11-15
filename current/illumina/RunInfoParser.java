package illumina;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Parses RunInfo.xml to extract read structure.
 * Minimal XML parsing - just finds Read elements.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class RunInfoParser{

	/**
	 * Parse RunInfo.xml and return read lengths.
	 * @param runFolder Path to run folder
	 * @return Array of read lengths in order [R1, I1, I2, R2, ...]
	 */
	public static int[] parseReadLengths(String runFolder) throws IOException{
		return parseReadStructure(runFolder).lengths;
	}

	/**
	 * Parse RunInfo.xml and return lengths plus indexed flags for each
	 * positive-length read term.
	 * Non-positive NumCycles terms are skipped to match
	 * {@link #parseReadLengths(String)}.
	 * @param runFolder path to a run folder containing RunInfo.xml
	 * @return a structure owning mutable arrays; lengths are in sequencing cycles
	 * and indexed flags are true only for {@code IsIndexedRead="Y"}
	 * @throws IOException if RunInfo.xml cannot be read
	 * @throws NumberFormatException if a present NumCycles attribute is not an
	 * integer
	 */
	public static ReadStructure parseReadStructure(String runFolder) throws IOException{
		String xmlPath=runFolder+"/RunInfo.xml";
		BufferedReader br=new BufferedReader(new FileReader(xmlPath));

		ArrayList<Integer> lengths=new ArrayList<>();
		ArrayList<Boolean> indexed=new ArrayList<>();
		String line;

		while((line=br.readLine())!=null){
			//Look for <Read Number="X" NumCycles="Y" IsIndexedRead="Z"/>
			if(line.contains("<Read ") && line.contains("NumCycles=")){
				int numCycles=extractNumCycles(line);
				if(numCycles>0){
					lengths.add(numCycles);
					indexed.add(extractIndexed(line));
				}
			}
		}

		br.close();

		//Convert to array
		int[] result=new int[lengths.size()];
		boolean[] indexedResult=new boolean[indexed.size()];
		for(int i=0; i<lengths.size(); i++){
			result[i]=lengths.get(i);
			indexedResult[i]=indexed.get(i);
		}

		return new ReadStructure(result, indexedResult);
	}

	/**
	 * Extract NumCycles value from Read XML element.
	 */
	private static int extractNumCycles(String line){
		int idx=line.indexOf("NumCycles=\"");
		if(idx<0){return -1;}

		idx+=11; //Skip to value
		int end=line.indexOf('"', idx);
		if(end<0){return -1;}

		String value=line.substring(idx, end);
		return Integer.parseInt(value);
	}

	private static boolean extractIndexed(String line){
		int idx=line.indexOf("IsIndexedRead=\"");
		if(idx<0){return false;}
		idx+=15; //Skip to value
		int end=line.indexOf('"', idx);
		if(end<=idx){return false;}
		return line.substring(idx, end).equalsIgnoreCase("Y");
	}

	/**
	 * Read lengths and indexed-read flags from one RunInfo.xml.
	 * The arrays are owned by this object but intentionally exposed for
	 * lightweight package callers; callers should treat them as read-only.
	 */
	public static class ReadStructure{
		private ReadStructure(int[] lengths_, boolean[] indexed_){
			lengths=lengths_;
			indexed=indexed_;
		}

		/** Per-read lengths in sequencing cycles, in RunInfo order. */
		public final int[] lengths;

		/** Per-read indexed flags in RunInfo order; true only for IsIndexedRead=Y. */
		public final boolean[] indexed;
	}

	/**
	 * Test RunInfo.xml parsing.
	 */
	public static void main(String[] args) throws Exception{
		if(args.length<1){
			System.err.println("Usage: RunInfoParser <run_folder>");
			System.exit(1);
		}

		ReadStructure structure=parseReadStructure(args[0]);
		int[] lengths=structure.lengths;
		System.out.println("Read structure:");
		for(int i=0; i<lengths.length; i++){
			System.out.println("  Read " + (i+1) + ": " + lengths[i] +
					" cycles, indexed="+structure.indexed[i]);
		}

		int total=0;
		for(int len : lengths){total+=len;}
		System.out.println("Total cycles: " + total);
	}
}
