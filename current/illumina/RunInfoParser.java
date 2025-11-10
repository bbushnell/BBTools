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
public class RunInfoParser {

	/**
	 * Parse RunInfo.xml and return read lengths.
	 * @param runFolder Path to run folder
	 * @return Array of read lengths in order [R1, I1, I2, R2, ...]
	 */
	public static int[] parseReadLengths(String runFolder) throws IOException {
		String xmlPath=runFolder+"/RunInfo.xml";
		BufferedReader br=new BufferedReader(new FileReader(xmlPath));

		ArrayList<Integer> lengths=new ArrayList<>();
		String line;

		while((line=br.readLine())!=null){
			//Look for <Read Number="X" NumCycles="Y" IsIndexedRead="Z"/>
			if(line.contains("<Read ") && line.contains("NumCycles=")){
				int numCycles=extractNumCycles(line);
				if(numCycles>0){
					lengths.add(numCycles);
				}
			}
		}

		br.close();

		//Convert to array
		int[] result=new int[lengths.size()];
		for(int i=0; i<lengths.size(); i++){
			result[i]=lengths.get(i);
		}

		return result;
	}

	/**
	 * Extract NumCycles value from Read XML element.
	 */
	private static int extractNumCycles(String line) {
		int idx=line.indexOf("NumCycles=\"");
		if(idx<0){return -1;}

		idx+=11; //Skip to value
		int end=line.indexOf('"', idx);
		if(end<0){return -1;}

		String value=line.substring(idx, end);
		return Integer.parseInt(value);
	}

	/**
	 * Test RunInfo.xml parsing.
	 */
	public static void main(String[] args) throws Exception {
		if(args.length<1){
			System.err.println("Usage: RunInfoParser <run_folder>");
			System.exit(1);
		}

		int[] lengths=parseReadLengths(args[0]);
		System.out.println("Read structure:");
		for(int i=0; i<lengths.length; i++){
			System.out.println("  Read " + (i+1) + ": " + lengths[i] + " cycles");
		}

		int total=0;
		for(int len : lengths){total+=len;}
		System.out.println("Total cycles: " + total);
	}
}
