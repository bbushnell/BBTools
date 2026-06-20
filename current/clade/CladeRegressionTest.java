package clade;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

/**
 * Regression test for QuickClade machine-format output.
 * Reads oneline output, checks self-matches, value ranges, and optionally
 * diffs against a baseline file.
 *
 * Usage:
 *   java clade.CladeRegressionTest results.txt              — semantic checks only
 *   java clade.CladeRegressionTest results.txt baseline.txt  — semantic + exact diff
 *
 * Exit code 0 = all passed, 1 = failures found.
 *
 * @author Chloe
 */
public class CladeRegressionTest {

	public static void main(String[] args) throws IOException {
		if(args.length<1){
			System.out.println("Usage: java clade.CladeRegressionTest results.txt [baseline.txt]");
			System.exit(1);
		}

		String resultsFile=args[0];
		String baselineFile=(args.length>1 ? args[1] : null);

		int pass=0, fail=0, warn=0;

		// Read all lines from results file
		ArrayList<String> resultLines=readLines(resultsFile);
		if(resultLines.isEmpty()){
			System.out.println("FAIL: Results file is empty");
			System.exit(1);
		}

		// Parse header
		String header=resultLines.get(0);
		if(!header.startsWith("#QueryName")){
			System.out.println("FAIL: First line is not a header (expected #QueryName...)");
			System.exit(1);
		}
		System.out.println("PASS: Header present"); pass++;

		//CLEVER [verified in-file]: parses columns BY NAME via findCol (not by hardcoded index), so the test is robust
		//to column reordering/insertion in the machine format; the parseInt/parseFloat helpers tolerate col==-1
		//(missing column) by returning -1, degrading gracefully rather than crashing. This is the pattern
		//CladeCalibrator#001 should adopt -- it hardcodes indices and would silently miscalibrate on a format shift.
		String[] colNames=header.split("\t");
		int colQueryName=findCol(colNames, "#QueryName");
		int colRTaxID=findCol(colNames, "R_TaxID");
		int colK5dif=findCol(colNames, "k5dif");
		int colK4dif=findCol(colNames, "k4dif");
		int colSSU=findCol(colNames, "ssuID");
		int colWKID=findCol(colNames, "WKID");
		int colMatches=findCol(colNames, "KmerMatches");
		int colANI=findCol(colNames, "ANI");

		System.out.println("  Columns: R_TaxID="+colRTaxID+" k5dif="+colK5dif
			+" WKID="+colWKID+" KmerMatches="+colMatches+" ssuID="+colSSU);

		// Group result lines by query name
		LinkedHashMap<String, ArrayList<String[]>> byQuery=new LinkedHashMap<>();
		for(int i=1; i<resultLines.size(); i++){
			String line=resultLines.get(i);
			if(line.isEmpty()) continue;
			String[] fields=line.split("\t");
			if(fields.length<colRTaxID+1) continue;
			String qname=fields[colQueryName];
			byQuery.computeIfAbsent(qname, k->new ArrayList<>()).add(fields);
		}

		System.out.println("PASS: "+byQuery.size()+" distinct queries, "
			+(resultLines.size()-1)+" result lines"); pass++;

		// Check each query
		int checked=0, selfFound=0;
		for(Entry<String,ArrayList<String[]>> entry : byQuery.entrySet()){
			String qname=entry.getKey();
			ArrayList<String[]> hits=entry.getValue();

			// Extract expected TID from filename like "tid_227321_Aspergillus_nidulans"
			int expectedTID=extractTID(qname);
			if(expectedTID<=0) continue;
			checked++;

			String shortName=shortName(qname);

			// Find self-match
			String[] selfMatch=null;
			for(String[] fields : hits){
				int rtid=parseInt(fields, colRTaxID);
				if(rtid==expectedTID){
					selfMatch=fields;
					break;
				}
			}

			if(selfMatch==null){
				System.out.println("WARN: "+shortName+" — self-match TID "+expectedTID+" not in top results"); warn++;
				continue;
			}
			selfFound++;

			// Check k5dif for self-match
			// Threshold 0.3: eukaryotes vary up to ~0.26 due to assembly differences
			if(colK5dif>=0){
				float k5=parseFloat(selfMatch, colK5dif);
				if(k5>=0 && k5<0.3f){
					System.out.println("PASS: "+shortName+" self k5dif="+k5); pass++;
				}else{
					System.out.println("FAIL: "+shortName+" self k5dif="+k5+" (expected < 0.3)"); fail++;
				}
			}

			// Check KmerMatches for self-match (> 0 means sketch attachment worked)
			if(colMatches>=0){
				int matches=parseInt(selfMatch, colMatches);
				if(matches>0){
					System.out.println("PASS: "+shortName+" self matches="+matches); pass++;
				}else{
					System.out.println("FAIL: "+shortName+" self matches="+matches+" (expected > 0)"); fail++;
				}
			}

			// Check WKID for self-match
			if(colWKID>=0){
				float wkid=parseFloat(selfMatch, colWKID);
				if(wkid>0.5f){
					System.out.println("PASS: "+shortName+" self WKID="+wkid); pass++;
				}else if(wkid>=0){
					System.out.println("WARN: "+shortName+" self WKID="+wkid+" (expected > 0.5)"); warn++;
				}
			}

			// Check SSU ANI for self-match (should be high when present)
			if(colSSU>=0){
				float ssu=parseFloat(selfMatch, colSSU);
				if(ssu>0.8f){
					System.out.println("PASS: "+shortName+" self ssuID="+ssu); pass++;
				}else if(ssu>0f){
					System.out.println("WARN: "+shortName+" self ssuID="+ssu+" (expected > 0.8)"); warn++;
				}
				// ssu==0 means no SSU data for this pair, not a failure
			}
		}

		if(checked>0){
			System.out.println("PASS: Checked "+checked+" queries, "+selfFound+"/"+checked+" had self-matches"); pass++;
		}
		if(checked>0 && selfFound==0){
			System.out.println("FAIL: No self-matches found for any query"); fail++;
		}

		// Global SSU check: count how many results have SSU data at all
		if(colSSU>=0){
			int ssuPresent=0, ssuTotal=0;
			for(Entry<String,ArrayList<String[]>> entry : byQuery.entrySet()){
				for(String[] fields : entry.getValue()){
					ssuTotal++;
					float ssu=parseFloat(fields, colSSU);
					if(ssu>0f) ssuPresent++;
				}
			}
			if(ssuPresent>0){
				System.out.println("PASS: SSU data present in "+ssuPresent+"/"+ssuTotal+" result lines"); pass++;
			}else{
				System.out.println("FAIL: No SSU data in any result line (SSU alignment broken?)"); fail++;
			}
		}

		// Baseline diff
		if(baselineFile!=null){
			System.out.println();
			System.out.println("=== Baseline Diff ===");
			ArrayList<String> baselineLines=readLines(baselineFile);
			if(resultLines.size()!=baselineLines.size()){
				System.out.println("FAIL: Line count differs — baseline="+baselineLines.size()
					+" current="+resultLines.size()); fail++;
			}else{
				int diffs=0;
				for(int i=0; i<resultLines.size(); i++){
					if(!resultLines.get(i).equals(baselineLines.get(i))){
						if(diffs<10){
							System.out.println("  Line "+(i+1)+" differs:");
							System.out.println("    baseline: "+truncate(baselineLines.get(i), 120));
							System.out.println("    current:  "+truncate(resultLines.get(i), 120));
						}
						diffs++;
					}
				}
				if(diffs==0){
					System.out.println("PASS: Output matches baseline exactly"); pass++;
				}else{
					System.out.println("FAIL: "+diffs+" lines differ from baseline"); fail++;
				}
			}
		}

		// Summary
		System.out.println();
		System.out.println("==============================");
		System.out.println("  PASS: "+pass+"  FAIL: "+fail+"  WARN: "+warn);
		System.out.println("==============================");
		System.out.println(fail>0 ? "RESULT: FAILED" : "RESULT: PASSED");
		System.exit(fail>0 ? 1 : 0);
	}

	static ArrayList<String> readLines(String path) throws IOException {
		ArrayList<String> lines=new ArrayList<>();
		try(BufferedReader br=new BufferedReader(new FileReader(path))){
			String line;
			while((line=br.readLine())!=null){
				lines.add(line);
			}
		}
		return lines;
	}

	static int findCol(String[] header, String name){
		for(int i=0; i<header.length; i++){
			if(header[i].equals(name)) return i;
		}
		return -1;
	}

	static int extractTID(String filename){
		int idx=filename.indexOf("tid_");
		if(idx<0) return -1;
		int start=idx+4;
		int end=start;
		while(end<filename.length() && Character.isDigit(filename.charAt(end))) end++;
		if(end==start) return -1;
		try{ return Integer.parseInt(filename.substring(start, end)); }
		catch(NumberFormatException e){ return -1; }
	}

	static String shortName(String filename){
		int idx=filename.indexOf("tid_");
		if(idx<0) return filename;
		String after=filename.substring(idx);
		int underscore2=after.indexOf('_', 4);
		if(underscore2<0) return after;
		String name=after.substring(underscore2+1);
		int dot=name.indexOf('.');
		if(dot>0) name=name.substring(0, dot);
		return name.replace('_', ' ');
	}

	static int parseInt(String[] fields, int col){
		if(col<0 || col>=fields.length) return -1;
		try{ return Integer.parseInt(fields[col].trim()); }
		catch(NumberFormatException e){ return -1; }
	}

	static float parseFloat(String[] fields, int col){
		if(col<0 || col>=fields.length) return -1;
		try{ return Float.parseFloat(fields[col].trim()); }
		catch(NumberFormatException e){ return -1; }
	}

	static String truncate(String s, int max){
		return s.length()<=max ? s : s.substring(0, max)+"...";
	}
}
