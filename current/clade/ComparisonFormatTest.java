package clade;

import structures.ByteBuilder;

/**
 * Tests Comparison output formatting to verify SSU tab separator fix.
 * Checks that columns are properly tab-separated in all output modes.
 *
 * @author Brian Bushnell, Ady
 * @date April 2026
 */
public class ComparisonFormatTest {

	public static void main(String[] args){
		int pass=0, fail=0;

		//Test all 4 combos: callSSU x MAKE_DDLS
		boolean[][] combos={{false,false},{true,false},{false,true},{true,true}};
		for(boolean[] combo : combos){
			Clade.callSSU=combo[0];
			Clade.MAKE_DDLS=combo[1];
			String label="callSSU="+combo[0]+", MAKE_DDLS="+combo[1];

			//Create test comparison with known values
			Clade q=new Clade(100, 5, "QueryOrganism");
			q.setBases(1000000);
			q.contigs=3;
			q.gc=0.45f;
			q.strandedness=0.5f;
			q.hh=0.01f;
			q.caga=0.005f;
			q.gcCompEntropy=0.9f;
			q.lineage="sk__Bacteria;k__Test";

			Clade r=new Clade(200, 5, "RefOrganism");
			r.setBases(2000000);
			r.contigs=5;
			r.gc=0.47f;
			r.strandedness=0.51f;
			r.hh=0.012f;
			r.caga=0.006f;
			r.gcCompEntropy=0.91f;
			r.lineage="sk__Bacteria;k__Test;p__Testphylum";

			Comparison comp=new Comparison();
			comp.query=q;
			comp.ref=r;
			comp.gcdif=0.02f;
			comp.strdif=0.01f;
			comp.hhdif=0.002f;
			comp.cagadif=0.001f;
			comp.entdif=0.01f;
			comp.k3dif=0.05f;
			comp.k4dif=0.08f;
			comp.k5dif=0.12f;
			comp.ssudif=0.03f; //ssu identity = 0.97
			if(Clade.MAKE_DDLS){
				comp.ani=0.9876f;
				comp.wkid=0.95f;
				comp.kid=0.90f;
				comp.completeness=1.0f;
				comp.kmerMatches=1500;
			}

			//Test 1: Machine header
			ByteBuilder hdr=Comparison.machineHeader(false);
			String hdrStr=hdr.toString();
			String[] hdrCols=hdrStr.split("\t");
			System.err.println("\n=== "+label+" ===");
			System.err.println("Header ("+hdrCols.length+" cols): "+hdrStr);

			//Verify no empty columns (sign of double-tab) and no concatenated columns
			boolean headerOK=true;
			for(int i=0; i<hdrCols.length; i++){
				if(hdrCols[i].isEmpty()){
					System.err.println("  FAIL: Empty column at index "+i+" (double tab?)");
					headerOK=false;
				}
				if(hdrCols[i].contains("k5dif") && hdrCols[i].length()>5){
					System.err.println("  FAIL: k5dif concatenated with something: '"+hdrCols[i]+"'");
					headerOK=false;
				}
			}
			if(combo[0] && !hdrStr.contains("\tssuID")){
				System.err.println("  FAIL: Missing ssuID column when callSSU=true");
				headerOK=false;
			}

			//Test 2: Machine data
			ByteBuilder data=new ByteBuilder();
			comp.appendResultMachine(false, data);
			String dataStr=data.toString();
			String[] dataCols=dataStr.split("\t");
			System.err.println("Data   ("+dataCols.length+" cols): "+dataStr);

			boolean dataOK=true;
			if(hdrCols.length!=dataCols.length){
				System.err.println("  FAIL: Header has "+hdrCols.length+" cols but data has "+dataCols.length);
				dataOK=false;
			}

			//Verify last column is lineage
			String lastCol=dataCols[dataCols.length-1];
			if(!lastCol.contains("Bacteria")){
				System.err.println("  FAIL: Last column should be lineage but got: '"+lastCol+"'");
				dataOK=false;
			}

			//Test 3: Human-readable
			ByteBuilder human=new ByteBuilder();
			comp.appendResultHuman(human, 0);
			System.err.println("Human:\n"+human.toString());

			//Test 4: Round-trip parse (machine data → constructor)
			boolean parseOK=true;
			try{
				Comparison parsed=new Comparison(dataStr, null);
				if(Math.abs(parsed.k5dif-0.12f)>0.001f){
					System.err.println("  FAIL: Parsed k5dif="+parsed.k5dif+" expected 0.12");
					parseOK=false;
				}
				if(parsed.ref.lineage==null || !parsed.ref.lineage.contains("Bacteria")){
					System.err.println("  FAIL: Parsed lineage='"+parsed.ref.lineage+"' expected *Bacteria*");
					parseOK=false;
				}
				if(combo[0]){
					if(Math.abs(parsed.ssudif-0.03f)>0.001f){
						System.err.println("  FAIL: Parsed ssudif="+parsed.ssudif+" expected 0.03");
						parseOK=false;
					}
				}
			}catch(Exception e){
				System.err.println("  FAIL: Parse threw "+e.getClass().getSimpleName()+": "+e.getMessage());
				parseOK=false;
			}

			if(headerOK){pass++;}else{fail++;}
			if(dataOK){pass++;}else{fail++;}
			if(parseOK){pass++;}else{fail++;}
			System.err.println("  Header: "+(headerOK?"PASS":"FAIL")
				+", Data: "+(dataOK?"PASS":"FAIL")
				+", Parse: "+(parseOK?"PASS":"FAIL"));
		}

		System.err.println("\n"+pass+" passed, "+fail+" failed.");
		if(fail>0){System.exit(1);}
	}
}
