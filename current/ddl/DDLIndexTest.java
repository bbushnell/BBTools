package ddl;

import java.util.ArrayList;
import java.util.Random;

import cardinality.DynamicDemiLog;
import stream.Read;

/**
 * Correctness gate for the CSR sketch index: builds a matrix {@link DDLIndex} and a
 * {@link DDLIndexCSR} from the same random records and verifies that populatedCells(),
 * query() (full per-clade count vector), and topHits() are IDENTICAL for every query.
 * CSR must be a drop-in replacement, so any difference is a bug.
 *
 * Usage: DDLIndexTest [nRecords=300] [buckets=512] [k=13] [maxHits=10]
 *
 * @author Noire
 * @date August 3, 2026
 */
public class DDLIndexTest {

	public static void main(String[] args){
		int nRecords=args.length>0 ? Integer.parseInt(args[0]) : 300;
		int buckets=args.length>1 ? Integer.parseInt(args[1]) : 512;
		int k=args.length>2 ? Integer.parseInt(args[2]) : 13;
		int maxHits=args.length>3 ? Integer.parseInt(args[3]) : 10;
		int threads=4;
		DynamicDemiLog.setExponent(5);

		Random rng=new Random(12345);
		byte[] acgt={'A','C','G','T'};
		ArrayList<DDLRecord> records=new ArrayList<>(nRecords);
		for(int i=0; i<nRecords; i++){
			int len=200+rng.nextInt(5000);
			byte[] bases=new byte[len];
			for(int j=0; j<len; j++){bases[j]=acgt[rng.nextInt(4)];}
			Read r=new Read(bases, null, "rec"+i, 0);
			DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
			ddl.hash(r);
			records.add(new DDLRecord(ddl, -1, -1, "rec"+i));
		}

		DDLIndex matrix=new DDLIndex(buckets);
		matrix.addAll(records, threads);
		DDLIndexCSR csr=new DDLIndexCSR(buckets);
		csr.addAll(records, threads);

		boolean ok=true;

		long pm=matrix.populatedCells(), pc=csr.populatedCells();
		if(pm!=pc){System.err.println("MISMATCH populatedCells: matrix="+pm+" csr="+pc); ok=false;}
		else{System.err.println("populatedCells match: "+pm);}

		int queryFails=0, topFails=0;
		for(int qi=0; qi<records.size(); qi++){
			DynamicDemiLog q=records.get(qi).ddl;
			int[] cm=matrix.query(q);
			int[] cc=csr.query(q);
			boolean qEq=(cm.length==cc.length);
			for(int j=0; qEq && j<cm.length; j++){if(cm[j]!=cc[j]){qEq=false;}}
			if(!qEq){queryFails++;}

			int[][] tm=matrix.topHits(q, maxHits);
			int[][] tc=csr.topHits(q, maxHits);
			boolean tEq=(tm.length==tc.length);
			for(int j=0; tEq && j<tm.length; j++){if(tm[j][0]!=tc[j][0] || tm[j][1]!=tc[j][1]){tEq=false;}}
			if(!tEq){topFails++;}
		}

		if(queryFails>0){System.err.println("query() MISMATCH on "+queryFails+"/"+records.size()); ok=false;}
		else{System.err.println("query() identical on all "+records.size()+" queries");}
		if(topFails>0){System.err.println("topHits() MISMATCH on "+topFails+"/"+records.size()); ok=false;}
		else{System.err.println("topHits() identical on all "+records.size()+" queries");}

		System.err.println(ok ? "PASS: CSR == matrix" : "FAIL: CSR != matrix");
		if(!ok){System.exit(1);}
	}
}
