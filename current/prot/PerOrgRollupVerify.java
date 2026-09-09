package prot;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/** Streaming, bounded-memory verifier for PerOrgRollup's paired dense (out1=) and sparse (out2=)
 * outputs. Confirms -- without trusting either file's own row count or the producer's "Orgs:" log
 * line -- that both outputs list exactly the same TIDs in the same strictly increasing order
 * (PerOrgRollup writes both from one shared TreeMap's iteration -- true by construction, but this
 * independently re-verifies it from the actual bytes rather than assuming it); every sparse rank
 * is in [0,topn), unique per row, with a nonnegative count; every dense family column agrees
 * EXACTLY with the sparse entry for that rank, including ranks the sparse row omits (an omission
 * must read as zero in dense); and the grand total of family copies in each output equals both
 * each other and a caller-supplied conservation anchor.
 *
 * Reads both files as one forward streaming pass with small per-row buffers reused across
 * iterations (topN-sized long[] and boolean[]) -- memory is bounded by topN, never by organism
 * count or file size, so it stays cheap even against the full ~36,760-organism v5 rollup.
 *
 * Usage: java -ea prot.PerOrgRollupVerify DENSE_TSV SPARSE_TSV TOPN EXPECTED_ORGS EXPECTED_TOTAL
 * Exits nonzero with a specific message on the first violation found; prints VERIFY_PASS on success.
 * Author: Eru. */
public class PerOrgRollupVerify {
	public static void main(String[] args) throws IOException {
		if(args.length!=5){
			throw new IllegalArgumentException("DENSE_TSV SPARSE_TSV TOPN EXPECTED_ORGS EXPECTED_TOTAL required");
		}
		final String densePath=args[0], sparsePath=args[1];
		final int topN=Integer.parseInt(args[2]);
		final long expectedOrgs=Long.parseLong(args[3]);
		final long expectedTotal=Long.parseLong(args[4]);
		if(topN<=0){throw new IllegalArgumentException("TOPN must be positive: "+topN);}

		try(BufferedReader dense=Files.newBufferedReader(Paths.get(densePath), StandardCharsets.UTF_8);
			BufferedReader sparse=Files.newBufferedReader(Paths.get(sparsePath), StandardCharsets.UTF_8)){

			final String dHeader=dense.readLine();
			final String sHeader=sparse.readLine();
			if(dHeader==null || !dHeader.startsWith("#")){
				throw new IllegalArgumentException("Dense output missing header: "+densePath);
			}
			if(sHeader==null || !sHeader.startsWith("#")){
				throw new IllegalArgumentException("Sparse output missing header: "+sparsePath);
			}

			final String[] dh=dHeader.split("\t", -1);
			final int expectedDenseCols=4+topN;
			if(dh.length!=expectedDenseCols){
				throw new IllegalArgumentException("Dense header has "+dh.length+" columns, expected "
					+expectedDenseCols+" (4 fixed + topn="+topN+")");
			}
			checkEq(dh[0], "#tid"); checkEq(dh[1], "length"); checkEq(dh[2], "gc"); checkEq(dh[3], "nshreds");
			for(int i=0; i<topN; i++){checkEq(dh[4+i], "fam_"+i);}

			final String[] sh=sHeader.split("\t", -1);
			if(sh.length!=16){
				throw new IllegalArgumentException("Sparse header has "+sh.length+" columns, expected 16");
			}
			checkEq(sh[0], "#tid"); checkEq(sh[1], "domain"); checkEq(sh[15], "families");

			long rows=0, prevTid=-1, denseSum=0, sparseSum=0;
			final long[] denseRow=new long[topN];
			final boolean[] rankSeen=new boolean[topN];
			String dLine, sLine;
			while(true){
				dLine=dense.readLine();
				sLine=sparse.readLine();
				if(dLine==null && sLine==null){break;}
				if(dLine==null || sLine==null){
					throw new IllegalArgumentException("Row-count mismatch between dense and sparse outputs "
						+"at row "+(rows+1)+": one file ended before the other.");
				}
				rows++;
				final String[] d=dLine.split("\t", -1);
				final String[] s=sLine.split("\t", -1);
				if(d.length!=expectedDenseCols){
					throw new IllegalArgumentException("Dense row "+rows+" has "+d.length+" columns, expected "
						+expectedDenseCols);
				}
				if(s.length!=16){
					throw new IllegalArgumentException("Sparse row "+rows+" has "+s.length+" columns, expected 16");
				}
				final long dTid=Long.parseLong(d[0]);
				final long sTid=Long.parseLong(s[0]);
				if(dTid!=sTid){
					throw new IllegalArgumentException("TID mismatch at row "+rows+": dense="+dTid+" sparse="+sTid);
				}
				if(dTid<=0 || dTid<=prevTid){
					throw new IllegalArgumentException("TID ordering violation at row "+rows+": tid="+dTid
						+" is not strictly greater than previous tid "+prevTid+" (expected strictly "
						+"increasing, unique).");
				}
				prevTid=dTid;

				for(int i=0; i<topN; i++){
					final long v=Long.parseLong(d[4+i]);
					if(v<0){
						throw new IllegalArgumentException("Negative dense family count at row "+rows
							+" (tid="+dTid+") rank "+i+": "+v);
					}
					denseRow[i]=v;
				}

				java.util.Arrays.fill(rankSeen, false);
				long rowSparseSum=0;
				final String famField=s[15];
				if(!famField.isEmpty()){
					for(String pair : famField.split(";")){
						if(pair.isEmpty()){continue;}
						final int colon=pair.indexOf(':');
						if(colon<=0){
							throw new IllegalArgumentException("Malformed famcounts pair at row "+rows
								+" (tid="+dTid+"): \""+pair+"\"");
						}
						final int rank=Integer.parseInt(pair.substring(0, colon));
						final long count=Long.parseLong(pair.substring(colon+1));
						if(rank<0 || rank>=topN){
							throw new IllegalArgumentException("Sparse rank out of bounds at row "+rows
								+" (tid="+dTid+"): rank="+rank+" topn="+topN);
						}
						if(rankSeen[rank]){
							throw new IllegalArgumentException("Duplicate sparse rank at row "+rows
								+" (tid="+dTid+"): rank="+rank);
						}
						rankSeen[rank]=true;
						if(count<0){
							throw new IllegalArgumentException("Negative sparse count at row "+rows
								+" (tid="+dTid+") rank "+rank+": "+count);
						}
						if(count!=denseRow[rank]){
							throw new IllegalArgumentException("Dense/sparse family-count mismatch at row "+rows
								+" (tid="+dTid+") rank "+rank+": dense="+denseRow[rank]+" sparse="+count);
						}
						rowSparseSum=Math.addExact(rowSparseSum,count);
					}
				}

				long rowDenseSum=0;
				for(int i=0; i<topN; i++){
					rowDenseSum=Math.addExact(rowDenseSum,denseRow[i]);
					if(!rankSeen[i] && denseRow[i]!=0){
						throw new IllegalArgumentException("Dense rank "+i+" at row "+rows+" (tid="+dTid
							+") is "+denseRow[i]+" but sparse omits it (an omission must mean zero).");
					}
				}
				if(rowDenseSum!=rowSparseSum){
					//structurally implied by the per-rank checks above, but stated as an explicit,
					//independently-computed invariant rather than only inferred from them.
					throw new IllegalArgumentException("Row-level dense/sparse sum mismatch at row "+rows
						+" (tid="+dTid+"): dense="+rowDenseSum+" sparse="+rowSparseSum);
				}
				denseSum=Math.addExact(denseSum,rowDenseSum);
				sparseSum=Math.addExact(sparseSum,rowSparseSum);
			}

			if(rows!=expectedOrgs){
				throw new IllegalArgumentException("Row count "+rows+" != expected organism count "+expectedOrgs);
			}
			if(denseSum!=sparseSum){
				throw new IllegalArgumentException("Grand total mismatch: dense="+denseSum+" sparse="+sparseSum);
			}
			if(denseSum!=expectedTotal){
				throw new IllegalArgumentException("Family-copy grand total "+denseSum
					+" != expected conservation anchor "+expectedTotal);
			}

			System.out.println("VERIFY_PASS rows="+rows+" denseSum="+denseSum+" sparseSum="+sparseSum);
		}
	}

	private static void checkEq(String actual, String expected){
		if(!actual.equals(expected)){
			throw new IllegalArgumentException("Header field mismatch: got \""+actual+"\", expected \""+expected+"\"");
		}
	}
}
