package ml;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

/**
 * Strict acceptance comparator for the BBNet export round-trip gate. Requirements,
 * each fatal: equal row counts; exact row IDs in exact order; exact, uniform value
 * widths in BOTH files; every value finite (NaN/Infinity rejected at parse); maximum
 * absolute difference within the given tolerance. The worst cell is reported either way.
 *
 * <p>Written for Yoimiya's review of bbtools-dev {@code 23a1245}: the earlier Python
 * comparator discarded row IDs, zip-truncated width mismatches, and let NaN evade
 * {@code max()}. Packaged copy of bbtools-dev {@code magqc_gpu_sweep/BBNetCompare.java}
 * (default package, SHA-256 {@code e4cbc429...} at bbtools-dev {@code c26965b}); behavior on
 * valid input is identical byte for byte; the only changes are the package, the usage
 * failure (throws instead of {@code System.exit}) and a reader {@code finally}.</p>
 *
 * <p>Usage: {@code java -ea ml.BBNetCompare <expected.tsv> <actual.tsv> <tolerance>}<br>
 * Files: {@code #}-prefixed lines ignored; data rows are {@code rowid<TAB>value[<TAB>value...]}.
 * Exit 0 with a PASS line, nonzero (uncaught exception) with a precise reason otherwise.</p>
 *
 * @author Ady (bbtools-dev, 2026-09), UMP45 (packaging, 2026-09-08)
 */
public class BBNetCompare {

	public static void main(String[] args) throws Exception{
		if(args.length!=3){
			throw new IllegalArgumentException("Usage: java -ea ml.BBNetCompare <expected.tsv> <actual.tsv> <tolerance>");
		}
		final double tol=Double.parseDouble(args[2]);
		if(!(tol>0) || !Double.isFinite(tol)){throw new RuntimeException("invalid tolerance "+args[2]);}
		final ArrayList<String[]> e=read(args[0]), a=read(args[1]);
		if(e.isEmpty()){throw new RuntimeException(args[0]+": zero data rows");}
		if(e.size()!=a.size()){
			throw new RuntimeException("row count mismatch: expected "+e.size()+" actual "+a.size());
		}
		final int width=e.get(0).length;
		if(width<2){throw new RuntimeException(args[0]+": no value columns (width "+width+")");}
		double worst=0; int worstRow=0, worstCol=1;
		for(int i=0; i<e.size(); i++){
			final String[] er=e.get(i), ar=a.get(i);
			if(er.length!=width){throw new RuntimeException(args[0]+" data row "+i+": width "+er.length+" != "+width);}
			if(ar.length!=width){throw new RuntimeException(args[1]+" data row "+i+": width "+ar.length+" != "+width);}
			if(!er[0].equals(ar[0])){
				throw new RuntimeException("row ID mismatch at data row "+i+": '"+er[0]+"' vs '"+ar[0]+"'");
			}
			for(int j=1; j<width; j++){
				final double x=Double.parseDouble(er[j]), y=Double.parseDouble(ar[j]);
				if(!Double.isFinite(x)){throw new RuntimeException(args[0]+" row "+er[0]+" col "+j+": non-finite value "+er[j]);}
				if(!Double.isFinite(y)){throw new RuntimeException(args[1]+" row "+ar[0]+" col "+j+": non-finite value "+ar[j]);}
				final double d=Math.abs(x-y);
				if(d>worst){worst=d; worstRow=i; worstCol=j;}
			}
		}
		// !(worst<=tol) rather than worst>tol: belt-and-braces so a NaN that somehow reached
		// here (it cannot, parse rejects) would FAIL, never pass.
		if(!(worst<=tol)){
			throw new RuntimeException("FAIL: max|diff|="+worst+" > tol="+tol
				+" at row "+e.get(worstRow)[0]+" col "+worstCol);
		}
		System.out.println("PASS rows="+e.size()+" outputs="+(width-1)
			+" max_abs_diff="+worst+" worst_row="+e.get(worstRow)[0]
			+" worst_col="+worstCol+" tol="+tol);
	}

	private static ArrayList<String[]> read(String path) throws Exception{
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		final BufferedReader br=new BufferedReader(new FileReader(path));
		try{
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				rows.add(line.split("\t", -1));
			}
		}finally{
			br.close();
		}
		return rows;
	}
}
