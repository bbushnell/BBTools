package cardinality;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

/**
 * Averages CF table columns over the last octave of cardinality.
 * <p>
 * Reads a v5 CF TSV (produced by {@code ddlcalibrate.sh cf=f out3=...}),
 * finds the maximum cardinality C in column 0, and for each remaining
 * column reports the arithmetic mean of rows where {@code card >= C/2}.
 * Used to extract terminal-bias ratios for {@code terminalMeanCF()} /
 * {@code terminalMeanPlusCF()} overrides in each estimator class.
 * <p>
 * Output: one line per column, {@code shortname\tavg}. Short name is the
 * column header with a trailing {@code _cf} stripped and lowercased.
 * <p>
 * Usage: {@code maketerminalcf.sh in=cffile.tsv.gz}
 *
 * @author Brian Bushnell, Chloe
 * @date April 13, 2026
 */
public class TerminalCfAverager {

	public static void main(String[] args) throws Exception {
		String in=null;
		for(String a : args){
			final int eq=a.indexOf('=');
			if(eq<0) continue;
			final String k=a.substring(0, eq).toLowerCase();
			final String v=a.substring(eq+1);
			if(k.equals("in")){in=v;}
		}
		if(in==null){
			System.err.println("Usage: maketerminalcf.sh in=cffile.tsv[.gz]");
			System.exit(1);
		}
		run(in);
	}

	static void run(final String path) throws Exception {
		final InputStream raw=new FileInputStream(path);
		final InputStream is=path.endsWith(".gz") ? new GZIPInputStream(raw) : raw;
		try(BufferedReader br=new BufferedReader(new InputStreamReader(is))){
			String line;
			String[] header=null;
			final ArrayList<double[]> rows=new ArrayList<>();
			while((line=br.readLine())!=null){
				if(line.isEmpty()) continue;
				if(line.startsWith("#")){
					// Header line begins with #TrueCard (or similar) followed by column names.
					// Last #-line before data is the column header.
					if(line.startsWith("#TrueCard") || line.startsWith("#Card")
							|| (line.contains("\t") && line.toLowerCase().contains("mean"))){
						final String[] parts=line.substring(1).split("\t");
						header=parts;
					}
					continue;
				}
				final String[] parts=line.split("\t");
				final double[] vals=new double[parts.length];
				for(int i=0; i<parts.length; i++){
					try{vals[i]=Double.parseDouble(parts[i]);}
					catch(NumberFormatException e){vals[i]=Double.NaN;}
				}
				rows.add(vals);
			}
			if(header==null){
				System.err.println("No header row found (expected line starting with '#TrueCard').");
				System.exit(2);
			}
			if(rows.isEmpty()){
				System.err.println("No data rows found.");
				System.exit(2);
			}

			// Find max cardinality (column 0).
			double maxCard=0;
			for(double[] r : rows){if(r[0]>maxCard){maxCard=r[0];}}
			final double threshold=maxCard/2.0;

			// Accumulate per-column sums for rows where card >= threshold.
			final int ncols=header.length;
			final double[] sums=new double[ncols];
			final int[] counts=new int[ncols];
			for(double[] r : rows){
				if(r[0]<threshold) continue;
				final int k=Math.min(r.length, ncols);
				for(int i=1; i<k; i++){
					if(!Double.isNaN(r[i])){
						sums[i]+=r[i];
						counts[i]++;
					}
				}
			}

			// Report: short name \t average for each column except the cardinality column.
			System.err.println("# MaxCard="+(long)maxCard+"  threshold="+(long)threshold
					+"  rowsInOctave="+(counts.length>1 ? counts[1] : 0));
			for(int i=1; i<ncols; i++){
				final String shortName=shorten(header[i]);
				final double avg=counts[i]>0 ? sums[i]/counts[i] : Double.NaN;
				System.out.println(shortName+"\t"+avg);
			}
		}
	}

	/** Strip trailing "_cf" and lowercase; map "meanh" to "mean+h" for clarity. */
	static String shorten(final String col){
		String s=col.trim();
		if(s.toLowerCase().endsWith("_cf")){s=s.substring(0, s.length()-3);}
		s=s.toLowerCase();
		if(s.equals("meanh")){s="mean+h";}
		return s;
	}

}
