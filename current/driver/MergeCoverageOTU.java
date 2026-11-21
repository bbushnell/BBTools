package driver;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Locale;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import jgi.CovStatsLine;
import shared.Tools;

/**
 * Merges coverage statistics for Operational Taxonomic Units (OTUs).
 * Groups rows by the text after the first space anywhere in the line and before
 * the first tab. The expected shape puts that space inside the ID; if the entire
 * line has no space, the whole ID is used. A space only in a later numeric term
 * violates that assumption. Groups are emitted in first-seen order with the OTU as ID.
 * Preserves the original header, known column order, optional omissions and
 * threshold/window labels; covered percent is recomputed from the retained counts.
 * Reads the complete table before opening the output, retaining one mutable row
 * per OTU. Aggregation uses {@link CovStatsLine#add(CovStatsLine)}; unresolved
 * medians and standard deviations serialize as NA; known singleton values remain numeric.
 * Header parsing uses global CovStatsLine state, so concurrent invocations with
 * different headers are unsafe.
 *
 * @author Brian Bushnell, Sucrose
 * @date Oct 13, 2014
 */
public class MergeCoverageOTU{
	
	/**
	 * Reads one coverage table and overwrites the output with its grouped rows.
	 * Only the first two arguments select files: either input/output positional names
	 * or in=/out= assignments in either order. This legacy parser does not preserve
	 * equals signs within assigned filenames and ignores later application arguments.
	 * The first input line must be the header; subsequent lines must be data rows.
	 * Closes successfully opened input after reading, including parsing/aggregation
	 * failures. A secondary close failure is suppressed beneath the original failure.
	 * A reported input-close error fails before output opens; normal completion starts
	 * and joins the output worker, and a reported output error fails after it exits.
	 * Failed input construction and exceptional output-worker cleanup are separate.
	 * @param args At least two arguments naming input and output
	 * @throws RuntimeException If input closure or the output worker reports an error
	 */
	public static void main(String[] args){
		String a=args[0];
		String b=args[1];
		String in=null, out=null;
		if(a.toLowerCase().startsWith("in=")){
			in=a.split("=")[1];
			out=b.split("=")[1];
		}else if(a.toLowerCase().startsWith("out=")){
			in=b.split("=")[1];
			out=a.split("=")[1];
		}else{
			in=a;
			out=b;
		}
		
		LinkedHashMap<String, CovStatsLine> map=new LinkedHashMap<String, CovStatsLine>();
		int count=0;
		ArrayList<String> headers=new ArrayList<String>();
		HeaderOrder order=null;
		final TextFile tf=new TextFile(in);
		Throwable primary=null;
		try{
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(count==0){
					assert(s.startsWith("#")) : "Expected a header line starting with #";
					//Validate before any output is opened; this object can represent only one value per meaning.
					order=new HeaderOrder(s);
					CovStatsLine.initializeHeader(s);
					headers.add(s);//FIXED (G11 2026-07-02, Brian-auth) [driver/MergeCoverageOTU#001] store the '#' header so the output loop actually writes it — headers was declared+iterated but never populated → merged output had no header line.
				}else{
					int space=s.indexOf(' ');
					String otu=s.substring(space+1, s.indexOf('\t'));
					CovStatsLine csl=new CovStatsLine(s);
					CovStatsLine old=map.get(otu);
					if(old==null){
						map.put(otu, csl);
					}else{
						old.add(csl);
					}
				}
				count++;
			}
		}catch(RuntimeException | Error failure){
			primary=failure;
			throw failure;
		}finally{
			closeInput(tf, primary);
		}
		
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		for(String s : headers){tsw.println(s);}
		for(String s : map.keySet()){
			CovStatsLine csl=map.get(s);
			csl.id=s;
			tsw.println(order.format(csl));
		}
		if(tsw.poisonAndWait()){
			throw new RuntimeException("driver.MergeCoverageOTU terminated in an error state; the output may be corrupt.");
		}
	}

	/** Closes the opened reader, preserving any earlier failure and reporting a close-only failure directly. */
	private static void closeInput(final TextFile input, final Throwable primary){
		try{
			if(input.close()){
				throw new RuntimeException("driver.MergeCoverageOTU terminated in an error state while closing the input; output was not written.");
			}
		}catch(RuntimeException | Error closeFailure){
			if(primary==null){throw closeFailure;}
			if(primary!=closeFailure){primary.addSuppressed(closeFailure);}
		}
	}

	/**
	 * Immutable output plan for the input header; labels/order and omissions stay intact.
	 * Does not change CovStatsLine's legacy aggregation semantics or its public serializer.
	 */
	private static final class HeaderOrder{
		/**
		 * Resolves each input column to exactly one supported field before output opens.
		 * Recognition ignores case and treats all Under_* labels as the same meaning.
		 * @param header First coverage line, with any leading hash characters
		 * @throws IllegalArgumentException For unknown columns or repeated meanings
		 */
		HeaderOrder(String header){
			assert(NAMES.length==FORMATS.length) : "HeaderOrder field lookup indexes FORMATS; column names and formats must align.";
			while(header.startsWith("#")){header=header.substring(1);}
			final String[] columns=header.split("\t", -1);
			fields=new int[columns.length];
			final boolean[] seen=new boolean[NAMES.length];
			final StringBuilder pattern=new StringBuilder();
			for(int i=0; i<columns.length; i++){
				final int field=field(columns[i]);
				if(seen[field]){throw new IllegalArgumentException("Duplicate coverage column meaning: "+columns[i]);}
				seen[field]=true;
				fields[i]=field;
				if(i>0){pattern.append('\t');}
				pattern.append(FORMATS[field]);
			}
			format=pattern.toString();
		}

		/**
		 * Selects values in input-column order and formats one row without a terminator.
		 * Finite fractional fields use four decimals; integral fields retain their long values.
		 * @param row Current aggregate; its ID is already replaced with the OTU name
		 * @return Tab-delimited row containing only the columns requested by the header
		 */
		String format(final CovStatsLine row){
			final Object[] values=new Object[fields.length];
			for(int i=0; i<fields.length; i++){
				switch(fields[i]){
					case ID: values[i]=row.id; break;
					case AVG_FOLD: values[i]=row.avgFold; break;
					case LENGTH: values[i]=row.length; break;
					case REF_GC: values[i]=row.refGC; break;
					case COVERED_PERCENT: values[i]=row.coveredPercent(); break;
					case COVERED_BASES: values[i]=row.coveredBases; break;
					case PLUS_READS: values[i]=row.plusReads; break;
					case MINUS_READS: values[i]=row.minusReads; break;
					case MEDIAN_FOLD: values[i]=row.medianText(); break;
					case UNDER_MIN: values[i]=row.underMin; break;
					case READ_GC: values[i]=row.readGC; break;
					case STD_DEV: values[i]=row.stdDevText(); break;
					default: throw new AssertionError("HeaderOrder constructed an unsupported field: "+fields[i]);
				}
			}
			return Tools.format(format, values);
		}

		/**
		 * Finds the semantic field index; only recognition normalizes Under_* labels.
		 * The original threshold/window spelling stays in the separately retained header.
		 * @throws IllegalArgumentException If the column cannot be represented
		 */
		private static int field(final String name){
			final String key=name.toLowerCase(Locale.ROOT);
			if(key.startsWith("under_")){return UNDER_MIN;}
			for(int i=0; i<NAMES.length; i++){
				if(key.equals(NAMES[i])){return i;}
			}
			throw new IllegalArgumentException("Unsupported coverage column: "+name);
		}

		/** Semantic field index for each input column, in the same order. */
		private final int[] fields;
		/** Complete printf pattern assembled once for this header. */
		private final String format;
		/** Shared indices into the name and format tables and the row-selection switch. */
		private static final int ID=0, AVG_FOLD=1, LENGTH=2, REF_GC=3, COVERED_PERCENT=4, COVERED_BASES=5,
				PLUS_READS=6, MINUS_READS=7, MEDIAN_FOLD=8, UNDER_MIN=9, READ_GC=10, STD_DEV=11;
		/** Lowercase recognized column meanings, including normalized Under_* labels. */
		private static final String[] NAMES={"id", "avg_fold", "length", "ref_gc", "covered_percent", "covered_bases",
				"plus_reads", "minus_reads", "median_fold", "under_min", "read_gc", "std_dev"};
		/** Formatting policy paired by index with NAMES. */
		private static final String[] FORMATS={"%s", "%.4f", "%d", "%.4f", "%.4f", "%d", "%d", "%d", "%s", "%d", "%.4f", "%s"};
	}
	
}
