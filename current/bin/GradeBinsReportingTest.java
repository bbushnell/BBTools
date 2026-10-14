package bin;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.ArrayList;

import map.IntLongHashMap;
import structures.LongList;

/** Verifies unknown contig truth and stable summary positions for shared callers.
 * @author Brian Bushnell, Nilou
 */
public final class GradeBinsReportingTest {

	public static void main(String[] args) throws Exception{
		require(args.length==1 && args[0].startsWith("dir="), "Expected dir=prepared_fixture");
		final String dir=args[0].substring(4);
		//GradeBins.process enables this before loading contigs without coverage (Cluster.add).
		BinObject.grading=true;
		BinObject.FILL_NORM_DEPTH=false;
		final Cluster a=GradeBins.loadCluster(dir+"/bins/a.fa"), b=GradeBins.loadCluster(dir+"/bins/b.fa");
		final IntLongHashMap sizes=new IntLongHashMap();
		for(Contig c : a){sizes.increment(c.labelTaxid, c.size());}
		for(Contig c : b){sizes.increment(c.labelTaxid, c.size());}
		a.calcContam(sizes); b.calcContam(sizes);
		final ArrayList<BinStats> bins=new ArrayList<BinStats>();
		bins.add(new BinStats(a, "a")); bins.add(new BinStats(b, "b"));
		require(bins.get(0).badContigs==1 && bins.get(1).badContigs==0,
			"The labeled fixture must contain one foreign contig among three");
		final Field field=GradeBins.class.getDeclaredField("outstream");
		field.setAccessible(true);
		final PrintStream previous=(PrintStream)field.get(null);
		try{
			final String labeled=capture(field, bins, 2);
			checkBad(labeled, "33.333");
			final String order=labels(labeled);
			for(int source : new int[]{BinObject.NONE, BinObject.CHECKM2, BinObject.EUKCC, BinObject.PROKCC}){
				// A mixed set cannot claim a complete per-contig contamination count.
				bins.get(1).truthSource=source;
				for(int taxa : new int[]{0, 2}){
					final String estimated=capture(field, bins, taxa);
					checkBad(estimated, "NA");
					require(order.equals(labels(estimated)), "Summary field order changed with truth availability");
				}
			}
			final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
			final PrintStream stream=new PrintStream(bytes, true, "UTF-8");
			field.set(null, stream);
			final LongList lengths=new LongList(); lengths.add(100); lengths.add(40);
			GradeBins.printL90(lengths, 140);
			require(bytes.toString("UTF-8").contains("L50:"),
				"L/N lines must use the same report stream as the adjacent summary sections");
			stream.close();
		}finally{field.set(null, previous);}
		System.out.println("GRADEBINS_REPORTING_PASS labeled=numeric estimates=NA mixed=NA order=identical redirected_LN=present");
	}

	private static String capture(Field field, ArrayList<BinStats> bins, int taxa) throws Exception{
		assert(bins!=null) : "The shared reporting methods require a completed bin list";
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final PrintStream stream=new PrintStream(bytes, true, "UTF-8");
		field.set(null, stream);
		GradeBins.printCleanDirty(bins);
		GradeBins.printScore(bins, 140, 3, taxa, true);
		stream.close();
		return bytes.toString("UTF-8");
	}

	private static void checkBad(String text, String expected){
		int count=0;
		for(String line : text.split("\n")){
			if(line.startsWith("Bad Contigs:")){
				require(line.substring(line.indexOf('\t')+1).trim().equals(expected), "Incorrect Bad Contigs value: "+line);
				count++;
			}
		}
		require(count==2, "Both existing summary sections must retain the Bad Contigs field");
	}

	private static String labels(String text){
		final StringBuilder result=new StringBuilder();
		for(String line : text.split("\n")){
			final int tab=line.indexOf('\t');
			require(tab>0, "Summary field lacks its value separator: "+line);
			result.append(line.substring(0, tab)).append('\n');
		}
		return result.toString();
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
