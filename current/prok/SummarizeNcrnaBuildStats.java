package prok;

import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Summarizes retained final-consensus identity histograms without alignment.
 * ByteFile ownership and flag parsing follow template.A_SampleByteFile.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaBuildStats{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaBuildStats.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String labels=null;
			String cutoffs="0.64,0.70,0.75";
			long expected=-1;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("labels")){labels=b;}
				else if(a.equals("cutoffs")){cutoffs=b;}
				else if(a.equals("expectedmembers")){expected=Parse.parseKMG(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown build-summary parameter: "+arg);}
			}
			if(parser.in1==null || parser.out1==null || labels==null || cutoffs==null || expected<1
					|| parser.overwrite || parser.append || parser.in2!=null || parser.out2!=null
					|| parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Use in=comma-separated build_stats.tsv files, labels=matching names, out=fresh.tsv and expectedmembers=N");
			}
			summarize(parser.in1.split(",", -1), labels.split(",", -1), parser.out1, expected, parseCuts(cutoffs));
			pp.outstream.println("NCRNA_FINAL_CONSENSUS_IDENTITY_HISTOGRAM_COMPLETE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Cutoffs must coincide with histogram bin centers to preserve strict bounds. */
	private static int[] parseCuts(final String text){
		final String[] terms=text.split(",", -1);
		final int[] cuts=new int[terms.length];
		int previous=-1;
		for(int i=0; i<terms.length; i++){
			final double value=Double.parseDouble(terms[i]);
			final int bin=(int)Math.round(value*1000);
			if(!Double.isFinite(value) || value<0 || value>1 || Math.abs(value*1000-bin)>0.000001 || bin<=previous){
				throw new IllegalArgumentException("Cutoffs need increasing distinct bin centers in [0,1], at resolution 0.001: "+text);
			}
			cuts[i]=bin;
			previous=bin;
		}
		return cuts;
	}

	/** No sequence or alignment inputs are used; each count represents one member. */
	public static void summarize(final String[] inputs, final String[] labels, final String output,
			final long expected, final int[] cuts){
		if(inputs.length==0 || inputs.length!=labels.length || expected<1 || cuts.length==0){
			throw new IllegalArgumentException("Each input histogram needs one label and a positive expected population");
		}
		final String[] paths=new String[inputs.length+1];
		System.arraycopy(inputs, 0, paths, 0, inputs.length);
		paths[inputs.length]=output;
		if(!Tools.testForDuplicateFiles(true, paths)){throw new IllegalArgumentException("Histogram inputs and output must not alias");}
		final ByteBuilder row=new ByteBuilder();
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(output)){
			out.writer.println("policy\tmembers\tq25_rounded\tmedian_rounded\tq75_rounded\tcut\tbelow_min\tbelow_max\tbelow_percent_min\tbelow_percent_max");
			for(int i=0; i<inputs.length; i++){
				if(labels[i].isEmpty() || labels[i].indexOf('\t')>=0 || labels[i].indexOf('\n')>=0){
					throw new IllegalArgumentException("Histogram labels must be nonempty TSV fields");
				}
				for(int j=0; j<i; j++){
					if(labels[i].equals(labels[j])){throw new IllegalArgumentException("Each histogram needs a distinct label: "+labels[i]);}
				}
				final long[] histogram=readHistogram(inputs[i], expected);
				final int q1=quantile(histogram, expected/4+(expected%4>0 ? 1 : 0));
				final int q2=quantile(histogram, expected/2+expected%2);
				final int q3=quantile(histogram, expected-expected/4);
				for(final int cut:cuts){
					final long lower=below(histogram, cut);
					final long upper=Math.addExact(lower, histogram[cut]);
					assert(lower<=upper && upper<=expected) : "Rounded threshold-bin uncertainty must stay inside the conserved member population";
					row.clear().append(labels[i]).tab().append(expected).tab()
						.append(q1/1000.0, 3).tab().append(q2/1000.0, 3).tab().append(q3/1000.0, 3).tab()
						.append(cut/1000.0, 3).tab().append(lower).tab().append(upper).tab()
						.append(100.0*lower/expected, 4).tab().append(100.0*upper/expected, 4).nl();
					out.writer.print(row);
				}
			}
		}
	}

	/** Reads the producer's rounded histogram and checks its complete denominator. */
	public static long[] readHistogram(final String path, final long expected){
		final long[] histogram=new long[1001];
		final boolean[] seen=new boolean[1001];
		final LineParser1 parser=new LineParser1('\t');
		long members=0;
		try(NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			RrnaResourceIO.header(input.reader, "kind\tvalue\tcount", path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				parser.set(line);
				if(parser.terms()!=3){throw new IllegalArgumentException("Build-stat row needs three fields: "+path);}
				if(!parser.termEquals("identity_nearest_thousandth", 0)){continue;}
				final int bin=parser.parseInt(1);
				final long count=parser.parseLong(2);
				if(bin<0 || bin>=histogram.length || count<1 || seen[bin]){
					throw new IllegalArgumentException("Identity histogram needs distinct bins 0..1000 and positive counts: "+path);
				}
				seen[bin]=true;
				histogram[bin]=count;
				members=Math.addExact(members, count);
			}
		}
		if(members!=expected){throw new IllegalArgumentException("Final-consensus identity population "+members+" differs from expected "+expected+": "+path);}
		return histogram;
	}

	/** Nearest-rank quantile of the stored bins, not an unrounded point estimate. */
	public static int quantile(final long[] histogram, final long rank){
		if(rank<1){throw new IllegalArgumentException("Quantile rank must be positive");}
		long count=0;
		for(int bin=0; bin<histogram.length; bin++){
			count=Math.addExact(count, histogram[bin]);
			if(count>=rank){return bin;}
		}
		throw new IllegalArgumentException("Quantile rank exceeds the complete histogram population: "+rank);
	}

	/** Only bins strictly below the cut are certain; the cut bin stays uncertain. */
	public static long below(final long[] histogram, final int cut){
		if(cut<0 || cut>=histogram.length){throw new IllegalArgumentException("Identity cut is outside histogram bins: "+cut);}
		long count=0;
		for(int bin=0; bin<cut; bin++){count=Math.addExact(count, histogram[bin]);}
		return count;
	}
}
