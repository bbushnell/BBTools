package assemble;

import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Reports exact read-kmer evidence along selected reference windows.
 * Coordinates describe forward-reference kmer ends, one-based and inclusive.
 * Reference sequence selects queries only; it never contributes counts.
 * @author Fischl
 */
public final class KmerWindowDepth {

	/** Loads one exact table and emits depths without assembling or mutating it. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, KmerWindowDepth.class, false);
		try{new KmerWindowDepth(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Separates reference-window options from the existing native table parser. */
	private KmerWindowDepth(final String[] args){
		assert(args!=null) : "PreParser supplies expanded diagnostic arguments.";
		final ArrayList<String> tableArgs=new ArrayList<String>();
		String ranges=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("ref")){reference=b;}
			else if(a.equals("scaffold")){scaffold=b;}
			else if(a.equals("regions")){ranges=b;}
			else if(a.equals("out")){output=b;}
			else{tableArgs.add(arg);}
		}
		if(reference==null || scaffold==null || ranges==null || output==null){
			throw new IllegalArgumentException("Required: in=reads ref=fasta scaffold=name regions=start-end[,start-end] out=tsv k=N");
		}
		for(String range : ranges.split(",", -1)){
			final String[] ends=range.split("-", -1);
			if(ends.length!=2){throw new IllegalArgumentException("Invalid inclusive region: "+range);}
			regions.add(new int[]{Integer.parseInt(ends[0]), Integer.parseInt(ends[1])});
		}
		// Match the assembly's unfiltered canonical counts, not an approximate sketch.
		tableArgs.add("prefilter=f");
		tableArgs.add("prealloc=f");
		tableArgs.add("hashonly=f");
		tableArgs.add("rcomp=t");
		tableArgs.add("minprob=0");
		tableArgs.add("minprobmain=f");
		tableArgs.add("packed=t");
		tables=new KmerTableSetU(tableArgs.toArray(new String[0]), 0);
		k=tables.kbig;
		final ArrayList<String> inputs=new ArrayList<String>(tables.in1);
		inputs.addAll(tables.in2);
		inputs.addAll(tables.extra);
		TadpoleGraph.checkPaths(reference, inputs, output);
		if(!Tools.testInputFiles(false, true, reference) ||
				!Tools.testOutputFiles(false, false, false, output)){
			throw new IllegalArgumentException("Missing reference or existing/unwritable output.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates every requested base before loading reads or opening output. */
	private byte[] referenceBases(){
		byte[] bases=null;
		for(Contig contig : TadpoleGraph.readContigs(reference)){
			final String name=contig.name;
			int length=0;
			while(length<name.length() && !Character.isWhitespace(name.charAt(length))){length++;}
			if(!scaffold.equals(name.substring(0, length))){continue;}
			if(bases!=null){throw new IllegalArgumentException("Duplicate scaffold identifier: "+scaffold);}
			bases=contig.bases;
		}
		if(bases==null){throw new IllegalArgumentException("Scaffold not found: "+scaffold);}
		for(int[] region : regions){
			if(region[0]<k || region[1]<region[0] || region[1]>bases.length){
				throw new IllegalArgumentException("Invalid kmer-end region "+region[0]+"-"+region[1]+" for K="+k+", length="+bases.length);
			}
			for(int pos=region[0]-k; pos<region[1]; pos++){
				if(AminoAcid.baseToNumber[bases[pos]]<0){
					throw new IllegalArgumentException("Ambiguous reference base at "+scaffold+":"+(pos+1));
				}
			}
		}
		return bases;
	}

	/** Owns writer and table cleanup even when parsing or output fails. */
	private void process(){
		final byte[] bases=referenceBases();
		ByteStreamWriter writer=null;
		boolean complete=false, writeError=false;
		long rows=0;
		try{
			tables.process(new Timer());
			writer=new ByteStreamWriter(output, false, false, false);
			writer.start();
			writer.print("region\tscaffold\tk\tstart1\tend1\tword\tdepth\tleft_A\tleft_C\tleft_G\tleft_T\tright_A\tright_C\tright_G\tright_T\tref_previous\tref_next\n");
			final Kmer word=new Kmer(k), neighbor=new Kmer(k);
			final ByteBuilder row=new ByteBuilder(k+200);
			for(int index=0; index<regions.size(); index++){
				final int[] region=regions.get(index);
				word.clear();
				for(int pos=region[0]-k; pos<region[0]-1; pos++){
					word.addRightNumeric(AminoAcid.baseToNumber[bases[pos]]);
				}
				for(int end=region[0]; end<=region[1]; end++){
					word.addRightNumeric(AminoAcid.baseToNumber[bases[end-1]]);
					assert(word.len>=k) : "Each output row must contain a complete query word.";
					row.clear().append(index+1).tab().append(scaffold).tab().append(k);
					row.tab().append(end-k+1).tab().append(end).tab().append(bases, end-k, k);
					row.tab().append(depth(word));
					for(int base=0; base<4; base++){
						neighbor.setFrom(word); neighbor.addLeftNumeric(base);
						row.tab().append(depth(neighbor));
					}
					for(int base=0; base<4; base++){
						neighbor.setFrom(word); neighbor.addRightNumeric(base);
						row.tab().append(depth(neighbor));
					}
					row.tab().append(end>k ? (char)bases[end-k-1] : '.');
					row.tab().append(end<bases.length ? (char)bases[end] : '.').nl();
					writer.print(row);
					rows++;
				}
			}
			complete=true;
		}finally{
			if(writer!=null){writeError=writer.poisonAndWait();}
			tables.clear();
			if(writeError && !complete){System.err.println("Window output also failed during diagnostic cleanup.");}
		}
		if(writeError){throw new RuntimeException("Window output failed: "+output);}
		System.err.println("KMER_WINDOW_DEPTH_PASS rows="+rows+" k="+k);
	}

	/** Normalizes the table's absent-key sentinel to a zero observation count. */
	private int depth(final Kmer word){
		assert(word.len>=k) : "Exact depth queries cannot use incomplete reference words.";
		return Math.max(0, tables.getCount(word));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String reference, scaffold, output;
	private final ArrayList<int[]> regions=new ArrayList<int[]>();
	private final KmerTableSetU tables;
	private final int k;
}
