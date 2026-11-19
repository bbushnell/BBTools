package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Exports selected original consensus/HBM pairs without fitting or scoring.
 * The keep= FASTA must be an unchanged subset of in=. The names= file contains
 * distinct native IDs to add, disjoint from keep=. Unknown or duplicate names
 * fail before the output directory is created. Original order, full headers,
 * upper-case bases and all HBM counts are preserved by NcrnaModelPairs.
 * Its existing writer may normalize FASTA wrapping/HBM encoding. Every written
 * pair is decoded and compared with the selected original before success.
 * This is a paired library, not a loadable family or caller verification.
 * Single-use and single-threaded; arrays/models are owned by this instance.
 * Adapted from template.A_SampleByteFile using existing native library I/O.
 * @author Raiden
 */
public final class ExportNcrnaModelSubset{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ExportNcrnaModelSubset.class, false);
		try{new ExportNcrnaModelSubset(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private ExportNcrnaModelSubset(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		String modelPath=null, keepPath=null, namePath=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("hbm")){modelPath=b;}
			else if(a.equals("keep")){keepPath=b;}
			else if(a.equals("names")){namePath=b;}
			else{require(parser.parse(arg, a, b), "Unknown paired subset parameter: "+arg);}
		}
		input=parser.in1;
		hbm=modelPath;
		keep=keepPath;
		names=namePath;
		require(input!=null && hbm!=null && keep!=null && names!=null && parser.out1!=null
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply in=original.fa hbm=original.hbm keep=unchanged_subset.fa names=unique_additions.txt out=fresh_directory");
		directory=new File(parser.out1);
		require(!directory.exists() && Tools.testInputFiles(false, true, input, hbm, keep, names),
			"Require readable original pairs/selection and a fresh output directory");
		Shared.setThreads(1);
		Shared.TRIM_READ_DESCRIPTION=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final Library original=NcrnaSupport.loadLibrary(input), retained=NcrnaSupport.loadLibrary(keep);
		final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		for(int i=0; i<original.names.length; i++){index.put(original.names[i], i);}
		final boolean[] selected=new boolean[original.names.length];
		for(int i=0; i<retained.names.length; i++){
			final int position=index.get(retained.names[i]);
			require(position>=0, "Retained model absent from original library: "+retained.names[i]);
			require(!selected[position] && original.headers[position].equals(retained.headers[i])
				&& Arrays.equals(original.bases[position], retained.bases[i]),
				"Retained FASTA must preserve a unique original header and bases: "+retained.names[i]);
			selected[position]=true;
		}
		int added=0;
		try(Input in=new Input(names)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				require(row.length>0, "Empty model name in addition list");
				for(final byte b:row){require(b>32 && b<127, "Addition list requires one native ASCII identifier per line, without descriptions");}
				final String name=new String(row, StandardCharsets.US_ASCII);
				final int position=index.get(name);
				require(position>=0, "Addition model absent from original library: "+name);
				require(!selected[position], "Duplicate addition or already-retained model: "+name);
				selected[position]=true;
				added++;
			}
		}
		// The existing selector validates EVERY original FASTA/HBM pair, including
		// unselected pairs. The complete record is retained; nothing is refitted.
		final NcrnaModelPairs expected=new NcrnaModelPairs(original, hbm, selected);
		require(directory.mkdir(), "Could not create fresh paired-library directory: "+directory);
		final String filename=ReduceNcrnaModels.modelFilename(hbm);
		expected.write(directory, filename);
		verify(expected, filename);
		try(Output summary=new Output(new File(directory, "summary.tsv").getPath())){
			summary.writer.println("original_models\tretained_models\tadded_models\texported_models\tqualification");
			summary.writer.print(new ByteBuilder().append(original.names.length).tab().append(retained.names.length).tab()
				.append(added).tab().append(expected.names.length).append("\tORIGINAL_PAIRED_RECORDS_VERIFIED_NOT_CALLER_VERIFIED\n"));
		}
	}

	/** Checks all full decoded records; deliberately makes no whole-file byte claim. */
	private void verify(final NcrnaModelPairs expected, final String filename){
		assert(expected.models!=null) : "Paired export always includes complete original HBM records";
		final Library exported=NcrnaSupport.loadLibrary(new File(directory, "consensus.fa").getPath());
		final boolean[] all=new boolean[exported.names.length];
		Arrays.fill(all, true);
		final NcrnaModelPairs actual=new NcrnaModelPairs(exported, new File(directory, filename).getPath(), all);
		require(Arrays.equals(expected.names, actual.names), "Exported model names/order differ from original selection");
		for(int i=0; i<expected.names.length; i++){
			require(expected.same(actual, i), "Exported full consensus/HBM record differs from original: "+expected.names[i]);
		}
	}

	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String input, hbm, keep, names;
	private final File directory;
}
