package test.prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Hand-stated graph/codec fixtures and a native one-versus-many worker scan.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSupportTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSupportTest.class, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args, pp.outstream);
		}
	}

	/** Main retains selected-log ownership even when malformed-input controls fail. */
	private static void run(final String[] args, final PrintStream outstream) throws Exception{
		final Parser parser=new Parser();
		String mode="test";
		String installation=null;
		for(final String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){
				mode=b;
			}else if(a.equals("installation")){
				installation=b;
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown test option: "+arg);
			}
		}
		if(mode.equals("badlibrary")){
			Shared.setBufferLen(4);
			NcrnaSupport.loadLibrary(parser.in1);
			throw new IllegalStateException("Malformed library unexpectedly returned normally");
		}
		if(parser.out1==null){throw new IllegalArgumentException("out= must name a fresh fixture directory");}
		final File directory=new File(parser.out1);
		if(!directory.mkdir()){throw new IllegalArgumentException("Fixture directory must be new: "+directory);}
		graphCases(directory);
		textOwnerCases(directory);
		scanCase(directory);
		malformedLibrary(directory, installation);
		outstream.println("Support graph checks passed: literal original/reduced graphs, malformed inputs, native worker scan and text-owner cleanup.");
	}

	private static void graphCases(final File directory){
		final String fixture=HEADER+
			"q0\t0\toP000\t0:oP000,1:o@000\n"+
			"q1\t1\to@000\t1:o@000,2:o:000\n"+
			"q2\t2\to0000\t.\n";
		final String original=path(directory, "original.tsv");
		write(original, fixture);
		final Graph graph=NcrnaSupportIO.read(original);
		assert(!graph.reduced && graph.queryNames.length==3 && graph.models.length==4) : "Literal original graph includes the below-cutoff third query";
		assert(graph.identities[3]==0.65625f) : "The first colon separates model ID; a later colon is an A48 digit";
		assert(graph.bestModels[2]==2 && graph.bestIdentities[2]==0.5f) : "Original unsupported genes retain their full-library assignment";
		final Graph reduced=graph.subset(new boolean[]{false, false, true});
		assert(reduced.reduced && reduced.queryNames.length==3 && reduced.modelNames[0].equals("m2")) : "Reduced graph retains query population and remaps its model axis";
		assert(reduced.bestModels[0]==-1 && Float.isNaN(reduced.bestIdentities[0])) : "Newly unsupported genes require NONE/NA";
		assert(reduced.bestModels[1]==0 && reduced.bestIdentities[1]==0.65625f) : "Best retained supporter uses the stored edge identity";
		assert(reduced.bestModels[2]==-1 && Float.isNaN(reduced.bestIdentities[2])) : "Reduced scope must not mislabel original below-cutoff best as retained support";
		final String output=path(directory, "reduced.tsv");
		NcrnaSupportIO.write(reduced, output);
		final Graph reread=NcrnaSupportIO.read(output);
		assert(Arrays.equals(reread.models, new int[]{0}) && reread.identities[0]==0.65625f) : "Consumer must recover the literal expected reduced edge";
		reject(directory, "duplicate.tsv", fixture.replace("q2\t", "q1\t"));
		reject(directory, "truncated.tsv", fixture.substring(0, fixture.lastIndexOf("q2\t")));
		reject(directory, "bad_float.tsv", fixture.replace("2:o:000", "2:zzzzz"));
		reject(directory, "duplicate_edge.tsv", fixture.replace("0:oP000,1:o@000", "0:oP000,0:o@000"));
		reject(directory, "bad_best.tsv", fixture.replace("q0\t0\toP000", "q0\t2\toP000"));
		reject(directory, "bad_scope.tsv", fixture.replace("ORIGINAL", "UNKNOWN"));
		reject(directory, "best_overflow.tsv", fixture.replace("q0\t0\t", "q0\t4294967296\t"));
		reject(directory, "edge_overflow.tsv", fixture.replace("0:oP000", "18446744073709551616:oP000"));
		reject(directory, "count_overflow.tsv", fixture.replace("#models\t3", "#models\t18446744073709551619"));
		reject(directory, "axis_overflow.tsv", fixture.replace("#model\t0", "#model\t4294967296"));
	}

	/** Native writer errorState exercises close failures without a fault framework. */
	private static void textOwnerCases(final File directory){
		final ByteStreamWriter successful;
		try(Output output=new Output(path(directory, "owner_success.txt"))){
			successful=output.writer;
			successful.print("complete\n");
		}
		assert(successful.getState()==Thread.State.TERMINATED && !successful.errorState)
			: "Successful owner close must finish and join its native writer";

		final IllegalArgumentException primary=new IllegalArgumentException("Known fixture body failure");
		ByteStreamWriter failed=null;
		try(Output output=new Output(path(directory, "owner_body_failure.txt"))){
			failed=output.writer;
			failed.errorState=true;
			failed.print("before failure\n");
			throw primary;
		}catch(IllegalArgumentException caught){
			assert(caught==primary && caught.getSuppressed().length==1
				&& caught.getSuppressed()[0] instanceof IllegalStateException)
				: "A native close error must remain suppressed behind the original body exception";
		}
		assert(failed!=null && failed.getState()==Thread.State.TERMINATED)
			: "Body and close failure together must still join the native writer";

		ByteStreamWriter closeOnly=null;
		boolean rejected=false;
		try(Output output=new Output(path(directory, "owner_close_failure.txt"))){
			closeOnly=output.writer;
			closeOnly.errorState=true;
			closeOnly.print("close failure\n");
		}catch(IllegalStateException caught){
			rejected=caught.getMessage().contains("Could not finish writing") && caught.getSuppressed().length==0;
		}
		assert(rejected && closeOnly!=null && closeOnly.getState()==Thread.State.TERMINATED)
			: "A close-only error must propagate after its native writer has terminated";
	}

	/** Both workers and the serial path must retain every identical gene exactly once. */
	private static void scanCase(final File directory){
		final String first="ACGTGCTAGTCAGATCGATGCTACGATCGTAGCATGCTAGTCAGCTAGCATGCTACGATCGTAC";
		final String second="GGATATGCAATCGGATACCGTATGCCGATGATACCTGACTGAGCTTCGATACGCGTAGATCAT";
		final String unsupported="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		final String ref=path(directory, "models.fa");
		final String genes=path(directory, "genes.fa");
		write(ref, ">m0 first model\n"+first+"\n>m1 second model\n"+second+"\n");
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<12000; i++){
			final int kind=i%3;
			bb.append('>').append(kind==0 ? 'a' : kind==1 ? 'b' : 'x').append(i).nl();
			bb.append(kind==0 ? first : kind==1 ? second : unsupported).nl();
		}
		write(genes, bb.toString());
		final Library library=NcrnaSupport.loadLibrary(ref);
		assert(library.headers[0].equals("m0 first model") && library.names[0].equals("m0")) : "Keep original FASTA descriptions separately from the graph's single-token model IDs";
		final FileFormat input=FileFormat.testInput(genes, FileFormat.FASTA, null, true, true);
		final Library retained=NcrnaSupport.loadLibrary(genes);
		final IntList members=new IntList();
		for(int i=retained.names.length-1; i>=0; i-=2){members.add(i);}
		for(final int threads : new int[]{1, 3}){
			final NcrnaSupport.ScanResult scan=NcrnaSupport.scanWithStats(input, library, 1, threads, 1);
			final Graph graph=scan.graph;
			assert(scan.batches>1 && (threads==1 || scan.activeWorkers>1)) : "The concurrency fixture must actually reach multiple batches and multiple consumers";
			assert(graph.queryNames.length==12000 && graph.models.length==8000) : "Multiple native batches must retain all genes and exactly the known identical edges";
			for(int query=0; query<12000; query++){
				final char kind=graph.queryNames[query].charAt(0);
				final int from=graph.starts[query], to=graph.starts[query+1];
				if(kind=='x'){
					assert(from==to) : "An all-T sequence is not identical to either diverse model";
				}else{
					final int expected=kind=='a' ? 0 : 1;
					assert(to==from+1 && graph.models[from]==expected && graph.identities[from]==1) : "Each exact gene must retain its corresponding model, independent of worker order";
				}
			}
			System.err.println("Scan fixture threads="+threads+" batches="+scan.batches+" active_workers="+scan.activeWorkers);
			final NcrnaSupport.ScanResult memory=NcrnaSupport.scanWithStats(retained, null, library, 1, threads);
			assert(memory.alignments==24000 && memory.graph.queryNames.length==12000)
				: "In-memory scanning must score every retained gene against both models: "+memory.alignments;
			compareScores(graph, memory.graph);
			final NcrnaSupport.ScanResult subset=NcrnaSupport.scanWithStats(retained, members, library, 1, threads);
			assert(subset.alignments==12000 && subset.graph.queryNames.length==6000)
				: "Member scanning must score only the selected half of the corpus: "+subset.alignments;
			compareScores(graph, subset.graph);
			for(final String name:subset.graph.queryNames){
				assert((Integer.parseInt(name.substring(1))&1)==1) : "Reversed member indices select odd-numbered genes only: "+name;
			}
		}
	}

	/** Compares named results, allowing legitimate worker-dependent row ordering. */
	private static void compareScores(final Graph full, final Graph selected){
		final ObjectIntMap<String> rows=new ObjectIntMap<String>(String.class);
		for(int i=0; i<full.queryNames.length; i++){rows.put(full.queryNames[i], i);}
		for(int q=0; q<selected.queryNames.length; q++){
			final int original=rows.get(selected.queryNames[q]);
			assert(original>=0) : "The selected scan may not invent a gene: "+selected.queryNames[q];
			assert(full.bestModels[original]==selected.bestModels[q]
				&& full.bestIdentities[original]==selected.bestIdentities[q])
				: "Stream and memory scans must agree on the overall winner and native score for "+selected.queryNames[q];
			final int from=full.starts[original], start=selected.starts[q];
			final int count=full.starts[original+1]-from;
			assert(count==selected.starts[q+1]-start) : "Every support edge must survive memory scanning for "+selected.queryNames[q];
			for(int edge=0; edge<count; edge++){
				assert(full.models[from+edge]==selected.models[start+edge]
					&& full.identities[from+edge]==selected.identities[start+edge])
					: "Retained support edges must preserve model and native identity for "+selected.queryNames[q]+", edge="+edge;
			}
		}
	}

	/** Failure on the first record must terminate despite queued unread FASTA data. */
	private static void malformedLibrary(final File directory, final String installation) throws Exception{
		if(installation==null){throw new IllegalArgumentException("Fixture launcher must supply installation=");}
		final String filename=path(directory, "bad_library.fa");
		final ByteBuilder bb=new ByteBuilder();
		bb.append(">bad").append((byte)1).append("name\nACGTACGT\n");
		for(int i=0; i<65536; i++){
			bb.append(">good").append(i).append("\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
		}
		try(Output output=new Output(filename)){output.writer.print(bb);}
		final Process child=new ProcessBuilder("bash", installation+"/dev/testncrnasupport.sh",
			"mode=badlibrary", "in="+filename).redirectOutput(new File(directory, "bad_library.out"))
			.redirectError(new File(directory, "bad_library.err")).start();
		try{
			final boolean completed=child.waitFor(20, TimeUnit.SECONDS);
			assert(completed) : "Malformed first record must not strand a native FASTA producer";
			assert(child.exitValue()!=0) : "Malformed library must terminate unsuccessfully";
			final String error=ReadWrite.readString(path(directory, "bad_library.err"));
			assert(error.contains("FASTA headers must be printable ASCII")) : "The child must fail on the intended malformed header, not an unrelated startup error";
		}finally{
			child.destroyForcibly();
		}
	}

	private static void reject(final File directory, final String name, final String text){
		final String filename=path(directory, name);
		write(filename, text);
		boolean rejected=false;
		try{
			NcrnaSupportIO.read(filename);
		}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "Malformed graph fixture must fail: "+name;
	}

	private static String path(final File directory, final String name){return new File(directory, name).getPath();}

	private static void write(final String path, final String text){
		try(Output output=new Output(path)){output.writer.print(text);}
	}

	private static final String HEADER="#ncrna_support_graph\t2\n#scope\tORIGINAL\n#aligner\tquantum\n"+
		"#cutoff_a48\to8000\n#models\t3\n#model\t0\tm0\n#model\t1\tm1\n#model\t2\tm2\n#queries\t3\n"+
		"query_id\tbest_model_index\tbest_identity_a48\tsupporters\n";
}
