package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec.Model;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Development experiment: regenerate models sequentially while retaining member support.
 * Native CLI and sequence I/O follow template.A_SampleStreamer. Consensus generation
 * remains the public trnaconsensus.sh tool; alignment remains NcrnaSupport.
 * Models are visited in input order. A candidate and its member assignments are
 * committed together only when every former member still has a supporting model.
 * A final whole-corpus scan settles other migrations; it never starts another round.
 * This driver is not a production reducer stage or a genomic accuracy evaluator.
 * Preserved standalone experiment: created to regenerate a supplied fixed library
 * without first reducing it. ReduceNcrnaModels now owns the normal reduction,
 * protection and regeneration workflow, but has no equivalent standalone mode.
 * Retained pending disposition for future direct-library regeneration; unused in
 * the shipping caller and not evidence of a performance or accuracy advantage.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSequentialRegenerationExperiment{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns diagnostics across both argument validation and the complete experiment. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSequentialRegenerationExperiment.class, false);
		try{
			boolean protect=false;
			for(String arg:pp.args){if(arg.equalsIgnoreCase("mode=protect")){protect=true;}}
			if(protect){
				NcrnaCallerReversionExperiment.run(pp.args, pp.outstream);
			}else{
				new NcrnaSequentialRegenerationExperiment(pp.args).process(pp.outstream);
			}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses a complete corpus, a fixed library and an explicit public consensus recipe. */
	private NcrnaSequentialRegenerationExperiment(final String[] args){
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("ref")){reference=b;}
			else if(a.equals("hbm")){hbm=b;}
			else if(a.equals("graph")){graphPath=b;}
			else if(a.equals("consensusconfig")){consensusConfig=b;}
			else if(a.equals("installation")){installation=b;}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(a.equals("cutoff")){cutoff=Float.parseFloat(b);}
			else if(a.equals("childheapmib")){childHeap=Parse.parseIntKMG(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		input=parser.in1;
		threads=parser.workers()>0 ? parser.workers() : Shared.threads();
		if(input==null || reference==null || hbm==null || consensusConfig==null || installation==null
				|| parser.out1==null || expectedGenes<1 || expectedModels<1 || threads<1 || childHeap<1
				|| !(cutoff>=0 && cutoff<=1)){
			throw new IllegalArgumentException("Supply in/ref/hbm/consensusconfig/installation/out, explicit positive gene/model counts, cutoff=[0,1] and positive threads/childheapmib");
		}
		if(parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1 || parser.append){
			throw new IllegalArgumentException("Sequential regeneration requires the entire unpaired corpus and fresh output");
		}
		if(!Tools.testInputFiles(true, false, input, reference, hbm, graphPath, consensusConfig)){
			throw new IllegalArgumentException("Cannot read an experiment input");
		}
		root=new File(parser.out1).getAbsoluteFile();
		if(root.exists() || !root.getParentFile().isDirectory()){
			throw new IllegalArgumentException("out= must be a fresh directory under an existing parent");
		}
		child=new NcrnaChildProcess(new File(installation));
		ReadWrite.setZipThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Experiment          ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains genes once; only the edited model's current members are rescored per step. */
	private void process(final PrintStream diagnostics) throws Exception{
		final Timer timer=new Timer();
		library=NcrnaSupport.loadLibrary(reference);
		if(library.names.length!=expectedModels){throw new IllegalArgumentException("Unexpected base model count");}
		models=loadModels(library, hbm);
		headers=library.headers.clone();
		genes=readFasta(input);
		if(genes.size()!=expectedGenes){throw new IllegalArgumentException("Expected complete corpus of "+expectedGenes+", found "+genes.size());}
		index=new ObjectIntMap<String>(String.class);
		for(int gene=0; gene<genes.size(); gene++){
			final Read read=genes.get(gene);
			read.id=NcrnaSupport.recordName(read.id);
			if(TrnaConsensusBuilder.parseAnticodon(read.id)!=null || index.contains(read.id)){
				throw new IllegalArgumentException("Experiment gene IDs must be unique and must not create anticodon subgroups: "+read.id);
			}
			index.put(read.id, gene);
		}
		final Graph initial=graphPath==null ? scan(input, library) : NcrnaSupportIO.read(graphPath);
		checkAxes(initial, library, expectedGenes);
		owners=assignments(initial, index);
		initialOwners=owners.clone();
		if(!root.mkdir()){throw new IllegalStateException("Could not create fresh experiment directory");}
		regenerateModels(diagnostics);
		assert(attempted==kept+reverted && attempted+empty==expectedModels)
			: "Every fixed model index must have exactly one kept, reverted or empty decision";
		final Graph finalGraph=scan(input, library);
		checkAxes(finalGraph, library, expectedGenes);
		final int[] finalOwners=assignments(finalGraph, index);
		NcrnaSupportIO.write(finalGraph, path("final.support.tsv"));
		int initiallySupported=0, finallySupported=0, lost=0, gained=0, finalMigrations=0;
		try(Output output=new Output(path("assignments.tsv"))){
			output.writer.print("gene\tinitial_model\tafter_series_model\tfinal_model\n");
			final ByteBuilder row=new ByteBuilder();
			for(int gene=0; gene<owners.length; gene++){
				if(initialOwners[gene]>=0){initiallySupported++;}
				if(finalOwners[gene]>=0){finallySupported++;}
				if(initialOwners[gene]>=0 && finalOwners[gene]<0){lost++;}
				if(initialOwners[gene]<0 && finalOwners[gene]>=0){gained++;}
				if(owners[gene]>=0 && owners[gene]!=finalOwners[gene]){finalMigrations++;}
				row.clear();
				row.append(genes.get(gene).id).tab().append(name(library, initialOwners[gene])).tab();
				row.append(name(library, owners[gene])).tab().append(name(library, finalOwners[gene])).nl();
				output.writer.print(row);
			}
		}
		if(lost!=0){throw new IllegalStateException("Sequential support protection failed for "+lost+" initially supported genes; inspect assignments.tsv");}
		final ArrayList<Read> consensus=new ArrayList<Read>(expectedModels);
		try(NucleotideHbmCodec.Writer output=new NucleotideHbmCodec.Writer(path("models.hbmc"), true, false)){
			for(int model=0; model<expectedModels; model++){
				assert(Arrays.equals(library.bases[model], models[model].bases)) : "Accepted consensus and HBM must retain one coordinate frame";
				consensus.add(new Read(library.bases[model], null, headers[model], model));
				output.add(models[model]);
			}
			output.finish();
		}
		writeFasta(path("consensus.fa"), consensus);
		timer.stop();
		final ByteBuilder report=new ByteBuilder("models\tgenes\tinitial_supported\tinitial_unsupported\tfinal_supported\tnewly_supported\tfinal_migrations\tattempted\tkept\treverted\tempty\tno_candidate\tmember_visits\tmember_migrations\talignments\tseconds\n");
		report.append(expectedModels).tab().append(expectedGenes).tab().append(initiallySupported).tab().append(expectedGenes-initiallySupported).tab();
		report.append(finallySupported).tab().append(gained).tab().append(finalMigrations).tab().append(attempted).tab().append(kept).tab();
		report.append(reverted).tab().append(empty).tab().append(noCandidate).tab().append(memberVisits).tab().append(memberMigrations).tab();
		report.append(alignments).tab().append(timer.elapsed*1e-9, 6).nl();
		try(Output output=new Output(path("summary.tsv"))){output.writer.print(report);}
		diagnostics.println("Sequential regeneration completed: "+kept+" kept, "+reverted+" reverted, "+empty+" empty unchanged; "+lost+" supported genes lost.");
	}

	/** Visits each model once in fixed order; newly arriving members join later visits. */
	private void regenerateModels(final PrintStream diagnostics) throws Exception{
		try(Output steps=new Output(path("steps.tsv")); Output moves=new Output(path("moves.tsv"));
				Output failures=new Output(path("unsupported_candidates.tsv"))){
			steps.writer.print("model_index\tmodel\tdecision\tmembers\tunsupported\tmoved\told_length\tnew_length\tseconds\n");
			moves.writer.print("step\tgene\tprevious_model\tnew_model\tidentity\n");
			failures.writer.print("step\tgene\tbest_model\tbest_identity\n");
			final ByteBuilder detail=new ByteBuilder();
			for(int model=0; model<expectedModels; model++){
				final Timer stepTimer=new Timer();
				final IntList members=new IntList();
				for(int gene=0; gene<owners.length; gene++){
					if(owners[gene]==model){members.add(gene);}
				}
				final int oldLength=library.bases[model].length;
				int newLength=oldLength, unsupported=0, moved=0;
				String decision="EMPTY_UNCHANGED";
				if(members.isEmpty()){empty++;}
				else{
					attempted++;
					memberVisits+=members.size;
					final File directory=new File(root, "model_"+model);
					if(!directory.mkdir()){throw new IllegalStateException("Cannot create fresh model directory");}
					final ArrayList<Read> selected=new ArrayList<Read>(members.size);
					for(int i=0; i<members.size; i++){selected.add(genes.get(members.get(i)));}
					final String memberPath=new File(directory, "members.fa").getPath();
					writeFasta(memberPath, selected);
					final Candidate candidate=regenerate(directory, library.names[model]);
					if(candidate==null){
						decision="NO_CANDIDATE_REVERTED";
						reverted++;
						noCandidate++;
					}else{
						newLength=candidate.model.bases.length;
						final byte[][] proposedBases=library.bases.clone();
						proposedBases[model]=candidate.model.bases;
						final Library proposed=new Library(library.names, proposedBases);
						final Graph scored=scan(memberPath, proposed);
						checkAxes(scored, proposed, members.size);
						final int[] geneIndices=new int[scored.queryNames.length];
						for(int row=0; row<geneIndices.length; row++){
							final int gene=index.get(scored.queryNames[row]);
							if(gene<0 || owners[gene]!=model){throw new IllegalStateException("Member scan returned a gene outside the edited model");}
							geneIndices[row]=gene;
							if(scored.bestIdentities[row]<cutoff){
								unsupported++;
								detail.clear();
								final ByteBuilder rowText=detail.append(model).tab().append(scored.queryNames[row]).tab();
								rowText.append(library.names[scored.bestModels[row]]).tab().append(scored.bestIdentities[row], 8).nl();
								failures.writer.print(rowText);
							}
						}
						if(unsupported>0){
							decision="UNSUPPORTED_REVERTED";
							reverted++;
						}else{
							// Nothing was mutated before acceptance: rollback retains the old pair and owners.
							library=proposed;
							models[model]=candidate.model;
							headers[model]=candidate.header;
							for(int row=0; row<geneIndices.length; row++){
								final int gene=geneIndices[row], next=scored.bestModels[row];
								if(next!=model){
									moved++;
									detail.clear();
									final ByteBuilder rowText=detail.append(model).tab().append(genes.get(gene).id).tab();
									rowText.append(library.names[model]).tab().append(library.names[next]).tab().append(scored.bestIdentities[row], 8).nl();
									moves.writer.print(rowText);
								}
								owners[gene]=next;
							}
							decision="KEPT";
							kept++;
							memberMigrations+=moved;
						}
					}
				}
				stepTimer.stop();
				final ByteBuilder row=new ByteBuilder().append(model).tab().append(library.names[model]).tab().append(decision).tab();
				row.append(members.size).tab().append(unsupported).tab().append(moved).tab().append(oldLength).tab().append(newLength).tab();
				row.append(stepTimer.elapsed*1e-9, 6).nl();
				steps.writer.print(row);
				diagnostics.println("Model "+(model+1)+"/"+expectedModels+": "+decision+", members="+members.size+", moved="+moved);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Native Operations      ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the public non-exec launcher with existing Unix-group ownership, preserving paired output. */
	private Candidate regenerate(final File directory, final String modelName) throws Exception{
		final String fasta=new File(directory, "candidate.fa").getPath();
		final String modelPath=new File(directory, "candidate.hbmc").getPath();
		child.run("trnaconsensus.sh", new String[]{"config="+consensusConfig,
			"in="+new File(directory, "members.fa").getPath(), "out="+fasta, "outmodel="+modelPath},
			null, new File(directory, "consensus"), 1, childHeap);
		final ArrayList<Read> records=readFasta(fasta);
		final Model candidate;
		try(NucleotideHbmCodec.Reader reader=new NucleotideHbmCodec.Reader(modelPath)){
			candidate=reader.next();
			if(reader.next()!=null){throw new IllegalStateException("One member group produced multiple HBMs");}
		}
		if(records.isEmpty() && candidate==null){return null;}
		if(records.size()!=1 || candidate==null){throw new IllegalStateException("Public regeneration must produce exactly one paired model or no candidate");}
		final Read record=records.get(0);
		Tools.toUpperCase(candidate.bases);
		if(!Arrays.equals(record.bases, candidate.bases)){throw new IllegalStateException("Regenerated consensus and HBM sequences differ");}
		int tab=1;
		while(tab<candidate.header.length && candidate.header[tab]!='\t'){tab++;}
		final String originalHeader=new String(candidate.header, 1, tab-1, StandardCharsets.US_ASCII);
		final String generatedName=NcrnaSupport.recordName(record.id);
		if(!NcrnaSupport.recordName(originalHeader).equals(generatedName)){throw new IllegalStateException("Regenerated consensus and HBM names differ");}
		final String header=modelName+record.id.substring(generatedName.length());
		final ByteBuilder hbmHeader=new ByteBuilder().append('>').append(header);
		hbmHeader.append(new String(candidate.header, tab, candidate.header.length-tab, StandardCharsets.US_ASCII));
		return new Candidate(header, new Model(hbmHeader.toBytes(), candidate.bases, candidate.counts));
	}

	/** Loads the base pairs once and checks their fixed ordered model axis. */
	static Model[] loadModels(final Library library, final String hbm){
		final Model[] result=new Model[library.names.length];
		try(NucleotideHbmCodec.Reader reader=new NucleotideHbmCodec.Reader(hbm)){
			for(int model=0; model<result.length; model++){
				final Model row=reader.next();
				if(row==null){throw new IllegalArgumentException("Fewer base HBMs than consensuses");}
				Tools.toUpperCase(row.bases);
				int tab=1;
				while(tab<row.header.length && row.header[tab]!='\t'){tab++;}
				final String name=NcrnaSupport.recordName(new String(row.header, 1, tab-1, StandardCharsets.US_ASCII));
				if(!name.equals(library.names[model]) || !Arrays.equals(row.bases, library.bases[model])){
					throw new IllegalArgumentException("Base consensus/HBM pair differs at model "+model);
				}
				result[model]=row;
			}
			if(reader.next()!=null){throw new IllegalArgumentException("More base HBMs than consensuses");}
		}
		return result;
	}

	/** Uses the existing exhaustive engine; fixed input model order resolves equal-score ties. */
	private Graph scan(final String fasta, final Library library){
		final FileFormat ff=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		final NcrnaSupport.ScanResult result=NcrnaSupport.scanWithStats(ff, library, cutoff, threads, 1);
		alignments+=result.alignments;
		return result.graph;
	}

	/** Converts named query rows to the retained corpus axis, including unsupported genes. */
	private int[] assignments(final Graph graph, final ObjectIntMap<String> index){
		final int[] result=new int[expectedGenes];
		Arrays.fill(result, -2);
		for(int row=0; row<graph.queryNames.length; row++){
			final int gene=index.get(graph.queryNames[row]);
			if(gene<0 || result[gene]!=-2){throw new IllegalArgumentException("Graph query axis differs from the complete corpus");}
			result[gene]=graph.starts[row]==graph.starts[row+1] ? -1 : graph.bestModels[row];
		}
		for(int owner:result){if(owner==-2){throw new IllegalArgumentException("Graph omitted a corpus gene");}}
		return result;
	}

	private void checkAxes(final Graph graph, final Library library, final int genes){
		if(graph.cutoff!=cutoff || graph.queryNames.length!=genes || !Arrays.equals(graph.modelNames, library.names)){
			throw new IllegalArgumentException("Support graph cutoff or named axes differ from the experiment");
		}
	}

	/** Retains native Reads; a failed consumer halts rather than stranding a queued producer. */
	private static ArrayList<Read> readFasta(final String path){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		if(!ff.fasta()){throw new IllegalArgumentException("Experiment sequences must be FASTA: "+path);}
		final Streamer input=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final ArrayList<Read> reads=new ArrayList<Read>();
		try{
			input.start();
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read:batch){
					if(read.mate!=null || read.bases==null || read.bases.length<1){throw new IllegalArgumentException("Expected nonempty unpaired gene sequences");}
					if(!read.validated()){read.validate(true);}
					Tools.toUpperCase(read.bases);
					reads.add(read);
				}
			}
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		if(ReadWrite.closeStream(input)){throw new IllegalStateException("FASTA input failed: "+path);}
		return reads;
	}

	/** Keeps the submitted batch stable until the native writer has finished using it. */
	static void writeFasta(final String path, final ArrayList<Read> records){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.FASTA, null, true, false, false, true);
		final Writer writer=WriterFactory.makeWriter(ff, null, 1, null, false);
		try{
			writer.start();
			writer.add(records, 0);
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		if(writer.poisonAndWait()){throw new IllegalStateException("FASTA output failed: "+path);}
		assert(writer.readsWritten()==records.size()) : "Every declared member must reach the public consensus tool";
	}

	private String path(final String name){return new File(root, name).getPath();}
	private static String name(final Library library, final int model){return model<0 ? "NONE" : library.names[model];}

	/** A renamed candidate whose consensus and native HBM counts share one coordinate frame. */
	private static final class Candidate{
		Candidate(final String header_, final Model model_){
			header=header_;
			model=model_;
		}
		final String header;
		final Model model;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, reference, hbm, graphPath, consensusConfig, installation;
	private int expectedGenes=-1, expectedModels=-1, threads=1, childHeap=1024;
	private float cutoff=Float.NaN;
	private final File root;
	private final NcrnaChildProcess child;
	private int attempted, kept, reverted, empty, noCandidate;
	private long memberVisits, memberMigrations, alignments;
	/** Coordinator-owned state; scans borrow an immutable candidate library. */
	private Library library;
	private Model[] models;
	private String[] headers;
	private ArrayList<Read> genes;
	private ObjectIntMap<String> index;
	private int[] owners, initialOwners;
}
