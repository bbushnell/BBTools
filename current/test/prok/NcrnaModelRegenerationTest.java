package test.prok;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaBoundaryComparison;
import prok.NcrnaCallerProtection;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaSupport;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec;
import prok.ReduceNcrnaModels;
import shared.Shared;
import structures.ByteBuilder;

/** Finite product-stage controls using the public builder and native support engine.
 * Reflection is confined to this fixture so private stage helpers need no public
 * testing API. Literal allele expectations are independent of stage decisions.
 * The complete CLI case uses one installed consensus as a synthetic genome.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelRegenerationTest{

	/** Runs member rollback/migration and one ordinary caller integration fixture. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaModelRegenerationTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			check(parser.out1!=null, "Fixture requires out=");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			write(new File(root, "consensus.config"), new ByteBuilder("passes=2\nminid=0.3\nendtrim=0\nstochastic=f\n"));
			for(final String mode:new String[]{"keep", "revert", "move", "short"}){serial(root, mode);}
			restoredSupportLoss(root);
			integration(root);
			fallback(root, 1, false);
			fallback(root, 3, false);
			fallback(root, 3, true);
			System.out.println("MODEL_REGENERATION_TEST_PASS keep revert move short restored-support-loss original-graph caller-baseline per-end CLI budget no-progress inconsistent-input fallback-outputs");
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** The ambiguous parent supports two alleles; only a supported alternative permits migration. */
	private static void serial(final File root, final String mode) throws Exception{
		final File work=new File(root, mode);
		check(work.mkdir(), "Each serial fixture needs fresh output");
		final int length=mode.equals("short") ? 40 : 200;
		final byte[] allele=new byte[length], parent=new byte[length], minority=new byte[length], other=new byte[length];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		final Random random=new Random(71);
		for(int i=0; i<length; i++){
			final int base=random.nextInt(4);
			allele[i]=alphabet[base];
			parent[i]=(i&1)==0 ? (byte)'N' : allele[i];
			minority[i]=(i&1)==0 ? alphabet[(base+1)&3] : allele[i];
			other[i]='T';
		}
		if(mode.equals("move")){
			System.arraycopy(minority, 0, other, 0, length);
			for(int i=0; i<length; i++){
				if(i%7==0 || i%7==3){other[i]=alphabet[(AminoAcid.baseToNumber[other[i]]+1)&3];}
			}
		}
		final File fasta=new File(work, "base.fa"), hbm=new File(work, "base.hbmc");
		paired(fasta, hbm, new String[]{"model0", "model1"}, new byte[][]{parent, other});
		final ByteBuilder corpus=new ByteBuilder();
		for(int i=0; i<3; i++){corpus.append(">member").append(i).nl().append(allele).nl();}
		if(mode.equals("revert") || mode.equals("move")){corpus.append(">minority\n").append(minority).nl();}
		final File genesPath=new File(work, "genes.fa");
		write(genesPath, corpus);
		final Library library=NcrnaSupport.loadLibrary(fasta.getPath());
		final Library genes=NcrnaSupport.loadLibrary(genesPath.getPath());
		final Graph graph=NcrnaSupport.scan(FileFormat.testInput(genesPath.getPath(), FileFormat.FASTA,
			null, true, true), library, .7f, 2, 1);
		final int[] owners=new int[genes.names.length];
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int i=0; i<genes.names.length; i++){indices.put(genes.names[i], i);}
		for(int i=0; i<graph.queryNames.length; i++){owners[indices.get(graph.queryNames[i])]=graph.bestModels[i];}
		for(final int owner:owners){check(owner==0, "All fixture genes must initially belong to the ambiguous parent");}
		final Class<?> pairsClass=Class.forName("prok.NcrnaModelPairs");
		final Constructor<?> pairsConstructor=pairsClass.getDeclaredConstructor(Library.class, String.class, boolean[].class);
		pairsConstructor.setAccessible(true);
		final Object pairs=pairsConstructor.newInstance(library, hbm.getPath(), new boolean[]{true, true});
		final Class<?> stageClass=Class.forName("prok.NcrnaModelRegeneration");
		final Constructor<?> stageConstructor=stageClass.getDeclaredConstructor(String.class, long.class, int.class, int.class);
		stageConstructor.setAccessible(true);
		final Object stage=stageConstructor.newInstance(new File(root, "consensus.config").getPath(), 512L, 2, 2);
		field(stage, "root").set(stage, work);
		final Method serial=stageClass.getDeclaredMethod("serial", pairsClass, Library.class, ObjectIntMap.class,
			int[].class, float.class, PrintStream.class);
		serial.setAccessible(true);
		serial.invoke(stage, pairs, genes, indices, owners, .7f, System.err);
		final Method view=pairsClass.getDeclaredMethod("library");
		view.setAccessible(true);
		final Library after=(Library)view.invoke(pairs);
		final boolean accepted=mode.equals("keep") || mode.equals("move");
		check(Arrays.equals(after.bases[0], accepted ? allele : parent), "Literal consensus allele must match the accepted or reverted case");
		check(Arrays.equals(after.bases[1], mode.equals("move") ? minority : other), "Empty model stays unchanged; migrated member rebuilds the later model");
		check((Long)field(stage, "memberMigrations").get(stage)==(mode.equals("move") ? 1L : 0L), "Only the minority with a supporting alternative may migrate");
		if(!accepted){
			for(final int owner:owners){check(owner==0, "Rollback must preserve assignments as well as sequence");}
			final Method write=pairsClass.getDeclaredMethod("write", File.class, String.class);
			write.setAccessible(true);
			final File result=new File(work, "restored");
			check(result.mkdir(), "Restored output must be fresh");
			write.invoke(pairs, result, "models.hbmc");
			compare(hbm, new File(result, "models.hbmc"));
		}
	}

	/** Real rebuilding and scoring on unchanged genes, followed by a controlled
	 * paired restoration. This tests production retention selection independently
	 * of the caller trigger; integration and fallback exercise the actual caller.
	 */
	private static void restoredSupportLoss(final File root) throws Exception{
		final File work=new File(root, "restored_support_loss");
		check(work.mkdir(), "Support-loss fixture requires fresh output");
		final byte[] q=new byte[200];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		final Random random=new Random(71);
		for(int i=0; i<q.length; i++){q[i]=alphabet[random.nextInt(4)];}
		final byte[] p=q.clone(), x=q.clone();
		// Distribute substitutions: contiguous mismatch blocks can align across
		// gaps and remain supported despite their low position-wise identity.
		for(int i=0; i<q.length; i++){
			if(i%20<3){p[i]=alphabet[(AminoAcid.baseToNumber[p[i]]+1)&3];}
			if(i%20>=3 && i%20<6){x[i]=alphabet[(AminoAcid.baseToNumber[x[i]]+1)&3];}
		}
		final byte[] r=x.clone(), b=x.clone();
		for(int i=0; i<q.length; i++){
			if(i%20>=6 && i%20<12){
				r[i]=alphabet[(AminoAcid.baseToNumber[r[i]]+1)&3];
				b[i]='N';
			}
		}
		final File fasta=new File(work, "base.fa"), hbm=new File(work, "base.hbmc");
		paired(fasta, hbm, new String[]{"modelA", "modelB"}, new byte[][]{p, b});
		final Library library=NcrnaSupport.loadLibrary(fasta.getPath());
		final Library genes=new Library(new String[]{"q0", "q1", "q2", "r0", "r1", "r2", "minority"},
			new byte[][]{q, q, q, r, r, r, x});
		final Graph initial=NcrnaSupport.scanWithStats(genes, null, library, .8f, 2).graph;
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int i=0; i<genes.names.length; i++){indices.put(genes.names[i], i);}
		final Class<?> stageClass=Class.forName("prok.NcrnaModelRegeneration");
		final Method assignments=stageClass.getDeclaredMethod("assignments", Graph.class, ObjectIntMap.class);
		assignments.setAccessible(true);
		final int[] before=(int[])assignments.invoke(null, initial, indices);
		NcrnaSupportIO.write(initial, new File(work, "initial.support.tsv").getPath());
		check(Arrays.equals(before, new int[]{0, 0, 0, 1, 1, 1, 1}),
			"All genes start supported; the minority starts in modelB: "+Arrays.toString(before));
		final Class<?> pairsClass=Class.forName("prok.NcrnaModelPairs");
		final Constructor<?> pairsConstructor=pairsClass.getDeclaredConstructor(Library.class, String.class, boolean[].class);
		pairsConstructor.setAccessible(true);
		final Object baseline=pairsConstructor.newInstance(library, hbm.getPath(), new boolean[]{true, true});
		final Method copy=pairsClass.getDeclaredMethod("copy"), view=pairsClass.getDeclaredMethod("library");
		copy.setAccessible(true);
		view.setAccessible(true);
		final Object current=copy.invoke(baseline);
		final Constructor<?> stageConstructor=stageClass.getDeclaredConstructor(String.class, long.class, int.class, int.class);
		stageConstructor.setAccessible(true);
		final Object stage=stageConstructor.newInstance(new File(root, "consensus.config").getPath(), 512L, 2, 2);
		field(stage, "root").set(stage, work);
		final Method serial=stageClass.getDeclaredMethod("serial", pairsClass, Library.class, ObjectIntMap.class,
			int[].class, float.class, PrintStream.class);
		serial.setAccessible(true);
		final int[] owners=before.clone();
		serial.invoke(stage, current, genes, indices, owners, .8f, System.err);
		final Library rebuilt=(Library)view.invoke(current);
		check(Arrays.equals(rebuilt.bases[0], q) && Arrays.equals(rebuilt.bases[1], r),
			"Majority votes must rebuild the two literal alleles");
		final Graph series=NcrnaSupport.scanWithStats(genes, null, rebuilt, .8f, 2).graph;
		final int[] afterSeries=(int[])assignments.invoke(null, series, indices);
		NcrnaSupportIO.write(series, new File(work, "rebuilt.support.tsv").getPath());
		check(Arrays.equals(afterSeries, new int[]{0, 0, 0, 1, 1, 1, 0}),
			"Rebuilding retains all genes and migrates the minority to modelA: "+Arrays.toString(afterSeries));
		final Method restore=pairsClass.getDeclaredMethod("restore", pairsClass, int.class);
		restore.setAccessible(true);
		check((Boolean)restore.invoke(current, baseline, 0), "The controlled restoration must replace a changed pair");
		final Graph candidate=NcrnaSupport.scanWithStats(genes, null, (Library)view.invoke(current), .8f, 2).graph;
		final int[] after=(int[])assignments.invoke(null, candidate, indices);
		NcrnaSupportIO.write(candidate, new File(work, "restored.support.tsv").getPath());
		check(Arrays.equals(after, new int[]{0, 0, 0, 1, 1, 1, -1}),
			"Only the minority loses support after restoring modelA: "+Arrays.toString(after)
			+", best identities="+Arrays.toString(candidate.bestIdentities));
		final ArrayList<Locus> loci=new ArrayList<Locus>();
		loci.add(new Locus("synthetic", 1, 200, (byte)0, "modelA"));
		final NcrnaBoundaryComparison.Result baselineCalls=NcrnaBoundaryComparison.grade(loci, loci, null);
		final NcrnaBoundaryComparison.Result candidateCalls=NcrnaBoundaryComparison.grade(loci, loci, null);
		field(stage, "selectedCalls").set(stage, candidateCalls);
		final Method select=stageClass.getDeclaredMethod("selectResult", pairsClass, pairsClass, Graph.class,
			Graph.class, NcrnaBoundaryComparison.Result.class, boolean.class, int[].class, int[].class);
		select.setAccessible(true);
		final Object retention=select.invoke(stage, baseline, current, initial, candidate, baselineCalls, true, before, after);
		check((Integer)field(retention, "lost").get(retention)==1 && (Integer)field(retention, "gained").get(retention)==0,
			"The measured restoration loses exactly one supported gene and gains none");
		check(field(stage, "pairs").get(stage)==baseline && field(stage, "graph").get(stage)==initial
			&& field(stage, "selectedCalls").get(stage)==baselineCalls,
			"Production selection must retain baseline pairs, assignments and caller results together despite caller acceptance");
		check(field(stage, "status").get(stage).toString().contains("regeneration lost corpus support"),
			"Support loss must be identified as the fallback reason");
	}

	/** One known installed consensus provides intrinsic complete-locus truth, not a biological benchmark. */
	private static void integration(final File root){
		final File work=new File(root, "cli");
		check(work.mkdir(), "CLI fixture requires fresh output");
		final String installed=new FamilyResourceResolver().resolve("sixs_rf01685", "consensus", null).path;
		final Library source=NcrnaSupport.loadLibrary(installed);
		final byte[] gene=source.bases[0];
		final File ref=new File(work, "base.fa"), hbm=new File(work, "base.hbmc"), corpus=new File(work, "genes.fa");
		paired(ref, hbm, new String[]{source.names[0]}, new byte[][]{gene});
		write(corpus, new ByteBuilder(">member0\n").append(gene).append("\n>member1\n").append(gene).nl());
		write(new File(work, "genomes.txt"), new ByteBuilder(ref.getPath()).nl());
		write(new File(work, "truth.tsv"), new ByteBuilder("seqid\tfamily\tstart\tstop\tstrand\n")
			.append(source.names[0]).append("\tfixture\t1\t").append(gene.length).append("\t+\n"));
		write(new File(work, "caller.config"), new ByteBuilder("taxonomy=f\ncallcds=f\nncrnaboundarynet=f\n"));
		final File out=new File(work, "result");
		ReduceNcrnaModels.main(new String[]{"in="+corpus, "ref="+ref, "hbm="+hbm, "cutoff=0.7", "target=1",
			"out="+out, "t=2", "genomes="+new File(work, "genomes.txt"), "truth="+new File(work, "truth.tsv"),
			"family=sixs_rf01685", "truthfamily=fixture", "callerconfig="+new File(work, "caller.config"),
			"callerheap=1g", "regenerate=t", "consensusconfig="+new File(root, "consensus.config"), "consensusheap=512m"});
		final Graph original=NcrnaSupportIO.read(new File(out, "original.support.tsv").getPath());
		final Graph reduced=NcrnaSupportIO.read(new File(out, "members.tsv").getPath());
		check(!original.reduced && reduced.reduced && reduced.queryNames.length==2, "Original and final graphs retain their separate scopes and full corpus");
		check(contains(new File(out, "summary.txt"), "Regeneration passed caller and corpus retention"), "The complete product stage must pass both retention conditions");
		final File calls=new File(out, "regeneration/caller_results.tsv");
		check(lines(calls)==3, "Regeneration must reuse its baseline and make exactly one new successful call");
		check(contains(calls, "exact_5prime\texact_3prime\tboth_exact"), "Per-end counts must be ordinary product output");
	}

	/** Injects literal saved baseline matches while the real caller rejects a homopolymer model.
	 * This is an orchestration control, not evidence that the injected baseline was called.
	 * One allowed new call exercises exhaustion; three allow restoration and no-progress.
	 * The changed-corpus case deliberately violates the immutable-input assumption:
	 * its retained graph describes the earlier input, not safe coverage of the changed input.
	 */
	private static void fallback(final File root, final int limit, final boolean changedCorpus) throws Exception{
		final File work=new File(root, changedCorpus ? "inconsistent_input" : "fallback_"+limit);
		check(work.mkdir(), "Fallback fixture needs a fresh directory");
		final String genome=new File(root, "cli/base.fa").getPath();
		final Library real=NcrnaSupport.loadLibrary(genome);
		final byte[] homopolymer=new byte[real.bases[0].length];
		Arrays.fill(homopolymer, (byte)'T');
		final boolean migrated=limit>1 && !changedCorpus;
		final File ref=migrated ? new File(root, "move/base.fa") : new File(work, "base.fa");
		final File hbm=migrated ? new File(root, "move/base.hbmc") : new File(work, "base.hbmc");
		if(!migrated){paired(ref, hbm, real.names, new byte[][]{changedCorpus ? real.bases[0] : homopolymer});}
		final Library library=NcrnaSupport.loadLibrary(ref.getPath());
		final File corpus=changedCorpus ? new File(work, "changed.fa") : migrated ? new File(root, "move/genes.fa") : ref;
		if(changedCorpus){
			// Deliberate fault injection: change bases after the retained entry snapshot,
			// keeping its gene ID. This forces the final named-support failure branch.
			write(corpus, new ByteBuilder(">").append(real.names[0]).nl().append(homopolymer).nl());
		}
		final boolean[] keptModels=new boolean[library.names.length];
		Arrays.fill(keptModels, true);
		final Graph initial=NcrnaSupport.scan(FileFormat.testInput((changedCorpus ? ref : corpus).getPath(), FileFormat.FASTA,
			null, true, true), library, .7f, 2, 1).subset(keptModels);
		final Constructor<?> callerConstructor=NcrnaCallerProtection.class.getDeclaredConstructor(String.class,
			String.class, String.class, String.class, String.class, int.class, String.class, int.class);
		callerConstructor.setAccessible(true);
		final NcrnaCallerProtection caller=(NcrnaCallerProtection)callerConstructor.newInstance(
			new File(root, "cli/caller.config").getPath(), new File(root, "cli/genomes.txt").getPath(),
			new File(root, "cli/truth.tsv").getPath(), "sixs_rf01685", "fixture", 2, "1g", 3);
		@SuppressWarnings("unchecked")
		final ArrayList<Locus> truth=(ArrayList<Locus>)field(caller, "references").get(caller);
		final ArrayList<Locus> saved=new ArrayList<Locus>();
		final Locus locus=truth.get(0);
		saved.add(new Locus(locus.seqid, locus.start, locus.stop, locus.strand, library.names[0]));
		field(caller, "selectedResult").set(caller, NcrnaBoundaryComparison.grade(truth, saved, null));
		final Class<?> pairsClass=Class.forName("prok.NcrnaModelPairs");
		final Constructor<?> pairsConstructor=pairsClass.getDeclaredConstructor(Library.class, String.class, boolean[].class);
		pairsConstructor.setAccessible(true);
		final Object baseline=pairsConstructor.newInstance(library, hbm.getPath(), keptModels);
		final Class<?> stageClass=Class.forName("prok.NcrnaModelRegeneration");
		final Constructor<?> stageConstructor=stageClass.getDeclaredConstructor(String.class, long.class, int.class, int.class);
		stageConstructor.setAccessible(true);
		final Object stage=stageConstructor.newInstance(new File(root, "consensus.config").getPath(), 512L, limit, 2);
		final Method process=stageClass.getDeclaredMethod("process", String.class, pairsClass, Graph.class,
			NcrnaCallerProtection.class, File.class, PrintStream.class);
		process.setAccessible(true);
		final File result=new File(work, "stage");
		process.invoke(stage, corpus.getPath(), baseline, initial, caller, result, System.err);
		check(field(stage, "pairs").get(stage)==baseline && field(stage, "graph").get(stage)==initial,
			"A failed stage must select the baseline pairs and its exact graph together");
		check(contains(new File(result, "summary.tsv"), changedCorpus ? "regeneration lost corpus support" : limit==1 ? "budget exhausted" : "no progress"),
			"Budget exhaustion, unchanged-winner stop and final corpus loss must be reported distinctly");
		check(contains(new File(result, "exported.caller.tsv"), "0\t1\t0\t0\t1\t1\t1\t0\t0"),
			"Fallback per-end counts must describe the literal baseline, not the rejected zero-call trial");
		if(migrated){
			check(contains(new File(result, "candidate.assignments.tsv"), "minority\tmodel0\tmodel1\tmodel1"),
				"The rejected candidate must actually assign the minority to another model");
			check(contains(new File(result, "final.assignments.tsv"), "minority\tmodel0\tmodel0\tmodel0"),
				"Fallback assignments must restore baseline ownership, not report the rejected migration");
		}else{
			check(contains(new File(result, "final.assignments.tsv"), real.names[0]+"\t"+real.names[0]+"\t"+real.names[0]+"\t"+real.names[0]),
				"Final assignments must name the exported baseline model");
			if(changedCorpus){
				check(contains(new File(result, "candidate.assignments.tsv"), real.names[0]+"\t"+real.names[0]+"\t"+real.names[0]+"\tNONE"),
					"The injected corpus change must produce actual measured support loss before fallback");
			}
		}
		final Method export=pairsClass.getDeclaredMethod("write", File.class, String.class);
		export.setAccessible(true);
		final File selected=new File(work, "selected");
		check(selected.mkdir(), "Selected output requires a fresh directory");
		export.invoke(baseline, selected, "models.hbmc");
		compare(ref, new File(selected, "consensus.fa"));
		compare(hbm, new File(selected, "models.hbmc"));
	}

	private static Field field(final Object object, final String name) throws Exception{
		final Field field=object.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Literal sequence/count pairs; rollback is compared against all native HBM bytes. */
	private static void paired(final File fasta, final File hbm, final String[] names, final byte[][] sequences){
		final ByteBuilder text=new ByteBuilder();
		final ArrayList<BaseGraph> models=new ArrayList<BaseGraph>();
		for(int i=0; i<names.length; i++){
			text.append('>').append(names[i]).nl().append(sequences[i]).nl();
			final BaseGraph graph=new BaseGraph(names[i], sequences[i], null, i, 0);
			BaseGraphHelper.initForScoring(graph);
			models.add(graph);
		}
		write(fasta, text);
		NucleotideHbmCodec.writeModels(models, hbm.getPath(), false);
	}

	private static void write(final File path, final ByteBuilder text){
		try(Output output=new Output(path.getPath())){output.writer.print(text);}
	}

	private static void compare(final File a, final File b){
		try(Input first=new Input(a.getPath()); Input second=new Input(b.getPath())){
			for(byte[] line=first.reader.nextLine(); line!=null; line=first.reader.nextLine()){
				check(Arrays.equals(line, second.reader.nextLine()), "Rollback must preserve every HBM field");
			}
			check(second.reader.nextLine()==null, "Rollback must not append models");
		}
	}

	private static boolean contains(final File path, final String text){
		try(Input input=new Input(path.getPath())){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				if(new String(row, StandardCharsets.US_ASCII).contains(text)){return true;}
			}
		}
		return false;
	}

	private static int lines(final File path){
		int count=0;
		try(Input input=new Input(path.getPath())){
			while(input.reader.nextLine()!=null){count++;}
		}
		return count;
	}

	private static void check(final boolean condition, final String reason){
		if(!condition){throw new AssertionError(reason);}
	}
}
