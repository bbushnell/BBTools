package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import fileIO.ReadWrite;
import map.LongHashSet;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Makes exact-assembly five-fold exclusions for endpoint-table and NN learning.
 * Outer test assemblies remain outside every fold's teaching population. Each
 * fold's ordinary exclusion audit is consumable by NcrnaEndpointSummary, so
 * signed-error histograms are rebuilt from the same genome exclusion union.
 * No sequences, alignments, tables or models are generated here.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaEndpointFoldInputs{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointFoldInputs.class, false);
		try{
			final NcrnaEndpointFoldInputs tool=new NcrnaEndpointFoldInputs(pp.args);
			tool.process();
			pp.outstream.println(tool.mode.equals("prepare") ? "NCRNA_ENDPOINT_FOLD_BINDINGS_COMPLETE" : "NCRNA_ENDPOINT_FOLD_SPLIT_COMPLETE_NO_OUTER_TRAINING");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaEndpointFoldInputs(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){mode=b;}
			else if(a.equals("mapping")){mapping=b;}
			else if(a.equals("fold")){selectedFold=Integer.parseInt(b);}
			else if(a.equals("samples")){samples=Integer.parseInt(b);}
			else if(a.equals("expected")){
				final String[] fields=b.split(",", -1);
				require(fields.length==6, "expected= binds six producer counts: train5,train3,monitor5,monitor3,evaluation5,evaluation3");
				expected=new long[6];
				for(int i=0; i<6; i++){
					expected[i]=Long.parseLong(fields[i]);
					require(expected[i]>0, "Each producer matrix must have a positive accepted row count");
				}
			}
			else if(a.equals("original")){original=b;}
			else if(a.equals("graph")){graph=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("exclude")){exclude=b;}
			else if(a.equals("audit")){audit=b;}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("allassemblies")){allAssemblies=Parse.parseBoolean(b);}
			else{require(parser.parse(arg, a, b), "Unknown fold-binding option: "+arg);}
		}
		manifest=parser.in1;
		output=parser.out1;
		require(parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1
			&& !parser.overwrite && !parser.append, "Bind the entire original population without sampling or overwrite");
		require(manifest!=null && output!=null && !new File(output).exists(), "Require in= and a fresh out= directory");
		if(mode.equals("prepare")){
			require(original!=null && graph!=null && roster!=null && loci!=null && exclude!=null && audit!=null && seed>=0
				&& mapping==null && selectedFold==-1 && samples==-1 && expected==null,
				"Preparation requires original= graph= roster= loci= exclude=outer_loci audit=outer_panel_audit seed=; split flags are unsupported");
			require(Tools.testInputFiles(false, true, manifest, original, graph, roster, loci, exclude, audit), "Unreadable preparation binding");
		}else{
			require(mode.equals("split") && mapping!=null && selectedFold>=0 && selectedFold<5 && samples>0 && expected!=null
				&& original==null && graph==null && roster==null && loci==null && exclude==null && audit==null && seed==-1 && !allAssemblies,
				"mode=split requires in=fold_vector_directory mapping=prepared_loci.tsv fold=0..4 samples= expected=six_producer_counts; no ignored prepare bindings");
			require(new File(manifest).isDirectory() && Tools.testInputFiles(false, true, mapping), "Unreadable fold matrices or locus mapping");
		}
	}

	public void process(){
		if(mode.equals("split")){split();return;}
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, exclude);
		final ObjectMap<String, String> fasta=readManifest();
		final ArrayList<byte[]> outerRows=new ArrayList<byte[]>();
		final ObjectSet<String> outer=readOuter(outerRows, fasta);
		final ObjectIntMap<String> assignments=new ObjectIntMap<String>(String.class);
		final ArrayList<String> training=new ArrayList<String>();
		for(int gene=0; gene<population.ids.length; gene++){
			final String assembly=population.assemblyName(gene);
			require(fasta.contains(assembly), "Original assembly is missing from the frozen corpus manifest: "+assembly);
			final boolean held=outer.contains(assembly);
			require(population.excluded[gene]==held, "Outer audit and locus exclusions disagree: "+population.ids[gene]);
			if(!held && !assignments.contains(assembly)){
				assignments.put(assembly, 0);
				training.add(assembly);
			}
		}
		// Witness support can come from an assembly without a teaching locus.
		// Assign those assemblies too when explicitly requested, so no support
		// row can escape the fold exclusion used by its later real-call labels.
		if(allAssemblies){
			for(final String assembly:fasta.keys()){
				if(assembly!=null && !outer.contains(assembly) && !assignments.contains(assembly)){
					assignments.put(assembly, 0);
					training.add(assembly);
				}
			}
		}
		require(training.size()>=5, "Five genome folds need at least five non-outer assemblies");
		Collections.sort(training);
		Collections.shuffle(training, new Random(seed));
		for(int i=0; i<training.size(); i++){assignments.put(training.get(i), i%5);}
		// Writing order is stable lexical order, independent of graph traversal.
		Collections.sort(training);
		require(new File(output).mkdir(), "Cannot create fresh fold-binding directory");
		if(allAssemblies){writeAssemblies(fasta, outer, assignments);}
		try(Output mapping=new Output(path(output, "loci.tsv")); Output summary=new Output(path(output, "folds.tsv"))){
			mapping.writer.println("locus_id\tsource\tassembly\trole\tfold\tpartial");
			final ByteBuilder row=new ByteBuilder();
			for(int gene=0; gene<population.ids.length; gene++){
				final String assembly=population.assemblyName(gene);
				final boolean held=outer.contains(assembly);
				mapping.writer.print(row.clear().append(population.ids[gene]).tab().append(population.sourceName(gene))
					.tab().append(assembly).tab().append(held ? "OUTER_TEST" : "INTERNAL_FOLD")
					.tab().append(held ? -1 : assignments.get(assembly)).tab().append(population.partial[gene] ? 1 : 0).nl());
			}
			summary.writer.println("fold\tseed\tinternal_assemblies\tfold_held_assemblies\touter_loci\tfold_held_loci\tfold_train_loci\tpartial_loci\texclude\taudit\tfold_paths");
			for(int fold=0; fold<5; fold++){
				writeFold(population, fasta, outerRows, assignments, training, fold, summary);
			}
		}
	}

	/** Stable complete corpus mapping, shared by support learning and real-call
	 * evaluation. The old locus-only mode retains its original bytes and folds. */
	private void writeAssemblies(final ObjectMap<String, String> fasta, final ObjectSet<String> outer,
			final ObjectIntMap<String> assignments){
		assert(allAssemblies && assignments.size()+outer.size()==fasta.size())
			: "Every canonical assembly must have exactly one outer or internal role";
		final ArrayList<String> names=new ArrayList<String>(fasta.size());
		for(final String assembly:fasta.keys()){if(assembly!=null){names.add(assembly);}}
		Collections.sort(names);
		final ByteBuilder row=new ByteBuilder();
		try(Output out=new Output(path(output, "assemblies.tsv"))){
			out.writer.println(ASSEMBLY_HEADER);
			for(final String name:names){
				final boolean held=outer.contains(name);
				final int fold=held ? -1 : assignments.get(name);
				require(held || fold>=0 && fold<5, "Support assemblies must have a known fold before any data is learned");
				out.writer.print(row.clear().append(name).tab().append(fasta.get(name)).tab()
					.append(held ? "OUTER_TEST" : "INTERNAL_FOLD").tab().append(fold).nl());
			}
		}
	}

	/** Checks every training/monitor/evaluation row, then emits only fold-held
	 * evaluation rows. OUTER_TEST rows are counted and never enter gate fitting.
	 */
	private void split(){
		final ObjectMap<String, Gene> genes=readMapping();
		require(new File(output).mkdir(), "Split output directory must be fresh");
		final boolean bgzf=ReadWrite.USE_BGZF;
		final boolean force=ReadWrite.FORCE_BGZIP;
		final boolean zipExplicit=ReadWrite.SET_ZIP_THREADS;
		final int zipThreads=ReadWrite.MAX_ZIP_THREADS();
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		ReadWrite.setZipThreads(1);
		try(Output summary=new Output(path(output, "split.tsv"))){
			summary.writer.println("fold\tinput_split\tend\tinput_rows\tfold_train\tfold_held\touter_test_dropped\tgate_rows");
			for(final String split:new String[]{"train", "monitor", "evaluation"}){
				for(int end=0; end<2; end++){splitMatrix(genes, split, end, summary);}
			}
		}finally{
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
			ReadWrite.setZipThreads(zipThreads);
			ReadWrite.SET_ZIP_THREADS=zipExplicit;
		}
	}

	private ObjectMap<String, Gene> readMapping(){
		final ObjectMap<String, Gene> genes=new ObjectMap<String, Gene>(String.class, Gene.class);
		final ObjectIntMap<String> assemblies=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(mapping)){
			RrnaResourceIO.header(input.reader, "locus_id\tsource\tassembly\trole\tfold\tpartial", mapping);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==6, "Fold mapping requires six explicit fields");
				for(int i=0; i<6; i++){require(p.length(i)>0, "Fold mapping fields cannot be empty");}
				final int fold=p.parseInt(4);
				final boolean outer=p.termEquals("OUTER_TEST", 3);
				require(outer ? fold==-1 : p.termEquals("INTERNAL_FOLD", 3) && fold>=0 && fold<5, "Unknown outer/internal fold role");
				require(p.termEquals("0", 5) || p.termEquals("1", 5), "Partial disposition must be literal0/1");
				final Gene gene=new Gene(p.parseString(0), p.parseString(1), p.parseString(2), fold, p.termEquals("1", 5), genes.size());
				require(!genes.contains(gene.id), "Duplicate original locus in fold mapping");
				final int state=fold+1;
				final int previous=assemblies.get(gene.assembly);
				require(previous<0 || previous==state, "One exact assembly cannot span fold roles");
				assemblies.put(gene.assembly, state);
				genes.put(gene.id, gene);
			}
		}
		require(genes.size()>0, "Fold mapping cannot be empty");
		return genes;
	}

	private void splitMatrix(final ObjectMap<String, Gene> genes, final String split, final int end, final Output summary){
		final String ending=end==0 ? "5prime" : "3prime";
		final String stem="r10."+split+'.'+ending;
		final boolean evaluation=split.equals("evaluation");
		final long expectedRows=expected[2*(evaluation ? 2 : split.equals("monitor") ? 1 : 0)+end];
		final LineParser1 values=new LineParser1('\t');
		final LineParser1 metadata=new LineParser1('\t');
		final LongHashSet seen=new LongHashSet();
		final ByteBuilder buffer=new ByteBuilder();
		long total=0;
		long held=0;
		long outer=0;
		Gene last=null;
		try(Input matrix=new Input(path(manifest, stem+".tsv")); Input rows=new Input(path(manifest, stem+".rows.tsv"));
				Output gate=evaluation ? new Output(path(output, "r10.fold_held."+ending+".tsv.gz")) : null;
				Output gateRows=evaluation ? new Output(path(output, "r10.fold_held."+ending+".rows.tsv")) : null){
			RrnaResourceIO.header(matrix.reader, "#dims\t69\t22\t0", stem);
			RrnaResourceIO.header(rows.reader, NcrnaPositionalFeatureVariants.ROW_HEADER, stem);
			if(evaluation){
				gate.writer.println("#dims\t69\t22\t0");
				gateRows.writer.println(NcrnaPositionalFeatureVariants.ROW_HEADER);
			}
			for(byte[] line=matrix.reader.nextLine(); line!=null; line=matrix.reader.nextLine()){
				final byte[] row=rows.reader.nextLine();
				require(row!=null, "Fold sidecar ended before matrix");
				values.set(line);
				metadata.set(row);
				require(values.terms()==91 && metadata.terms()==11, "Fold matrix/sidecar width changed");
				for(int i=0; i<91; i++){require(values.length(i)>0, "Empty fold matrix field");}
				for(int i=0; i<11; i++){require(metadata.length(i)>0, "Empty fold sidecar field");}
				final Gene gene=last!=null && metadata.termEquals(last.id, 0) ? last : genes.get(metadata.parseString(0));
				require(gene!=null && !gene.partial && metadata.termEquals(gene.source, 1) && metadata.termEquals(gene.assembly, 2)
					&& metadata.termEquals(split, 3), "Fold row does not match its exact original mapping");
				last=gene;
				final int sample=metadata.parseInt(6);
				require(sample>=0 && sample<samples && seen.add(((long)gene.ordinal<<32)|sample), "Invalid or duplicate fold trial");
				final boolean isOuter=gene.fold==-1;
				final boolean isHeld=gene.fold==selectedFold;
				require(evaluation ? isOuter || isHeld : !isOuter && !isHeld, "Outer/fold-held genomes must never teach; fold-training rows must never enter evaluation");
				final int rawStart=metadata.parseInt(7);
				final int rawStop=metadata.parseInt(8);
				final int trueStart=metadata.parseInt(9);
				final int trueStop=metadata.parseInt(10);
				require(rawStart>=0 && rawStop>=rawStart && trueStart>=0 && trueStop>=trueStart && metadata.parseInt(4)>=0,
					"Fold row coordinates/model must remain valid");
				final long correction=end==0 ? (long)trueStart-rawStart : (long)trueStop-rawStop;
				final int target=correction<-10 || correction>10 ? 21 : (int)correction+10;
				for(int i=0; i<22; i++){require(values.termEquals(i==target ? "1" : "0", 69+i), "Fold target and sidecar correction disagree");}
				total++;
				if(isOuter){outer++;}
				if(isHeld){
					held++;
					gate.writer.print(buffer.clear().append(line).nl());
					metadata.setBounds(3);
					final int from=metadata.a();
					final int to=metadata.b();
					gateRows.writer.print(buffer.clear().append(row, 0, from).append("FOLD_HELD").append(row, to, row.length-to).nl());
				}
			}
			require(rows.reader.nextLine()==null && total==expectedRows,
				"Fold matrix and sidecar must retain the complete producer count: "+stem+" expected="+expectedRows+" actual="+total);
			require(!evaluation || held>0 && outer>0 && held+outer==total, "Evaluation must explicitly account for held and outer rows");
		}
		summary.writer.print(buffer.clear().append(selectedFold).tab().append(split).tab().append(ending).tab().append(total)
			.tab().append(evaluation ? 0 : total).tab().append(held).tab().append(outer).tab().append(held).nl());
	}

	private static final class Gene{
		Gene(final String id_, final String source_, final String assembly_, final int fold_, final boolean partial_, final int ordinal_){
			id=id_;
			source=source_;
			assembly=assembly_;
			fold=fold_;
			partial=partial_;
			ordinal=ordinal_;
		}
		final String id, source, assembly;
		final int fold, ordinal;
		final boolean partial;
	}

	private ObjectMap<String, String> readManifest(){
		final ObjectMap<String, String> result=new ObjectMap<String, String>(String.class, String.class);
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(manifest)){
			RrnaResourceIO.header(input.reader, MANIFEST_HEADER, manifest);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==12, "Frozen corpus manifest requires twelve fields");
				final String assembly=p.parseString(0);
				final String sequence=p.parseString(3);
				require(!assembly.isEmpty() && !assembly.equals(".") && !assembly.equals("NA") && !result.contains(assembly)
					&& sequence.startsWith("/") && paths.add(sequence), "Corpus assemblies and canonical FASTA paths must be explicit and unique");
				result.put(assembly, sequence);
				alias(sequence, assembly);
				alias(p.parseString(7), assembly);
			}
		}
		require(result.size()>0, "The frozen corpus manifest must not be empty");
		return result;
	}

	private void alias(final String path, final String assembly){
		final String previous=byPath.get(path);
		require(path.startsWith("/") && (previous==null || previous.equals(assembly)),
			"A frozen FASTA alias must identify exactly one assembly");
		byPath.put(path, assembly);
	}

	/** Preserve the original outer audit rows; require exact mapped assemblies. */
	private ObjectSet<String> readOuter(final ArrayList<byte[]> retained, final ObjectMap<String, String> fasta){
		final ObjectSet<String> result=new ObjectSet<String>(String.class);
		final ObjectIntMap<String> ordinals=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(audit)){
			RrnaResourceIO.header(input.reader, AUDIT_LABEL, audit);
			RrnaResourceIO.header(input.reader, AUDIT_HEADER, audit);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==5, "Outer audit requires five fields");
				final String panel=p.parseString(0);
				final String assembly=p.parseString(3);
				final int previous=ordinals.get(panel);
				require(panel.startsWith("/") && p.termStartsWith("/", 2) && p.parseInt(1)==(previous<0 ? 1 : previous+1),
					"Outer panel rows must be explicit, ordered and dense");
				require(p.termEquals("EXCLUDED_ASSEMBLY", 4) && fasta.contains(assembly),
					"The accepted outer panels have no unmatched aliases; unknown outer assemblies cannot enter folds");
				require(assembly.equals(byPath.get(p.parseString(2))), "Outer audit FASTA path disagrees with the frozen assembly identity");
				ordinals.put(panel, previous<0 ? 1 : previous+1);
				result.add(assembly);
				retained.add(line);
			}
		}
		require(ordinals.size()==2 && result.size()>0, "Bind both existing outer panels before generating folds");
		return result;
	}

	private void writeFold(final NcrnaLearningPopulation population, final ObjectMap<String, String> fasta,
			final ArrayList<byte[]> outerRows, final ObjectIntMap<String> assignments,
			final ArrayList<String> training, final int fold, final Output summary){
		assert(fold>=0 && fold<5) : "Every output describes exactly one of the five assembly folds";
		final String root=path(output, "fold"+fold);
		require(new File(root).mkdir(), "Fold output directory must be fresh");
		final String exclusions=path(root, "exclude.txt");
		final String panelAudit=path(root, "panel_sources.tsv");
		final String panelPaths=path(root, "held_paths.txt");
		long outerLoci=0;
		long heldLoci=0;
		long trainLoci=0;
		long partialLoci=0;
		int assemblies=0;
		try(Output excluded=new Output(exclusions); Output auditOut=new Output(panelAudit); Output paths=new Output(panelPaths)){
			auditOut.writer.println(AUDIT_LABEL);
			auditOut.writer.println(AUDIT_HEADER);
			for(final byte[] row:outerRows){auditOut.writer.print(new ByteBuilder().append(row).nl());}
			final ByteBuilder row=new ByteBuilder();
			for(final String assembly:training){
				if(assignments.get(assembly)!=fold){continue;}
				final String sequence=fasta.get(assembly);
				paths.writer.print(row.clear().append(sequence).nl());
				auditOut.writer.print(row.clear().append(panelPaths).tab().append(++assemblies).tab().append(sequence)
					.tab().append(assembly).append("\tEXCLUDED_ASSEMBLY\n"));
			}
			for(int gene=0; gene<population.ids.length; gene++){
				final boolean outer=population.excluded[gene];
				final boolean held=!outer && assignments.get(population.assemblyName(gene))==fold;
				if(outer || held){excluded.writer.print(row.clear().append(population.ids[gene]).nl());}
				if(population.partial[gene]){partialLoci++;}
				else if(outer){outerLoci++;}
				else if(held){heldLoci++;}
				else{trainLoci++;}
			}
		}
		require(assemblies>0 && outerLoci>0 && heldLoci>0 && trainLoci>0
			&& outerLoci+heldLoci+trainLoci+partialLoci==population.ids.length,
			"Each fold must retain nonempty disjoint outer, held and teaching cohorts and account for all originals");
		summary.writer.print(new ByteBuilder().append(fold).tab().append(seed).tab().append(training.size()).tab().append(assemblies)
			.tab().append(outerLoci).tab().append(heldLoci).tab().append(trainLoci).tab().append(partialLoci)
			.tab().append(exclusions).tab().append(panelAudit).tab().append(panelPaths).nl());
	}

	private static String path(final String directory, final String name){return new File(directory, name).getAbsolutePath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	public static final String AUDIT_LABEL="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result";
	public static final String AUDIT_HEADER="panel\trow1\tfasta\tassembly\tstatus";
	public static final String MANIFEST_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	public static final String ASSEMBLY_HEADER="assembly\tfasta\trole\tfold";
	private String manifest, output, original, graph, roster, loci, exclude, audit;
	private String mode="prepare", mapping;
	private int selectedFold=-1, samples=-1;
	private long[] expected;
	private final ObjectMap<String, String> byPath=new ObjectMap<String, String>(String.class, String.class);
	private long seed=-1;
	private boolean allAssemblies;
}
