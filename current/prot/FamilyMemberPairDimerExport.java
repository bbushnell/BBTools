package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Counts actual adjacent substitutions from retained member-pair traces and derives two dimer alphabets.
 * LONG70/SHORT28, identity below40, overlapping adjacency and the existing fixed-odds reducer are explicit.
 * This command does no alignment, family construction or evaluation.
 * @author Brian, Fischl
 */
public final class FamilyMemberPairDimerExport{

	private FamilyMemberPairDimerExport(){}

	/** Normal BBTools entry; input/output failure invalidates the complete derivation. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilyMemberPairDimerExport.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	/** Replay both accepted views, validate the historical partition, then publish counts and two cuts. */
	public static void run(String[] args) throws Exception{
		final HashMap<String, String> options=parse(args);
		require(required(options, "verification").equals("structural") && required(options, "population").equals("span_gated"),
			"Retained stream byte identity is inherited explicitly; LONG70/SHORT28 is the selected population");
		require(required(options, "t").equals("1"), "Retained-stream derivation uses one bounded reader and reducer");
		Shared.setThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final Inputs inputs=new Inputs();
		final Path matrixPath=inputs.file(Paths.get(required(options, "matrixfile")));
		final String matrixPin=required(options, "matrixfilesha80");
		checkPin(matrixPath, matrixPin);
		require(matrixPin.equals(ProteinSubstitutionMatrix.CUSTOM40_FILE_SHA80), "Retained traces must replay against the original POSSUM40 bytes");
		final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.load("custom40", matrixPath.toString());
		require(matrix.scale==4 && matrix.linearGap==16 && matrix.sha80.equals(ProteinSubstitutionMatrix.CUSTOM40_BINDING_SHA80),
			"Replay must preserve the original scale4/gap16 scoring and identity selection");
		final Path historical=inputs.file(Paths.get(required(options, "consensus")));
		final Path baseline=inputs.file(Paths.get(required(options, "baseline")));
		final String historicalPin=required(options, "consensussha80"), baselinePin=required(options, "baselinesha80");
		checkPin(historical, historicalPin); checkPin(baseline, baselinePin);
		final int[] oldMap=readBaseline(baseline);
		final Path output=Paths.get(required(options, "out"));
		require(output.isAbsolute() && output.equals(output.normalize()) && !Files.exists(output) && !Files.isSymbolicLink(output)
			&& output.getParent().equals(output.getParent().toRealPath()), "Output must be a fresh directory beneath a canonical existing parent");
		final Counter longCounts=readView(Paths.get(required(options, "long")), "long", buckets(options, "longbuckets"), matrix, inputs);
		final Counter shortCounts=readView(Paths.get(required(options, "short")), "short", buckets(options, "shortbuckets"), matrix, inputs);
		final Counter pooled=new Counter();pooled.add(longCounts);pooled.add(shortCounts);pooled.verify();
		final long[] consensus=readHistorical(historical);
		final FamilyDimerSubstitutionTable.Table oldTable=new FamilyDimerSubstitutionTable.Table(consensus, DIMERS);
		final FamilyDimerSubstitutionTable.Table pairTable=new FamilyDimerSubstitutionTable.Table(pooled.ordered, DIMERS);
		require(oldTable.total>0 && pairTable.total>0, "Both sources must contain actual adjacent observations before normalization");
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		// Preserve the historical integer-endpoint arithmetic exactly for the control and pure pairwise cut.
		final FamilyDimerSubstitutionTable.Cut control=cut(oldTable.odds, labels);
		require(samePartition(control.map, oldMap), "Historical consensus-only reduction must reproduce the current120 membership groups");
		final double[] oldFrequency=observed(consensus), pairFrequency=observed(pooled.ordered);
		final double[] mixture=blend(oldFrequency, pairFrequency, .5), mixedOdds=odds(mixture);
		final FamilyDimerSubstitutionTable.Cut pairCut=cut(pairTable.odds, labels), mixedCut=cut(mixedOdds, labels);
		inputs.verify();
		checkPin(matrixPath, matrixPin);checkPin(historical, historicalPin);checkPin(baseline, baselinePin);
		Files.createDirectory(output);
		writeOrdered(output.resolve("ordered_counts.tsv"), longCounts, shortCounts, pooled, labels);
		writeFrequencies(output.resolve("frequencies.tsv"), oldFrequency, pairFrequency, mixture, mixedOdds, labels);
		final ByteBuilder summary=new ByteBuilder("view\tband\talignments\tpath_columns\teligible_columns\tvalid_runs\tdimer_observations\n");
		appendSummary(summary, "LONG", longCounts);appendSummary(summary, "SHORT", shortCounts);appendSummary(summary, "pooled", pooled);
		write(output.resolve("summary.tsv"), summary);
		final ByteBuilder comparison=new ByteBuilder("candidate\tshared_groups\tunchanged_dimers\tchanged_dimers\n");
		writeCut(output, "consensus_control", control, oldMap, labels, comparison);
		writeCut(output, "pairwise", pairCut, oldMap, labels, comparison);
		writeCut(output, "blend50", mixedCut, oldMap, labels, comparison);
		write(output.resolve("comparison.tsv"), comparison);
		final ByteBuilder recipe=new ByteBuilder("key\tvalue\npopulation\tspan_gated_LONG70_SHORT28\nidentity\twhole_path_including_gaps_lt40\n")
			.append("adjacency\toverlapping_consecutive_standard_paired_columns;gap_or_X_breaks\n")
			.append("blend\tseparately_normalized_unordered_frequencies;recomputed_endpoint_marginals\n")
			.append("reducer\tFamilyDimerSubstitutionTable.reduce;fixed_original_odds;minimum_updates;exact_partition_ties\n")
			.append("zero_odds\tretained_without_smoothing\ncut\t120\ngeometry\tnonoverlap_K6_pairs_1_3_5\n")
			.append("stream_byte_identity\tinherited_structural_not_rehashed\nmonomer_control\texact_both_gates_and_all_bands\n")
			.append("historical_control\texact_current_partition_membership_not_class_numbers\n")
			.append("bias\tfamily_cap;POSSUM40_alignment_selection;identity_selection;historical_population_and_aligner_difference\n")
			.append("evaluation\tnot_performed\n");
		for(String name:OPTIONS){if(options.containsKey(name)){recipe.append("argument_"+name).tab().append(options.get(name)).nl();}}
		write(output.resolve("recipe.tsv"), recipe);
		write(output.resolve("PASS"), new ByteBuilder("FAMILY_MEMBER_PAIR_DIMER_EXPORT_PASS candidates=2 control=exact_current120\n"));
	}

	/*--------------------------------------------------------------*/
	/*----------------       Validated Streams      ----------------*/
	/*--------------------------------------------------------------*/

	/** The callback receives only replay-validated positive traces; standard-column runs never cross records. */
	public static final class Counter implements FamilyMemberPairSubstitutions.TraceConsumer{
		/** Select primary below40 traces and walk adjacency without retaining the reader's reusable buffers. */
		@Override
		public void accept(byte[] query, byte[] reference, HbmPositionModel.TailResult trace, int band, boolean spanPass){
			require(query!=null && reference!=null && trace!=null && trace.path!=null && band>=0 && band<3,
				"Dimer observation requires a validated positive trace and its original identity band");
			if(band!=0 || !spanPass){return;}
			countSelectedTrace(query,reference,trace);
		}

		/** Count one externally selected validated path; callers own the explicit identity/span filter. */
		public void countSelectedTrace(byte[] query,byte[] reference,HbmPositionModel.TailResult trace){
			require(query!=null && reference!=null && trace!=null && trace.path!=null && trace.path.length>0,
				"Selected adjacency counting requires the replay-validated nonempty trace and sequence arrays");
			int q=trace.queryStart, r=trace.start, previousQ=-1, previousR=-1;
			for(byte op:trace.path){
				if(op=='m'){
					require(q>=0 && q<query.length && r>=0 && r<reference.length, "Paired columns must remain within validated endpoints");
					final int a=query[q++], b=reference[r++];
					if(a>=0 && a<20 && b>=0 && b<20){
						monomer[20*a+b]=Math.addExact(monomer[20*a+b], 1);
						columns=Math.addExact(columns, 1);
						if(previousQ>=0){
							final int cell=DIMERS*(20*previousQ+a)+20*previousR+b;
							ordered[cell]=Math.addExact(ordered[cell], 1);dimers=Math.addExact(dimers, 1);
						}else{runs=Math.addExact(runs, 1);}
						previousQ=a;previousR=b;
					}else{previousQ=previousR=-1;}
				}else if(op=='I'){q++;previousQ=previousR=-1;}
				else if(op=='D'){r++;previousQ=previousR=-1;}
				else{throw new IllegalArgumentException("Unsupported retained trace operation: "+(char)op);}
			}
			require(q==trace.queryEnd+1 && r==trace.end+1, "Dimer counting must consume exactly the validated suffix endpoints");
			alignments=Math.addExact(alignments, 1);pathColumns=Math.addExact(pathColumns, trace.path.length);
		}

		/** Pool raw LONG/SHORT observations before separately normalizing the pairwise source. */
		public void add(Counter other){
			require(other!=null && other!=this, "Pooling requires a distinct view counter");other.verify();
			alignments=Math.addExact(alignments, other.alignments);pathColumns=Math.addExact(pathColumns, other.pathColumns);
			columns=Math.addExact(columns, other.columns);runs=Math.addExact(runs, other.runs);dimers=Math.addExact(dimers, other.dimers);
			for(int i=0; i<ordered.length; i++){ordered[i]=Math.addExact(ordered[i], other.ordered[i]);}
			for(int i=0; i<monomer.length; i++){monomer[i]=Math.addExact(monomer[i], other.monomer[i]);}
		}

		/** Every eligible run of L paired columns supplies L-1 actual adjacent observations. */
		public void verify(){
			require(alignments>=0 && pathColumns>=columns && columns>=runs && runs>=0 && dimers==columns-runs,
				"Overlapping adjacency must conserve eligible_columns-valid_runs, including singleton runs");
			require(sum(ordered)==dimers && sum(monomer)==columns, "Directed dimer and monomer tables must preserve all observed mass");
		}

		public long alignments, pathColumns, columns, runs, dimers;
		public final long[] ordered=new long[DIMERS*DIMERS], monomer=new long[400];
	}

	/** Validate a complete counter output, replay every declared bucket once and compare both monomer populations. */
	private static Counter readView(Path root, String view, int buckets, ProteinSubstitutionMatrix matrix, Inputs inputs) throws Exception{
		require(root.isAbsolute() && Files.isDirectory(root) && root.equals(root.toRealPath()), "Pair output root must be canonical");
		final Path pass=inputs.file(root.resolve("PASS")), recipe=inputs.file(root.resolve("recipe.tsv"));
		require(new String(Files.readAllBytes(pass), StandardCharsets.US_ASCII).equals("FAMILY_MEMBER_PAIR_SUBSTITUTIONS_PASS\n"),
			"Only a completed step1 stream population may supply derived alphabets");
		final HashMap<String, String> settings=readSettings(recipe);
		final int span=view.equals("long") ? 70 : 28;
		require(view.equals(settings.get("view")) && Integer.toString(span).equals(settings.get("span"))
			&& "span_gated".equals(settings.get("primary_population")) && "100".equals(settings.get("pair_cap"))
			&& "lower_member_ordinal_is_query".equals(settings.get("orientation")) && matrix.binding().equals(settings.get("matrix_binding")),
			"Stream recipe must bind the selected view, span, cap, orientation and original matrix");
		final FamilyMemberPairMatrixExport.CountTable ungated=FamilyMemberPairMatrixExport.read(inputs.file(root.resolve("ungated.ordered_counts.tsv")));
		final FamilyMemberPairMatrixExport.CountTable gated=FamilyMemberPairMatrixExport.read(inputs.file(root.resolve("gated.ordered_counts.tsv")));
		final FamilyMemberPairMatrixExport.PairSummary summary=FamilyMemberPairMatrixExport.readSummary(
			inputs.file(root.resolve("summary.tsv")), view, ungated, gated);
		final Counter counter=new Counter();
		final FamilyMemberPairSubstitutions.Counts replay=new FamilyMemberPairSubstitutions.Counts();
		final ByteFile reader=ByteFile.makeByteFile1(inputs.file(root.resolve("streams.tsv")).toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(Arrays.equals(reader.nextLine(), ascii("bucket\tstream")), "Require the complete retained bucket manifest");
			for(int bucket=0; bucket<buckets; bucket++){
				final byte[] line=reader.nextLine();require(line!=null, "Retained stream manifest omitted a declared bucket");parser.set(line);
				final String leaf=String.format(Locale.ROOT, "pairs_%03d.tsv.gz", bucket);
				require(parser.terms()==2 && parser.parseInt(0)==bucket && parser.termEquals(leaf, 1),
					"Every bucket must appear once in canonical order with its original retained-stream filename");
				replay.add(FamilyMemberPairSubstitutions.readStream(inputs.file(root.resolve(leaf)), matrix, span, counter));
			}
			require(reader.nextLine()==null, "Retained stream manifest has more buckets than the bound population");
		}finally{require(!reader.close(), "Retained stream manifest read failed");}
		replay.verify();counter.verify();
		require(replay.selected==summary.census[3] && replay.failed==summary.census[4], "Replayed success/failure outcomes must equal step1 sampling census");
		for(int gate=0; gate<2; gate++){
			final FamilyMemberPairMatrixExport.CountTable table=gate==0 ? ungated : gated;
			for(int band=0; band<3; band++){
				require(Arrays.equals(replay.ordered[gate][band], table.count[band]) && replay.alignments[gate][band]==table.alignments[band]
					&& replay.eligible[gate][band]==table.columns[band] && replay.columns[gate][band]==summary.pathColumns[gate][band],
					"All retained traces must reproduce the exact step1 monomer cells, paths and populations");
			}
		}
		require(counter.alignments==replay.alignments[1][0] && counter.pathColumns==replay.columns[1][0]
			&& counter.columns==replay.eligible[1][0] && Arrays.equals(counter.monomer, replay.ordered[1][0]),
			"Adjacent-dimer walker must select exactly the primary below40 traces and their standard paired columns");
		return counter;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Frequencies and Cuts   ----------------*/
	/*--------------------------------------------------------------*/

	/** Read the complete historical three-band table; retain only the already selected below40 observations. */
	public static long[] readHistorical(Path path){
		require(path!=null, "Historical dimer counts are required");
		final long[] out=new long[DIMERS*DIMERS];
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(Arrays.equals(reader.nextLine(), ascii("band\tquery_dimer\tconsensus_dimer\tcount")), "Historical directed-dimer schema differs");
			for(String band:FamilyDimerSubstitutionTable.BANDS){for(int a=0; a<DIMERS; a++){for(int b=0; b<DIMERS; b++){
				final byte[] line=reader.nextLine();require(line!=null, "Historical dimer table is incomplete");parser.set(line);
				require(parser.terms()==4 && parser.termEquals(band, 0) && parser.termEquals(labels[a], 1) && parser.termEquals(labels[b], 2),
					"Historical cells must occur exactly once in original dictionary/band order");
				final long value=parser.parseLong(3);require(value>=0, "Historical directed counts cannot be negative");
				if(band.equals("<40")){out[DIMERS*a+b]=value;}
			}}}
			require(reader.nextLine()==null, "Unexpected trailing historical dimer cell");
		}finally{require(!reader.close(), "Historical dimer input failed");}
		return out;
	}

	/** Fold directed endpoints once; preserve the internal order of each two-residue symbol. */
	public static double[] observed(long[] directed){
		require(directed!=null && directed.length==DIMERS*DIMERS, "Dimer observations require the complete400-by400 dictionary");
		final long total=sum(directed);require(total>0, "An empty source cannot define a normalized dimer distribution");
		final double[] out=new double[directed.length];
		for(int a=0; a<DIMERS; a++){for(int b=a; b<DIMERS; b++){
			final long value=a==b ? directed[DIMERS*a+b] : Math.addExact(directed[DIMERS*a+b], directed[DIMERS*b+a]);
			out[DIMERS*a+b]=value/(double)total;
		}}
		return out;
	}

	/** Mix source probabilities, not unequal sample counts or already derived odds. */
	public static double[] blend(double[] first, double[] second, double fraction){
		validate(first);validate(second);
		require(Double.isFinite(fraction) && fraction>=0 && fraction<=1, "Blend weight must lie in0..1");
		final double[] out=new double[first.length];
		for(int i=0; i<out.length; i++){out[i]=(1-fraction)*first[i]+fraction*second[i];}
		return out;
	}

	/** Recompute endpoint marginals from mixed observed mass; absent observations retain zero odds. */
	public static double[] odds(double[] frequency){
		validate(frequency);
		final double[] endpoint=new double[DIMERS], out=new double[DIMERS*DIMERS];
		for(int a=0; a<DIMERS; a++){for(int b=a; b<DIMERS; b++){
			final double half=.5*frequency[DIMERS*a+b];endpoint[a]+=half;endpoint[b]+=half;
		}}
		for(int a=0; a<DIMERS; a++){for(int b=a; b<DIMERS; b++){
			final double expected=endpoint[a]*endpoint[b]*(a==b ? 1 : 2);
			final double value=expected==0 ? 0 : frequency[DIMERS*a+b]/expected;
			require(Double.isFinite(value) && value>=0 && (expected>0 || frequency[DIMERS*a+b]==0),
				"Nonnegative observed mass must have valid endpoint expectations; no smoothing is introduced");
			out[DIMERS*a+b]=out[DIMERS*b+a]=value;
		}}
		return out;
	}

	/** Reject double-counted lower triangles and implicit rescaling before any hierarchy is constructed. */
	private static void validate(double[] frequency){
		require(frequency!=null && frequency.length==DIMERS*DIMERS, "Dimer frequencies require the complete upper triangle");
		double total=0;
		for(int a=0; a<DIMERS; a++){for(int b=0; b<DIMERS; b++){
			final double value=frequency[DIMERS*a+b];
			require(Double.isFinite(value) && value>=0 && (a<=b || value==0), "Unordered mass must be finite, nonnegative and upper triangular");
			total+=value;
		}}
		require(Math.abs(total-1)<1e-10, "Separately normalized dimer probabilities must sum to one");
	}

	/** Reuse the adopted complete-link implementation without changing its arithmetic or tie policy. */
	private static FamilyDimerSubstitutionTable.Cut cut(double[] values, String[] labels){
		final FamilyDimerSubstitutionTable.Hierarchy hierarchy=FamilyDimerSubstitutionTable.reduce(values, labels, 120);
		require(hierarchy.cuts.size()==1 && hierarchy.cuts.get(0).symbols==120, "Each derivation must publish exactly its requested120-group partition");
		return hierarchy.cuts.get(0);
	}

	/** Membership-set comparison ignores harmless numeric class renumbering. */
	public static boolean samePartition(int[] first, int[] second){
		require(first!=null && second!=null && first.length==second.length, "Partition comparison requires the same complete dictionary");
		for(int a=0; a<first.length; a++){
			require(first[a]>=0 && second[a]>=0, "Every dictionary symbol must belong to a group");
			for(int b=0; b<a; b++){if((first[a]==first[b])!=(second[a]==second[b])){return false;}}
		}
		return true;
	}

	/** Parse the literal adopted nonoverlapping alphabet; never infer group membership from a class count. */
	public static int[] readBaseline(Path path) throws Exception{
		String value=null;
		for(String line:Files.readAllLines(path, StandardCharsets.US_ASCII)){
			if(line.startsWith("alphabet=")){require(value==null, "Baseline must contain one alphabet definition");value=line.substring(9);}
		}
		require(value!=null && value.startsWith("dimer:1,3,5:"), "Baseline must be the adopted nonoverlapping d120 K6 alphabet");
		final String[] groups=value.substring("dimer:1,3,5:".length()).split("/", -1);
		require(groups.length==120, "Baseline must retain all120 current groups");
		final int[] map=new int[DIMERS];Arrays.fill(map, -1);
		for(int group=0; group<groups.length; group++){
			for(String token:groups[group].split(",", -1)){
				require(token.length()==2, "Every partition member must be a two-residue token");
				final int a=residue(token.charAt(0)), b=residue(token.charAt(1)), id=20*a+b;
				require(map[id]<0, "A dimer cannot appear in two baseline groups");map[id]=group;
			}
		}
		for(int group:map){require(group>=0, "Baseline partition must contain all400 dimers");}
		return map;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Output and Binding    ----------------*/
	/*--------------------------------------------------------------*/

	/** Publish actual directed counts, not a product of monomer frequencies. */
	private static void writeOrdered(Path path, Counter first, Counter second, Counter pooled, String[] labels){
		writeOrdered(path,first,second,pooled,labels,"<40");
	}

	/** Preserve the selected interval in every directed-count row. */
	private static void writeOrdered(Path path,Counter first,Counter second,Counter pooled,String[] labels,String band){
		final ByteStreamWriter writer=writer(path);final ByteBuilder row=new ByteBuilder();
		final String[] names={"LONG", "SHORT", "pooled"};final Counter[] counters={first, second, pooled};
		try{
			writer.println("view\tband\tquery_dimer\treference_dimer\tcount");
			for(int view=0; view<counters.length; view++){for(int a=0; a<DIMERS; a++){for(int b=0; b<DIMERS; b++){
				row.clear().append(names[view]).tab().append(band).tab().append(labels[a]).tab().append(labels[b]).tab().append(counters[view].ordered[DIMERS*a+b]);
				writer.println(row);
			}}}
		}finally{require(!writer.poisonAndWait(), "Directed adjacent-count publication failed");}
	}

	/** Publish the released30–50 interval without reducing an extra alphabet or synthesizing counts. */
	public static void writeBandCounts(Path output,Counter first,Counter second,Counter pooled){
		writeBandCounts(output,first,second,pooled,"30to50");
	}

	/** Publish the selected interval without reducing an extra alphabet or synthesizing counts. */
	public static void writeBandCounts(Path output,Counter first,Counter second,Counter pooled,String band){
		require(Files.isDirectory(output),"Band output must belong to the caller's fresh derivation directory");
		first.verify();second.verify();pooled.verify();
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		writeOrdered(output.resolve("dimer_ordered_counts.tsv.gz"),first,second,pooled,labels,band);
		final ByteBuilder summary=new ByteBuilder("view\tband\talignments\tpath_columns\teligible_columns\tvalid_runs\tdimer_observations\n");
		final String[] names={"LONG","SHORT","pooled"};final Counter[] counters={first,second,pooled};
		for(int i=0;i<3;i++){
			final Counter c=counters[i];
			summary.append(names[i]).tab().append(band).tab().append(c.alignments).tab().append(c.pathColumns).tab()
				.append(c.columns).tab().append(c.runs).tab().append(c.dimers).nl();
		}
		write(output.resolve("dimer_summary.tsv"),summary);
		final FamilyDimerSubstitutionTable.Table table=new FamilyDimerSubstitutionTable.Table(pooled.ordered,DIMERS);
		require(table.total>0,"The selected band must contain actual adjacent observations before defining odds");
		final ByteStreamWriter writer=writer(output.resolve("dimer_symmetric.tsv.gz"));final ByteBuilder row=new ByteBuilder();
		try{
			writer.println("band\tdimer_i\tdimer_j\tcount\tendpoint_i\tendpoint_j\tobserved_frequency\texpected_frequency\todds");
			for(int i=0;i<DIMERS;i++){for(int j=i;j<DIMERS;j++){
				final double pi=table.endpoints[i]/(2.0*table.total),pj=table.endpoints[j]/(2.0*table.total);
				row.clear().append(band).tab().append(labels[i]).tab().append(labels[j]).tab().append(table.counts[i*DIMERS+j]).tab()
					.append(table.endpoints[i]).tab().append(table.endpoints[j]).tab()
					.append(Double.toString(table.counts[i*DIMERS+j]/(double)table.total)).tab()
					.append(Double.toString((i==j ? 1 : 2)*pi*pj)).tab().append(Double.toString(table.odds[i*DIMERS+j]));
				writer.println(row);
			}}
		}finally{require(!writer.poisonAndWait(),"Band dimer odds publication failed");}
	}

	/** Retain enough normalized evidence to independently reproduce the mixture and its odds. */
	private static void writeFrequencies(Path path, double[] first, double[] second, double[] mixed, double[] odds, String[] labels){
		final ByteStreamWriter writer=writer(path);final ByteBuilder row=new ByteBuilder();
		try{
			writer.println("dimer_i\tdimer_j\tconsensus_frequency\tpairwise_frequency\tblend50_frequency\tblend50_odds");
			for(int a=0; a<DIMERS; a++){for(int b=a; b<DIMERS; b++){
				final int cell=DIMERS*a+b;
				row.clear().append(labels[a]).tab().append(labels[b]).tab().append(Double.toString(first[cell])).tab()
					.append(Double.toString(second[cell])).tab().append(Double.toString(mixed[cell])).tab().append(Double.toString(odds[cell]));
				writer.println(row);
			}}
		}finally{require(!writer.poisonAndWait(), "Normalized dimer evidence publication failed");}
	}

	/** Report complete group memberships and exact group/dimer retention, independent of numeric IDs. */
	private static void writeCut(Path output, String name, FamilyDimerSubstitutionTable.Cut cut, int[] baseline,
			String[] labels, ByteBuilder summary){
		final String[] oldGroups=memberships(baseline, labels), newGroups=memberships(cut.map, labels);
		final HashSet<String> shared=new HashSet<String>();int unchanged=0;
		final ByteBuilder rows=new ByteBuilder("dimer_id\tdimer\tclass_id\told_membership\tnew_membership\tunchanged\n");
		for(int i=0; i<DIMERS; i++){
			final boolean same=oldGroups[i].equals(newGroups[i]);
			if(same){unchanged++;shared.add(newGroups[i]);}
			rows.append(i).tab().append(labels[i]).tab().append(cut.map[i]).tab().append(oldGroups[i]).tab().append(newGroups[i]).tab()
				.append(same ? "true" : "false").nl();
		}
		write(output.resolve(name+".map.tsv"), rows);
		write(output.resolve(name+".alphabet.txt"), new ByteBuilder("dimer:1,3,5:").append(String.join("/", cut.groups)).nl());
		summary.append(name).tab().append(shared.size()).tab().append(unchanged).tab().append(DIMERS-unchanged).nl();
	}

	/** Canonical numeric-dictionary membership strings remain stable when classes or tokens are reordered. */
	private static String[] memberships(int[] map, String[] labels){
		require(map.length==labels.length, "Membership labels must cover the complete partition");
		final StringBuilder[] grouped=new StringBuilder[map.length];
		for(int i=0; i<map.length; i++){
			require(map[i]>=0 && map[i]<map.length, "Class index must fit the complete dictionary");
			StringBuilder group=grouped[map[i]];
			if(group==null){group=new StringBuilder();grouped[map[i]]=group;}
			if(group.length()>0){group.append(',');}group.append(labels[i]);
		}
		final String[] text=new String[map.length];
		for(int i=0; i<text.length; i++){if(grouped[i]!=null){text[i]=grouped[i].toString();}}
		final String[] out=new String[map.length];for(int i=0; i<out.length; i++){out[i]=text[map[i]];}return out;
	}

	/** Keep path columns, eligible monomers, run boundaries and actual dimers separately visible. */
	private static void appendSummary(ByteBuilder out, String view, Counter counts){
		counts.verify();out.append(view).append("\t<40\t").append(counts.alignments).tab().append(counts.pathColumns).tab()
			.append(counts.columns).tab().append(counts.runs).tab().append(counts.dimers).nl();
	}

	/** Metadata stability is a structural bracket, not a new claim of byte identity for large streams. */
	private static final class Inputs{
		/** Remember each distinct canonical input before its first read. */
		Path file(Path path) throws Exception{
			require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path) && path.equals(path.toRealPath()),
				"Every derived input must be an explicit canonical regular file: "+path);
			final BasicFileAttributes attributes=Files.readAttributes(path, BasicFileAttributes.class);
			require(before.put(path, attributes)==null, "An input path cannot silently supply two different derivation roles: "+path);
			return path;
		}
		/** Detect changed or replaced inputs before any derived output is published. */
		void verify() throws Exception{
			for(Map.Entry<Path, BasicFileAttributes> entry:before.entrySet()){
				final Path path=entry.getKey();final BasicFileAttributes old=entry.getValue();
				require(!Files.isSymbolicLink(path) && path.equals(path.toRealPath()), "Input canonical location changed during derivation");
				final BasicFileAttributes current=Files.readAttributes(path, BasicFileAttributes.class);
				require(current.isRegularFile() && current.size()==old.size() && current.lastModifiedTime().equals(old.lastModifiedTime())
					&& Objects.equals(current.fileKey(), old.fileKey()), "Input metadata changed during derivation: "+path);
			}
		}
		final HashMap<Path, BasicFileAttributes> before=new HashMap<Path, BasicFileAttributes>();
	}

	/** A control pin is a measured sha80; 'unmeasured' is permitted only for the explicitly structural stream population. */
	private static void checkPin(Path path, String pin){
		DigestSuffix.requireSuffix(pin, "Dimer derivation control");
		require(pin.equals(DigestSuffix.file(path.toString())), "Control input differs from its declared bytes: "+path);
	}
	/** Read the small provenance recipe without tolerating duplicate or missing values. */
	private static HashMap<String, String> readSettings(Path path) throws Exception{
		final HashMap<String, String> settings=new HashMap<String, String>();
		final java.util.List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		require(!lines.isEmpty() && lines.get(0).equals("key\tvalue"), "Counter recipe requires its key/value header");
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t", -1);
			require(fields.length==2 && !fields[1].isEmpty() && settings.put(fields[0], fields[1])==null, "Counter recipe fields must be nonempty and unique");
		}
		return settings;
	}
	/** Reject unrecognized or repeated options before opening inputs. */
	private static HashMap<String, String> parse(String[] args){
		final HashMap<String, String> values=new HashMap<String, String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');require(pair.length==2 && !pair[1].isEmpty(), "Require nonempty flag=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.contains(key) && values.put(key, pair[1])==null, "Unknown or duplicate derivation option: "+key);
		}
		return values;
	}
	private static int buckets(HashMap<String, String> values, String key){
		final int count=Integer.parseInt(required(values, key));require(count>0 && count<=10000, "Bucket count must bind a bounded nonempty stream manifest");return count;
	}
	private static String required(HashMap<String, String> values, String key){
		final String value=values.get(key);require(value!=null, "Missing option: "+key);return value;
	}
	/** Partition coordinates use only the same canonical residue dictionary as the trace counter. */
	private static int residue(char letter){
		require(letter>='A' && letter<='Z', "Dimer tokens require canonical uppercase residues");
		final int code=AminoAcid.acidToNumber[letter];
		require(code>=0 && code<20 && AminoAcid.numberToAcid[code]==letter, "Dimer tokens cannot contain unknown or alternate residues");return code;
	}
	private static long sum(long[] values){
		long total=0;for(long value:values){require(value>=0, "Observation counts cannot be negative or overflowed");total=Math.addExact(total, value);}return total;
	}
	private static byte[] ascii(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static ByteStreamWriter writer(Path path){
		require(!Files.exists(path) && !Files.isSymbolicLink(path), "Refusing to replace derivation evidence");
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);writer.start();return writer;
	}
	private static void write(Path path, ByteBuilder data){
		final ByteStreamWriter writer=writer(path);try{writer.print(data);}finally{require(!writer.poisonAndWait(), "Dimer output failed: "+path);}
	}
	private static void require(boolean valid, String message){if(!valid){throw new IllegalArgumentException(message);}}

	private static final int DIMERS=400;
	private static final java.util.List<String> OPTIONS=Arrays.asList("long", "short", "longbuckets", "shortbuckets", "consensus", "consensussha80",
		"baseline", "baselinesha80", "matrixfile", "matrixfilesha80", "verification", "population", "t", "out");
}
