package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Derives pairwise and blended monomer matrices from explicitly bound primary count tables.
 * Raw LONG/SHORT counts are pooled first. Consensus and pairwise observed frequencies
 * are normalized separately, blended, and only then converted to log odds.
 * @author Brian, Fischl
 */
public final class FamilyMemberPairMatrixExport{

	private FamilyMemberPairMatrixExport(){}

	/** BBTools command entry; no alignment, clustering, calibration or caller execution occurs here. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilyMemberPairMatrixExport.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	/** Derive four candidates and require the original consensus calculation to reproduce POSSUM40. */
	public static void run(String[] args) throws Exception{
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String argument:args){
			final String[] pair=Parse.splitOnFirst(argument, '=');
			require(pair.length==2 && !pair[1].isEmpty(), "Matrix derivation requires explicit nonempty flag=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.contains(key) && options.put(key, pair[1])==null, "Unknown or duplicate argument: "+key);
		}
		require(required(options, "population").equals("span_gated"), "Brian selected LONG70/SHORT28 span-gated counts as primary");
		require(required(options, "t").equals("1"), "Small aggregate matrix derivation uses one worker");
		final Path[] inputs=new Path[INPUT_NAMES.length];
		final String[] pins=new String[INPUT_NAMES.length];
		for(int i=0; i<INPUT_NAMES.length; i++){
			inputs[i]=Paths.get(required(options, INPUT_NAMES[i]));
			pins[i]=required(options, INPUT_NAMES[i]+"sha80");
			verifyInput(inputs[i], pins[i]);
		}
		require(pins[2].equals("9dfc9b7b521f352f55af"), "The consensus blend source must be the accepted historical POSSUM40 count table");
		final Path output=Paths.get(required(options, "out"));
		require(output.isAbsolute() && output.equals(output.normalize()) && !Files.exists(output)
			&& !Files.isSymbolicLink(output) && Files.isDirectory(output.getParent()), "Require a fresh absolute output directory");
		final CountTable longCounts=read(inputs[0]), shortCounts=read(inputs[1]), consensus=read(inputs[2]);
		final CountTable pooled=pool(longCounts, shortCounts);
		final CountTable longUngated=read(inputs[4]), shortUngated=read(inputs[5]);
		final CountTable pooledUngated=pool(longUngated, shortUngated);
		final PairSummary longSummary=readSummary(inputs[6], "long", longUngated, longCounts);
		final PairSummary shortSummary=readSummary(inputs[7], "short", shortUngated, shortCounts);
		final PairSummary pooledSummary=poolSummaries(longSummary, shortSummary);
		final ProteinSubstitutionMatrix baseline=ProteinSubstitutionMatrix.read(inputs[3].toString());
		require(baseline.sha80.equals(ProteinSubstitutionMatrix.CUSTOM40_BINDING_SHA80)
			&& pins[3].equals(ProteinSubstitutionMatrix.CUSTOM40_FILE_SHA80), "Cell differences must use the accepted POSSUM40 binding and bytes");
		final double[] oldFrequency=observed(consensus.count[0]), newFrequency=observed(pooled.count[0]);
		final Derived control=derive(oldFrequency, "custom40");
		require(control.matrix.same(baseline), "Consensus-only normalization must exactly reproduce accepted POSSUM40 scores");
		final ByteBuilder controlBytes=matrixText(control.matrix);
		require(Arrays.equals(controlBytes.toBytes(), Files.readAllBytes(inputs[3])), "Consensus-only control must reproduce the original matrix artifact bytes");
		final Derived[] candidates=new Derived[FRACTIONS.length];
		for(int i=0; i<FRACTIONS.length; i++){
			candidates[i]=derive(blend(oldFrequency, newFrequency, FRACTIONS[i]), NAMES[i]);
		}
		for(int i=0; i<inputs.length; i++){verifyInput(inputs[i], pins[i]);}
		Files.createDirectory(output);
		write(output.resolve("consensus_control.tsv"), controlBytes);
		write(output.resolve("pooled.gated.ordered_counts.tsv"), countText(pooled));
		write(output.resolve("pooled.ungated.ordered_counts.tsv"), countText(pooledUngated));
		final ByteBuilder cells=new ByteBuilder("matrix\tpairwise_fraction\taa1\taa2\tobserved_unordered\texpected_unordered\tlog_odds_bits\tscore\tpossum40_score\tdelta\n");
		final ByteBuilder ranges=new ByteBuilder("matrix\tpairwise_fraction\tminimum\tmaximum\tdiagonal_minimum\tdiagonal_maximum\tscale\ti0\tgap\tsha80\tscore_representative\n");
		for(int k=0; k<candidates.length; k++){
			final Derived candidate=candidates[k];
			final Path path=output.resolve(NAMES[k]+".tsv");
			write(path, matrixText(candidate.matrix));
			require(ProteinSubstitutionMatrix.read(path.toString()).same(candidate.matrix), "Exported matrix must reload with the exact portable binding");
			int minimum=Integer.MAX_VALUE, maximum=Integer.MIN_VALUE, diagonalMin=Integer.MAX_VALUE, diagonalMax=Integer.MIN_VALUE;
			for(int a=0; a<20; a++){for(int b=0; b<20; b++){
				final int score=candidate.matrix.score((byte)a, (byte)b), old=baseline.score((byte)a, (byte)b);
				minimum=Math.min(minimum, score); maximum=Math.max(maximum, score);
				if(a==b){diagonalMin=Math.min(diagonalMin, score);diagonalMax=Math.max(diagonalMax, score);}
				final int cell=20*Math.min(a, b)+Math.max(a, b);
				cells.append(NAMES[k]).tab().append(Double.toString(FRACTIONS[k])).tab().append(letter(a)).tab().append(letter(b)).tab()
					.append(Double.toString(candidate.observed[cell])).tab().append(Double.toString(candidate.expected[cell])).tab()
					.append(Double.toString(candidate.bits[cell])).tab().append(score).tab().append(old).tab().append(score-old).nl();
			}}
			String representative=NAMES[k];
			for(int previous=0; previous<k; previous++){
				if(sameScores(candidate.matrix, candidates[previous].matrix)){representative=NAMES[previous];break;}
			}
			if(sameScores(candidate.matrix, baseline)){representative="POSSUM40";}
			ranges.append(NAMES[k]).tab().append(Double.toString(FRACTIONS[k])).tab().append(minimum).tab().append(maximum).tab()
				.append(diagonalMin).tab().append(diagonalMax).append("\t4\t4\t16\t").append(DigestSuffix.file(path.toString())).tab().append(representative).nl();
		}
		write(output.resolve("cells.tsv"), cells);
		write(output.resolve("ranges.tsv"), ranges);
		final ByteBuilder totals=new ByteBuilder(SUMMARY_HEADER).nl();
		appendTotals(totals, "LONG", longSummary, longUngated, longCounts);
		appendTotals(totals, "SHORT", shortSummary, shortUngated, shortCounts);
		appendTotals(totals, "pooled", pooledSummary, pooledUngated, pooled);
		write(output.resolve("populations.tsv"), totals);
		final ByteBuilder sampling=new ByteBuilder("view");
		for(String key:CENSUS_NAMES){sampling.tab().append(key);}
		sampling.nl();
		appendSampling(sampling, "LONG", longSummary); appendSampling(sampling, "SHORT", shortSummary);
		appendSampling(sampling, "pooled", pooledSummary);
		write(output.resolve("sampling.tsv"), sampling);
		final ByteBuilder differences=new ByteBuilder("band\tremoved_alignments\tremoved_columns\tremoved_eligible_columns\n");
		for(int band=0; band<3; band++){
			differences.append(BANDS[band]).tab().append(pooledUngated.alignments[band]-pooled.alignments[band]).tab()
				.append(pooledSummary.pathColumns[0][band]-pooledSummary.pathColumns[1][band]).tab()
				.append(pooledUngated.columns[band]-pooled.columns[band]).nl();
		}
		write(output.resolve("span_effect.tsv"), differences);
		final ByteBuilder recipe=new ByteBuilder("key\tvalue\npopulation\tspan_gated_LONG70_SHORT28\ntraining_band\t<40\n")
			.append("formula\tround_half_away(8*log2(observed/expected))\nscale\t4\ni0\t4\ngap\t16\n")
			.append("blend\tnormalized_unordered_observed_frequencies_before_marginals_and_log_odds\n")
			.append("weighting\traw_eligible_columns_after_family_cap\nsmoothing\tnone;zero_cells_fail\nclipping\tnone\n")
			.append("nonstandard\tBLOSUM62_X_times4\nevaluation\tnot_performed\n");
		for(int i=0; i<inputs.length; i++){
			recipe.append(INPUT_NAMES[i]).tab().append(inputs[i].toString()).nl()
				.append(INPUT_NAMES[i]+"_sha80").tab().append(pins[i]).nl();
		}
		write(output.resolve("recipe.tsv"), recipe);
		write(output.resolve("PASS"), new ByteBuilder("FAMILY_MEMBER_PAIR_MATRIX_EXPORT_PASS candidates=4 control=exact_possum40\n"));
	}

	/*--------------------------------------------------------------*/
	/*----------------      Counts and Frequencies   ----------------*/
	/*--------------------------------------------------------------*/

	/** Read every band/cell and validate repeated marginals and denominators, including empty bands. */
	public static CountTable read(Path path){
		require(path!=null, "Count input path is required");
		final CountTable table=new CountTable();
		final long[][] query=new long[3][400], reference=new long[3][400];
		final boolean[][] seen=new boolean[3][400];
		final int[] rows=new int[3];
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(Arrays.equals(reader.nextLine(), HEADER.getBytes(StandardCharsets.US_ASCII)), "Expected complete eight-field directed-count schema");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				require(parser.terms()==8, "Directed-count records have exactly eight fields");
				final int band=parser.termEquals("<40", 0) ? 0 : parser.termEquals("40-60", 0) ? 1 : parser.termEquals(">=60", 0) ? 2 : -1;
				final int a=residue(parser, 1), b=residue(parser, 2), cell=20*a+b;
				require(band>=0 && !seen[band][cell], "Every identity-band directed cell must occur exactly once");
				seen[band][cell]=true;
				final long count=parser.parseLong(3), q=parser.parseLong(4), r=parser.parseLong(5);
				final long alignments=parser.parseLong(6), columns=parser.parseLong(7);
				require(count>=0 && q>=0 && r>=0 && alignments>=0 && columns>=0 && (alignments>0 || columns==0),
					"Count populations cannot be negative or contain columns without an alignment");
				if(rows[band]++==0){table.alignments[band]=alignments;table.columns[band]=columns;}
				require(alignments==table.alignments[band] && columns==table.columns[band], "Repeated population denominators must agree");
				table.count[band][cell]=count; query[band][cell]=q; reference[band][cell]=r;
			}
		}finally{require(!reader.close(), "Directed-count input failed");}
		for(int band=0; band<3; band++){
			require(rows[band]==400, "All400 ordered cells are required even for an empty identity band");
			for(int a=0; a<20; a++){
				long q=0, r=0;
				for(int b=0; b<20; b++){q=Math.addExact(q, table.count[band][20*a+b]);r=Math.addExact(r, table.count[band][20*b+a]);}
				for(int b=0; b<20; b++){
					require(query[band][20*a+b]==q && reference[band][20*b+a]==r, "Stored endpoint counts must equal directed row/column sums");
				}
			}
		}
		table.verify();
		return table;
	}

	/** Bind the counter's full-path totals and sampling census to both corresponding directed tables. */
	public static PairSummary readSummary(Path path, String view, CountTable ungated, CountTable gated){
		require(path!=null && ("long".equals(view) || "short".equals(view)), "Summary requires an explicit LONG or SHORT view");
		ungated.verify(); gated.verify();
		final PairSummary summary=new PairSummary();
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(Arrays.equals(reader.nextLine(), SUMMARY_HEADER.getBytes(StandardCharsets.US_ASCII)), "Expected the counter's full-path summary schema");
			for(int gate=0; gate<2; gate++){for(int band=0; band<3; band++){
				final byte[] line=reader.nextLine();require(line!=null, "Summary must contain both gates and all three bands");
				parser.set(line);
				final CountTable table=gate==0 ? ungated : gated;
				require(parser.terms()==6 && parser.termEquals(view, 0)
					&& parser.termEquals(gate==0 ? "ungated" : "span_gated", 1) && parser.termEquals(BANDS[band], 2),
					"Summary view, gate and band must match their canonical count-table coordinates");
				final long alignments=parser.parseLong(3), columns=parser.parseLong(4), eligible=parser.parseLong(5);
				require(alignments==table.alignments[band] && eligible==table.columns[band]
					&& columns>=eligible && columns>=alignments && columns<=Math.multiplyExact(200, alignments),
					"Summary must preserve count-table denominators and nonempty paths of at most200 columns");
				summary.pathColumns[gate][band]=columns;
			}}
			for(int i=0; i<CENSUS_NAMES.length; i++){
				final byte[] line=reader.nextLine();require(line!=null, "Sampling census must retain every counter field");
				parser.set(line);
				require(parser.terms()==2 && parser.termEquals("#"+CENSUS_NAMES[i], 0), "Sampling census fields must retain the counter's canonical order");
				summary.census[i]=parser.parseLong(1);
				require(summary.census[i]>=0, "Sampling census values must be nonnegative");
			}
			require(reader.nextLine()==null, "Unexpected trailing summary record");
		}finally{require(!reader.close(), "Pair summary input failed");}
		long aligned=0;
		for(int band=0; band<3; band++){
			aligned=Math.addExact(aligned, ungated.alignments[band]);
			require(gated.alignments[band]<=ungated.alignments[band]
				&& summary.pathColumns[1][band]<=summary.pathColumns[0][band], "Span gating must retain a subset of the same paths");
			for(int cell=0; cell<400; cell++){
				require(gated.count[band][cell]<=ungated.count[band][cell], "Span gating cannot add a directed observation");
			}
		}
		final long[] c=summary.census;
		require(c[3]==Math.addExact(aligned, c[4]) && c[3]<=c[2], "Selected pairs must equal returned paths plus explicit failed attempts");
		require(Math.addExact(Math.addExact(c[5], c[6]), c[7])<=c[0] && c[1]>=c[7],
			"Capped, empty and singleton families are disjoint subsets of the family census");
		return summary;
	}

	/** Sum both view summaries without reconstructing path length from standard paired residues. */
	public static PairSummary poolSummaries(PairSummary first, PairSummary second){
		require(first!=null && second!=null, "Pooled summary requires both validated views");
		final PairSummary pooled=new PairSummary();
		for(int gate=0; gate<2; gate++){for(int band=0; band<3; band++){
			pooled.pathColumns[gate][band]=Math.addExact(first.pathColumns[gate][band], second.pathColumns[gate][band]);
		}}
		for(int i=0; i<CENSUS_NAMES.length; i++){pooled.census[i]=Math.addExact(first.census[i], second.census[i]);}
		return pooled;
	}

	/** Pool raw view counts before normalization; a view with more columns supplies proportionate mass. */
	public static CountTable pool(CountTable first, CountTable second){
		require(first!=null && second!=null, "Pooling requires both declared view populations");
		first.verify(); second.verify();
		final CountTable out=new CountTable();
		for(int band=0; band<3; band++){
			out.alignments[band]=Math.addExact(first.alignments[band], second.alignments[band]);
			out.columns[band]=Math.addExact(first.columns[band], second.columns[band]);
			for(int cell=0; cell<400; cell++){out.count[band][cell]=Math.addExact(first.count[band][cell], second.count[band][cell]);}
		}
		out.verify();
		return out;
	}

	/** Fold two directed off-diagonals once; a diagonal contributes once to observed mass. */
	public static double[] observed(long[] directed){
		require(directed!=null && directed.length==400, "Monomer frequencies require the complete20-by20 directed matrix");
		long total=0;
		for(long value:directed){require(value>=0, "Directed counts cannot be negative");total=Math.addExact(total, value);}
		require(total>0, "The selected <40 population has no eligible columns; no matrix can be derived");
		final double[] out=new double[400];
		for(int a=0; a<20; a++){for(int b=a; b<20; b++){
			final long count=a==b ? directed[20*a+b] : Math.addExact(directed[20*a+b], directed[20*b+a]);
			out[20*a+b]=count/(double)total;
		}}
		return out;
	}

	/** Mix separately normalized probabilities; unequal source sample sizes must not redefine the blend. */
	public static double[] blend(double[] consensus, double[] pairwise, double fraction){
		validateFrequency(consensus); validateFrequency(pairwise);
		require(Double.isFinite(fraction) && fraction>=0 && fraction<=1, "Pairwise blend fraction must lie in0..1");
		final double[] out=new double[400];
		for(int i=0; i<400; i++){out[i]=(1-fraction)*consensus[i]+fraction*pairwise[i];}
		return out;
	}

	/** Recompute marginals from the selected observed distribution before calculating independence odds. */
	public static Derived derive(double[] frequency, String name){
		validateFrequency(frequency);
		final Derived out=new Derived(frequency);
		final double[] endpoint=new double[20];
		for(int a=0; a<20; a++){for(int b=a; b<20; b++){
			final double half=0.5*frequency[20*a+b];
			endpoint[a]+=half; endpoint[b]+=half;
		}}
		final int[][] scores=new int[20][20];
		for(int a=0; a<20; a++){for(int b=a; b<20; b++){
			final int cell=20*a+b;
			require(frequency[cell]>0 && endpoint[a]>0 && endpoint[b]>0, "Zero probability needs an explicit smoothing decision; no pseudocount is inferred");
			out.expected[cell]=endpoint[a]*endpoint[b]*(a==b ? 1 : 2);
			out.bits[cell]=StrictMath.log(frequency[cell]/out.expected[cell])/StrictMath.log(2);
			scores[a][b]=scores[b][a]=FamilySubstitutionMatrixExport.roundHalfAway(8*out.bits[cell]);
		}}
		out.matrix=ProteinSubstitutionMatrix.fromScores(name, 4, 16, scores);
		return out;
	}

	/** Probability validation rejects incomplete or already symmetrically duplicated tables. */
	private static void validateFrequency(double[] values){
		require(values!=null && values.length==400, "Frequency array must retain the20-by20 upper triangle");
		double total=0;
		for(int a=0; a<20; a++){for(int b=0; b<20; b++){
			final double value=values[20*a+b];
			require(Double.isFinite(value) && value>=0 && (a<=b || value==0), "Observed probabilities must be finite, nonnegative and upper-triangular");
			total+=value;
		}}
		require(Math.abs(total-1)<1e-12, "Observed unordered probabilities must sum to one without implicit renormalization");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Output and State       ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable-by-convention derived evidence; lower-triangle probabilities remain zero. */
	public static final class Derived{
		private Derived(double[] frequency){observed=frequency.clone();}
		public ProteinSubstitutionMatrix matrix;
		public final double[] observed, expected=new double[400], bits=new double[400];
	}

	/** Complete directed counts for three disjoint identity bands. */
	public static final class CountTable{
		public final long[][] count=new long[3][400];
		public final long[] alignments=new long[3], columns=new long[3];
		/** Check mass conservation before any source table participates in pooling. */
		public void verify(){
			for(int band=0; band<3; band++){
				long sum=0;
				for(long value:count[band]){require(value>=0, "Count cells cannot be negative");sum=Math.addExact(sum, value);}
				require(sum==columns[band] && alignments[band]>=0 && (alignments[band]>0 || sum==0), "Population counts must conserve eligible columns");
			}
		}
	}

	/** Counter summary data; path columns include gaps and unknowns excluded from matrix counts. */
	public static final class PairSummary{
		public final long[][] pathColumns=new long[2][3];
		public final long[] census=new long[CENSUS_NAMES.length];
	}

	/** Preserve the historical residue order as well as scores for the byte-exact POSSUM40 control. */
	public static ByteBuilder matrixText(ProteinSubstitutionMatrix matrix){
		require(matrix!=null, "Matrix text requires a validated portable matrix");
		require(matrix.scale==4 && matrix.linearGap==16, "This exporter writes only the declared scale4/gap16 study matrices");
		final ByteBuilder out=new ByteBuilder("#bbgap_substitution_matrix\t1\nname\t").append(matrix.name).append("\nscale\t4\ngap\t16\n");
		for(int i=0; i<20; i++){for(int j=0; j<20; j++){
			final byte a=AminoAcid.acidToNumber[OUTPUT_ALPHABET.charAt(i)], b=AminoAcid.acidToNumber[OUTPUT_ALPHABET.charAt(j)];
			out.append("score\t").append(OUTPUT_ALPHABET.charAt(i)).tab().append(OUTPUT_ALPHABET.charAt(j)).tab().append(matrix.score(a, b)).nl();
		}}
		return out;
	}

	/** Distinct scientific treatments depend on all scores and units, not different artifact names. */
	public static boolean sameScores(ProteinSubstitutionMatrix first, ProteinSubstitutionMatrix second){
		require(first!=null && second!=null, "Matrix comparison requires two validated artifacts");
		if(first.scale!=second.scale || first.linearGap!=second.linearGap){return false;}
		for(byte a=0; a<20; a++){for(byte b=0; b<20; b++){
			if(first.score(a, b)!=second.score(a, b)){return false;}
		}}
		return true;
	}

	/** Serialize pooled counts with exact recomputed directed marginals. */
	public static ByteBuilder countText(CountTable table){
		table.verify();
		final ByteBuilder out=new ByteBuilder(HEADER).nl();
		for(int band=0; band<3; band++){
			final long[] row=new long[20], column=new long[20];
			for(int a=0; a<20; a++){for(int b=0; b<20; b++){row[a]+=table.count[band][20*a+b];column[b]+=table.count[band][20*a+b];}}
			for(int a=0; a<20; a++){for(int b=0; b<20; b++){
				out.append(BANDS[band]).tab().append(letter(a)).tab().append(letter(b)).tab().append(table.count[band][20*a+b]).tab()
					.append(row[a]).tab().append(column[b]).tab().append(table.alignments[band]).tab().append(table.columns[band]).nl();
			}}
		}
		return out;
	}
	/** Keep full-path and eligible-column populations distinct for each gate and view. */
	private static void appendTotals(ByteBuilder out, String name, PairSummary summary, CountTable ungated, CountTable gated){
		for(int gate=0; gate<2; gate++){
			final CountTable counts=gate==0 ? ungated : gated;
			for(int band=0; band<3; band++){
				out.append(name).tab().append(gate==0 ? "ungated" : "span_gated").tab().append(BANDS[band]).tab()
					.append(counts.alignments[band]).tab().append(summary.pathColumns[gate][band]).tab().append(counts.columns[band]).nl();
			}
		}
	}
	/** Publish pooled failure/cap counters alongside their separate LONG and SHORT populations. */
	private static void appendSampling(ByteBuilder out, String name, PairSummary summary){
		out.append(name);
		for(long value:summary.census){out.tab().append(value);}
		out.nl();
	}
	private static char letter(int code){return (char)AminoAcid.numberToAcid[code];}
	private static int residue(LineParser1 parser, int field){
		require(parser.length(field)==1, "A residue column must contain exactly one standard uppercase letter");
		final byte letter=parser.parseByte(field, 0);
		require(letter>='A' && letter<='Z', "Residue columns must use standard uppercase amino acids");
		final int code=AminoAcid.acidToNumber[letter];
		require(code>=0 && code<20, "Unknown/ambiguous residues cannot be matrix coordinates");
		return code;
	}
	private static String required(HashMap<String, String> options, String key){
		final String value=options.get(key);require(value!=null, "Missing option: "+key);return value;
	}
	private static void verifyInput(Path path, String pin){
		DigestSuffix.requireSuffix(pin, "Matrix derivation input");
		require(path.isAbsolute() && path.equals(path.normalize()) && Files.isRegularFile(path) && !Files.isSymbolicLink(path), "Require a regular absolute count/control input");
		require(pin.equals(DigestSuffix.file(path.toString())), "Derivation input differs from its recorded byte identity: "+path);
	}
	private static void write(Path path, ByteBuilder data){
		require(!Files.exists(path) && !Files.isSymbolicLink(path), "Refusing to replace derived matrix evidence");
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);writer.start();
		try{writer.print(data);}finally{require(!writer.poisonAndWait(), "Matrix output failed");}
	}
	private static void require(boolean valid, String message){if(!valid){throw new IllegalArgumentException(message);}}

	private static final String HEADER="band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs";
	private static final String SUMMARY_HEADER="view\tpopulation\tband\talignments\tcolumns\teligible_columns";
	private static final String[] CENSUS_NAMES={"families", "members", "possible_pairs", "selected_pairs", "failed_alignments",
		"cap_hit_families", "empty_families", "singleton_families"};
	private static final String OUTPUT_ALPHABET="ARNDCEQGHILKMFPSTWYV";
	private static final String[] BANDS={"<40", "40-60", ">=60"};
	private static final String[] INPUT_NAMES={"long", "short", "consensus", "baseline", "longungated", "shortungated", "longsummary", "shortsummary"};
	private static final String[] NAMES={"pairwise40", "blend25", "blend50", "blend75"};
	private static final double[] FRACTIONS={1, .25, .5, .75};
	private static final java.util.List<String> OPTIONS=Arrays.asList("long", "longsha80", "short", "shortsha80", "consensus", "consensussha80",
		"baseline", "baselinesha80", "longungated", "longungatedsha80", "shortungated", "shortungatedsha80",
		"longsummary", "longsummarysha80", "shortsummary", "shortsummarysha80", "population", "t", "out");
}
