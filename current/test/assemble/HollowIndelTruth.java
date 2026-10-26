package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.function.Consumer;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.GlocalAligner;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Builds controlled indel inputs and grades entire extensions at their known origins.
 * Scrabble uses the whole query and only an origin-local reference window.
 * @author Brian Bushnell, Fischl */
public final class HollowIndelTruth {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses native BBTools options and runs one benchmark operation. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowIndelTruth.class, false);
		try{new HollowIndelTruth(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires fresh outputs and explicit truth paths for grading. */
	private HollowIndelTruth(final String[] args, final PrintStream log_){
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){mode=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("shreds")){shreds=b;}
			else if(a.equals("details")){details=b;}
			else if(a.equals("hist")){hist=b;}
			else if(a.equals("label")){label=b;}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("deletions")){expectedD=Integer.parseInt(b);}
			else if(a.equals("insertions")){expectedI=Integer.parseInt(b);}
			else if(a.equals("substitutions")){expectedS=Integer.parseInt(b);}
			else if(a.equals("insrate")){insRate=Double.parseDouble(b);}
			else if(a.equals("delrate")){delRate=Double.parseDouble(b);}
			else if(a.equals("maxextension")){maxExtension=Integer.parseInt(b);}
			else if(a.equals("seedchanges")){seedChanges=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		Parser.processQuality();
		input=parser.in1;
		output=parser.out1;
		comparisonInput=parser.in2;
		comparisonOutput=parser.out2;
		require(!parser.overwrite && !parser.append, "Benchmark outputs must be new.");
		require(mode.equals("test") || (input!=null && output!=null && details!=null), "Specify in, out and details.");
		require(expectedD==0 || expectedD==1, "Expected baseline deletion count must be zero or one.");
		require(expectedI>=0 && expectedS>=0 && expectedD+expectedI+expectedS<=1,
				"This panel permits at most one injected seed edit.");
		require(insRate>=0 && delRate>=0 && insRate+delRate<=1, "Invalid read indel rates.");
		require(!seedChanges || mode.equals("audit"), "seedchanges is supported only by audit mode.");
		FastaReadInputStream.SPLIT_READS=false;
		if(mode.equals("compare")){
			require(shreds!=null && comparisonInput!=null && comparisonOutput!=null && maxExtension>0,
					"Compare requires shreds=original.fa in=baseline.fa in2=treatment.fa out=selected_seeds.fa out2=selected_baseline.fa details=changes.tsv.");
			Tools.testInputFiles(true, true, input, comparisonInput, shreds);
			for(String path : new String[]{input, comparisonInput, shreds}){
				require(FileFormat.testInput(path, FileFormat.FASTA, null, true, true).fasta(), "Compare requires FASTA inputs: "+path);
			}
			require(Tools.testOutputFiles(false, false, false, output, comparisonOutput, details), "Compare outputs must be new writable files.");
			Tools.testForDuplicateFiles(true, output, comparisonOutput, details);
			for(String path : new String[]{output, comparisonOutput}){
				require(FileFormat.testOutput(path, FileFormat.FASTA, null, true, false, false, true).fasta(), "Selected comparison outputs must be FASTA: "+path);
			}
		}
	}

	/** Keeps input construction and scoring independent of assembly behavior. */
	private void process(){
		if(mode.equals("test")){selfTest(); return;}
		if(mode.equals("compare")){compareExtensions(); return;}
		if(mode.equals("delete") || mode.equals("insert") || mode.equals("substitute") || mode.equals("readindels")){mutate(); return;}
		if(mode.equals("diagnose")){
			require(reference!=null, "Diagnosis requires the exact reference window.");
			diagnoseAlignment();
			return;
		}
		if(mode.equals("audit")){
			require(reference!=null && shreds!=null && maxExtension>0,
					"Suffix audit requires ref, shreds and a positive maxextension.");
			auditSuffixes();
			return;
		}
		require(mode.equals("grade") && reference!=null && shreds!=null && hist!=null && maxExtension>0,
				"Grade requires ref, shreds, hist and a positive maxextension.");
		grade();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Input generation       ----------------*/
	/*--------------------------------------------------------------*/

	/** Injects one central deletion per shred, or independent 0.1% I and D per read base.
	 * Read evidence starts from the previous Q20 corpus, preserving its existing SNPs. */
	private void mutate(){
		final Random random=new Random(seed);
		final long[] counts=new long[5];
		final ByteBuilder bb=new ByteBuilder();
		try(ReportWriter resource=new ReportWriter(details)){
			final ByteStreamWriter report=resource.writer;
			if(!mode.equals("readindels")){report.print("id\tedit\tzero_based_position\toriginal_base\tnew_base\n");}
			walk(input, output, r->{
				counts[0]++;
				counts[1]+=r.length();
				if(!mode.equals("readindels")){
					require(r.length()==250 && r.quality==null, "Delete mode expects 250bp FASTA shreds.");
					for(byte base : r.bases){require(AminoAcid.isFullyDefined(base), "Seed edit requires defined ACGT bases.");}
					final int position=25+random.nextInt(200);
					final byte old=r.bases[position];
					final byte replacement=mode.equals("delete") ? (byte)'-' : mode.equals("insert") ?
							ALPHABET[random.nextInt(4)] : ALPHABET[(AminoAcid.baseToNumber[old]+1+random.nextInt(3))&3];
					bb.clear().append(r.id).tab().append(mode).tab().append(position).tab().append((char)old)
							.tab().append((char)replacement).nl();
					report.print(bb);
					if(mode.equals("delete")){r.bases=delete(r.bases, position); counts[4]++;}
					else if(mode.equals("insert")){r.bases=insert(r.bases, position, replacement); counts[3]++;}
					else{r.bases[position]=replacement;}
				}else{
					require(r.length()==250 && r.quality!=null, "Evidence must be the original 250bp Q20 FASTQ.");
					final ByteBuilder bases=new ByteBuilder(260), quals=new ByteBuilder(260);
					for(int i=0; i<r.length(); i++){
						require(r.quality[i]==20, "All original qualities must be Q20.");
						final double draw=random.nextDouble();
						if(draw<insRate){
							bases.append(r.bases[i]).append(ALPHABET[random.nextInt(4)]);
							quals.append(r.quality[i]).append((byte)20);
							counts[3]++;
						}else if(draw<insRate+delRate){counts[4]++;}
						else{bases.append(r.bases[i]); quals.append(r.quality[i]);}
					}
					r.bases=bases.toBytes();
					r.quality=quals.toBytes();
				}
				counts[2]+=r.length();
			});
			if(mode.equals("readindels")){
				bb.clear().append("reads\tinput_bases\toutput_bases\tinserted_bases\tdeleted_bases\n");
				for(int i=0; i<counts.length; i++){if(i>0){bb.tab();} bb.append(counts[i]);}
				report.print(bb.nl());
			}
		}
		require(counts[0]>0 && counts[2]==counts[1]+counts[3]-counts[4], "Mutation length accounting failed.");
		log.println("HOLLOW_INDEL_INPUT_PASS reads="+counts[0]+" inserted="+counts[3]+" deleted="+counts[4]);
	}

	/** Removes one base without changing any remaining base or its order. */
	private static byte[] delete(final byte[] bases, final int position){
		assert(position>=0 && position<bases.length) : "Deletion position must be inside the source sequence.";
		final byte[] result=new byte[bases.length-1];
		System.arraycopy(bases, 0, result, 0, position);
		System.arraycopy(bases, position+1, result, position, result.length-position);
		return result;
	}

	/** Inserts one base before the chosen original position without changing the source suffix. */
	private static byte[] insert(final byte[] bases, final int position, final byte base){
		assert(position>=0 && position<bases.length) : "The injected insertion must lie inside the original seed.";
		final byte[] result=new byte[bases.length+1];
		System.arraycopy(bases, 0, result, 0, position);
		result[position]=base;
		System.arraycopy(bases, position, result, position+1, bases.length-position);
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Grading            ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects every treatment truncation, including correct suffixes absent from error reports.
	 * Synchronous ordered inputs avoid whole-library maps and blocked producer cleanup.
	 * Both selected outputs retain original records; no sequence is synthesized or corrected. */
	private void compareExtensions(){
		assert(mode.equals("compare") && shreds!=null && comparisonInput!=null) :
				"Compare-mode validation supplies three ordered, unpaired FASTA libraries.";
		long reads=0, changed=0, removed=0, batchID=0;
		final ByteBuilder row=new ByteBuilder();
		ArrayList<Read> selectedSeeds=new ArrayList<Read>(256), selectedBaseline=new ArrayList<Read>(256);
		try(SequenceCursor seeds=new SequenceCursor(shreds);
				SequenceCursor baseline=new SequenceCursor(input);
				SequenceCursor treatment=new SequenceCursor(comparisonInput);
				ReportWriter report=new ReportWriter(details);
				SequenceOutput seedOutput=new SequenceOutput(output);
				SequenceOutput baselineOutput=new SequenceOutput(comparisonOutput)){
			report.writer.print("id\tseed_length\tbaseline_length\ttreatment_length\tremoved_bases\n");
			for(Read seedRead=seeds.next(); seedRead!=null; seedRead=seeds.next()){
				final Read base=baseline.next(), treated=treatment.next();
				final int lost=compareExtension(seedRead, base, treated, maxExtension);
				reads++;
				if(lost==0){continue;}
				changed++; removed+=lost;
				row.clear().append(seedRead.id).tab().append(seedRead.length()).tab()
						.append(base.length()).tab().append(treated.length()).tab().append(lost).nl();
				report.writer.print(row);
				selectedSeeds.add(seedRead); selectedBaseline.add(base);
				if(selectedSeeds.size()==256){
					seedOutput.writer.add(selectedSeeds, batchID);
					baselineOutput.writer.add(selectedBaseline, batchID++);
					selectedSeeds=new ArrayList<Read>(256); selectedBaseline=new ArrayList<Read>(256);
				}
			}
			require(baseline.next()==null && treatment.next()==null, "Extension libraries contain more records than the original seeds.");
			require(reads>0, "Comparison requires a nonempty original seed library.");
			if(!selectedSeeds.isEmpty()){
				seedOutput.writer.add(selectedSeeds, batchID);
				baselineOutput.writer.add(selectedBaseline, batchID);
			}
		}
		log.println("HOLLOW_EXTENSION_COMPARE_PASS sequences="+reads+" changed="+changed+" removed_bases="+removed);
	}

	/** Validates an ordered pair and returns the exact number of withheld suffix bases. */
	private static int compareExtension(final Read seedRead, final Read baseline, final Read treatment, final int cap){
		require(seedRead!=null && baseline!=null && treatment!=null, "Missing seed or extension record in ordered comparison.");
		require(seedRead.id!=null && seedRead.id.indexOf('\t')<0 && seedRead.id.equals(baseline.id) && seedRead.id.equals(treatment.id),
				"Comparison ID/order mismatch: "+seedRead.id);
		require(seedRead.mate==null && baseline.mate==null && treatment.mate==null && seedRead.length()>0,
				"Compare requires nonempty unpaired seeds and unpaired extensions.");
		require(treatment.length()>=seedRead.length() && baseline.length()>=treatment.length() &&
				(long)baseline.length()-seedRead.length()<=cap, "Treatment grew or extension length bounds failed: "+seedRead.id);
		for(int i=0; i<seedRead.length(); i++){
			require(seedRead.bases[i]==baseline.bases[i], "Baseline changed an original seed: "+seedRead.id);
		}
		for(int i=0; i<treatment.length(); i++){
			require(baseline.bases[i]==treatment.bases[i], "Treatment is not a prefix of the baseline: "+seedRead.id);
		}
		assert(baseline.length()>=treatment.length()) : "The validated prefix relation must yield a nonnegative removed-base count.";
		return baseline.length()-treatment.length();
	}

	/** Reproduces one seed alignment with both Scrabble APIs and a band visualization.
	 * The known one-edit source construction supplies a score independent of Scrabble. */
	private void diagnoseAlignment(){
		final Read[] pair=new Read[2];
		walk(input, null, r->{require(pair[0]==null, "Diagnosis needs exactly one query."); pair[0]=r;});
		walk(reference, null, r->{require(pair[1]==null, "Diagnosis needs exactly one reference window."); pair[1]=r;});
		require(pair[0]!=null && pair[1]!=null, "Missing diagnosis input.");
		final byte[] query=pair[0].bases, ref=pair[1].bases;
		validateSeed(query, ref, 0, expectedS, expectedI, expectedD);
		int queryUnknown=0, refUnknown=0;
		for(byte b : query){if(!AminoAcid.isFullyDefined(b)){queryUnknown++;}}
		for(byte b : ref){if(!AminoAcid.isFullyDefined(b)){refUnknown++;}}
		final AlignmentStats stats=new AlignmentStats(true);
		ScrabbleAligner.alignAndTraceStatic(query, ref, stats);
		final AlignmentStats glocal=new AlignmentStats(true), quantum=new AlignmentStats(true);
		GlocalAligner.alignAndTraceStatic(query, ref, glocal);
		QuantumAligner.alignAndTraceStatic(query, ref, quantum);
		final int[] pos=new int[7];
		require(!new java.io.File(details+".matrix").exists(), "Refuse to overwrite the band visualization.");
		ScrabbleAligner.output=details+".matrix";
		try{ScrabbleAligner.alignStatic(query, ref, pos);}
		finally{ScrabbleAligner.output=null;}
		final int knownScore=250-2*expectedS-expectedI-2*expectedD;
		final ByteBuilder bb=new ByteBuilder();
		bb.append("query_length\treference_length\tquery_unknown\treference_unknown\tknown_path_score")
				.append("\ttrace_score\tplain_score\tM\tS\tI\tD\trstart\trstop")
				.append("\tglocal_score\tglocal_S\tglocal_I\tglocal_D\tquantum_score\tquantum_S\tquantum_I\tquantum_D\n")
				.append(query.length).tab().append(ref.length).tab().append(queryUnknown).tab().append(refUnknown)
				.tab().append(knownScore).tab().append(stats.score).tab().append(pos[2]).tab().append(stats.matches)
				.tab().append(stats.subs).tab().append(stats.ins).tab().append(stats.dels)
				.tab().append(stats.rStart).tab().append(stats.rStop)
				.tab().append(glocal.score).tab().append(glocal.subs).tab().append(glocal.ins).tab().append(glocal.dels)
				.tab().append(quantum.score).tab().append(quantum.subs).tab().append(quantum.ins).tab().append(quantum.dels).nl();
		write(output, bb);
		bb.clear().append(">query\n").append(query).nl().append(">reference_window\n").append(ref).nl()
				.append(">scrabble_operations\n").append(stats.matchString).nl()
				.append(">glocal_operations\n").append(glocal.matchString).nl()
				.append(">quantum_operations\n").append(quantum.matchString).nl();
		write(details, bb);
		log.println("HOLLOW_SCRABBLE_DIAGNOSIS_PASS known_score="+knownScore+" trace_score="+stats.score);
	}

	/** Independently compares only new bases against their fixed source coordinates.
	 * Seed edits are internal, so the first new base always follows original position249.
	 * Wrapped comparisons establish sequence compatibility, not chromosome topology. */
	private void auditSuffixes(){
		final Map<String,Read> refs=new HashMap<String,Read>();
		walk(reference, null, r->{require(refs.put(r.id, r)==null, "Duplicate reference ID.");});
		final Map<String,Seed> seeds=new HashMap<String,Seed>();
		walk(shreds, null, r->{
			final int split=r.id.lastIndexOf('_'), dash=r.id.lastIndexOf('-');
			require(split>0 && dash>split, "Not a Shred coordinate label.");
			final Read ref=refs.get(r.id.substring(0, split));
			final int start=Parse.parseInt(r.id, split+1, dash);
			final int stop=Parse.parseInt(r.id, dash+1, r.id.length());
			require(ref!=null && start>=0 && stop<ref.length() && stop-start==249
					&& r.length()==250-expectedD+expectedI, "Seed coordinates or length invalid.");
			validateSeed(r.bases, ref.bases, start, expectedS, expectedI, expectedD);
			final Seed s=new Seed(r, ref, start, seeds.size());
			require(seeds.put(r.id, s)==null, "Duplicate seed ID.");
		});
		require(!seeds.isEmpty(), "No seeds.");
		final BitSet seen=new BitSet(seeds.size());
		// Inputs, extended, added, linear exact, internal disagreement, boundary
		// sequences, circular-compatible boundaries, pre-boundary disagreeing bases.
		final long[] totals=new long[8];
		// Changed seeds/bases, improved/worsened seeds, residual seed S/I/D, shifted alignments.
		final long[] seedTotals=new long[9];
		final AlignmentStats seedStats=new AlignmentStats(true);
		final int[] counts=new int[4];
		final ByteBuilder row=new ByteBuilder();
		try(ReportWriter resource=new ReportWriter(details)){
			resource.writer.print("id\tadded\tlinear_disagreements\tbeyond_reference\twrapped_disagreements\tfirst_disagreement_offset\tclassification\tobserved_suffix\tcoordinate_suffix"+
					(seedChanges ? "\tseed_changes\tseed_S\tseed_I\tseed_D\tseed_N\tseed_shifted\tnew_seed_mismatches_no_indels\n" : "\n"));
			walk(input, null, r->{
				final Seed s=seeds.get(r.id);
				require(s!=null && !seen.get(s.index), "Unknown or duplicate output ID: "+r.id);
				seen.set(s.index);
				final int added=r.length()-s.read.length();
				require(added>=0 && added<=maxExtension, "Output violates right-extension length limit.");
				int changed=0, newWrong=0;
				for(int i=0; i<s.read.length(); i++){
					if(r.bases[i]==s.read.bases[i]){continue;}
					changed++;
					// Equal total edit counts can hide fixing one SNP while introducing
					// another. Exact coordinates expose that exchange for non-indel seeds.
					if(expectedI+expectedD==0 && r.bases[i]!=s.ref.bases[s.start+i]){newWrong++;}
				}
				require(seedChanges || changed==0, "Original seed was changed.");
				boolean shifted=false;
				if(seedChanges){
					seedStats.clear();
					seedStats.subs=expectedS; seedStats.ins=expectedI; seedStats.dels=expectedD;
					if(changed>0){
						alignChangedSeed(r.bases, s.read.length(), s.ref.bases, s.start, seedStats);
						shifted=seedStats.rStart!=0 || seedStats.rStop!=249;
						seedTotals[0]++; seedTotals[1]+=changed;
						final boolean sameIndels=seedStats.ins==expectedI && seedStats.dels==expectedD;
						if(!shifted && sameIndels && seedStats.ns==0 && seedStats.subs<expectedS){seedTotals[2]++;}
						if(shifted || !sameIndels || seedStats.ns>0 || seedStats.subs>expectedS || newWrong>0){seedTotals[3]++;}
					}
					seedTotals[4]+=seedStats.subs;
					seedTotals[5]+=seedStats.ins;
					seedTotals[6]+=seedStats.dels;
					if(shifted){seedTotals[7]++;}
					seedTotals[8]+=newWrong;
				}
				totals[0]++;
				if(added==0 && changed==0){return;}
				if(added>0){totals[1]++;}
				totals[2]+=added;
				compareSuffix(r.bases, s.read.length(), s.ref.bases, s.start+250, counts);
				final boolean boundary=counts[1]>0;
				final String classification;
				if(boundary){
					totals[5]++;
					if(counts[0]+counts[2]==0){totals[6]++; classification="circular_compatible";}
					else{classification="boundary_disagreement";}
				}else if(counts[0]>0){totals[4]++; classification="internal_disagreement";}
				else{if(added>0){totals[3]++;} classification=added>0 ? "linear_exact" : "unextended";}
				totals[7]+=counts[0];
				if(!boundary && counts[0]==0 && changed==0){return;}
				row.clear().append(r.id).tab().append(added);
				for(int count : counts){row.tab().append(count);}
				row.tab().append(classification).tab();
				for(int i=s.read.length(); i<r.length(); i++){row.append((char)r.bases[i]);}
				row.tab();
				for(int i=0; i<added; i++){row.append((char)s.ref.bases[(s.start+250+i)%s.ref.length()]);}
				if(seedChanges){
					row.tab().append(changed).tab().append(seedStats.subs).tab().append(seedStats.ins)
							.tab().append(seedStats.dels).tab().append(seedStats.ns).tab().append(shifted ? 1 : 0).tab().append(newWrong);
				}
				resource.writer.print(row.nl());
			});
		}
		require(seen.cardinality()==seeds.size(), "Missing output seeds.");
		require(totals[3]+totals[4]+totals[5]==totals[1], "Suffix classification must partition all extended seeds.");
		row.clear().append("arm\tsequences\textended_sequences\tadded_bases\tlinear_exact\tinternal_disagreement")
				.append("\tboundary_sequences\tcircular_compatible\tlinear_disagreeing_bases");
		if(seedChanges){row.append("\tchanged_seeds\tchanged_seed_bases\timproved_seeds\tworsened_seeds\tseed_substitutions\tseed_insertions\tseed_deletions\tshifted_seed_alignments\tnew_seed_mismatches_no_indels");}
		row.nl().append(label);
		for(long total : totals){row.tab().append(total);}
		if(seedChanges){for(long total : seedTotals){row.tab().append(total);}}
		write(output, row.nl());
		log.println("HOLLOW_SUFFIX_AUDIT_PASS sequences="+totals[0]+" extended="+totals[1]);
	}

	/** Aligns only a changed original-length prefix against its exact250bp origin.
	 * Glocal has no heuristic path pruning. Every query base must be accounted
	 * for; shifted reference endpoints are reported as adverse rather than being
	 * hidden by selecting a different genomic locus. Unchanged seeds already
	 * passed validateSeed and do not pay for alignment. */
	private static void alignChangedSeed(final byte[] query, final int seedLength,
			final byte[] reference, final int start, final AlignmentStats stats){
		require(seedLength>=249 && seedLength<=251 && query.length>=seedLength
				&& start>=0 && start+250<=reference.length, "Changed-seed grading requires the frozen250bp origin.");
		stats.clear();
		GlocalAligner.alignAndTraceStatic(Arrays.copyOf(query, seedLength),
				Arrays.copyOfRange(reference, start, start+250), stats);
		require(stats.matches+stats.subs+stats.ins+stats.ns==seedLength,
				"Changed-seed traceback must account for all original-prefix bases; clipping cannot hide a false correction.");
		require(stats.matches+stats.subs+stats.dels+stats.ns==stats.rStop-stats.rStart+1,
				"Changed-seed edit counts must explain the entire consumed origin interval.");
	}

	/** Verifies the exact injected edit independently of any alignment heuristic.
	 * Equal prefix/suffix lengths bound every legal placement of a single indel,
	 * including equivalent placements within a homopolymer. */
	private static void validateSeed(final byte[] seed, final byte[] ref, final int start,
			final int subs, final int ins, final int dels){
		require(start>=0 && start+250<=ref.length && seed.length==250+ins-dels
				&& subs+ins+dels<=1, "Seed validation requires a known 250bp source and at most one edit.");
		for(byte base : seed){require(AminoAcid.isFullyDefined(base), "Seed contains an unknown base.");}
		if(ins+dels==0){
			int differences=0;
			for(int i=0; i<250; i++){
				if(seed[i]!=ref[start+i]){
					require(i>=25 && i<=224, "Injected substitution must be central.");
					differences++;
				}
			}
			require(differences==subs, "Seed does not contain exactly its declared substitution count.");
			return;
		}
		int prefix=0, suffix=0;
		final int minLength=Math.min(250, seed.length);
		while(prefix<minLength && seed[prefix]==ref[start+prefix]){prefix++;}
		while(suffix<minLength && seed[seed.length-1-suffix]==ref[start+249-suffix]){suffix++;}
		final int minimumPosition=Math.max(25, 250-dels-suffix);
		final int maximumPosition=Math.min(224, prefix);
		require(minimumPosition<=maximumPosition, "Seed is not exactly one central insertion/deletion.");
	}

	/** Counts fixed-coordinate disagreements; unknown bases are not evidence of a match.
	 * Result entries are linear disagreement, overrun length, wrapped disagreement,
	 * and first disagreement offset (-1 for none). No edit distance is inferred. */
	private static void compareSuffix(final byte[] query, final int seedLength,
			final byte[] ref, final int sourceEnd, final int[] counts){
		require(seedLength>=0 && seedLength<=query.length && sourceEnd>=0 && sourceEnd<=ref.length
				&& ref.length>0 && counts.length==4, "Suffix coordinates must follow the frozen seed.");
		Arrays.fill(counts, 0);
		counts[3]=-1;
		for(int i=seedLength; i<query.length; i++){
			final int offset=i-seedLength;
			final long position=(long)sourceEnd+offset;
			final boolean beyond=position>=ref.length;
			if(beyond){counts[1]++;}
			final byte truth=ref[(int)(position%ref.length)];
			if(!AminoAcid.isFullyDefined(truth) || !AminoAcid.isFullyDefined(query[i]) || truth!=query[i]){
				counts[beyond ? 2 : 0]++;
				if(counts[3]<0){counts[3]=offset;}
			}
		}
	}

	/** Audits every seed, then measures MSDI for all complete, prefix-preserving outputs. */
	private void grade(){
		final Map<String,Read> refs=new HashMap<String,Read>();
		walk(reference, null, r->{require(refs.put(r.id, r)==null, "Duplicate reference ID.");});
		final Map<String,Seed> seeds=new HashMap<String,Seed>();
		final AlignmentStats stats=new AlignmentStats(true);
		walk(shreds, null, r->{
			final int split=r.id.lastIndexOf('_'), dash=r.id.lastIndexOf('-');
			require(split>0 && dash>split, "Not a Shred coordinate label.");
			final Read ref=refs.get(r.id.substring(0, split));
			final int start=Parse.parseInt(r.id, split+1, dash);
			final int stop=Parse.parseInt(r.id, dash+1, r.id.length());
			require(ref!=null && start>=0 && stop<ref.length() && stop-start==249
					&& r.length()==250-expectedD+expectedI, "Seed coordinates or length invalid.");
			//TODO: Probable bug - Scrabble gives 3S/12I/11D for the known single-I
			//g0600 seed at3156000 (D270). Keep this failure loud; audit mode checks
			//the injected seed and extension directly without trusting that alignment.
			align(r.bases, ref.bases, start, stats);
			require(goodBaseline(stats), "Seed lacks exactly the expected baseline edits: "+r.id+" "+stats);
			final Seed s=new Seed(r, ref, start, seeds.size());
			require(seeds.put(r.id, s)==null, "Duplicate seed ID.");
		});
		require(!seeds.isEmpty(), "No seeds.");
		final BitSet seen=new BitSet(seeds.size());
		final long[] totals=new long[11];
		final long[] substitutions=new long[251+maxExtension], excess=new long[2*(251+maxExtension)];
		final ByteBuilder row=new ByteBuilder();
		try(ReportWriter resource=new ReportWriter(details)){
			final ByteStreamWriter errors=resource.writer;
			errors.print("id\tadded\tM\tS\tD\tI\tN\trstart\trstop\texcess_edits\n");
			walk(input, null, r->{
				final Seed s=seeds.get(r.id);
				require(s!=null && !seen.get(s.index), "Unknown or duplicate output ID: "+r.id);
				seen.set(s.index);
				final int added=r.length()-s.read.length();
				require(added>=0 && added<=maxExtension, "Output violates right-extension length limit.");
				for(int i=0; i<s.read.length(); i++){require(r.bases[i]==s.read.bases[i], "Original seed was changed.");}
				align(r.bases, s.ref.bases, s.start, stats);
				totals[0]++;
				if(added==0){require(goodBaseline(stats), "Unextended baseline changed."); return;}
				final int subDeviation=Math.abs(stats.subs-expectedS);
				final int extra=subDeviation+Math.abs(stats.ins-expectedI)+Math.abs(stats.dels-expectedD)+stats.ns;
				totals[1]++;
				totals[2]+=added;
				totals[3]+=stats.subs;
				totals[4]+=stats.ins;
				totals[5]+=stats.dels;
				totals[6]+=stats.ns;
				if(subDeviation>0){totals[7]++;}
				if(!goodBaseline(stats)){totals[8]++;}
				if(stats.rStart!=0){totals[9]++;}
				totals[10]=Math.max(totals[10], subDeviation);
				require(subDeviation<substitutions.length && extra<excess.length, "Histogram range exceeded.");
				substitutions[subDeviation]++;
				excess[extra]++;
				if(!goodBaseline(stats)){
					row.clear().append(r.id).tab().append(added).tab().append(stats.matches).tab().append(stats.subs)
							.tab().append(stats.dels).tab().append(stats.ins).tab().append(stats.ns)
							.tab().append(stats.rStart).tab().append(stats.rStop).tab().append(extra).nl();
					errors.print(row);
				}
			});
		}
		require(seen.cardinality()==seeds.size(), "Missing output seeds.");
		row.clear().append("arm\tsequences\textended_sequences\tadded_bases\tsubstitutions\tinsertions\tdeletions")
				.append("\tNs\tsequences_with_substitution_deviation\tsequences_with_excess_edits\tshifted_origins\tmax_substitution_deviation")
				.append("\tfraction_with_substitution_deviation\tfraction_with_excess_edits\n").append(label);
		for(long total : totals){row.tab().append(total);}
		row.tab().append(totals[1]==0 ? 0 : totals[7]/(double)totals[1], 8)
				.tab().append(totals[1]==0 ? 0 : totals[8]/(double)totals[1], 8).nl();
		write(output, row);
		row.clear().append("edit_count\textended_with_substitution_deviation\textended_with_excess_edits\n");
		for(int i=0; i<excess.length; i++){
			final long subs=i<substitutions.length ? substitutions[i] : 0;
			if(subs+excess[i]>0){row.append(i).tab().append(subs).tab().append(excess[i]).nl();}
		}
		write(hist, row);
		log.println("HOLLOW_INDEL_GRADE_PASS seeds="+seeds.size()+" extended="+totals[1]);
	}

	/** Uses full-query traceback with no genome-wide best-hit selection.
	 * Extra reference suffix allows true output insertions/deletions without a forced endpoint. */
	private static void align(final byte[] query, final byte[] ref, final int start, final AlignmentStats stats){
		require(start>=0 && start<ref.length, "Origin outside reference.");
		final byte[] window=Arrays.copyOfRange(ref, start, Math.min(ref.length, start+query.length+33));
		stats.clear();
		ScrabbleAligner.alignAndTraceStatic(query, window, stats);
		require(stats.matches+stats.subs+stats.ins+stats.ns==query.length,
				"Scrabble must account for every query base; no clipping may hide extension errors.");
		require(stats.matches+stats.subs+stats.dels+stats.ns==stats.rStop-stats.rStart+1,
				"Trace counts must match the consumed reference interval.");
	}

	/** A shifted origin or any edit beyond the known deletion is not an exact continuation. */
	private static boolean good(final AlignmentStats stats, final int deletions){
		return stats.rStart==0 && stats.subs==0 && stats.ins==0 && stats.dels==deletions && stats.ns==0;
	}

	/** Requires exactly the known baseline edit vector, not merely a low total edit count. */
	private boolean goodBaseline(final AlignmentStats stats){
		return stats.rStart==0 && stats.subs==expectedS && stats.ins==expectedI && stats.dels==expectedD && stats.ns==0;
	}

	/** Checks the scorer with independently constructed edits, including an error near the end. */
	private static void selfTest(){
		final Random random=new Random(429001);
		final byte[] ref=new byte[450];
		for(int i=0; i<ref.length; i++){ref[i]=ALPHABET[random.nextInt(4)];}
		final AlignmentStats stats=new AlignmentStats(true);
		int checks=0;
		for(int position : new int[]{25, 125, 224}){
			final byte[] inserted=insert(Arrays.copyOf(ref, 350), position, (byte)'A');
			validateSeed(insert(Arrays.copyOf(ref, 250), position, (byte)'A'), ref, 0, 0, 1, 0);
			validateSeed(delete(Arrays.copyOf(ref, 250), position), ref, 0, 0, 0, 1);
			align(inserted, ref, 0, stats);
			require(stats.ins==1 && stats.subs==0 && stats.dels==0 && stats.rStart==0, "Insertion fixture not exactly one I.");
			checks++;
			final byte[] substituted=Arrays.copyOf(ref, 350);
			substituted[position]=substituted[position]=='A' ? (byte)'C' : (byte)'A';
			align(substituted, ref, 0, stats);
			require(stats.subs==1 && stats.ins==0 && stats.dels==0 && stats.rStart==0, "Substitution fixture not exactly one S.");
			checks++;
			final byte[] q=delete(Arrays.copyOf(ref, 350), position);
			align(q, ref, 0, stats);
			require(good(stats, 1), "Single-deletion fixture was not exactly one D.");
			checks++;
			q[300]=q[300]=='A' ? (byte)'C' : (byte)'A';
			align(q, ref, 0, stats);
			require(stats.subs==1 && stats.dels==1 && stats.ins==0 && !good(stats, 1), "Added substitution was hidden.");
			checks++;
		}
		align(Arrays.copyOf(ref, 350), ref, 0, stats);
		require(good(stats, 0), "Clean fixture failed.");
		final byte[] twice=delete(delete(Arrays.copyOf(ref, 350), 220), 100);
		align(twice, ref, 0, stats);
		require(stats.dels==2 && !good(stats, 1), "Second deletion was hidden.");
		System.out.println("HOLLOW_INDEL_TRUTH_TEST_PASS checks="+(checks+2));
		final int[] suffix=new int[4];
		compareSuffix(new byte[]{'A', 'C', 'G'}, 1, new byte[]{'A', 'C', 'G'}, 1, suffix);
		require(Arrays.equals(suffix, new int[]{0, 0, 0, -1}), "Exact linear suffix failed.");
		compareSuffix(new byte[]{'G', 'A', 'C'}, 1, new byte[]{'A', 'C', 'G'}, 3, suffix);
		require(Arrays.equals(suffix, new int[]{0, 2, 0, -1}), "Circular-compatible suffix failed.");
		compareSuffix(new byte[]{'G', 'T', 'C'}, 1, new byte[]{'A', 'C', 'G'}, 3, suffix);
		require(Arrays.equals(suffix, new int[]{0, 2, 1, 0}), "Wrong wrapped suffix was hidden.");
		compareSuffix(new byte[]{'A', 'T', 'G'}, 1, new byte[]{'A', 'C', 'G'}, 1, suffix);
		require(Arrays.equals(suffix, new int[]{1, 0, 0, 0}), "Internal suffix error was hidden.");
		compareSuffix(new byte[]{'A', 'N'}, 1, new byte[]{'A', 'N'}, 1, suffix);
		require(Arrays.equals(suffix, new int[]{1, 0, 0, 0}), "Unknown bases cannot validate a suffix.");
		System.out.println("HOLLOW_SUFFIX_TEST_PASS checks=5");
		comparisonTest();
		alignChangedSeed(Arrays.copyOf(ref, 350), 250, ref, 0, stats);
		require(good(stats, 0) && stats.rStop==249, "A repaired substitution must leave a clean origin-anchored seed.");
		for(int position : new int[]{25, 125, 224}){
			final byte[] deletion=delete(Arrays.copyOf(ref, 250), position);
			alignChangedSeed(deletion, deletion.length, ref, 0, stats);
			require(good(stats, 1) && stats.rStop==249, "Changed-seed grading lost the known deletion.");
			final byte[] insertion=insert(Arrays.copyOf(ref, 250), position, (byte)'A');
			alignChangedSeed(insertion, insertion.length, ref, 0, stats);
			require(stats.subs==0 && stats.ins==1 && stats.dels==0 && stats.rStart==0 && stats.rStop==249,
					"Changed-seed grading lost the known insertion.");
		}
		final byte[] wrong=Arrays.copyOf(ref, 250);
		wrong[100]=wrong[100]=='A' ? (byte)'C' : (byte)'A';
		alignChangedSeed(wrong, wrong.length, ref, 0, stats);
		require(stats.subs==1 && stats.ins==0 && stats.dels==0 && stats.rStart==0 && stats.rStop==249,
				"A false seed correction must remain visible in changed-seed grading.");
		System.out.println("HOLLOW_CHANGED_SEED_TEST_PASS checks=8");
	}

	/** Exercises the comparison independently of file I/O and genomic labels. */
	private static void comparisonTest(){
		final Read seed=new Read(new byte[]{'A', 'C'}, null, "seed", 0);
		final Read baseline=new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "seed", 0);
		final Read shorter=new Read(new byte[]{'A', 'C', 'G'}, null, "seed", 0);
		require(compareExtension(seed, baseline, baseline, 2)==0, "Equal extensions must not be selected.");
		require(compareExtension(seed, baseline, shorter, 2)==1, "Partial truncation lost-base accounting failed.");
		require(compareExtension(seed, baseline, seed, 2)==2, "Full suffix truncation must retain the entire seed.");
		final Read wrongID=new Read(new byte[]{'A', 'C'}, null, "other", 0);
		final Read wrongSeed=new Read(new byte[]{'A', 'T'}, null, "seed", 0);
		final Read divergent=new Read(new byte[]{'A', 'C', 'T'}, null, "seed", 0);
		for(Read bad : new Read[]{null, wrongID, wrongSeed, divergent}){
			boolean rejected=false;
			try{compareExtension(seed, baseline, bad, 2);}
			catch(IllegalArgumentException expected){rejected=true;}
			require(rejected, "Missing, reordered, changed-seed or non-prefix comparison was accepted.");
		}
		boolean rejected=false;
		try{compareExtension(seed, shorter, baseline, 2);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Treatment growth must not be silently treated as truncation.");
		rejected=false;
		try{compareExtension(wrongSeed, baseline, shorter, 2);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "A changed original seed must not be hidden by matching output prefixes.");
		rejected=false;
		try{compareExtension(shorter, baseline, seed, 2);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "A treatment shorter than its seed must be rejected.");
		rejected=false;
		try{compareExtension(seed, baseline, shorter, 1);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "The extension cap must be checked independently of the removed length.");
		System.out.println("HOLLOW_EXTENSION_COMPARE_TEST_PASS checks=11");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Native I/O           ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains one native input batch, with no background record producer to strand on failure. */
	private static final class SequenceCursor implements AutoCloseable {
		/** Starts a host-driven single-end stream for an ordered comparison input. */
		SequenceCursor(final String path){
			assert(path!=null) : "Compare-mode input validation supplies a named original or extension library.";
			stream=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, true),
					null, true, -1, false, true, 0);
			stream.start();
		}
		/** Returns each record once and remains at EOF after exhausting the source. */
		Read next(){
			if(finished){return null;}
			while(batch==null || index>=batch.list.size()){
				batch=stream.nextList(); index=0;
				if(batch==null){finished=true; return null;}
			}
			assert(index<batch.list.size()) : "A loaded nonempty batch must contain the next ordered record.";
			return batch.list.get(index++);
		}
		@Override
		public void close(){require(!ReadWrite.closeStream(stream), "Comparison sequence input failed.");}
		private final Streamer stream;
		private ListNum<Read> batch;
		private int index;
		private boolean finished;
	}

	/** Owns a native FASTA writer; only selected batches are submitted with dense IDs. */
	private static final class SequenceOutput implements AutoCloseable {
		/** Uses the normal sequence writer rather than formatting FASTA in a tabular helper. */
		SequenceOutput(final String path){
			assert(path!=null) : "Compare-mode output validation supplies a fresh selected-record destination.";
			writer=WriterFactory.makeWriter(FileFormat.testOutput(path, FileFormat.FASTA, null, true, false, false, true),
					null, 1, null, false);
			writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Comparison sequence output failed.");}
		final Writer writer;
	}

	/** Streams single-end batches; optional sequence output uses the native Writer. */
	private static void walk(final String path, final String out, final Consumer<Read> action){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		final Streamer stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final Writer writer=out==null ? null : WriterFactory.makeWriter(
				FileFormat.testOutput(out, FileFormat.FASTA, null, true, false, false, true), null, 1, null, false);
		stream.start();
		if(writer!=null){writer.start();}
		try{
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read r : batch){require(r.mate==null, "Single-end inputs required."); action.accept(r);}
				if(writer!=null){writer.addReads(batch);}
			}
		}finally{
			final boolean writeError=writer!=null && writer.poisonAndWait();
			final boolean readError=ReadWrite.closeStream(stream);
			require(!writeError && !readError, "Native sequence I/O failed.");
		}
	}

	/** Starts a fresh native tabular output. */
	private static ByteStreamWriter writer(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, true);
		writer.start();
		return writer;
	}

	/** Writes and checks a small report. */
	private static void write(final String path, final ByteBuilder bytes){
		try(ReportWriter resource=new ReportWriter(path)){
			resource.writer.print(bytes);
		}
	}

	/** Always sends the terminal and joins the report thread, including on validation failure.
	 * Try-with-resources preserves the primary exception and suppresses any cleanup failure. */
	private static final class ReportWriter implements AutoCloseable {
		/** Opens and starts the writer owned by this resource. */
		ReportWriter(final String path){writer=HollowIndelTruth.writer(path);}

		/** ByteStreamWriter is a non-daemon worker and must receive poison before scope exit. */
		@Override
		public void close(){require(!writer.poisonAndWait(), "Report write failed.");}

		final ByteStreamWriter writer;
	}

	/** Checks scientific invariants even without JVM assertions. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Keeps immutable origin and input sequence for each output identity. */
	private static final class Seed {
		Seed(final Read read_, final Read ref_, final int start_, final int index_){
			assert(read_!=null && ref_!=null && start_>=0) : "Scoring requires input and reference coordinates.";
			read=read_; ref=ref_; start=start_; index=index_;
		}
		final Read read, ref;
		final int start, index;
	}

	private String mode="grade", input, output, details, hist, reference, shreds;
	private String comparisonInput, comparisonOutput;
	private String label="arm";
	private long seed=429002;
	private int expectedD=0, expectedI=0, expectedS=0, maxExtension=100;
	private double insRate=0.001, delRate=0.001;
	private boolean seedChanges=false;
	private final PrintStream log;
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
