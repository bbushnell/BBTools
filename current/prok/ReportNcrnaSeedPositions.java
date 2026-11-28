package prok;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.LongHashSet;
import map.LongIntMap;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import shared.KillSwitch;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;

/** Reference-gene lengths and empirical seed-edge caps.
 * Length mode streams once and retains primitive lengths only. Position mode
 * streams sense-oriented genes twice and retains distances and reference identities,
 * never the gene sequences. Reference-only windows have no genomic clipping,
 * neighboring-gene hits, claims, length gate, alignments or caller decisions.
 * @author Brian, Diona, Raiden
 */
public final class ReportNcrnaSeedPositions{

	/** Owns the command-line lifetime; fatal consumer errors must not strand native I/O workers. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaSeedPositions.class, false);
		try{new ReportNcrnaSeedPositions(pp.args, pp.outstream).process();}catch(OutOfMemoryError failure){
			KillSwitch.memKill(failure);
			throw failure;
		}catch(RuntimeException | AssertionError failure){
			failure.printStackTrace(pp.outstream);
			KillSwitch.kill("Seed-position report failed; an input worker must not remain blocked after consumer failure");
			throw failure;
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepts a complete FASTA cohort and fresh output prefix; position mode also needs seeds. */
	private ReportNcrnaSeedPositions(final String[] args, final PrintStream log_){
		log=log_;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		String orientation=null;
		boolean percentileSpecified=false, paddingSpecified=false, populationSpecified=false, metadataSpecified=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("seeds")){
				seedPath=b;
			}else if(a.equals("mode")){
				require("positions".equals(b) || "lengths".equals(b) || "bounds".equals(b), "mode must be positions, lengths or bounds");
				lengthsOnly="lengths".equals(b);
				boundsOnly="bounds".equals(b);
			}else if(a.equals("member")){
				member=b;
			}else if(a.equals("population")){
				require("original_genes".equals(b) || "scored_reference_intervals".equals(b),
					"population must explicitly name original_genes or scored_reference_intervals");
				constructionPopulation="original_genes".equals(b) ? NcrnaOverlapWindowBuilder.POPULATION_ORIGINALS
					: NcrnaOverlapWindowBuilder.POPULATION_SCORED;
				populationSpecified=true;
			}else if(a.equals("metadata")){
				writeMetadata=Parse.parseBoolean(b);
				metadataSpecified=true;
			}else if(a.equals("boundpercentile")){
				boundPercentile=Double.parseDouble(b);
				percentileSpecified=true;
			}else if(a.equals("boundpadding")){
				boundPadding=Parse.parseIntKMG(b);
				paddingSpecified=true;
			}else if(a.equals("boundcombine")){
				require("union".equals(b) || "intersection".equals(b), "boundcombine must be union or intersection");
				boundCombination="union".equals(b) ? NcrnaOverlapWindowBuilder.UNION : NcrnaOverlapWindowBuilder.INTERSECTION;
			}else if(a.equals("expectedgenes")){
				expectedGenes=Parse.parseIntKMG(b);
			}else if(a.equals("expectedkeys")){
				expectedKeys=Parse.parseIntKMG(b);
			}else if(a.equals("orientation")){
				orientation=b;
			}else{
				require(parser.parse(arg, a, b), "Unknown seed-position argument: "+arg);
			}
		}
		input=parser.in1;
		out=parser.out1;
		require(!percentileSpecified || boundsOnly, "boundpercentile applies only to mode=bounds, never historical position or length diagnostics");
		require(!paddingSpecified || boundsOnly, "boundpadding applies only to mode=bounds, never historical position or length diagnostics");
		NcrnaOverlapWindowBuilder.capRule(boundPercentile, boundPadding);
		require(!populationSpecified || boundsOnly, "population is declared only by mode=bounds; it is not inferred from reference names");
		require(!metadataSpecified || boundsOnly, "metadata applies only to mode=bounds");
		require(boundCombination==null || boundsOnly && writeMetadata, "boundcombine requires bounds mode and metadata=t; legacy tables keep their existing external policy");
		require(!populationSpecified || writeMetadata, "An explicit construction population requires metadata=t; legacy bounds cannot record it");
		require(input!=null && out!=null && member!=null && !member.isEmpty(),
			"Supply one member, its complete reference-gene FASTA and a fresh output prefix");
		if(lengthsOnly){
			require(seedPath==null && expectedKeys==-1 && (expectedGenes==-1 || expectedGenes>0),
				"Length-only mode measures complete sequence records; omit seed options and optionally supply expectedgenes>0");
		}else{
			require(seedPath!=null && expectedGenes>0 && expectedKeys>0 && "sense".equals(orientation),
				"Position mode requires sense-oriented genes, seeds and expected genes/keys");
		}
		require(parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1
			&& !parser.overwrite && !parser.append && Shared.threads()==1,
			"Exact reference reporting uses t=1, complete unpaired genes and fresh outputs; no sampling or overwriting");
		Shared.TRIM_READ_DESCRIPTION=false;
		if(lengthsOnly){
			RrnaResourceIO.fresh(out, ".lengths.tsv");
		}else if(boundsOnly){
			RrnaResourceIO.fresh(out, ".bounds.tsv", ".lengths.tsv");
		}else{
			RrnaResourceIO.fresh(out, ".keys.tsv", ".lengths.tsv", ".windows.tsv", ".window_histogram.tsv", ".references.tsv");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Measures lengths alone, or builds occurrence-weighted caps and reference-only windows. */
	private void process(){
		if(lengthsOnly){
			measureLengths();
			return;
		}
		final LongHashSet seeds=ProkObject.loadLongKmers(seedPath, K);
		require(seeds!=null && seeds.size()==expectedKeys, "Loaded forward seed set must equal its declared complete axis");
		keys=seeds.toArray();
		Arrays.sort(keys);
		index=new LongIntMap(keys.length*2);
		stats=new Key[keys.length];
		for(int i=0; i<keys.length; i++){
			index.put(keys[i], i);
			stats[i]=new Key();
		}
		scratch=new Scratch(keys.length);
		firstPass();
		lengths.sort();
		blanket=Math.max(0, percentile(lengths.array, lengths.size, 999, 1000)-K);
		long checkedOccurrences=0, checkedReferences=0;
		for(final Key key:stats){
			key.finish();
			checkedOccurrences+=key.left.size;
			checkedReferences+=key.references;
			require(key.references<=expectedGenes, "A key contributes at most once per distinct reference");
		}
		require(checkedOccurrences==occurrences && checkedReferences==referenceKeys, "Per-key and per-reference populations must conserve each other");
		if(boundsOnly){
			writeBounds();
			writeRecordLengths();
			log.println("NCRNA_SEED_BOUNDS_REPORTED member="+member+" genes="+expectedGenes+" keys="+keys.length
				+" occurrences="+occurrences+" distinct_reference_key_pairs="+referenceKeys
				+" gene_max="+lengths.get(lengths.size-1)+" estimator="+NcrnaOverlapWindowBuilder.capRule(boundPercentile, boundPadding));
			return;
		}
		for(int method=0; method<METHODS.length; method++){
			for(int rule=0; rule<RULES.length; rule++){
				for(int stratum=0; stratum<STRATA.length; stratum++){windows[method][rule][stratum]=new Windows();}
			}
		}
		secondPass();
		writeKeys();
		writeLengths();
		writeWindows();
		log.println("NCRNA_SEED_POSITIONS_REPORTED member="+member+" genes="+expectedGenes+" keys="+keys.length
			+" occurrences="+occurrences+" distinct_reference_key_pairs="+referenceKeys+"; "+QUALIFICATION);
	}
	/** Returns an unstarted native FASTA owner; each pass starts and drains its own reader. */
	private SequenceInput reader(){
		final FileFormat format=FileFormat.testInput(input, FileFormat.FA, null, true, true);
		require(format.fasta(), "Reference measurements require whole reference-gene FASTA records");
		return new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1));
	}
	/** Reads each FASTA record once, including stdin; no sequence, ID or seed state is retained.
	 * Record multiplicity is preserved. The caller names the reference population;
	 * this mode neither deduplicates sequences nor judges annotation validity. */
	private void measureLengths(){
		try(SequenceInput in=reader()){
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					require(read!=null && read.mate==null && read.bases!=null && read.length()>0,
						"Length mode requires one complete nonempty unpaired gene per record; next record="+(lengths.size+1));
					lengths.add(read.length());
				}
			}
		}
		require(lengths.size>0 && (expectedGenes<0 || lengths.size==expectedGenes),
			"Length statistics need a nonempty complete reference cohort: observed="+lengths.size+", expected="+expectedGenes);
		lengths.sort();
		writeRecordLengths();
		log.println("NCRNA_REFERENCE_LENGTHS_REPORTED member="+member+" records="+lengths.size);
	}
	/** Shared record-length output; bounds mode must not claim unmeasured second-pass window strata. */
	private void writeRecordLengths(){
		require(lengths.size>0, "Record-length output requires an actual nonempty sorted full-member cohort");
		final int p999=percentile(lengths.array, lengths.size, 999, 1000);
		try(Output output=new Output(out+".lengths.tsv")){
			output.writer.println("member\trecords\tgene_min\tgene_median\tgene_p99\tgene_p999\tgene_max\tp999_minus17_cap\tmax_minus17_cap\tqualification");
			row.clear().append(member).tab().append(lengths.size);
			appendQuantiles(row, lengths);
			output.writer.print(row.tab().append(Math.max(0, p999-K))
				.tab().append(Math.max(0, lengths.get(lengths.size-1)-K))
				.append("\tREFERENCE_RECORD_LENGTHS_ONLY;NO_SEED_WINDOW_OR_CALLER_MEASUREMENT\n"));
		}
	}
	/** Stores endpoint distances in bases and unique IDs/lengths; does not retain sequences. */
	private void firstPass(){
		try(SequenceInput in=reader()){
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					require(read!=null && read.mate==null && read.bases!=null && read.length()>0, "One nonempty gene per unpaired FASTA record");
					final String id=RrnaResourceIO.first(read.id);
					require(!id.isEmpty() && genes.put(id, new Gene(read.length()))==null, "Reference IDs must be unique within this member: "+id);
					lengths.add(read.length());
					scratch.scan(read.bases, index);
					occurrences+=scratch.hitKeys.size;
					referenceKeys+=scratch.distinct.size;
					for(int i=0; i<scratch.hitKeys.size; i++){
						final Key key=stats[scratch.hitKeys.get(i)];
						key.left.add(scratch.starts.get(i));
						key.right.add(read.length()-K-scratch.starts.get(i));
					}
					for(int i=0; i<scratch.distinct.size; i++){stats[scratch.distinct.get(i)].references++;}
				}
			}
		}
		require(genes.size()==expectedGenes && lengths.size==expectedGenes, "The declared full reference population must survive the first pass");
	}
	/** Reuses the same genes to measure caps; all reference windows are zero-based inclusive. */
	private void secondPass(){
		long observedHits=0, observedPairs=0;
		int observedGenes=0;
		try(Output rows=new Output(out+".references.tsv"); SequenceInput in=reader()){
			rows.writer.println("member\treference_id\tgene_length\toccurrences\tdistinct_keys\tmethod\trule\thas_window\tstart0\tstop0\twindow_length\tfully_contains_gene\tqualification");
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					final String id=RrnaResourceIO.first(read.id);
					final Gene gene=genes.get(id);
					require(gene!=null && !gene.seen && read.mate==null && read.length()==gene.length,
						"Both streaming passes need identical reference identities and lengths: "+id);
					gene.seen=true;
					observedGenes++;
					scratch.scan(read.bases, index);
					observedHits+=scratch.hitKeys.size;
					observedPairs+=scratch.distinct.size;
					final int stratum=scratch.distinct.size==0 ? 0 : scratch.distinct.size==1 ? 1 : 2;
					geneStrata[stratum]++;
					geneStrata[3]++;
					for(int method=0; method<METHODS.length; method++){
						scratch.proposals(stats, method, blanket);
						containedHits[method]+=scratch.containedHits;
						for(int i=0; i<scratch.hitLengths.size; i++){hitLengths[method].add(scratch.hitLengths.get(i));}
						for(int rule=0; rule<RULES.length; rule++){
							final boolean present=scratch.distinct.size>0;
							final long start=present ? scratch.combineStart(rule) : 0, stop=present ? scratch.combineStop(rule) : -1;
							windows[method][rule][stratum].add(present, start, stop, gene.length);
							windows[method][rule][3].add(present, start, stop, gene.length);
							row.clear().append(member).tab();
							BuildNcrnaReference.field(row, id).tab().append(gene.length)
								.tab().append(scratch.hitKeys.size).tab().append(scratch.distinct.size).tab().append(METHODS[method]).tab()
								.append(RULES[rule]).tab().append(present);
							if(present){
								row.tab().append(start).tab().append(stop).tab().append(Math.max(0L, stop-start+1))
									.tab().append(start<=0 && stop>=gene.length-1L);
							}else{
								row.append("\tNA\tNA\tNA\tNA");
							}
							rows.writer.print(row.tab().append(QUALIFICATION).nl());
						}
					}
				}
			}
		}
		require(observedGenes==expectedGenes && observedHits==occurrences && observedPairs==referenceKeys,
			"Second-pass gene, occurrence and distinct-key populations must conserve the first pass");
	}

	/** Emits every selected forward key, including unavailable caps for unobserved keys. */
	private void writeKeys(){
		try(Output output=new Output(out+".keys.tsv")){
			output.writer.println("member\tforward_17mer\toccurrences\tdistinct_references\tleft_cap_bp\tright_cap_bp\tper_hit_window_bp\testimator");
			for(int i=0; i<keys.length; i++){
				final Key key=stats[i];
				row.clear().append(member).tab().appendKmer(keys[i], K).tab().append(key.left.size).tab().append(key.references);
				if(key.left.size==0){
					row.append("\tNA\tNA\tNA\tNO_OCCURRENCES");
				}else{
					row.tab().append(key.leftCap).tab().append(key.rightCap).tab().append((long)K+key.leftCap+key.rightCap)
						.tab().append(key.left.size<1000 ? "MAX_N_LT_1000" : "NEAREST_RANK_P999");
				}
				output.writer.print(row.nl());
			}
		}
	}
	/** Makes serving margins from the complete admitted originals in one counting
	 * pass. Each end uses its own occurrence-weighted percentile plus explicit padding; the cap
	 * is the FULL member maximum minus seventeen, never the key's observed maximum. */
	private void writeBounds(){
		final int maximum=lengths.get(lengths.size-1);
		final String estimator=NcrnaOverlapWindowBuilder.capRule(boundPercentile, boundPadding);
		try(Output output=new Output(out+".bounds.tsv")){
			output.writer.println(NcrnaOverlapWindowBuilder.CAP_HEADER);
			if(writeMetadata){
				final int ceiling=chainCeiling(lengths.array, lengths.size);
				output.writer.print(new ByteBuilder().append(NcrnaOverlapWindowBuilder.CEILING_TAG).tab().append(ceiling).tab()
					.append(lengths.size>=1000 ? NcrnaOverlapWindowBuilder.CEILING_P999 : NcrnaOverlapWindowBuilder.CEILING_MAX).nl());
				output.writer.print(new ByteBuilder().append(NcrnaOverlapWindowBuilder.POPULATION_TAG).tab().append(constructionPopulation).nl());
				if(boundCombination!=null){
					output.writer.print(new ByteBuilder().append(NcrnaOverlapWindowBuilder.COMBINATION_TAG).tab().append(boundCombination).nl());
				}
			}
			for(int i=0; i<keys.length; i++){
				final Key key=stats[i];
				row.clear().append(member).tab().appendKmer(keys[i], K).tab().append(key.left.size).tab().append(key.references);
				if(key.left.size==0){row.append("\tNA\tNA\tNA\tNO_OCCURRENCES");}
				else{
					final int left=boundedMargin(key.left.array, key.left.size, maximum, boundPercentile, boundPadding);
					final int right=boundedMargin(key.right.array, key.right.size, maximum, boundPercentile, boundPadding);
					row.tab().append(left).tab().append(right).tab().append((long)K+left+right)
						.tab().append(estimator);
				}
				output.writer.print(row.nl());
			}
		}
	}
	/** Input lengths are sorted complete originals. The sample-size rule excludes one giant at n=1000. */
	public static int chainCeiling(final int[] sorted, final int size){
		require(sorted!=null && size>0 && size<=sorted.length,
			"Chain ceilings require the complete nonempty original-gene length population");
		final int length=size>=1000 ? percentile(sorted, size, 999, 1000) : sorted[size-1];
		require(length>K, "Protected K17 windows require a positive gene-length-minus17 chain ceiling");
		return length-K;
	}
	/** Reports complete member lengths and the distinct-key strata, without dropping no-hit genes. */
	private void writeLengths(){
		int zero=0, below100=0, below1000=0;
		for(final Key key:stats){
			if(key.left.size==0){zero++;}
			if(key.left.size<100){below100++;}
			if(key.left.size<1000){below1000++;}
		}
		try(Output output=new Output(out+".lengths.tsv")){
			output.writer.println("member\treferences\tgene_min\tgene_median\tgene_p99\tgene_p999\tgene_max\tblanket_p999_cap\tblanket_max_cap\tkeys\tzero_occurrence_keys\tkeys_lt100_occurrences\tkeys_lt1000_occurrences\tgenes_zero_keys\tgenes_one_key\tgenes_two_plus_keys\tqualification");
			row.clear().append(member).tab().append(expectedGenes);
			appendQuantiles(row, lengths);
			row.tab().append(blanket).tab().append(Math.max(0, lengths.get(lengths.size-1)-K)).tab().append(keys.length)
				.tab().append(zero).tab().append(below100).tab().append(below1000);
			for(int stratum=0; stratum<3; stratum++){row.tab().append(geneStrata[stratum]);}
			output.writer.print(row.tab().append(QUALIFICATION).nl());
		}
	}
	/** Writes occurrence and reference distributions, with each method/rule/stratum kept separate. */
	private void writeWindows(){
		try(Output summary=new Output(out+".windows.tsv"); Output histogram=new Output(out+".window_histogram.tsv")){
			summary.writer.println("member\tmethod\trule\tstratum\treferences_or_hits\tno_window\tempty_windows\tfully_contained\tnot_contained_with_window\tlength_min\tlength_median\tlength_p99\tlength_p999\tlength_max\tfailure_fraction_with_window\tfailure_fraction_all\tqualification");
			histogram.writer.println("member\tmethod\trule\tstratum\twindow_length\tcount");
			for(int method=0; method<METHODS.length; method++){
				final LongList hits=hitLengths[method];
				hits.sort();
				require(hits.size==occurrences, "Every observed occurrence contributes one per-hit window per method");
				windowRow(summary, histogram, method, "PER_HIT", "ALL", hits.size, 0, 0, containedHits[method], hits);
				for(int rule=0; rule<RULES.length; rule++){
					for(int stratum=0; stratum<STRATA.length; stratum++){
						final Windows w=windows[method][rule][stratum];
						w.lengths.sort();
						require(w.references==geneStrata[stratum] && w.lengths.size+w.absent==w.references,
							"Every reference must remain in its explicit seed-count stratum");
						windowRow(summary, histogram, method, RULES[rule], STRATA[stratum], w.references, w.absent, w.empty, w.contained, w.lengths);
					}
				}
			}
		}
	}
	/** Conditional failure excludes absent windows; all-reference failure counts them as failures. */
	private void windowRow(final Output output, final Output histogram, final int method, final String rule,
		final String stratum, final long count, final long absent, final long empty, final long contained, final LongList values){
		final long present=count-absent;
		assert(count>=0 && absent>=0 && absent<=count && empty>=0 && empty<=present
			&& contained>=0 && contained<=present-empty && values.size==present)
			: "Window populations must partition into absent, empty, contained and other present windows: "
			+count+"/"+absent+"/"+empty+"/"+contained+" lengths="+values.size;
		row.clear().append(member).tab().append(METHODS[method]).tab().append(rule).tab().append(stratum).tab().append(count).tab().append(absent).tab().append(empty);
		row.tab().append(contained).tab().append(present-contained);
		appendQuantiles(row, values);
		if(present==0){row.append("\tNA");}else{row.tab().append((double)(present-contained)/present, 8);}
		if(count==0){row.append("\tNA");}else{row.tab().append((double)(count-contained)/count, 8);}
		output.writer.print(row.tab().append(QUALIFICATION).nl());
		for(int i=0; i<values.size; ){
			int end=i+1;
			while(end<values.size && values.get(end)==values.get(i)){end++;}
			histogram.writer.print(row.clear().append(member).tab().append(METHODS[method]).tab().append(rule).tab().append(stratum)
				.tab().append(values.get(i)).tab().append(end-i).nl());
			i=end;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends min, nearest-rank median, p99, p99.9 and max from sorted gene lengths in bases. */
	private static void appendQuantiles(final ByteBuilder out, final IntList values){
		require(values.size>0, "Gene-length quantiles require the complete nonempty reference cohort");
		out.tab().append(values.get(0)).tab().append(percentile(values.array, values.size, 1, 2)).tab().append(percentile(values.array, values.size, 99, 100))
			.tab().append(percentile(values.array, values.size, 999, 1000)).tab().append(values.get(values.size-1));
	}
	/** Uses the same ranks for sorted window lengths; an absent population has NA quantiles. */
	private static void appendQuantiles(final ByteBuilder out, final LongList values){
		if(values.size==0){
			out.append("\tNA\tNA\tNA\tNA\tNA");
			return;
		}
		out.tab().append(values.get(0)).tab().append(values.get(rank(values.size, 1, 2))).tab().append(values.get(rank(values.size, 99, 100)))
			.tab().append(values.get(rank(values.size, 999, 1000))).tab().append(values.get(values.size-1));
	}
	/** Returns the zero-based index for ceil(n*numerator/denominator), with no interpolation. */
	public static int rank(final int n, final int numerator, final int denominator){
		require(n>0 && numerator>0 && numerator<=denominator, "Nearest rank requires a nonempty sample and probability in (0,1]");
		return (int)(((long)numerator*n+denominator-1)/denominator)-1;
	}
	/** Borrows a sorted populated prefix; values and returned caps have the same base-pair units. */
	public static int percentile(final int[] sorted, final int n, final int numerator, final int denominator){
		require(sorted!=null && n<=sorted.length, "Quantiles use the populated sorted prefix");
		return sorted[rank(n, numerator, denominator)];
	}
	/** Independent end-distance nearest rank; the caller supplies the maximum
	 * length of the complete member cohort, including genes without this key. */
	public static int boundedMargin(final int[] sorted, final int n, final int memberMaximum){
		return boundedMargin(sorted, n, memberMaximum, 99);
	}
	/** The explicit serving percentile uses the same occurrence weights and
	 * full-member cap as p99; only its nearest rank changes. */
	public static int boundedMargin(final int[] sorted, final int n, final int memberMaximum, final int percentile){
		return boundedMargin(sorted, n, memberMaximum, percentile, 10);
	}
	/** Exact tenth-percent ranks also express p99.9; padding never exceeds the full-member cap. */
	public static int boundedMargin(final int[] sorted, final int n, final int memberMaximum, final double percentile, final int padding){
		NcrnaOverlapWindowBuilder.capRule(percentile, padding);
		require(memberMaximum>0, "Serving padding is bounded by the measured full member maximum, not a per-key gene subset");
		final int numerator=(int)Math.round(percentile*10);
		final int margin=percentile(sorted, n, numerator, 1000);
		require(margin>=0 && (long)margin+K<=memberMaximum,
			"An observed key-end distance plus its seed footprint must fit within the complete member maximum");
		return (int)Math.min(Math.max(0L, memberMaximum-(long)K), margin+(long)padding);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns two primitive occurrence arrays; distinct-reference counts do not weight the quantiles. */
	public static final class Key{
		public final IntList left=new IntList(1), right=new IntList(1);
		public long references;
		public int leftCap=-1, rightCap=-1;
		/** Sorts each edge independently; p99.9 equals the maximum when n is below 1000. */
		public void finish(){
			require(left.size==right.size && references<=left.size, "Each occurrence owns two distances and each reference at least one occurrence");
			left.sort();
			right.sort();
			if(left.size>0){
				leftCap=percentile(left.array, left.size, 999, 1000);
				rightCap=percentile(right.array, right.size, 999, 1000);
			}
		}
	}
	/** Reusable single-consumer state: scan resets active hits, proposals replace their key envelopes. */
	public static final class Scratch{
		/** Allocates one envelope slot for every key in the fixed shared axis. */
		public Scratch(final int keys){
			require(keys>0, "Scratch belongs to a complete seed axis");
			present=new boolean[keys];
			left=new long[keys];
			right=new long[keys];
		}
		/** Records every forward 17-mer occurrence; non-ACGT bases break the rolling word. */
		public void scan(final byte[] bases, final LongIntMap index){
			require(bases!=null && index!=null, "Scan only actual sense-oriented bases and a fixed forward seed index");
			geneLength=bases.length;
			for(int i=0; i<distinct.size; i++){present[distinct.get(i)]=false;}
			distinct.clear();
			hitKeys.clear();
			starts.clear();
			long word=0;
			int length=0;
			final byte[] numbers=AminoAcid.baseToNumber;
			for(int i=0; i<bases.length; i++){
				final int symbol=bases[i]&255;
				require(symbol<numbers.length, "Reference base is outside the native nucleotide alphabet");
				final int value=numbers[symbol];
				if(value<0){
					length=0;
					word=0;
					continue;
				}
				word=((word<<2)|value)&MASK;
				length++;
				if(length<K){continue;}
				final int key=index.get(word);
				if(key<0){continue;}
				require(key<present.length, "Loaded seed index must stay inside its fixed axis");
				hitKeys.add(key);
				starts.add(i-K+1);
				if(!present[key]){
					present[key]=true;
					distinct.add(key);
				}
			}
		}
		/** Builds inclusive unclipped proposals, then envelopes repeated occurrences of each key. */
		public void proposals(final Key[] stats, final int method, final int blanket){
			require(stats!=null && stats.length==present.length && method>=0 && method<3 && blanket>=0, "Bounded seed/recipe axes required");
			hitLengths.clear();
			containedHits=0;
			for(int i=0; i<distinct.size; i++){
				final int key=distinct.get(i);
				left[key]=Long.MAX_VALUE;
				right[key]=Long.MIN_VALUE;
			}
			for(int i=0; i<hitKeys.size; i++){
				final int key=hitKeys.get(i), start=starts.get(i);
				final long a, z;
				if(method==0){
					require(stats[key].leftCap>=0 && stats[key].rightCap>=0, "Observed keys need measured caps");
					a=(long)start-stats[key].leftCap;
					z=(long)start+K-1+stats[key].rightCap;
				}else{
					final int pad=method==1 ? blanket : 450;
					final long center=(long)start+K/2;
					a=center-pad;
					z=center+pad+K;
				}
				left[key]=Math.min(left[key], a);
				right[key]=Math.max(right[key], z);
				hitLengths.add(z-a+1);
				if(a<=0 && z>=geneLength-1L){containedHits++;}
			}
		}
		/** Uses minimum left edge for union, maximum for intersection across distinct-key envelopes. */
		public long combineStart(final int rule){
			require(distinct.size>0 && (rule==0 || rule==1), "Combine one envelope per distinct key, with explicit union or intersection");
			long result=rule==0?Long.MAX_VALUE:Long.MIN_VALUE;
			for(int i=0; i<distinct.size; i++){
				final long value=left[distinct.get(i)];
				result=rule==0?Math.min(result, value):Math.max(result, value);
			}
			return result;
		}
		/** Uses maximum right edge for union, minimum for intersection; stop below start means empty. */
		public long combineStop(final int rule){
			require(distinct.size>0 && (rule==0 || rule==1), "Combine one envelope per distinct key, with explicit union or intersection");
			long result=rule==0?Long.MIN_VALUE:Long.MAX_VALUE;
			for(int i=0; i<distinct.size; i++){
				final long value=right[distinct.get(i)];
				result=rule==0?Math.max(result, value):Math.min(result, value);
			}
			return result;
		}
		public final IntList hitKeys=new IntList(), starts=new IntList(), distinct=new IntList();
		public final LongList hitLengths=new LongList();
		public int containedHits;
		private int geneLength;
		private final boolean[] present;
		private final long[] left, right;
	}
	/** Retains only the first-pass length and second-pass identity-visit state. */
	private static final class Gene{
		Gene(final int length_){length=length_;}
		final int length;
		boolean seen;
	}
	/** Owns one reference population and its present-window lengths, including empty intersections. */
	private static final class Windows{
		/** Adds one full gene: absence, empty interval, and incomplete containment remain distinct. */
		void add(final boolean present, final long start, final long stop, final int length){
			assert(length>0):"Window containment uses the validated positive full-gene length";
			references++;
			if(!present){
				absent++;
				return;
			}
			final long size=Math.max(0L, stop-start+1);
			lengths.add(size);
			if(size==0){empty++;}
			if(start<=0 && stop>=length-1L){contained++;}
		}
		long references, absent, empty, contained;
		final LongList lengths=new LongList();
	}
	/** Rejects malformed inputs or violated population/geometry contracts without fallback. */
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream log;
	private String input, out, seedPath, member;
	private String constructionPopulation=NcrnaOverlapWindowBuilder.POPULATION_ORIGINALS;
	private boolean lengthsOnly, boundsOnly, writeMetadata;
	private int expectedGenes=-1, expectedKeys=-1, blanket;
	private double boundPercentile=99;
	private int boundPadding=10;
	private String boundCombination;
	private long[] keys;
	private LongIntMap index;
	private Key[] stats;
	private Scratch scratch;
	private long occurrences, referenceKeys;
	private final IntList lengths=new IntList();
	private final ObjectMap<String, Gene> genes=new ObjectMap<String, Gene>(String.class, Gene.class);
	private final LongList[] hitLengths={new LongList(), new LongList(), new LongList()};
	private final long[] containedHits=new long[3];
	private final Windows[][][] windows=new Windows[3][2][4];
	private final long[] geneStrata=new long[4];
	private final ByteBuilder row=new ByteBuilder();

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int K=17;
	private static final long MASK=(1L<<(2*K))-1;
	private static final String[] METHODS={"PER_KEY_P999", "MEMBER_LENGTH_P999_CAP", "LEGACY_PAD450"}, RULES={"UNION", "INTERSECTION"}, STRATA={"ZERO_KEYS", "ONE_KEY", "TWO_PLUS_KEYS", "ALL"};
	private static final String QUALIFICATION="IN_CORPUS_REFERENCE_ONLY;SENSE_BASES;UNCLIPPED;REPEATS_ENVELOPED_PER_KEY;NO_GENOMIC_FLANKS_CLAIMS_OR_CALLER_REPLAY";
}
