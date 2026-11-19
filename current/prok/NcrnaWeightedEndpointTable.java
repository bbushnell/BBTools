package prok;

import prok.NcrnaSignedEndpointErrors.Distribution;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RrnaPositionalKmerTable.Table;
import parse.LineParser1;
import parse.Parse;
import structures.ByteBuilder;
import structures.DoubleList;
import structures.IntList;
import structures.LongList;

/** Optional exact signed-error background for an existing positional geometry.
 * The borrowed table still owns the correction radius and in-window counts.
 * Additional observed error offsets own dense primitive rows; overlapping rows
 * alias the borrowed table and are never counted twice. Count on one worker,
 * finish after the borrowed table, then share immutable scores/counts. Training
 * workers own separate leave-one-locus-out views. This changes no legacy table.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWeightedEndpointTable{

	/** Binds one member/end's complete signed histogram to one unfinalized table. */
	public NcrnaWeightedEndpointTable(final Table geometry_, final NcrnaSignedEndpointErrors errors){
		require(geometry_!=null && geometry_.counts!=null && !geometry_.finished
			&& errors!=null && geometry_.model.equals(errors.member),
			"Weighted background requires an unfinalized counted table and its exact member's signed errors");
		geometry=geometry_;
		final Distribution distribution=errors.distribution(geometry.end.equals("5prime") ? 0 : 1);
		require(distribution.incorrect()>0, "A wrong-end background needs observed nonzero errors, never a uniform fallback");
		final int size=distribution.size()+(distribution.zero==0 ? 1 : 0);
		offsets=new int[size];
		observed=new long[size];
		weights=new double[size];
		counts=new long[size][];
		valid=new long[size];
		clipped=new long[size];
		ambiguous=new long[size];
		int next=0, zero=-1;
		for(int bin=0; bin<distribution.size(); bin++){
			final int offset=distribution.offset(bin);
			if(offset>0 && zero<0){zero=next++;}
			if(offset==0){zero=next;}
			offsets[next]=offset;
			observed[next]=distribution.count(bin);
			weights[next]=distribution.backgroundWeight(bin);
			next++;
		}
		if(zero<0){zero=next++;}
		anchor=zero;
		require(next==size && offsets[anchor]==0 && weights[anchor]==0,
			"Measured background support must include one distinct positive anchor with zero wrong-end weight");
		for(int bin=0; bin<size; bin++){
			final int offset=offsets[bin];
			counts[bin]=inCorrectionWindow(offset) ? geometry.counts[offset+geometry.radius] : new long[geometry.space];
		}
	}

	/** Takes validated saved metadata and scores, without any training-count rows. */
	private NcrnaWeightedEndpointTable(final Table geometry_, final int[] offsets_, final long[] observed_,
			final double[] weights_, final long[] valid_, final long[] clipped_, final long[] ambiguous_,
			final int anchor_, final double[] scores_){
		assert(geometry_.counts==null && geometry_.records>0 && scores_.length==geometry_.space)
			: "The signed-score reader publishes a complete inference-only key space, never reconstructed training counts";
		geometry=geometry_;
		offsets=offsets_;
		observed=observed_;
		weights=weights_;
		valid=valid_;
		clipped=clipped_;
		ambiguous=ambiguous_;
		anchor=anchor_;
		counts=null;
		records=geometry.records;
		scores=scores_;
		finished=true;
		geometry.scores=scores;
		geometry.finished=true;
	}

	/** Reads the exact dense signed-background format written by writeScores().
	 * Keeps one double score per key and small signed-bin metadata. Histogram counts
	 * describe caller matches; table_loci describes teaching loci, so those totals
	 * deliberately need not agree. Each offset's observation accounting must agree
	 * with table_loci. Saved doubles use an exact cold-path decimal round trip.
	 */
	public static NcrnaWeightedEndpointTable loadScores(final String path, final String member){
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "#NCRNA_SIGNED_BACKGROUND_SCORES_V1", path);
			RrnaResourceIO.header(input.reader, "member\tend\tk\tradius_bp\tprobe_offset_bp\talpha_per_key\ttable_loci", path);
			fields(parser, input.reader.nextLine(), 7);
			require(member!=null && parser.termEquals(member, 0), "Signed scores must name the requested exact member");
			final Table table=RrnaPositionalKmerTable.inferenceGeometry(parser.parseString(0), parser.parseString(1),
				parser.parseInt(2), parser.parseInt(4), exactDouble(parser, 5), parser.parseInt(3));
			table.records=parser.parseLong(6);
			require(table.records>1, "Signed scores require the producer's nonempty teaching population");
			RrnaResourceIO.header(input.reader, "#offset_bp\tobserved_errors\tbackground_weight\tvalid\tclipped\tambiguous", path);
			final IntList offsets=new IntList();
			final DoubleList weights=new DoubleList();
			final LongList observed=new LongList(), valid=new LongList(), clipped=new LongList(), ambiguous=new LongList();
			int anchor=-1, previous=0;
			long incorrect=0;
			byte[] row=input.reader.nextLine();
			while(row!=null && row.length>0 && row[0]=='#'){
				fields(parser, row, 6);
				parser.setBounds(0);
				require(parser.length(0)>1, "A signed-background bin needs its integer offset after '#'");
				final int offset=Parse.parseInt(row, parser.a()+1, parser.b());
				final long count=parser.parseLong(1), usable=parser.parseLong(3), cut=parser.parseLong(4), bad=parser.parseLong(5);
				final double weight=exactDouble(parser, 2);
				require((offsets.size==0 || offset>previous) && count>=0 && usable>0 && cut>=0 && bad>=0
					&& Double.isFinite(weight) && weight>=0
					&& Math.addExact(Math.addExact(usable, cut), bad)==table.records,
					"Signed bins must be ordered, unique, nonnegative and account for every teaching locus");
				if(offset==0){
					anchor=offsets.size;
					require(weight==0, "The positive anchor has zero wrong-end background weight");
				}else{
					require(count>0, "Every nonzero signed offset must have measured error support");
					incorrect=Math.addExact(incorrect, count);
				}
				offsets.add(offset);
				observed.add(count);
				weights.add(weight);
				valid.add(usable);
				clipped.add(cut);
				ambiguous.add(bad);
				previous=offset;
				row=input.reader.nextLine();
			}
			require(anchor>=0 && incorrect>0, "Signed scores require a positive anchor and a measured nonzero-error distribution");
			final double[] savedWeights=weights.toArray();
			final long[] savedObserved=observed.toArray();
			for(int bin=0; bin<savedWeights.length; bin++){
				final double expected=bin==anchor ? 0 : savedObserved[bin]/(double)incorrect;
				require(savedWeights[bin]==expected, "Saved background weights must match Distribution.backgroundWeight's nonzero-conditioned counts");
			}
			fields(parser, row, 2);
			require(parser.termEquals("kmer", 0) && parser.termEquals("log_score", 1), "Missing signed-score key header");
			final double[] scores=new double[table.space];
			final boolean[] seen=new boolean[table.space];
			int keys=0;
			for(row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				fields(parser, row, 2);
				parser.setBounds(0);
				require(parser.length(0)==table.k, "Signed-score keys must have the declared word length");
				final int word=RrnaPositionalKmerTable.encode(row, parser.a(), table.k);
				require(word>=0 && !seen[word], "Invalid or duplicate signed-score key");
				final double value=exactDouble(parser, 1);
				require(Double.isFinite(value), "Nonfinite signed score would corrupt caller inputs");
				seen[word]=true;
				scores[word]=value;
				keys++;
			}
			require(keys==table.space, "Signed-score format writes every key; no uniform or neutral omitted-key fallback is valid");
			return new NcrnaWeightedEndpointTable(table, offsets.toArray(), savedObserved, savedWeights,
				valid.toArray(), clipped.toArray(), ambiguous.toArray(), anchor, scores);
		}
	}

	/** Empty fields must not become numeric zeros in saved model metadata or scores. */
	private static void fields(final LineParser1 parser, final byte[] row, final int count){
		require(row!=null, "Signed-score file ended before a required record");
		parser.set(row);
		require(parser.terms()==count, "Signed-score record has the wrong field count");
		for(int field=0; field<count; field++){
			require(parser.length(field)>0, "Signed-score fields must be nonempty");
		}
	}

	/** Exact parsing is confined to resource loading, never the per-call feature path. */
	private static double exactDouble(final LineParser1 parser, final int term){
		parser.setBounds(term);
		return Parse.parseDoubleSlow(parser.line(), parser.a(), parser.b());
	}

	/** Call once per same locus AFTER the borrowed correction table counts it. */
	public void add(final byte[] bases, final int trueEnd){
		require(!finished && bases!=null && trueEnd>=0 && trueEnd<bases.length,
			"Count only real sense-oriented endpoints before weighted finalization");
		records++;
		for(int bin=0; bin<offsets.length; bin++){
			if(inCorrectionWindow(offsets[bin])){continue;}
			final long position=(long)trueEnd+geometry.offset+offsets[bin];
			if(position<0 || position+geometry.k>bases.length){
				clipped[bin]++;
				continue;
			}
			final int word=RrnaPositionalKmerTable.encode(bases, (int)position, geometry.k);
			if(word<0){ambiguous[bin]++;}
			else{
				counts[bin][word]++;
				valid[bin]++;
			}
		}
	}

	/** Counts/denominators are checked before exact nonzero-conditioned scores open. */
	public void finish(){
		require(!finished && geometry.finished && records==geometry.records && records>1,
			"Finalize the same nonempty correction population before its weighted background");
		for(int bin=0; bin<offsets.length; bin++){
			if(inCorrectionWindow(offsets[bin])){
				final int site=offsets[bin]+geometry.radius;
				valid[bin]=geometry.valid[site];
				clipped[bin]=geometry.clipped[site];
				ambiguous[bin]=geometry.ambiguous[site];
			}
			long sum=0;
			for(final long count:counts[bin]){
				require(count>=0, "Count overflow cannot define a weighted endpoint frequency");
				sum=Math.addExact(sum, count);
			}
			require(sum==valid[bin] && valid[bin]+clipped[bin]+ambiguous[bin]==records,
				"Each measured offset must account for every locus as valid, clipped or ambiguous");
			require(valid[bin]>0, "A positively weighted signed offset or anchor has no observed words; do not omit or renormalize it");
		}
		scores=new double[geometry.space];
		for(int word=0; word<scores.length; word++){scores[word]=computeScore(word, null, null);}
		finished=true;
	}

	/** Scores use each error offset's own denominator; zero never enters background. */
	private double computeScore(final int word, final int[] own, final double[] remaining){
		assert(word>=0 && word<geometry.space) : "Weighted table lookup must name a valid packed key";
		final double prior=geometry.alpha*geometry.space;
		final double numerator=(counts[anchor][word]-(own!=null && own[anchor]==word ? 1 : 0)+geometry.alpha)
			/(remaining==null ? valid[anchor]+prior : remaining[anchor]);
		double background=0;
		for(int bin=0; bin<offsets.length; bin++){
			if(bin==anchor){continue;}
			final double frequency=(counts[bin][word]-(own!=null && own[bin]==word ? 1 : 0)+geometry.alpha)
				/(remaining==null ? valid[bin]+prior : remaining[bin]);
			background+=weights[bin]*frequency;
		}
		require(numerator>0 && background>0 && Double.isFinite(numerator) && Double.isFinite(background),
			"Separate smoothed anchor/background frequencies must remain finite and positive");
		final double score=Math.log(numerator/background);
		require(Double.isFinite(score), "Finite positive measured frequencies must yield a finite log enrichment");
		return score;
	}

	public double scoreKey(final int word){
		require(finished && word>=0 && word<scores.length, "Use a finalized weighted table and valid packed key");
		return scores[word];
	}

	/** Correction width is unchanged even when background errors extend far away.
	 * False means a clipped/ambiguous word; the incomplete output must not be used.
	 * A computed nonfinite score is an error, never missing-word state.
	 */
	public boolean scoreWindow(final byte[] bases, final int candidateEnd, final float[] output, final int from){
		require(finished && bases!=null && output!=null && from>=0 && from<=output.length-geometry.sites,
			"Provide the complete correction window for finalized weighted scores");
		for(int site=0; site<geometry.sites; site++){
			final int word=wordAt(bases, (long)candidateEnd+geometry.offset+site-geometry.radius);
			if(word<0){return false;}
			final float value=(float)scores[word];
			require(Float.isFinite(value), "A measured weighted score must fit the finite neural input domain");
			output[from+site]=value;
		}
		return true;
	}

	/** Experimental score-only schema explicitly identifies the changed background.
	 * Training counts remain in memory for LOO; the recorded producer regenerates them.
	 * This format is not accepted by the legacy uniform-background resource loader.
	 */
	public void writeScores(final String path){
		require(finished, "Write weighted scores only after all counted offsets have been validated");
		try(Output output=new Output(path)){
			final ByteBuilder buffer=new ByteBuilder();
			output.writer.println("#NCRNA_SIGNED_BACKGROUND_SCORES_V1");
			output.writer.println("member\tend\tk\tradius_bp\tprobe_offset_bp\talpha_per_key\ttable_loci");
			output.writer.print(buffer.append(geometry.model).tab().append(geometry.end).tab().append(geometry.k)
				.tab().append(geometry.radius).tab().append(geometry.offset).tab().append(Double.toString(geometry.alpha)).tab().append(records).nl());
			output.writer.println("#offset_bp\tobserved_errors\tbackground_weight\tvalid\tclipped\tambiguous");
			for(int bin=0; bin<offsets.length; bin++){
				output.writer.print(buffer.clear().append('#').append(offsets[bin]).tab().append(observed[bin])
					.tab().append(Double.toString(weights[bin])).tab().append(valid[bin]).tab().append(clipped[bin]).tab().append(ambiguous[bin]).nl());
			}
			output.writer.println("kmer\tlog_score");
			for(int word=0; word<scores.length; word++){
				buffer.clear();
				RrnaPositionalKmerTable.appendKmer(buffer, word, geometry.k);
				output.writer.print(buffer.tab().append(Double.toString(scores[word])).nl());
			}
		}
	}

	/** Independent worker state for exact literal omitted-locus normalization. */
	public final class View{
		private View(){
			require(finished && counts!=null, "Training views require immutable training counts; score-only serving tables cannot omit a locus");
			own=new int[offsets.length];
			remaining=new double[offsets.length];
		}
		/** Binds only counted words; a failed binding invalidates the prior locus. */
		public void set(final byte[] bases, final int trueEnd){
			ready=false;
			require(bases!=null && trueEnd>=0 && trueEnd<bases.length, "Bind a real counted sense endpoint before weighted LOO");
			final double prior=geometry.alpha*geometry.space;
			for(int bin=0; bin<offsets.length; bin++){
				final int word=wordAt(bases, (long)trueEnd+geometry.offset+offsets[bin]);
				own[bin]=word;
				require(word<0 || counts[bin][word]>0, "An omitted locus must have contributed each of its valid weighted words");
				final long observations=valid[bin]-(word>=0 ? 1 : 0);
				require(observations>0, "Removing this locus leaves an unsupported measured offset; do not renormalize background weights");
				remaining[bin]=observations+prior;
			}
			ready=true;
		}
		public double scoreKey(final int word){
			require(ready && word>=0 && word<geometry.space, "Bind a counted locus before weighted leave-one-out key lookup");
			return computeScore(word, own, remaining);
		}
		private final int[] own;
		private final double[] remaining;
		private boolean ready;
	}

	public View leaveOneOut(){return new View();}
	private boolean inCorrectionWindow(final int offset){return offset>=-geometry.radius && offset<=geometry.radius;}
	private int wordAt(final byte[] bases, final long position){
		return position<0 || position+geometry.k>bases.length ? -1
			: RrnaPositionalKmerTable.encode(bases, (int)position, geometry.k);
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public final Table geometry;
	private final int[] offsets;
	private final long[] observed, valid, clipped, ambiguous;
	private final double[] weights;
	private final long[][] counts;
	private final int anchor;
	private long records;
	private double[] scores;
	private boolean finished;
}
