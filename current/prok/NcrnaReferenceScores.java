package prok;

import java.util.Arrays;
import java.math.BigDecimal;
import java.math.RoundingMode;

import map.IntLongHashMap;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.Pair;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Describes captured CM scores without selecting an admission or masking rule.
 * Exact model spans and native truncation remain separate from biological length.
 * The histogram bins printed score/GA at one-percent intervals, with an open
 * upper bin at two times GA. Floor curves describe 80 through 90 percent of GA.
 * Native inclusion is authoritative only at original GA, not at a new floor.
 * Infernal1.1.5 prints one decimal (cm_tophits_TabularTargets1); +/-0.05 bits
 * therefore bounds rounding uncertainty. Equality at either bound is uncertain.
 * A native-included event is definitely above every displayed sub-GA floor.
 * Competitors are direct valid same-strand other-family edges, never components.
 * Primitive arrays/maps own counters; only the existing retained event loader
 * allocates per-event records. All mutable state belongs to one reporting thread.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceScores{

	static void process(final String input, final String output, final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceGeometry geometry=NcrnaReferenceGeometry.load(input, catalog, roster, true);
		require(geometry.profiles().batchEvidence(), "Score curves require original-GA-bound native inclusion, not legacy printed GA_pass");
		final NcrnaReferenceScores report=new NcrnaReferenceScores(geometry, roster);
		report.count();
		report.writeHistogram(output+SUFFIXES[0]);
		report.writeCurves(output+SUFFIXES[1]);
		report.writeSummary(output+SUFFIXES[2]);
	}

	private NcrnaReferenceScores(final NcrnaReferenceGeometry geometry_, final Roster roster_){
		geometry=geometry_;
		roster=roster_;
		final int models=roster.queried(), events=geometry.hits().length;
		histograms=new IntLongHashMap[models];
		competitors=new IntLongHashMap[models];
		curves=new long[models][2][11][COLUMNS];
		totals=new long[models][3];
		printed=new int[events];
		definite=new int[events];
		possible=new int[events];
		otherPrinted=new int[events];
		otherDefinite=new int[events];
		otherPossible=new int[events];
		anyOther=new boolean[events];
		binThresholds=new long[models][BINS];
		floorLower=new long[models][11];
		floorUpper=new long[models][11];
		for(int i=0; i<models; i++){
			final Model model=roster.query(i);
			require(model.ga>0, "Score/GA ratios require positive model GA: "+model.accession);
			Math.multiplyExact(Math.addExact(model.length, 1), BINS*16);
			histograms[i]=new IntLongHashMap();
			competitors[i]=new IntLongHashMap();
			// Decimal arithmetic occurs only here, never per event or histogram cell.
			final BigDecimal ga=BigDecimal.valueOf(model.ga);
			for(int bin=0; bin<BINS; bin++){
				binThresholds[i][bin]=ga.multiply(BigDecimal.valueOf(bin)).divide(BigDecimal.TEN, 0, RoundingMode.CEILING).longValueExact();
			}
			for(int step=0; step<11; step++){
				final BigDecimal scaled=ga.multiply(BigDecimal.valueOf(80+step));
				floorLower[i][step]=scaled.setScale(0, RoundingMode.FLOOR).longValueExact();
				floorUpper[i][step]=scaled.setScale(0, RoundingMode.CEILING).longValueExact();
			}
		}
	}

	/** Combines independent floor-qualification masks over actual direct neighbors. */
	private void competitor(final Hit own, final Hit other){
		assert(own.input==other.input && own.target==other.target && own.strand==other.strand)
			: "Only validated same-input/target/strand overlap edges supply competing evidence";
		final int i=own.index;
		anyOther[i]=true;
		otherPrinted[i]|=printed[other.index];
		otherDefinite[i]|=definite[other.index];
		otherPossible[i]|=possible[other.index];
	}

	private void count(){
		for(final Hit hit:geometry.hits()){classify(hit);}
		for(final Pair pair:geometry.pairs()){
			final Hit a=pair.a, b=pair.b;
			if(a.valid() && b.valid() && a.strand==b.strand
				&& !roster.query(a.query).scope.family.equals(roster.query(b.query).scope.family)){
				competitor(a, b);
				competitor(b, a);
			}
		}
		for(final Hit hit:geometry.hits()){
			final int model=hit.query, included=hit.included ? 1 : 0, span=hit.to-hit.from+1;
			final double ga=roster.query(model).ga;
			final int bin=scoreBin(tenths(hit.bits), binThresholds[model]);
			final int key=(((span*BINS+bin)*4+hit.truncation)*2+included)*2+(hit.left || hit.right ? 1 : 0);
			histograms[model].add(key);
			final boolean any=anyOther[hit.index];
			if(any){competitors[model].add(key);}
			totals[model][included]++;
			if(hit.included!=(hit.bits>=ga)){totals[model][2]++;}
			for(int step=0; step<11; step++){
				final int bit=1<<step;
				final long[] row=curves[model][included][step];
				row[0]++;
				if((definite[hit.index]&bit)!=0){row[2]++;}
				else if((possible[hit.index]&bit)!=0){row[3]++;}
				else{row[4]++;}
				if((printed[hit.index]&bit)==0){continue;}
				row[1]++;
				if(span==roster.query(model).length){row[5]++;}
				if(hit.truncation==0){row[6]++;}else{row[7]++;}
				if(hit.left || hit.right){row[8]++;}
				if(any){row[9]++;}
				if((otherPrinted[hit.index]&bit)!=0){row[10]++;}
				if((otherDefinite[hit.index]&bit)!=0){row[11]++;}
				else if((otherPossible[hit.index]&bit)!=0){row[12]++;}
			}
		}
		for(int i=0; i<roster.queried(); i++){
			require(histograms[i].sum()==geometry.profiles().modelEvents(i), "Score histogram must conserve every model's saved events");
		}
	}

	/** Integer bounds compare 100*bits +/-5 with exact decimal GA*percent.
	 * n>T iff n>floor(T); n>=T iff n>=ceil(T). Bound equality stays uncertain.
	 * Native inclusion independently proves every sub-GA pass.
	 */
	private void classify(final Hit hit){
		final long score=tenths(hit.bits), low=score*10-5, high=score*10+5;
		for(int step=0; step<11; step++){
			final int bit=1<<step;
			if(score>=binThresholds[hit.query][80+step]){printed[hit.index]|=bit;}
			if(hit.included || low>floorLower[hit.query][step]){definite[hit.index]|=bit;}
			if(hit.included || high>=floorUpper[hit.query][step]){possible[hit.index]|=bit;}
		}
		assert((definite[hit.index]&~possible[hit.index])==0) : "Every definite floor pass must also be possible";
	}

	/** Restricts conversion to one-decimal scores with exact integer arithmetic headroom. */
	private static long tenths(final double bits){
		final long score=Math.round(bits*10);
		require(bits>=0 && score<=(1L<<50) && bits==score/10.0,
			"Exact rounding bounds require a nonnegative one-decimal score within precise integer conversion range");
		return score;
	}

	/** Largest exact percent boundary reached; the last bin is open above 200%. */
	private static int scoreBin(final long score, final long[] thresholds){
		assert(score>=0 && thresholds.length==BINS && thresholds[0]==0) : "Histogram boundaries start at zero printed score";
		int low=0, high=BINS;
		while(low<high){
			final int mid=(low+high)>>>1;
			if(score>=thresholds[mid]){low=mid+1;}else{high=mid;}
		}
		return low-1;
	}

	private void writeHistogram(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("accession\tmodel_length\tmodel_span\tnative_truncation\tnative_included\tcontig_edge\tprinted_score_GA_percent_bin\traw_events\tdirect_other_family_overlap_events");
			for(int model=0; model<roster.queried(); model++){
				final int[] keys=histograms[model].toArray();
				Arrays.sort(keys);
				for(final int key:keys){
					int rest=key;
					final int edge=rest%2;
					rest/=2;
					final int included=rest%2;
					rest/=2;
					final int truncation=rest%4;
					rest/=4;
					final int bin=rest%BINS, span=rest/BINS;
					b.clear().append(roster.query(model).accession).tab().append(roster.query(model).length).tab();
					b.append(span).tab().append(truncation).tab().append(included==1).tab().append(edge==1).tab();
					if(bin==BINS-1){b.append(">=200");}else{b.append('[').append(bin).append(',').append(bin+1).append(')');}
					b.tab().append(histograms[model].get(key)).tab().append(Math.max(0, competitors[model].get(key))).nl();
					out.writer.print(b);
				}
			}
		}
	}

	private void writeCurves(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("accession\tGA_bits\tfloor_percent_GA\tnative_included\tevents\tprinted_pass\tdefinite_pass\trounding_sensitive\tdefinite_fail\tprinted_pass_full_model_span\tprinted_pass_untruncated\tprinted_pass_native_truncated\tprinted_pass_contig_edge\tprinted_pass_any_direct_other_family\tprinted_pass_other_family_printed_pass\tprinted_pass_other_family_definite_pass\tprinted_pass_other_family_threshold_uncertain");
			for(int model=0; model<roster.queried(); model++){
				for(int included=0; included<2; included++){
					for(int step=0; step<11; step++){
						final long[] row=curves[model][included][step];
						require(row[0]==totals[model][included] && row[0]==row[2]+row[3]+row[4]
							&& row[1]==row[6]+row[7], "Every floor row must conserve native inclusion, uncertainty and truncation counts");
						b.clear().append(roster.query(model).accession).tab().append(roster.query(model).ga, 6).tab();
						b.append(80+step).tab().append(included==1);
						for(final long count:row){b.tab().append(count);}
						out.writer.print(b.nl());
					}
				}
			}
		}
	}

	/** Written last; zero-event models still have explicit denominator/inclusion rows. */
	private void writeSummary(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("accession\tinputs\traw_events\tnative_included\tnative_excluded\tprinted_GA_disagreements\tstatus");
			for(int i=0; i<roster.queried(); i++){
				b.clear().append(roster.query(i).accession).tab().append(geometry.profiles().observed()).tab();
				b.append(totals[i][0]+totals[i][1]).tab().append(totals[i][1]).tab().append(totals[i][0]).tab();
				out.writer.print(b.append(totals[i][2]).append("\tSCORE_DISTRIBUTIONS_NOT_ADJUDICATED\n"));
			}
		}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final NcrnaReferenceGeometry geometry;
	private final Roster roster;
	private final IntLongHashMap[] histograms, competitors;
	private final long[][][][] curves;
	private final long[][] totals;
	private final int[] printed, definite, possible, otherPrinted, otherDefinite, otherPossible;
	private final boolean[] anyOther;
	private final long[][] binThresholds, floorLower, floorUpper;
	private static final int BINS=201, COLUMNS=13;
	static final String[] SUFFIXES={".score_histogram.tsv", ".floor_curve.tsv", ".score_summary.tsv"};
}
