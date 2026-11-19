package prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.LongHashSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import structures.IntList;

/** Immutable original evaluation trials and training-supported model offsets.
 * This class never draws examples, chooses models or changes truth. Original
 * sidecar bytes remain available for an all-trial ledger, including held spans.
 * Coordinates are zero-based inclusive in the sense-oriented context sequence.
 * @author Brian, Ganyu
 */
public final class NcrnaFixedEndpointTrials{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Load both frozen end populations and the same ordered library's offsets. */
	public static NcrnaFixedEndpointTrials load(final NcrnaLearningPopulation population,
			final NcrnaFamily family, final String member, final String startRows, final String stopRows,
			final String offsets, final int minimumSupport, final int samples, final long[] expected){
		require(population!=null && family!=null && member!=null && offsets!=null && startRows!=null && stopRows!=null
			&& family.library!=null && family.modelNames!=null && family.library.length==family.modelNames.length
			&& minimumSupport>0 && samples>0 && expected!=null && expected.length==2 && expected[0]>0 && expected[1]>0,
			"Fixed trials require the joined population, exact library, both sidecars, support floor and accepted counts");
		final NcrnaFixedEndpointTrials result=new NcrnaFixedEndpointTrials(population, family, member, minimumSupport);
		result.readOffsets(offsets, minimumSupport);
		result.ends[0]=result.readEnd(startRows, samples, expected[0]);
		result.ends[1]=result.readEnd(stopRows, samples, expected[1]);
		result.validateCommonTrials();
		return result;
	}

	private NcrnaFixedEndpointTrials(final NcrnaLearningPopulation population_, final NcrnaFamily family, final String member_, final int minimumSupport){
		population=population_;
		library=family;
		member=member_;
		supportMinimum=minimumSupport;
		modelNames=new String[family.modelNames.length];
		for(int model=0; model<modelNames.length; model++){modelNames[model]=NcrnaSupport.recordName(family.modelNames[model]);}
		median5=new int[modelNames.length];
		median3=new int[modelNames.length];
		trainingSupport=new long[modelNames.length];
		supported=new boolean[modelNames.length];
		minimumLength=family.minLen;
		require(minimumLength>0, "Centered geometry must use the loaded family's positive minimum length");
		modelByGene=new int[population.ids.length];
		Arrays.fill(modelByGene, -1);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Input Validation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Other members may share the resource, but this member must match library order exactly. */
	private void readOffsets(final String path, final int minimumSupport){
		final LineParser1 p=new LineParser1('\t');
		int model=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, OFFSET_HEADER, path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==6, "Every supported-offset resource row has six fields");
				if(!p.termEquals(member, 0)){continue;}
				require(model<modelNames.length && p.termEquals(modelNames[model], 1),
					"Offset resource must preserve the selected member's complete ordered model axis");
				final long support=p.parseLong(2);
				final int floor=p.parseInt(3);
				final long five=p.parseLong(4), three=p.parseLong(5);
				require(support>=0 && floor==minimumSupport && five>=Integer.MIN_VALUE && five<=Integer.MAX_VALUE
					&& three>=Integer.MIN_VALUE && three<=Integer.MAX_VALUE && (support>=floor || five==0 && three==0),
					"Offset support floor must match the declared experiment; unsupported models must have zero corrections");
				median5[model]=(int)five;
				median3[model]=(int)three;
				trainingSupport[model]=support;
				supported[model]=support>=floor;
				model++;
			}
		}
		require(model==modelNames.length, "Offset resource cannot omit models or silently invent zero offsets");
	}

	/** Retain one end's original order, with primitive linked indices grouped by locus. */
	private End readEnd(final String path, final int samples, final long expected){
		require(expected<=Integer.MAX_VALUE, "Fixed evaluation populations must fit Java's indexed arrays");
		final int[] first=new int[population.ids.length];
		final int[] last=new int[first.length];
		Arrays.fill(first, -1);
		Arrays.fill(last, -1);
		final IntList genes=new IntList(), models=new IntList(), draws=new IntList(), next=new IntList();
		final IntList rawStart=new IntList(), rawStop=new IntList(), trueStart=new IntList(), trueStop=new IntList();
		final ArrayList<byte[]> originals=new ArrayList<byte[]>();
		final LongHashSet seen=new LongHashSet();
		final LineParser1 p=new LineParser1('\t');
		String locus=null;
		int gene=-1;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, ROW_HEADER, path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==11, "Frozen outer sidecars require the original eleven fields");
				for(int i=0; i<11; i++){require(p.length(i)>0, "A frozen trial cannot contain empty identity fields");}
				if(locus==null || !p.termEquals(locus, 0)){
					locus=p.parseString(0);
					gene=population.index(locus);
					require(gene>=0 && first[gene]<0, "Frozen trials must contain known contiguous original-locus blocks");
				}
				require(population.excluded[gene] && !population.partial[gene] && p.termEquals("evaluation", 3)
					&& p.termEquals(population.sourceName(gene), 1) && p.termEquals(population.assemblyName(gene), 2),
					"Only the exact original excluded assembly/source population may enter fixed outer trials");
				final int model=p.parseInt(4), sample=p.parseInt(6);
				final int start=p.parseInt(7), stop=p.parseInt(8), truthStart=p.parseInt(9), truthStop=p.parseInt(10);
				require(model>=0 && model<modelNames.length && p.termEquals(modelNames[model], 5)
					&& sample>=0 && sample<samples && seen.add(((long)gene<<32)|(sample&0xffffffffL)),
					"Frozen model/sample identities must be valid and unique; model reselection is forbidden");
				require(start>=0 && stop>=start && truthStart>=0 && truthStop>=truthStart
					&& truthStart==population.leftFlank[gene] && (long)truthStop-truthStart+1==population.length[gene],
					"Frozen true bounds must retain the joined original context geometry");
				require(modelByGene[gene]<0 || modelByGene[gene]==model, "Every original locus retains its existing single selected model");
				modelByGene[gene]=model;
				final int row=originals.size();
				require(row<expected, "Frozen sidecar extends beyond its accepted original population");
				if(first[gene]<0){first[gene]=row;}else{next.set(last[gene], row);}
				last[gene]=row;
				genes.add(gene);
				models.add(model);
				draws.add(sample);
				next.add(-1);
				rawStart.add(start);
				rawStop.add(stop);
				trueStart.add(truthStart);
				trueStop.add(truthStop);
				originals.add(line.clone());
			}
		}
		require(originals.size()==expected, "Frozen outer sidecar must close at its exact accepted trial count");
		return new End(path, first, next.toArray(), genes.toArray(), models.toArray(), draws.toArray(),
			rawStart.toArray(), rawStop.toArray(), trueStart.toArray(), trueStop.toArray(), originals.toArray(new byte[originals.size()][]));
	}

	/** The two ends may have different availability, but a shared trial has one original span. */
	private void validateCommonTrials(){
		for(int gene=0; gene<modelByGene.length; gene++){
			for(int five=ends[0].first[gene]; five>=0; five=ends[0].next[five]){
				for(int three=ends[1].first[gene]; three>=0; three=ends[1].next[three]){
					if(ends[0].sample[five]!=ends[1].sample[three]){continue;}
					require(ends[0].rawStart[five]==ends[1].rawStart[three] && ends[0].rawStop[five]==ends[1].rawStop[three]
						&& ends[0].trueStart[five]==ends[1].trueStart[three] && ends[0].trueStop[five]==ends[1].trueStop[three],
						"Shared start/stop sample identities cannot describe different original spans or truth");
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Access        ----------------*/
	/*--------------------------------------------------------------*/

	/** Bind a frozen trial to its unchanged context and subtract medians in sense coordinates. */
	public void center(final int end, final int row, final NcrnaLearningContext context, final Span target){
		require(end>=0 && end<2 && row>=0 && row<ends[end].size() && context!=null && context.hasFeatures() && target!=null,
			"Centered features need a known frozen trial, usable original context and worker-owned scratch");
		final End trials=ends[end];
		require(context.gene==trials.gene[row] && context.start==trials.trueStart[row] && context.stop==trials.trueStop[row]
			&& trials.rawStop[row]<context.read.length(), "Frozen raw and true bounds must still belong to this exact original context");
		final int model=trials.model[row];
		target.start=(long)trials.rawStart[row]-median5[model];
		target.stop=(long)trials.rawStop[row]-median3[model];
		target.crossed=target.stop<target.start;
		target.outside=!target.crossed && (target.start<0 || target.stop>=context.read.length());
		target.tooShort=!target.crossed && !target.outside && target.stop-target.start+1<minimumLength;
		target.valid=!target.crossed && !target.outside && !target.tooShort;
	}

	/** True only for the exact population instance used to validate this resource. */
	public boolean belongsTo(final NcrnaLearningPopulation other, final NcrnaFamily family){return population==other && library==family;}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------        Immutable Tables      ----------------*/
	/*--------------------------------------------------------------*/

	/** Worker-owned geometric scratch; invalid spans become holds, never replacement draws. */
	public static final class Span{
		public long start, stop;
		public boolean valid, crossed, outside, tooShort;
	}

	/** Parallel primitive arrays retain the complete fixed denominator and its original bytes. */
	public static final class End{
		private End(final String source_, final int[] first_, final int[] next_, final int[] gene_, final int[] model_, final int[] sample_,
				final int[] rawStart_, final int[] rawStop_, final int[] trueStart_, final int[] trueStop_, final byte[][] original_){
			source=source_;
			first=first_;
			next=next_;
			gene=gene_;
			model=model_;
			sample=sample_;
			rawStart=rawStart_;
			rawStop=rawStop_;
			trueStart=trueStart_;
			trueStop=trueStop_;
			original=original_;
		}
		public int size(){return original.length;}
		public final String source;
		public final int[] first, next, gene, model, sample, rawStart, rawStop, trueStart, trueStop;
		public final byte[][] original;
	}

	private final NcrnaLearningPopulation population;
	/** Borrowed immutable family used for both frozen model validation and feature construction. */
	private final NcrnaFamily library;
	private final String member;
	private final String[] modelNames;
	public final End[] ends=new End[2];
	public final int[] median5, median3, modelByGene;
	public final long[] trainingSupport;
	public final boolean[] supported;
	public final int minimumLength;
	public final int supportMinimum;
	public static final String ROW_HEADER="locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0";
	public static final String AUDIT_HEADER=ROW_HEADER+"\toriginal_row1\tcenter_start0\tcenter_stop0\tmedian5_bp\tmedian3_bp\tdisposition"
		+"\ttraining_n\tminimum_n\toffset_supported\tminimum_length";
	public static final String OFFSET_HEADER="member\tmodel\ttraining_n\tminimum_n\tapplied_median5_bp\tapplied_median3_bp";
}
