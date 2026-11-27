package prok;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Counts pooled-member positional tables for several radii in one context pass.
 * Adapted from A_SampleStreamerMT and NcrnaLearningTables. One counting worker
 * owns the dense positional arrays, avoiding a copy of every table per CPU.
 * All original contexts join their source-qualified population before eligibility
 * is applied. Missing geometry is counted explicitly, never fabricated as bases.
 * Tables become immutable to consumers after the worker joins and finish runs.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalTables{

	/** Counts both inside-K9 ends at every supplied radius, without withholding.
	 * @param alpha Positive per-word/per-position pseudocount, not total prior mass
	 * @return Finalized tables and original-locus availability on the population axis
	 */
	public static NcrnaPositionalTables build(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final double alpha){
		return build(population, fasta, member, radii, 9, 0, 9, 0, alpha, alpha);
	}

	/** Counts one explicitly placed word per end, retaining total prior mass for
	 * each word length. This overload never changes the legacy per-key-alpha API.
	 */
	public static NcrnaPositionalTables build(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double prior){
		final int fiveK=wordLength(startInside, startOutside), threeK=wordLength(stopInside, stopOutside);
		if(!(prior>0 && Double.isFinite(prior))){throw new IllegalArgumentException("Total positional prior must be finite and positive");}
		return build(population, fasta, member, radii, startInside, startOutside, stopInside, stopOutside,
			prior/(1<<(2*fiveK)), prior/(1<<(2*threeK)));
	}

	private static NcrnaPositionalTables build(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double fiveAlpha, final double threeAlpha){
		return build(population, fasta, member, radii, startInside, startOutside,
			stopInside, stopOutside, fiveAlpha, threeAlpha, false);
	}

	/** Optional K7..9 probes at the declared placement and the exclusive left/right
	 * windows. All six tables count the same eligible loci in one context pass.
	 */
	public static NcrnaPositionalTables buildTriple(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double prior){
		final int k=tripleWordLength(startInside, startOutside, stopInside, stopOutside);
		if(!(prior>0 && Double.isFinite(prior))){
			throw new IllegalArgumentException("The optional three-probe layout requires equal K7..9 at both ends and positive total prior");
		}
		final double alpha=prior/(1<<(2*k));
		return build(population, fasta, member, radii, startInside, startOutside,
			stopInside, stopOutside, alpha, alpha, true);
	}

	/** Optional measured full-error background; correction windows remain radius ten. */
	public static NcrnaPositionalTables buildFull(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double prior, final NcrnaSignedEndpointErrors errors){
		final int k=tripleWordLength(startInside, startOutside, stopInside, stopOutside);
		if(errors==null || !member.equals(errors.member) || radii==null || radii.length!=1 || radii[0]!=10
				|| !(prior>0 && Double.isFinite(prior))){
			throw new IllegalArgumentException("Full-error inputs require this member's measured errors, equal K7..9 probes, radius ten and positive prior");
		}
		final double alpha=prior/(1<<(2*k));
		return build(population, fasta, member, radii, startInside, startOutside,
			stopInside, stopOutside, alpha, alpha, true, errors);
	}

	private static NcrnaPositionalTables build(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double fiveAlpha, final double threeAlpha,
			final boolean triple){
		return build(population, fasta, member, radii, startInside, startOutside,
			stopInside, stopOutside, fiveAlpha, threeAlpha, triple, null);
	}

	private static NcrnaPositionalTables build(final NcrnaLearningPopulation population, final String fasta,
			final String member, final int[] radii, final int startInside, final int startOutside,
			final int stopInside, final int stopOutside, final double fiveAlpha, final double threeAlpha,
			final boolean triple, final NcrnaSignedEndpointErrors errors){
		if(population==null || population.split().selected!=0 || radii==null || radii.length==0){
			throw new IllegalArgumentException("Pooled positional refit needs all-locus population and explicit positive radii");
		}
		final NcrnaPositionalTables result=new NcrnaPositionalTables(population, member, radii,
			startInside, startOutside, stopInside, stopOutside, fiveAlpha, threeAlpha, triple, errors);
		result.process(fasta);
		return result;
	}

	/** Dense positional storage is deliberately bounded; larger K needs a separate design. */
	public static int wordLength(final int inside, final int outside){
		final long k=(long)inside+outside;
		if(inside<0 || outside<0 || k<5 || k>9){
			throw new IllegalArgumentException("Single-table positional geometry needs nonnegative inside/outside counts summing to5..9");
		}
		return (int)k;
	}

	/** Both ends share one word length so the three-probe neural layout is explicit. */
	public static int tripleWordLength(final int startInside, final int startOutside,
			final int stopInside, final int stopOutside){
		final int k=wordLength(startInside, startOutside);
		if(k<7 || k!=wordLength(stopInside, stopOutside)){
			throw new IllegalArgumentException("Three-probe ends require the same K7..9; each probe's full word must fit its declared geometry");
		}
		return k;
	}

	private NcrnaPositionalTables(final NcrnaLearningPopulation population_, final String member_, final int[] radii_,
			final int startInside, final int startOutside, final int stopInside, final int stopOutside,
			final double fiveAlpha, final double threeAlpha, final boolean triple, final NcrnaSignedEndpointErrors errors){
		population=population_;
		member=member_;
		radii=radii_.clone();
		five=new RrnaPositionalKmerTable.Table[radii.length];
		three=new RrnaPositionalKmerTable.Table[radii.length];
		fiveProbes=triple ? new RrnaPositionalKmerTable.Table[radii.length][3] : null;
		threeProbes=triple ? new RrnaPositionalKmerTable.Table[radii.length][3] : null;
		fullErrors=errors;
		fiveWeighted=errors==null ? null : new NcrnaWeightedEndpointTable[radii.length][3];
		threeWeighted=errors==null ? null : new NcrnaWeightedEndpointTable[radii.length][3];
		for(int r=0; r<radii.length; r++){
			if(radii[r]<1 || (r>0 && radii[r]<=radii[r-1])){
				throw new IllegalArgumentException("Positional radii must be positive, unique and increasing");
			}
			five[r]=new RrnaPositionalKmerTable.Table(member, "5prime", wordLength(startInside, startOutside), -startOutside, fiveAlpha, radii[r]);
			three[r]=new RrnaPositionalKmerTable.Table(member, "3prime", wordLength(stopInside, stopOutside), 1-stopInside, threeAlpha, radii[r]);
			if(triple){
				fiveProbes[r][0]=five[r];
				threeProbes[r][0]=three[r];
				for(int probe=1; probe<3; probe++){
					final int k=five[r].k, offset=probe==1 ? -k : 1;
					fiveProbes[r][probe]=new RrnaPositionalKmerTable.Table(member, "5prime", k, offset, fiveAlpha, radii[r]);
					threeProbes[r][probe]=new RrnaPositionalKmerTable.Table(member, "3prime", k, offset, threeAlpha, radii[r]);
				}
				if(errors!=null){
					for(int probe=0; probe<3; probe++){
						fiveWeighted[r][probe]=new NcrnaWeightedEndpointTable(fiveProbes[r][probe], errors);
						threeWeighted[r][probe]=new NcrnaWeightedEndpointTable(threeProbes[r][probe], errors);
					}
				}
			}
		}
		statusByGene=new byte[population.ids.length];
		Arrays.fill(statusByGene, (byte)-1);
	}

	/** The sole counter consumes native batches; the owner joins it before closing input. */
	private void process(final String fasta){
		final FileFormat format=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		if(!format.fasta()){throw new IllegalArgumentException("Positional contexts must use native FASTA input");}
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		try(NcrnaWorkers workers=new NcrnaWorkers();
				SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(format, null, false, -1, false, true, 1))){
			final Streamer input=owner.stream;
			final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(1);
			final ProcessThread counter=new ProcessThread(input, this);
			threads.add(counter);
			input.start();
			workers.startAndJoin(threads, input);
			if(!counter.success){throw new IllegalStateException("Positional context counting did not reach terminal input");}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		if(records!=population.ids.length){throw new IllegalArgumentException("Every original locus must have exactly one positional context");}
		for(final byte status:statusByGene){
			if(status<0){throw new IllegalArgumentException("An original locus is missing its positional context");}
		}
		assert(records==filtered+learning && learning==counted+noGeometry)
			: "All original contexts partition into filtered or eligible; eligible contexts have geometry or an explicit limitation";
		for(int r=0; r<radii.length; r++){
			five[r].finish();
			three[r].finish();
			if(triple()){
				for(int probe=1; probe<3; probe++){
					fiveProbes[r][probe].finish();
					threeProbes[r][probe].finish();
					assert(fiveProbes[r][probe].records==counted && threeProbes[r][probe].records==counted)
						: "Every extra probe must count the same eligible geometry-complete loci as the primary table";
				}
			}
			if(fullErrors!=null){
				for(int probe=0; probe<3; probe++){
					fiveWeighted[r][probe].finish();
					threeWeighted[r][probe].finish();
				}
			}
			assert(five[r].records==counted && three[r].records==counted)
				: "Every radius/end must be counted on the same eligible geometry-complete source loci";
		}
	}

	/** Native input is borrowed until join; all table mutation belongs to this one worker. */
	private static final class ProcessThread extends Thread{
		ProcessThread(final Streamer input_, final NcrnaPositionalTables owner_){
			input=input_;
			owner=owner_;
		}
		@Override
		public void run(){
			try{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(final Read read:batch){processRead(read);}
				}
				success=true;
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		}
		private void processRead(final Read read){
			context.set(read, owner.population);
			final int gene=context.gene;
			if(owner.statusByGene[gene]>=0){throw new IllegalArgumentException("Duplicate positional context: "+owner.population.ids[gene]);}
			owner.statusByGene[gene]=(byte)context.status;
			owner.availability[context.status]++;
			owner.records++;
			if(!owner.population.eligible(gene)){
				owner.filtered++;
				return;
			}
			owner.learning++;
			if(!context.hasGeometry()){
				owner.noGeometry++;
				return;
			}
			owner.counted++;
			for(int r=0; r<owner.radii.length; r++){
				owner.five[r].add(read.bases, context.start);
				owner.three[r].add(read.bases, context.stop);
				if(owner.triple()){
					for(int probe=1; probe<3; probe++){
						owner.fiveProbes[r][probe].add(read.bases, context.start);
						owner.threeProbes[r][probe].add(read.bases, context.stop);
					}
				}
				if(owner.fullErrors!=null){
					for(int probe=0; probe<3; probe++){
						owner.fiveWeighted[r][probe].add(read.bases, context.start);
						owner.threeWeighted[r][probe].add(read.bases, context.stop);
					}
				}
			}
		}
		private final Streamer input;
		private final NcrnaPositionalTables owner;
		private final NcrnaLearningContext context=new NcrnaLearningContext();
		boolean success;
	}

	/** Tables and vectors must share the exact original-locus population axis. */
	boolean belongsTo(final NcrnaLearningPopulation candidate){return population==candidate;}
	/** Whether the optional three-table, LOO and HBM-feature layout was requested. */
	public boolean triple(){return fiveProbes!=null;}
	/** Feature width declared by the native matrix header, excluding output classes. */
	public int inputs(final int r){return triple() ? Math.addExact(Math.multiplyExact(3, five[r].sites), 6) : five[r].inputs;}

	/** Finalized tables use the same radius axis and complete geometry-locus population. */
	public final RrnaPositionalKmerTable.Table[] five, three;
	/** Optional radius/probe axes; primary entries alias five/three, never recount them. */
	public final RrnaPositionalKmerTable.Table[][] fiveProbes, threeProbes;
	/** Null preserves legacy uniform-background counting and score semantics. */
	public final NcrnaSignedEndpointErrors fullErrors;
	public final NcrnaWeightedEndpointTable[][] fiveWeighted, threeWeighted;
	public final int[] radii;
	public final String member;
	/** All original-locus context statuses, including excluded/unavailable records. */
	public final byte[] statusByGene;
	public final long[] availability=new long[4];
	public long records, filtered, learning, counted, noGeometry;
	private final NcrnaLearningPopulation population;
}
