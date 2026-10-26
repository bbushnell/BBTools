package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import ukmer.HollowKmerLookup;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Reconstructs the16 ingress/next-base counts for retained right-extension seed edits.
 * Original seeds and original reads supply evidence; saved extensions never do.
 * This is a post-hoc diagnostic, not an execution trace or a correction policy.
 * @author Brian Bushnell, Fischl */
public final class HollowRescueChoices {

	/** Expands native configuration arguments and closes redirected diagnostics. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowRescueChoices.class, false);
		try{new HollowRescueChoices(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps evidence and ordered comparison inputs explicit and separate. */
	private HollowRescueChoices(final String[] args, final PrintStream log_){
		assert(args!=null) : "PreParser supplies expanded diagnostic arguments.";
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("seeds")){seeds=b;}
			else if(a.equals("queries")){queries=b;}
			else if(a.equals("stats")){stats=b;}
			else if(a.equals("k")){k=Parse.parseIntKMG(b);}
			else if(a.equals("support")){support=Parse.parseIntKMG(b);}
			else if(a.equals("test")){test=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		Parser.processQuality();
		reads=parser.in1; output=parser.out1;
		require(reads!=null && k>64 && support>0, "Specify in=original_reads, Hollow k>64 and positive support.");
		require(parser.in2==null && parser.out2==null && parser.maxReads<0 && !parser.overwrite && !parser.append,
				"Use complete unpaired inputs and new reports, without sampling or append.");
		Tools.testInputFiles(true, true, reads);
		if(!test){
			require(seeds!=null && queries!=null && output!=null && stats!=null, "Specify seeds, queries, out and stats.");
			Tools.testInputFiles(true, true, seeds, queries);
			Tools.testForDuplicateFiles(true, reads, seeds, queries, output, stats);
			require(Tools.testOutputFiles(false, false, false, output, stats), "Reports must be new writable files.");
			require(FileFormat.testInput(reads, FileFormat.FASTQ, null, true, true).fastq(), "Read evidence must be original FASTQ.");
		}
		FastaReadInputStream.SPLIT_READS=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads evidence once and streams each original/proposed pair without retaining libraries. */
	private void process(){
		if(test){selfTest(); return;}
		final KmerTableSetU tables=load(k, false);
		final HollowKmerLookup lookup=new HollowKmerLookup(tables);
		final int[] counts=new int[16];
		final ByteBuilder row=new ByteBuilder();
		long examined=0, changed=0, competitors=0, notMaximum=0;
		try(Cursor originals=new Cursor(seeds); Cursor proposals=new Cursor(queries); Report report=new Report(output)){
			row.append("id\tseed_length\tposition0\tcheckpoint\toriginal_base\tproposed_base\tnext_base\tcontext\tchosen_count\toriginal_max\tother_max\tproposed_is_max\tsupported_competitor");
			for(int i=0; i<16; i++){
				row.tab().append("count_").append((char)AminoAcid.numberToBase[i>>>2]).append((char)AminoAcid.numberToBase[i&3]);
			}
			report.writer.print(row.nl());
			for(Read seed=originals.next(); seed!=null; seed=originals.next()){
				final Read candidate=proposals.next();
				require(candidate!=null && seed.mate==null && candidate.mate==null && seed.id.equals(candidate.id),
						"Missing, paired or reordered proposed record for "+seed.id);
				require(seed.id.indexOf('\t')<0 && seed.length()>=k, "Use tab-free IDs and seeds at least k bases long.");
				final int position=changedPosition(seed.bases, candidate.bases);
				examined++;
				if(position<0){continue;}
				final int end=checkpoint(seed.length(), candidate.bases, position, k);
				final long left=arm(candidate.bases, end-k), right=arm(candidate.bases, end-32);
				lookup.fillRightCounts16(left, right, counts);
				final int original=AminoAcid.baseToNumber[seed.bases[position]];
				final int proposed=AminoAcid.baseToNumber[candidate.bases[position]];
				final int next=AminoAcid.baseToNumber[candidate.bases[end]];
				final int chosen=counts[(proposed<<2)|next], other=otherMaximum(counts, proposed);
				int maximum=0, originalMaximum=0;
				for(int i=0; i<16; i++){maximum=Math.max(maximum, counts[i]);}
				for(int i=0; i<4; i++){originalMaximum=Math.max(originalMaximum, counts[(original<<2)|i]);}
				final boolean isMaximum=chosen==maximum, competitor=other>=support;
				changed++; if(competitor){competitors++;} if(!isMaximum){notMaximum++;}
				row.clear().append(seed.id).tab().append(seed.length()).tab().append(position).tab().append(end)
						.tab().append((char)seed.bases[position]).tab().append((char)candidate.bases[position])
						.tab().append((char)candidate.bases[end]).tab().append(lookup.count(left, right))
						.tab().append(chosen).tab().append(originalMaximum).tab().append(other)
						.tab().append(isMaximum ? 1 : 0).tab().append(competitor ? 1 : 0);
				for(int count : counts){row.tab().append(count);}
				report.writer.print(row.nl());
			}
			require(proposals.next()==null, "Extra proposed records have no matching original seed.");
		}finally{tables.clear();}
		require(examined>0, "No original seeds were examined.");
		try(Report report=new Report(stats)){
			row.clear().append("k\tsupport\tseeds\tchanges\tsupported_competitors\tchosen_not_maximum\n")
					.append(k).tab().append(support).tab().append(examined).tab().append(changed)
					.tab().append(competitors).tab().append(notMaximum).nl();
			report.writer.print(row);
		}
		log.println("HOLLOW_RESCUE_CHOICES_PASS seeds="+examined+" changes="+changed+" k="+k);
	}

	/** Matches unfiltered extension evidence: every original seed and read counted once. */
	private KmerTableSetU load(final int span, final boolean tiny){
		assert(span>64) : "Endpoint-pair counts require a non-solid Hollow span.";
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList(
				"in="+reads, "k="+span, "hollow=t", "hashonly=t", "packed=t", "rcomp=t",
				"prefilter=f", "prealloc=f", "minprob=0", "minprobmain=f", "qtrim=f", "maq=0", "ecco=f",
				"merge=f", "reads=-1", "onepass=f", "interleaved=f", "t="+Shared.threads()));
		if(!tiny){options.add("extra="+seeds);}
		if(tiny){options.add("ways=7"); options.add("initialsize=1009");}
		final KmerTableSetU tables=new KmerTableSetU(options.toArray(new String[0]), 0);
		tables.process(new Timer());
		require(tables.hollow && tables.kbig==span, "Evidence must retain the exact requested Hollow geometry.");
		return tables;
	}

	/** Finds the sole retained seed substitution; appended bases are not seed edits. */
	private static int changedPosition(final byte[] original, final byte[] proposed){
		require(original!=null && proposed!=null && proposed.length>=original.length, "Proposed output must retain the original seed.");
		int position=-1;
		for(int i=0; i<original.length; i++){
			final int a=AminoAcid.baseToNumber[original[i]], b=AminoAcid.baseToNumber[proposed[i]];
			require(a>=0 && b>=0, "Seed comparison requires defined bases.");
			if(a==b){continue;}
			require(position<0, "More than one seed edit cannot identify a single recovery transaction.");
			position=i;
		}
		return position;
	}

	/** Inverts Tadpole2's position=end-k+32; a committed trial appended at least one base. */
	private static int checkpoint(final int seedLength, final byte[] bases, final int position, final int span){
		final long end=(long)position+span-32;
		require(position>=0 && position<seedLength && end>=seedLength && end>=span && end<bases.length,
				"Seed edit is incompatible with one right-extension rescue: position="+position+", checkpoint="+end);
		for(int i=0; i<=end; i++){require(AminoAcid.baseToNumber[bases[i]]>=0, "Rescue context and next base must be defined.");}
		return (int)end;
	}

	/** Packs one counted anchor; the edited ingress is outside both pre-step anchors. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start+32<=bases.length) : "Checkpoint geometry must provide32existing bases per anchor.";
		long value=0;
		for(int i=start; i<start+32; i++){value=(value<<2)|AminoAcid.baseToNumber[bases[i]];}
		return value;
	}

	/** Excludes all four successors of the proposed ingress, not just its chosen next base. */
	private static int otherMaximum(final int[] counts, final int proposed){
		assert(counts.length==16 && proposed>=0 && proposed<4) : "Matrix layout is (ingress<<2)|next, in query orientation.";
		int maximum=0;
		for(int i=0; i<16; i++){if((i>>>2)!=proposed){maximum=Math.max(maximum, counts[i]);}}
		return maximum;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Tests             ----------------*/
	/*--------------------------------------------------------------*/

	/** Cross-checks all16 direct packed queries, both orientations and checkpoint geometry. */
	private void selfTest(){
		final Random random=new Random(421130);
		int checks=0;
		for(int span : new int[]{96, 128, 192}){
			final KmerTableSetU tables=load(span, true);
			try{
				final HollowKmerLookup lookup=new HollowKmerLookup(tables);
				final Kmer word=new Kmer(span, true);
				for(int strand=0; strand<2; strand++){
					final int seedLength=span+40, end=span+50, position=82;
					final byte[] proposed=new byte[end+8];
					for(int i=0; i<proposed.length; i++){proposed[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					if(strand==1){AminoAcid.reverseComplementBasesInPlace(proposed);}
					final byte[] original=Arrays.copyOf(proposed, seedLength);
					original[position]=AminoAcid.numberToBase[(AminoAcid.baseToNumber[proposed[position]]+1)&3];
					final byte[] before=proposed.clone();
					require(changedPosition(original, proposed)==position && checkpoint(seedLength, proposed, position, span)==end,
							"Right-rescue geometry did not recover its original edit/checkpoint.");
					for(int i=0; i<16; i++){
						word.clear();
						for(int p=end-span; p<end; p++){word.addRight(proposed[p]);}
						word.substituteBase(32, i>>>2); word.addRightNumeric(i&3);
						tables.tables()[word.mod(tables.ways)].set(word, i+1);
					}
					final int[] counts=new int[16];
					lookup.fillRightCounts16(arm(proposed, end-span), arm(proposed, end-32), counts);
					for(int i=0; i<16; i++){require(counts[i]==i+1, "Endpoint matrix disagrees with packed substitution/shift at "+i); checks++;}
					for(int allele=0; allele<4; allele++){
						require(otherMaximum(counts, allele)==(allele==3 ? 12 : 16), "Competitor maximum included the proposed ingress."); checks++;
					}
					require(Arrays.equals(before, proposed) && changedPosition(original, original)==-1,
							"Read-only profiling mutated sequence or reported an unchanged seed."); checks+=2;
					final byte[] twice=original.clone(); twice[0]=AminoAcid.numberToBase[(AminoAcid.baseToNumber[twice[0]]+1)&3];
					boolean rejected=false;
					try{changedPosition(twice, proposed);}catch(IllegalArgumentException expected){rejected=true;}
					require(rejected, "Two seed changes must not be treated as one transaction."); checks++;
					rejected=false;
					try{checkpoint(seedLength, Arrays.copyOf(proposed, end), position, span);}catch(IllegalArgumentException expected){rejected=true;}
					require(rejected, "A retained correction must have a saved successor after its checkpoint."); checks++;
				}
			}finally{tables.clear();}
		}
		log.println("HOLLOW_RESCUE_CHOICES_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------            I/O               ----------------*/
	/*--------------------------------------------------------------*/

	/** Native synchronous FASTA iteration avoids a producer blocked after an input error. */
	private static final class Cursor implements AutoCloseable {
		/** Opens one unpaired comparison stream without background queue dependencies. */
		Cursor(final String path){
			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			require(ff.fasta(), "Ordered comparison inputs must be FASTA.");
			stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 0); stream.start();
		}
		/** Returns each record once and retains EOF across calls. */
		Read next(){
			if(finished){return null;}
			while(batch==null || index>=batch.list.size()){
				batch=stream.nextList(); index=0;
				if(batch==null){finished=true; return null;}
			}
			return batch.list.get(index++);
		}
		@Override
		public void close(){require(!ReadWrite.closeStream(stream), "Comparison input failed.");}
		private final Streamer stream;
		private ListNum<Read> batch;
		private int index;
		private boolean finished;
	}

	/** Flushes the native writer on success and failure without closing other outputs. */
	private static final class Report implements AutoCloseable {
		/** Starts one new diagnostic report. */
		Report(final String path){
			writer=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TXT, null, true, false, false, false)); writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Choice report output failed.");}
		final ByteStreamWriter writer;
	}

	/** Keeps scientific input checks active even without assertions. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String reads, seeds, queries, output, stats;
	private int k=192, support=4;
	private boolean test;
	private final PrintStream log;
}
