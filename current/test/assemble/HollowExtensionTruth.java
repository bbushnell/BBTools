package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.function.Consumer;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Audits synthetic inputs and scores only appended bases at Shred's exact source coordinates.
 * No alignment, best-hit selection or mutation truth is used by the assembler.
 * @author Brian Bushnell, Fischl */
public final class HollowExtensionTruth {

	/** Runs audit or grade mode using native sequence readers and tabular output. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowExtensionTruth.class, false);
		try{new HollowExtensionTruth(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses flags without accepting a mistyped truth/input path. */
	private HollowExtensionTruth(final String[] args, final PrintStream log_){
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase(Locale.ROOT), value=split.length>1 ? split[1] : null;
			if(key.equals("ref")){reference=value;}
			else if(key.equals("mutant")){mutant=value;}
			else if(key.equals("shreds")){shreds=value;}
			else if(key.equals("mode")){mode=value;}
			else if(key.equals("label")){label=value;}
			else if(key.equals("details")){details=value;}
			else if(key.equals("maxextension")){maxExtension=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		Parser.processQuality();
		input=parser.in1;
		output=parser.out1;
		require(reference!=null && input!=null && output!=null && shreds!=null, "Specify ref, in, out and shreds.");
		require(!parser.overwrite && !parser.append, "Truth reports require new output files.");
		require(mode.equals("audit") || mode.equals("grade"), "Use mode=audit or grade.");
		require(!mode.equals("audit") || mutant!=null, "Audit requires the substitution-mutated genome.");
		require(!mode.equals("grade") || maxExtension>0, "Grade requires a positive maxextension.");
		FastaReadInputStream.SPLIT_READS=false;
	}

	/** Resolves references and shreds once, then streams the much larger synthetic read input. */
	private void process(){
		final ArrayList<Read> refs=new ArrayList<Read>();
		final Map<String,Read> byName=new HashMap<String,Read>();
		walk(reference, r->{
			require(byName.put(r.id, r)==null, "Duplicate reference ID: "+r.id);
			refs.add(r);
		});
		require(!refs.isEmpty(), "Reference is empty.");
		final Map<String,Seed> seeds=new HashMap<String,Seed>();
		walk(shreds, r->{
			final int separator=r.id.lastIndexOf('_'), dash=r.id.lastIndexOf('-');
			require(separator>0 && dash>separator, "Not a Shred origin label: "+r.id);
			final Read ref=byName.get(r.id.substring(0, separator));
			require(ref!=null, "Shred reference ID is missing: "+r.id);
			final int start=Parse.parseInt(r.id, separator+1, dash);
			final int stop=Parse.parseInt(r.id, dash+1, r.id.length());
			require(r.length()==200 && stop-start+1==200 && start>=0 && stop<ref.length(),
					"Expected an in-bounds 200bp shred: "+r.id);
			final Seed seed=new Seed(r, ref, start, seeds.size());
			require(seeds.put(r.id, seed)==null, "Duplicate shred ID: "+r.id);
		});
		require(!seeds.isEmpty(), "No full-length shreds.");
		if(mode.equals("audit")){audit(refs, byName, seeds);}
		else{grade(seeds);}
	}

	/** Verifies mutation-only coordinates, shred sequence, Q20 and empirical read error rate. */
	private void audit(final ArrayList<Read> refs, final Map<String,Read> byName, final Map<String,Seed> seeds){
		final Map<String,Read> mutated=new HashMap<String,Read>();
		final long[] genome={0, 0, 0};
		walk(mutant, r->{
			final Read ref=byName.get(r.id);
			require(ref!=null && r.length()==ref.length(), "Mutation changed reference IDs or lengths.");
			require(mutated.put(r.id, r)==null, "Duplicate mutant ID.");
			genome[0]+=r.length();
			for(int i=0; i<r.length(); i++){
				if(!defined(ref.bases[i]) || !defined(r.bases[i])){genome[2]++;}
				else if(ref.bases[i]!=r.bases[i]){genome[1]++;}
			}
		});
		require(mutated.size()==refs.size(), "Mutant omitted a contig.");
		long expectedSeeds=0;
		for(Read ref : refs){expectedSeeds+=ref.length()/200;}
		require(seeds.size()==expectedSeeds, "Nonoverlapping shredding lost or duplicated a full-length piece.");
		for(Seed seed : seeds.values()){
			final byte[] source=mutated.get(seed.reference.id).bases;
			for(int i=0; i<200; i++){require(seed.read.bases[i]==source[seed.start+i], "Shred sequence does not match its mutant coordinates.");}
		}
		final long[] reads={0, 0, 0, 0, 0};
		walk(input, r->{
			require(r.length()==250 && r.quality!=null, "Synthetic input must be 250bp FASTQ.");
			final int file=header(r.id, "f_"), contig=header(r.id, "_c_");
			final int strand=header(r.id, "_s_"), start=header(r.id, "_p_");
			require(file==0 && contig>=0 && contig<refs.size() && (strand==0 || strand==1)
					&& header(r.id, "_r_")==250, "Unexpected RandomReadsMG origin.");
			final byte[] ref=refs.get(contig).bases;
			require(start>=0 && start+250<=ref.length, "Read origin exceeds its reference.");
			reads[0]++;
			reads[1]+=250;
			for(int i=0; i<250; i++){
				final byte truth=strand==0 ? ref[start+i] : AminoAcid.baseToComplementExtended[ref[start+249-i]];
				if(!defined(truth) || !defined(r.bases[i])){reads[3]++; continue;}
				require(r.quality[i]==20, "Defined synthetic bases must all have Q20.");
				reads[4]++;
				if(r.bases[i]!=truth){reads[2]++;}
			}
		});
		require(reads[0]>0 && genome[0]>0, "Audit must see reads and reference bases.");
		final ByteBuilder bb=new ByteBuilder();
		bb.append("metric\tvalue\n");
		row(bb, "reference_bases", genome[0]);
		row(bb, "genome_substitutions", genome[1]);
		row(bb, "genome_ambiguous_bases", genome[2]);
		row(bb, "shreds", seeds.size());
		row(bb, "reads", reads[0]);
		row(bb, "read_bases", reads[1]);
		row(bb, "read_substitutions", reads[2]);
		row(bb, "read_unscorable_bases", reads[3]);
		row(bb, "read_defined_q20_bases", reads[4]);
		bb.append("observed_depth\t").append(reads[1]/(double)genome[0], 8).nl();
		bb.append("observed_substitution_rate\t").append(reads[2]/(double)reads[4], 8).nl();
		write(output, bb);
		log.println("HOLLOW_INPUT_AUDIT_PASS reads="+reads[0]+" shreds="+seeds.size());
	}

	/** Scores suffixes only; rejects missing/duplicate IDs and any alteration to the seed. */
	private void grade(final Map<String,Seed> seeds){
		final Counts[] counts={new Counts(), new Counts(), new Counts(), new Counts(), new Counts()};
		final BitSet seen=new BitSet(seeds.size());
		final Score score=new Score();
		final ByteStreamWriter errors=details==null ? null : writer(details);
		final ByteBuilder detail=new ByteBuilder();
		if(errors!=null){errors.print("id\tseed_class\textended_bases\tmismatches\tbeyond_reference\tunscorable\n");}
		try{
			walk(input, r->{
				final Seed seed=seeds.get(r.id);
				require(seed!=null, "Unknown extended read ID: "+r.id);
				require(!seen.get(seed.index), "Duplicate extended read ID: "+r.id);
				seen.set(seed.index);
				score(r.bases, seed.read.bases, seed.reference.bases, seed.start, maxExtension, score);
				counts[0].add(score);
				counts[seed.group].add(score);
				if(errors!=null && score.bad+score.beyond+score.unknown>0){
					detail.clear().append(r.id).tab().append(GROUPS[seed.group]).tab().append(score.added).tab()
							.append(score.bad).tab().append(score.beyond).tab().append(score.unknown).nl();
					errors.print(detail);
				}
			});
		}finally{
			if(errors!=null){require(!errors.poisonAndWait(), "Could not finish error details.");}
		}
		require(seen.cardinality()==seeds.size(), "An output omitted shreds.");
		final ByteBuilder bb=new ByteBuilder();
		bb.append("arm\tgroup\tsequences\textended_sequences\textended_bases\tcorrect_bases\tmismatches")
				.append("\tbeyond_reference\tunscorable_bases\terror_free_extended_sequences\tincorrect_extended_sequences\tunscorable_only_sequences\n");
		for(int i=0; i<counts.length; i++){
			final Counts c=counts[i];
			bb.append(label).tab().append(GROUPS[i]).tab().append(c.sequences).tab().append(c.extended).tab()
					.append(c.added).tab().append(c.correct).tab().append(c.bad).tab().append(c.beyond).tab()
					.append(c.unknown).tab().append(c.exact).tab().append(c.incorrect).tab().append(c.uncertain).nl();
		}
		write(output, bb);
		log.println("HOLLOW_EXTENSION_GRADE_PASS sequences="+seen.cardinality());
	}

	/** No alignment is needed: mutations are substitutions and extensions are right-only. */
	static void score(final byte[] sequence, final byte[] seed, final byte[] ref,
			final int start, final int maximum, final Score score){
		require(start>=0 && start+seed.length<=ref.length, "Seed origin is outside reference.");
		require(sequence.length>=seed.length && sequence.length<=seed.length+maximum, "Output length violates right-extension limits.");
		for(int i=0; i<seed.length; i++){require(sequence[i]==seed[i], "Original shred bases were edited.");}
		score.added=sequence.length-seed.length;
		score.bad=score.beyond=score.unknown=0;
		for(int i=seed.length; i<sequence.length; i++){
			if(start+i>=ref.length){score.beyond++;}
			else if(!defined(ref[start+i]) || !defined(sequence[i])){score.unknown++;}
			else if(sequence[i]!=ref[start+i]){score.bad++;}
		}
		assert(score.bad+score.beyond+score.unknown<=score.added) : "Every appended base belongs to exactly one scoring category.";
	}

	/** Locates numeric RandomReadsMG fields without per-read split arrays or substrings. */
	private static int header(final String id, final String marker){
		final int position=id.indexOf(marker);
		if(position<0){throw new IllegalArgumentException("Missing RandomReadsMG field: "+marker);}
		final int start=position+marker.length(), stop=id.indexOf('_', start);
		if(stop<=start){throw new IllegalArgumentException("Missing field terminator: "+marker);}
		return Parse.parseInt(id, start, stop);
	}

	/** Shares the native streaming lifecycle; small references/shreds may be retained by consumers. */
	private static void walk(final String path, final Consumer<Read> consumer){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		final Streamer stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		stream.start();
		try{
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read r : batch){
					require(r.mate==null, "This benchmark uses single-end evidence and forward shreds.");
					consumer.accept(r);
				}
			}
		}finally{require(!ReadWrite.closeStream(stream), "Input stream failed: "+path);}
	}

	/** Opens a non-overwriting native tabular writer. */
	private static ByteStreamWriter writer(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, true);
		writer.start();
		return writer;
	}

	/** Writes a small complete report and checks writer completion. */
	private static void write(final String path, final ByteBuilder data){
		final ByteStreamWriter writer=writer(path);
		writer.print(data);
		require(!writer.poisonAndWait(), "Output stream failed: "+path);
	}

	/** Appends an exact integer metric rather than rounding a count through floating point. */
	private static void row(final ByteBuilder bb, final String name, final long value){bb.append(name).tab().append(value).nl();}
	private static boolean defined(final byte b){return AminoAcid.isFullyDefined(b);}
	static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Coordinates and error class are derived from the original shred, not its extension. */
	private static final class Seed {
		Seed(final Read read_, final Read reference_, final int start_, final int index_){
			read=read_;
			reference=reference_;
			start=start_;
			index=index_;
			boolean middle=false, anchor=false, ambiguous=false;
			for(int i=72; i<200; i++){
				final byte a=read.bases[i], b=reference.bases[start+i];
				if(!defined(a) || !defined(b)){ambiguous=true;}
				else if(a!=b){
					if(i<104 || i>=168){anchor=true;}else{middle=true;}
				}
			}
			group=ambiguous ? 4 : anchor ? 3 : middle ? 2 : 1;
		}
		final Read read, reference;
		final int start, index, group;
	}

	/** Reused per-sequence score, keeping the hot loop allocation-free. */
	static final class Score {
		int added, bad, beyond, unknown;
	}

	/** Non-overlapping suffix categories and complete-sequence accuracy totals. */
	private static final class Counts {
		void add(final Score s){
			assert(s.added>=s.bad+s.beyond+s.unknown) : "Suffix categories must partition appended bases.";
			sequences++;
			added+=s.added;
			bad+=s.bad;
			beyond+=s.beyond;
			unknown+=s.unknown;
			correct+=s.added-s.bad-s.beyond-s.unknown;
			if(s.added>0){
				extended++;
				if(s.bad+s.beyond>0){incorrect++;}
				else if(s.unknown>0){uncertain++;}
				else{exact++;}
			}
		}
		long sequences, extended, added, correct, bad, beyond, unknown, exact, incorrect, uncertain;
	}

	private String reference, mutant, shreds, input, output, details;
	private String mode="grade", label="arm";
	private int maxExtension=-1;
	private final PrintStream log;
	private static final String[] GROUPS={"all", "clean_seed", "middle_only_seed_error", "anchor_seed_error", "ambiguous_seed"};
}
