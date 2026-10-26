package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import assemble.HollowSeedGuard;
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
import ukmer.AbstractKmerTableU;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Measures original-read Solid64/96 evidence around saved single-base seed edits.
 * Seeds, proposed extensions and reference truth never contribute counts. This
 * diagnostic does not accept, reject or apply a correction.
 * @author Brian Bushnell, Fischl */
public final class HollowSeedSupport {

	/** Loads one exact count table and streams paired original/proposed records. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowSeedSupport.class, false);
		try{new HollowSeedSupport(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses native common flags; comparison inputs are explicit, not paired reads. */
	private HollowSeedSupport(final String[] args, final PrintStream log_){
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
		reads=parser.in1;
		output=parser.out1;
		require(reads!=null && (k==64 || k==96) && support>=1, "Specify in=original_reads, k=64 or96, and positive support.");
		require(parser.in2==null && parser.out2==null && parser.maxReads<0 && !parser.overwrite && !parser.append,
				"Use one complete original read input and fresh reports, without paired or sampled inputs.");
		Tools.testInputFiles(true, true, reads);
		if(!test){
			require(seeds!=null && queries!=null && output!=null && stats!=null, "Specify seeds, queries, out and stats.");
			Tools.testInputFiles(true, true, seeds, queries);
			Tools.testForDuplicateFiles(true, reads, seeds, queries, output, stats);
			require(Tools.testOutputFiles(false, false, false, output, stats), "Reports must be new writable files.");
			require(FileFormat.testInput(reads, FileFormat.FASTQ, null, true, true).fastq(), "Evidence must be the original FASTQ.");
		}
		FastaReadInputStream.SPLIT_READS=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Does not retain seed libraries or write per-window sequence profiles. */
	private void process(){
		if(test){selfTest(); return;}
		final KmerTableSetU tables=load(k, false);
		final Kmer originalWord=new Kmer(k, false), proposedWord=new Kmer(k, false);
		final int[] before=new int[3], after=new int[3];
		final ByteBuilder row=new ByteBuilder();
		long examined=0, changed=0, originalSupported=0, proposedSupported=0, bothSupported=0;
		try(Cursor originals=new Cursor(seeds); Cursor proposals=new Cursor(queries); Report report=new Report(output)){
			report.writer.print("id\tseed_length\tposition0\toriginal_base\tproposed_base\twindows\toriginal_min\toriginal_max\toriginal_supported_windows\tproposed_min\tproposed_max\tproposed_supported_windows\n");
			for(Read seed=originals.next(); seed!=null; seed=originals.next()){
				final Read candidate=proposals.next();
				require(candidate!=null && seed.mate==null && candidate.mate==null && seed.id.equals(candidate.id),
						"Missing, paired or reordered proposed record for "+seed.id);
				require(seed.id.indexOf('\t')<0 && seed.length()>=k, "Reports require tab-free IDs and seeds at least k bases long.");
				final int position=changedPosition(seed.bases, candidate.bases);
				examined++;
				if(position<0){continue;}
				final int windows=scan(seed.bases, seed.length(), position, tables, originalWord, support, before);
				require(scan(candidate.bases, seed.length(), position, tables, proposedWord, support, after)==windows,
						"Original and proposed support must cover the same seed-only window coordinates.");
				changed++;
				if(before[2]>0){originalSupported++;}
				if(after[2]>0){proposedSupported++;}
				if(before[2]>0 && after[2]>0){bothSupported++;}
				row.clear().append(seed.id).tab().append(seed.length()).tab().append(position)
						.tab().append((char)seed.bases[position]).tab().append((char)candidate.bases[position]).tab().append(windows);
				for(int n : before){row.tab().append(n);}
				for(int n : after){row.tab().append(n);}
				report.writer.print(row.nl());
			}
			require(proposals.next()==null, "Proposed output has records missing from the original seed stream.");
		}finally{tables.clear();}
		require(examined>0, "No seed records were examined.");
		try(Report report=new Report(stats)){
			row.clear().append("k\tsupport\tseeds\tchanged_seeds\toriginal_supported\tproposed_supported\tboth_supported\n")
					.append(k).tab().append(support).tab().append(examined).tab().append(changed)
					.tab().append(originalSupported).tab().append(proposedSupported).tab().append(bothSupported).nl();
			report.writer.print(row);
		}
		log.println("HOLLOW_SEED_SUPPORT_PASS seeds="+examined+" changes="+changed+" k="+k);
	}

	/** Matches the supplemental table's raw-read-only, unfiltered exact counting recipe. */
	private KmerTableSetU load(final int length, final boolean tiny){
		assert(length==64 || length==96) : "The comparison concerns the existing short-support table geometries.";
		final String[] args={"in="+reads, "k="+length, "hollow=f", "hashonly=f", "packed=t", "rcomp=t",
				"prefilter=f", "prealloc=f", "minprob=0", "minprobmain=f", "qtrim=f", "maq=0", "ecco=f",
				"merge=f", "reads=-1", "onepass=f", "interleaved=f", "t="+Shared.threads()};
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList(args));
		if(tiny){options.add("ways=7"); options.add("initialsize=1009");}
		final KmerTableSetU tables=new KmerTableSetU(options.toArray(new String[0]), 0);
		tables.process(new Timer());
		require(!tables.hollow && tables.kbig==length, "Evidence table must retain the requested solid geometry.");
		return tables;
	}

	/** Finds one substitution in the retained seed, ignoring appended bases. */
	private static int changedPosition(final byte[] original, final byte[] proposed){
		require(original!=null && proposed!=null && proposed.length>=original.length, "Proposed output must preserve seed length.");
		int position=-1;
		for(int i=0; i<original.length; i++){
			final int a=AminoAcid.baseToNumber[original[i]], b=AminoAcid.baseToNumber[proposed[i]];
			require(a>=0 && b>=0, "Seed support requires defined bases; unknowns cannot support a correction.");
			if(a==b){continue;}
			require(position<0, "This diagnostic requires at most one seed substitution per record.");
			position=i;
		}
		return position;
	}

	/** Counts every solid window containing the edit and wholly inside the original seed.
	 * Rolling words and the three-slot result are reused; no sequence is edited. */
	private static int scan(final byte[] bases, final int seedLength, final int position,
			final KmerTableSetU tables, final Kmer word, final int support, final int[] result){
		final int k=tables.kbig;
		assert(!word.hollow && !tables.hollow && word.kbig==k) :
				"Query words must use the same full solid key and span as the read evidence table; Kmer's boolean constructor flag selects Hollow geometry.";
		assert(position>=0 && position<seedLength && seedLength>=k && bases.length>=seedLength && result.length==3) :
				"Support windows must contain the single seed edit and exclude every appended base.";
		final int first=Math.max(0, position-k+1), last=Math.min(position, seedLength-k);
		word.clear();
		for(int i=first; i<first+k; i++){word.addRight(bases[i]);}
		result[0]=Integer.MAX_VALUE; result[1]=0; result[2]=0;
		for(int start=first; start<=last; start++){
			if(start>first){word.addRight(bases[start+k-1]);}
			final int depth=Math.max(0, tables.getCount(word));
			result[0]=Math.min(result[0], depth);
			result[1]=Math.max(result[1], depth);
			if(depth>=support){result[2]++;}
		}
		assert(result[0]<=result[1] && result[2]<=last-first+1) : "Window extrema and supported counts must conserve the inspected coordinates.";
		return last-first+1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Tests             ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests window bounds, both orientations and seed-only counts using in-memory graphs. */
	private void selfTest(){
		final Random random=new Random(421051);
		expectInvalid(new byte[]{'A', 'C'}, new byte[]{'T', 'G'});
		expectInvalid(new byte[]{'A', 'C'}, new byte[]{'A'});
		expectInvalid(new byte[]{'A', 'N'}, new byte[]{'A', 'N'});
		int checks=3;
		for(int length : new int[]{64, 96}){
			final KmerTableSetU tables=load(length, true);
			final Kmer word=new Kmer(length, false);
			final HollowSeedGuard guard=new HollowSeedGuard(tables, 2);
			try{
				for(int position : new int[]{0, 1, 31, 32, 63, 64, 95, 96, 127, 128, 249}){
					for(int strand=0; strand<2; strand++){
						final byte[] original=new byte[250];
						for(int i=0; i<original.length; i++){original[i]=AminoAcid.numberToBase[random.nextInt(4)];}
						if(strand!=0){AminoAcid.reverseComplementBasesInPlace(original);}
						final byte[] proposed=Arrays.copyOf(original, 350);
						Arrays.fill(proposed, 250, 350, (byte)'N');//Appended bases must never be queried.
						proposed[position]=AminoAcid.numberToBase[(AminoAcid.baseToNumber[original[position]]+1)&3];
						final int depth=position%4+1;
						int windows=0;
						final AbstractKmerTableU[] lanes=tables.tables();
						for(int start=0; start+length<=250; start++){
							if(start>position || start+length<=position){continue;}
							windows++;
							for(int kind=0; kind<2; kind++){
								final byte[] bases=kind==0 ? original : proposed;
								word.clear();
								for(int i=start; i<start+length; i++){word.addRight(bases[i]);}
								lanes[word.mod(lanes.length)].set(word, kind==0 ? depth : 5);
							}
						}
						final int[] before=new int[3], after=new int[3];
						require(changedPosition(original, proposed)==position, "Single-edit coordinate mismatch.");
						require(scan(original, 250, position, tables, word, 2, before)==windows &&
								scan(proposed, 250, position, tables, word, 2, after)==windows, "Rolling window bounds differ from exhaustive enumeration.");
						require(Arrays.equals(before, new int[]{depth, depth, depth>=2 ? windows : 0}) &&
								Arrays.equals(after, new int[]{5, 5, windows}), "Read-count extrema or support thresholds differ.");
						require(changedPosition(original, original)==-1, "Unchanged seeds must not produce support rows.");
						checks+=4;
						final byte[] extended=Arrays.copyOf(original, 350);
						final byte[] originalCopy=original.clone(), proposedCopy=proposed.clone();
						Arrays.fill(extended, 250, 350, (byte)'N');
						require(guard.protects(extended, 250, position)==(depth>=2), "Guard must match the all-window threshold and ignore appended unknowns.");
						final byte[] reverse=AminoAcid.reverseComplementBases(original);
						require(guard.protects(reverse, 250, 249-position)==(depth>=2), "Original-read canonical counts must protect both orientations identically.");
						require(!guard.protects(extended, 250, 250) && !guard.protects(original, length-1, 0),
								"Appended positions and seeds without a full solid window cannot receive original-seed protection.");
						// Weaken only the last overlapping window; all earlier supported
						// windows must not hide the final failure of the ALL predicate.
						final int last=Math.min(position, 250-length);
						word.clear();
						for(int i=last; i<last+length; i++){word.addRight(original[i]);}
						lanes[word.mod(lanes.length)].set(word, 1);
						require(!guard.protects(original, 250, position), "One weak overlapping window must prevent seed protection.");
						lanes[word.mod(lanes.length)].set(word, depth);
						extended[position]='N';
						require(!guard.protects(extended, 250, position), "An undefined seed window cannot supply exact protection.");
						require(Arrays.equals(original, originalCopy) && Arrays.equals(proposed, proposedCopy),
								"Protection checks must never mutate original or proposed bases.");
						checks+=6;
						checks+=exerciseRelativeGuard(tables, guard, original, proposed, position, word);
						checks+=exerciseProposedGuard(tables, guard, position, random, word);
						require(Arrays.equals(original, originalCopy) && Arrays.equals(proposed, proposedCopy),
								"Relative support must leave both original and proposed sequence arrays unchanged.");
						checks++;
					}
				}
			}finally{tables.clear();}
		}
		log.println("HOLLOW_SEED_SUPPORT_TEST_PASS checks="+checks);
	}

	/** Tests the fixed relative exception at equality, extrema and integer limits.
	 * All windows are original-seed-contained; appended Ns must remain irrelevant. */
	private static int exerciseRelativeGuard(final KmerTableSetU tables, final HollowSeedGuard guard,
			final byte[] original, final byte[] proposed, final int position, final Kmer word){
		final int length=original.length, replacement=AminoAcid.baseToNumber[proposed[position]];
		assert(replacement>=0 && original[position]!=proposed[position]) :
				"Relative fixtures compare a defined single substitution, not an unchanged allele.";
		setWindows(tables, original, length, position, 2, word);
		setWindows(tables, proposed, length, position, 7, word);
		require(guard.protects(original, length, position, replacement), "Below fourfold support must retain the original veto.");
		setWindows(tables, proposed, length, position, 8, word);
		require(!guard.protects(original, length, position, replacement), "Exactly fourfold support must release the veto.");
		final byte[] reverse=AminoAcid.reverseComplementBases(original);
		require(!guard.protects(reverse, length, length-1-position, 3-replacement), "Relative support must be reverse-complement invariant.");
		// setWindows leaves the last overlapping word available for an edge test.
		setWindows(tables, proposed, length, position, 8, word);
		tables.tables()[word.mod(tables.ways)].set(word, 7);
		require(guard.protects(original, length, position, replacement), "One weak proposed window must retain the veto.");
		setWindows(tables, proposed, length, position, 8, word);
		setWindows(tables, original, length, position, 2, word);
		tables.tables()[word.mod(tables.ways)].set(word, 3);
		require(guard.protects(original, length, position, replacement), "Use the strongest original window, not its minimum or last proposed depth.");
		setWindows(tables, original, length, position, Integer.MAX_VALUE, word);
		setWindows(tables, proposed, length, position, Integer.MAX_VALUE, word);
		require(guard.protects(original, length, position, replacement), "Fourfold multiplication must not overflow into a false release.");
		setWindows(tables, original, length, position, 1, word);
		require(!guard.protects(original, length, position, replacement), "Weak original evidence must still abstain from protection.");
		return 7;
	}

	/** Tests missing first/middle/last windows without deleting occupied table cells.
	 * Fresh random sequences avoid shared probe-chain effects; depth one must pass
	 * even when original evidence is too weak for the old protection predicate. */
	private static int exerciseProposedGuard(final KmerTableSetU tables, final HollowSeedGuard guard,
			final int position, final Random random, final Kmer word){
		assert(tables!=null && guard!=null && random!=null) : "Proposed-support fixtures need an independent populated evidence table.";
		int checks=0;
		for(int hole=0; hole<3; hole++){
			final int length=249+hole, site=Math.min(position, length-1), k=tables.kbig;
			final byte[] original=new byte[length];
			for(int i=0; i<length; i++){original[i]=AminoAcid.numberToBase[random.nextInt(4)];}
			final byte[] before=original.clone(), proposed=original.clone();
			final int replacement=(AminoAcid.baseToNumber[original[site]]+1)&3;
			proposed[site]=AminoAcid.numberToBase[replacement];
			setWindows(tables, original, length, site, 1, word);
			final int first=Math.max(0, site-k+1), last=Math.min(site, length-k);
			final int omitted=hole==0 ? first : hole==1 ? (first+last)/2 : last;
			for(int start=first; start<=last; start++){
				if(start==omitted){continue;}
				word.clear();
				for(int i=start; i<start+k; i++){word.addRight(proposed[i]);}
				tables.tables()[word.mod(tables.ways)].set(word, 1);
			}
			word.clear();
			for(int i=omitted; i<omitted+k; i++){word.addRight(proposed[i]);}
			require(tables.getCount(word)<1, "The omitted test word must not collide with an earlier controlled fixture.");
			require(!guard.protects(original, length, site, replacement), "Weak originals must preserve the existing relative-veto abstention.");
			require(guard.rejectsProposal(original, length, site, replacement), "One missing proposed window must veto the correction.");
			tables.tables()[word.mod(tables.ways)].set(word, 1);
			require(!guard.rejectsProposal(original, length, site, replacement), "Exactly depth one in every proposed window must suffice.");
			final byte[] reverse=AminoAcid.reverseComplementBases(original);
			require(!guard.rejectsProposal(reverse, length, length-1-site, 3-replacement), "Proposed evidence must be reverse-complement invariant.");
			final byte[] extended=Arrays.copyOf(original, length+100);
			Arrays.fill(extended, length, extended.length, (byte)'N');
			require(!guard.rejectsProposal(extended, length, site, replacement), "Appended unknown bases must be excluded from proposed evidence.");
			require(!guard.rejectsProposal(extended, length, length, replacement), "Appended correction positions are outside the seed-only veto.");
			require(!guard.rejectsProposal(original, k-1, Math.min(site, k-2), replacement), "Short seeds cannot impose a solid-window requirement.");
			extended[first==site ? first+1 : first]='N';
			require(guard.rejectsProposal(extended, length, site, replacement), "An undefined proposed seed window cannot provide exact support.");
			require(Arrays.equals(original, before), "Proposed-support inspection must never change the original seed.");
			require(tables.getCount(word)==1, "Proposed-support inspection must never alter the count table.");
			checks+=11;
		}
		return checks;
	}

	/** Installs exact test counts without generating reads or touching seed bytes. */
	private static void setWindows(final KmerTableSetU tables, final byte[] bases, final int length,
			final int position, final int count, final Kmer word){
		assert(length>=tables.kbig && length<=bases.length && position>=0 && position<length && count>0) :
				"Controlled support fixtures need defined seed bounds and occupied positive-count entries.";
		final int first=Math.max(0, position-tables.kbig+1), last=Math.min(position, length-tables.kbig);
		for(int start=first; start<=last; start++){
			word.clear();
			for(int i=start; i<start+tables.kbig; i++){word.addRight(bases[i]);}
			tables.tables()[word.mod(tables.ways)].set(word, count);
		}
	}

	/** Invalid comparisons must fail before an evidence row can be emitted. */
	private static void expectInvalid(final byte[] original, final byte[] proposed){
		boolean rejected=false;
		try{changedPosition(original, proposed);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Malformed seed comparison was silently accepted.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            I/O               ----------------*/
	/*--------------------------------------------------------------*/

	/** Synchronous FASTA iteration avoids producer deadlocks on comparison errors. */
	private static final class Cursor implements AutoCloseable {
		Cursor(final String path){
			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			require(ff.fasta(), "Seed and proposed inputs must be ordered FASTA files.");
			stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 0);
			stream.start();
		}
		/** Returns each record once, retaining EOF after the final batch. */
		Read next(){
			if(finished){return null;}
			while(batch==null || index>=batch.list.size()){
				batch=stream.nextList(); index=0;
				if(batch==null){finished=true; return null;}
			}
			return batch.list.get(index++);
		}
		@Override
		public void close(){require(!ReadWrite.closeStream(stream), "Seed comparison input failed.");}
		private final Streamer stream;
		private ListNum<Read> batch;
		private int index;
		private boolean finished;
	}

	/** Closes the report writer even if a malformed record aborts comparison. */
	private static final class Report implements AutoCloseable {
		Report(final String path){
			writer=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TXT, null, true, false, false, false));
			writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Support report output failed.");}
		final ByteStreamWriter writer;
	}

	/** Scientific input and result checks remain active without -ea. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String reads, seeds, queries, output, stats;
	private int k=96, support=2;
	private boolean test;
	private final PrintStream log;
}
