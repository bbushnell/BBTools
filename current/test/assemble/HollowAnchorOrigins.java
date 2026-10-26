package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Finds literal sequence occurrences of two 32bp anchors, without hashing.
 * Reports both orientations, interior disagreement and the query-oriented next base.
 * @author Brian Bushnell, Fischl */
public final class HollowAnchorOrigins {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Supports native config expansion and standard file arguments. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowAnchorOrigins.class, false);
		try{new HollowAnchorOrigins(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Selects explicit query ends from a saved KmerWindowDepth table. */
	private HollowAnchorOrigins(final String[] args, final PrintStream log_){
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("ref")){reference=b;}
			else if(a.equals("queryfasta")){queryFasta=Parse.parseBoolean(b);}
			else if(a.equals("scaffold")){scaffold=b;}
			else if(a.equals("k")){queryLength=Parse.parseIntKMG(b);}
			else if(a.equals("details")){summary=b;}
			else if(a.equals("label")){label=b;}
			else if(a.equals("paths")){paths=b;}
			else if(a.equals("flanks")){flanks=b;}
			else if(a.equals("successors")){successors=Parse.parseBoolean(b);}
			else if(a.equals("ends")){
				final String[] terms=b.split(",");
				ends=new int[terms.length];
				for(int i=0; i<terms.length; i++){ends[i]=Integer.parseInt(terms[i]);}
			}
			else if(a.equals("test")){test=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(!parser.overwrite && !parser.append, "Anchor-origin reports must be fresh.");
		require((paths==null)==(flanks==null), "paths and flanks must be supplied together.");
		require(test || (input!=null && output!=null && reference!=null && summary!=null && ends!=null),
				"Specify in=depth.tsv ref=genome.fa ends=end1[,end1] out=hits.tsv details=summary.tsv.");
		require(!queryFasta || (scaffold!=null && queryLength>64), "queryfasta=t requires scaffold=ID and k>64.");
		FastaReadInputStream.SPLIT_READS=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads each reference record once and probes both orientations of each retained key. */
	private void process(){
		if(test){selfTest(); return;}
		final ArrayList<Target> targets=loadTargets();
		if(paths!=null){attachPaths(targets);}
		final FileFormat ff=FileFormat.testInput(reference, FileFormat.FASTA, null, true, true);
		final Streamer stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final ByteBuilder row=new ByteBuilder();
		final int[] flankCounts=flanks==null ? null : new int[4];
		int references=0;
		stream.start();
		try(Report hits=new Report(output); Report flankReport=flanks==null ? null : new Report(flanks)){
			if(flankReport!=null){
				flankReport.writer.print("label\tread_ordinal\tread_id\tquery_id\tquery_end1\tsource_end1\tquery_role\tcandidate_base\tstrand\tread_word_start0\tleft_bases\tleft_differences\tright_bases\tright_differences\tquery_left\tread_left\tquery_right\tread_right\n");
			}
			hits.writer.print("label\tquery_end1\tspan\treference\tstart0\tend0\tstrand\tinterior_disagreements\tnext_base\tgenomic_word_query_orientation");
			if(successors){hits.writer.print("\tsource_end1\tquery_role\tcandidate_base");}
			hits.writer.print("\n");
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read r : batch){
					require(r.mate==null, "Reference must be unpaired.");
					Tools.toUpperCase(r.bases);
					references++;
					for(Target target : targets){
						for(int strand=0; strand<2; strand++){
							final byte[] query=strand==0 ? target.word : target.reverse;
							for(int start=0; start<=r.length()-query.length; start++){
								if(!matches(query, r.bases, start)){continue;}
								final int different=interiorDifferences(query, r.bases, start);
								final int nextPos=strand==0 ? start+query.length : start-1;
								final byte next=nextPos<0 || nextPos>=r.length() ? (byte)'.' : strand==0 ?
										r.bases[nextPos] : AminoAcid.baseToComplementExtended[r.bases[nextPos]];
								target.counts[0]++;
								if(different==0){target.counts[1]++;}
								target.counts[2+strand]++;
								final int base=AminoAcid.baseToNumber[next];
								target.counts[base<0 ? 8 : 4+base]++;
								row.clear().append(label).tab().append(target.end).tab().append(query.length)
										.tab().append(r.id).tab().append(start).tab().append(start+query.length-1)
										.tab().append(strand==0 ? '+' : '-').tab().append(different).tab().append((char)next).tab();
								for(int i=0; i<query.length; i++){
									row.append((char)(strand==0 ? r.bases[start+i] :
											AminoAcid.baseToComplementExtended[r.bases[start+query.length-1-i]]));
								}
								if(successors){appendRole(row, target);}
								hits.writer.print(row.nl());
								if(flankReport!=null){
									appendFlanks(row, target, r, references, start, strand, flankCounts);
									flankReport.writer.print(row.nl());
								}
							}
						}
					}
				}
			}
		}finally{require(!ReadWrite.closeStream(stream), "Reference input failed.");}
		require(references>0, "No reference records were examined.");
		try(Report totals=new Report(summary)){
			totals.writer.print("label\tquery_end1\tspan\toriented_hits\tfull_word_hits\tforward_hits\treverse_hits\tnext_A\tnext_C\tnext_G\tnext_T\tedge_or_ambiguous");
			if(successors){totals.writer.print("\tsource_end1\tquery_role\tcandidate_base");}
			totals.writer.print("\n");
			for(Target target : targets){
				final long[] c=target.counts;
				require(c[0]==c[2]+c[3] && c[0]==c[4]+c[5]+c[6]+c[7]+c[8] && c[1]<=c[0],
						"Orientation and neighbor totals must conserve every genomic anchor hit.");
				row.clear().append(label).tab().append(target.end).tab().append(target.word.length);
				for(long count : c){row.tab().append(count);}
				if(successors){appendRole(row, target);}
				totals.writer.print(row.nl());
			}
		}
		log.println("HOLLOW_ANCHOR_ORIGINS_PASS targets="+targets.size()+" references="+references);
	}

	/** Requires each requested end exactly once; the saved full word supplies the actual anchors. */
	private ArrayList<Target> loadTargets(){
		final ArrayList<Target> result=new ArrayList<Target>();
		final boolean[] seen=new boolean[ends.length];
		for(int i=0; i<ends.length; i++){
			require(ends[i]>64, "Query ends must include a complete long kmer.");
			for(int j=0; j<i; j++){require(ends[i]!=ends[j], "Duplicate requested query end.");}
		}
		if(queryFasta){return loadFastaTargets();}
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[] header=file.nextLine();
			require(header!=null, "Missing depth-table header.");
			lp.set(header);
			require(lp.terms()==17 && lp.termEquals("end1", 4) && lp.termEquals("word", 5), "Not a depth table.");
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				require(lp.terms()==17, "Malformed depth row.");
				final int end=lp.parseInt(4);
				for(int i=0; i<ends.length; i++){
					if(end!=ends[i]){continue;}
					require(!seen[i], "Duplicate depth row for requested end: "+end);
					seen[i]=true;
					final byte[] word=lp.parseByteArray(5);
					require(word.length==lp.parseInt(2) && word.length>64, "Word must match its long-kmer span.");
					for(byte base : word){require(AminoAcid.isFullyDefined(base), "Query word must be ACGT.");}
					result.add(new Target(end, word, (byte)'.', lp.parseString(0)));
				}
			}
		}finally{require(!file.close(), "Depth-table input failed.");}
		require(result.size()==ends.length, "Some requested query ends were missing.");
		return successors ? withSuccessors(result) : result;
	}

	/** Extracts literal query words from one saved assembly record, without building a count table. */
	private ArrayList<Target> loadFastaTargets(){
		assert(queryFasta && scaffold!=null && queryLength>64) : "Validated FASTA mode identifies one named path and its sampled span.";
		Read selected=null;
		for(Read read : ReadInputStream.toReads(input, FileFormat.FASTA, -1)){
			int length=0;
			while(length<read.id.length() && !Character.isWhitespace(read.id.charAt(length))){length++;}
			if(!scaffold.equals(read.id.substring(0, length))){continue;}
			require(selected==null, "Duplicate FASTA scaffold: "+scaffold);
			Tools.toUpperCase(read.bases);
			selected=read;
		}
		require(selected!=null, "FASTA scaffold not found: "+scaffold);
		final ArrayList<Target> result=new ArrayList<Target>();
		for(int end : ends){result.add(fromPath(selected, end, queryLength));}
		return successors ? withSuccessors(result) : result;
	}

	/** The one-based inclusive end is also the exclusive zero-based array bound. */
	private static Target fromPath(final Read read, final int end, final int span){
		require(read!=null && read.bases!=null && read.id!=null && span>64 && end>=span && end<=read.length(),
				"A saved FASTA query must contain the complete requested long-kmer interval.");
		final byte[] word=Arrays.copyOfRange(read.bases, end-span, end);
		for(byte base : word){require(AminoAcid.isFullyDefined(base), "FASTA query word must be ACGT.");}
		return new Target(end, word, (byte)'.', read.id);
	}

	/** Attaches only matching saved paths; verifies their overlap with each original query word. */
	private void attachPaths(final ArrayList<Target> targets){
		assert(paths!=null && targets!=null) : "Optional flank reports require saved paths and validated targets.";
		final HashMap<String, Read> byName=new HashMap<String, Read>();
		for(Read read : ReadInputStream.toReads(paths, FileFormat.FASTA, -1)){
			require(read.mate==null && byName.put(read.id, read)==null, "Saved paths must have unique unpaired IDs.");
			Tools.toUpperCase(read.bases);
		}
		for(Target target : targets){
			final Read read=byName.get(target.queryId);
			require(read!=null && read.length()>=target.end, "Missing or short saved query path: "+target.queryId);
			final int begin=target.end-target.word.length;
			require(begin>=0, "Saved target end must include the complete query span.");
			final int checked=target.word.length-(target.candidate=='.' ? 0 : 1);
			for(int i=0; i<checked; i++){
				require(read.bases[begin+i]==target.word[i], "Saved path differs from the original query word.");
			}
			target.path=read.bases;
		}
	}

	/** Measures literal flanks in query orientation, clipping both sides to their covered overlap. */
	private static void countFlanks(final byte[] path, final int end, final int span,
			final byte[] read, final int start, final int strand, final int[] counts){
		assert(path!=null && read!=null && counts.length==4 && end>=span && end<=path.length &&
				start>=0 && start+span<=read.length && (strand==0 || strand==1)) :
				"Validated anchor hits and saved path bounds define the overlap; uncovered tails are not mismatches.";
		final int orientedStart=strand==0 ? start : read.length-start-span;
		final int left=Math.min(end-span, orientedStart);
		final int right=Math.min(path.length-end, read.length-orientedStart-span);
		counts[0]=left; counts[1]=0; counts[2]=right; counts[3]=0;
		for(int i=1; i<=left; i++){
			final byte base=path[end-span-i];
			if(!AminoAcid.isFullyDefined(base) || base!=orientedBase(read, orientedStart-i, strand)){counts[1]++;}
		}
		for(int i=0; i<right; i++){
			final byte base=path[end+i];
			if(!AminoAcid.isFullyDefined(base) || base!=orientedBase(read, orientedStart+span+i, strand)){counts[3]++;}
		}
	}

	/** Converts a query-oriented coordinate without allocating a reversed read. */
	private static byte orientedBase(final byte[] read, final int position, final int strand){
		assert(position>=0 && position<read.length && (strand==0 || strand==1)) :
				"Flank overlap clipping must keep every oriented coordinate inside its read.";
		return strand==0 ? read[position] : AminoAcid.baseToComplementExtended[read[read.length-1-position]];
	}

	/** Emits covered bases as well as counts so an indel shift is not mistaken for many edit events. */
	private void appendFlanks(final ByteBuilder row, final Target target, final Read read,
			final int ordinal, final int start, final int strand, final int[] counts){
		assert(target.path!=null && ordinal>0) : "Attached path and input-record ordinal preserve query and duplicate-read identity.";
		countFlanks(target.path, target.end, target.word.length, read.bases, start, strand, counts);
		final int orientedStart=strand==0 ? start : read.length()-start-target.word.length;
		row.clear().append(label).tab().append(ordinal).tab().append(read.id).tab().append(target.queryId)
				.tab().append(target.end);
		appendRole(row, target);
		row.tab().append(strand==0 ? '+' : '-').tab().append(start);
		for(int count : counts){row.tab().append(count);}
		row.tab();
		for(int i=counts[0]; i>0; i--){row.append((char)target.path[target.end-target.word.length-i]);}
		row.tab();
		for(int i=counts[0]; i>0; i--){row.append((char)orientedBase(read.bases, orientedStart-i, strand));}
		row.tab();
		for(int i=0; i<counts[2]; i++){row.append((char)target.path[target.end+i]);}
		row.tab();
		for(int i=0; i<counts[2]; i++){row.append((char)orientedBase(read.bases, orientedStart+target.word.length+i, strand));}
	}

	/** Constructs each shifted A/C/G/T key from the saved word, without reference truth. */
	private static ArrayList<Target> withSuccessors(final ArrayList<Target> originals){
		assert(originals!=null) : "Successor expansion needs the validated saved-query targets from loadTargets.";
		final ArrayList<Target> expanded=new ArrayList<Target>();
		for(Target source : originals){
			require(source.candidate=='.' && source.end<Integer.MAX_VALUE,
					"Only original contexts can expand; successor coordinates must not overflow.");
			expanded.add(source);
			for(byte base : ALPHABET){
				final byte[] word=new byte[source.word.length];
				System.arraycopy(source.word, 1, word, 0, word.length-1);
				word[word.length-1]=base;
				expanded.add(new Target(source.end+1, word, base, source.queryId));
			}
		}
		return expanded;
	}

	/** Disambiguates overlapping requested contexts and successor coordinates in opt-in output. */
	private static void appendRole(final ByteBuilder row, final Target target){
		assert(row!=null && target!=null) : "Role columns must describe the same target as the preceding report row.";
		final boolean context=target.candidate=='.';
		row.tab().append(context ? target.end : target.end-1)
				.tab().append(context ? "context" : "successor").tab().append((char)target.candidate);
	}

	/** Deliberately uses literal bases, independent of Kmer hashing or canonicalization. */
	private static boolean matches(final byte[] query, final byte[] ref, final int start){
		assert(query.length>64 && start>=0 && start+query.length<=ref.length) : "Both full 32bp anchors must lie inside the reference.";
		for(int i=0; i<32; i++){if(query[i]!=ref[start+i]){return false;}}
		for(int i=query.length-32; i<query.length; i++){if(query[i]!=ref[start+i]){return false;}}
		return true;
	}

	/** Counts literal middle differences; zero means an identical solid word. */
	private static int interiorDifferences(final byte[] query, final byte[] ref, final int start){
		assert(matches(query, ref, start)) : "Interior comparison is defined only for exact anchor hits.";
		int count=0;
		for(int i=32; i<query.length-32; i++){if(query[i]!=ref[start+i]){count++;}}
		return count;
	}

	/** Tests exact matching, ignored interiors, sampled endpoints and reverse complements. */
	private static void selfTest(){
		final byte[] query=new byte[128];
		final Random random=new Random(4286901);
		for(int i=0; i<query.length; i++){query[i]=ALPHABET[random.nextInt(4)];}
		final byte[] ref=Arrays.copyOf(query, query.length);
		require(matches(query, ref, 0) && interiorDifferences(query, ref, 0)==0, "Exact fixture failed.");
		ref[64]=ref[64]=='A' ? (byte)'C' : (byte)'A';
		require(matches(query, ref, 0) && interiorDifferences(query, ref, 0)==1, "Interior mutation changed the anchors.");
		ref[31]=ref[31]=='A' ? (byte)'C' : (byte)'A';
		require(!matches(query, ref, 0), "Left endpoint mutation was ignored.");
		System.arraycopy(query, 0, ref, 0, query.length);
		ref[96]='N';
		require(!matches(query, ref, 0), "Unknown right-anchor base was accepted.");
		final Target target=new Target(128, query);
		final byte[] twice=Arrays.copyOf(target.reverse, target.reverse.length);
		AminoAcid.reverseComplementBasesInPlace(twice);
		require(Arrays.equals(query, twice), "Reverse-complement matching lost query orientation.");
		int checks=5;
		final ArrayList<Target> originals=new ArrayList<Target>();
		originals.add(target);
		final ArrayList<Target> expanded=withSuccessors(originals);
		require(expanded.size()==5 && expanded.get(0)==target, "Expansion must retain the original context followed by four successors.");
		checks++;
		for(int i=0; i<ALPHABET.length; i++){
			final Target child=expanded.get(i+1);
			final byte[] expected=Arrays.copyOfRange(query, 1, query.length+1);
			expected[expected.length-1]=ALPHABET[i];
			require(child.end==129 && child.candidate==ALPHABET[i] && Arrays.equals(child.word, expected),
					"Successor must drop exactly the first base and append the designated candidate.");
			final byte[] restored=Arrays.copyOf(child.reverse, child.reverse.length);
			AminoAcid.reverseComplementBasesInPlace(restored);
			require(Arrays.equals(restored, expected), "Successor reverse complement must preserve the shifted key.");
			final ByteBuilder role=new ByteBuilder();
			appendRole(role, child);
			require(role.toString().equals("\t128\tsuccessor\t"+(char)ALPHABET[i]),
					"Successor role must report the unshifted source coordinate and query-oriented candidate.");
			checks+=3;
		}
		require(Arrays.equals(query, twice) && target.candidate=='.', "Successor construction must not mutate the original word.");
		checks++;
		final ByteBuilder role=new ByteBuilder();
		appendRole(role, target);
		require(role.toString().equals("\t128\tcontext\t."), "Context role must remain distinct from every candidate.");
		checks++;
		final byte[] path=new byte[260];
		for(int i=0; i<path.length; i++){path[i]=ALPHABET[random.nextInt(4)];}
		final byte[] covered=Arrays.copyOfRange(path, 52, 230);
		final int[] flankCounts=new int[4];
		countFlanks(path, 200, 128, covered, 20, 0, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{20, 0, 30, 0}), "Covered flank lengths or exact agreement failed.");
		covered[19]=covered[19]=='A' ? (byte)'C' : (byte)'A';
		covered[148]=covered[148]=='A' ? (byte)'C' : (byte)'A';
		countFlanks(path, 200, 128, covered, 20, 0, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{20, 1, 30, 1}), "Flank mismatches must be separate from anchor/interior differences.");
		AminoAcid.reverseComplementBasesInPlace(covered);
		countFlanks(path, 200, 128, covered, 30, 1, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{20, 1, 30, 1}), "Reverse-oriented flanks must preserve coverage and counts.");
		final byte[] shorter=Arrays.copyOfRange(path, 0, 205);
		countFlanks(shorter, 200, 128, covered, 30, 1, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{20, 1, 5, 1}), "Uncovered query tail must not create apparent differences.");
		countFlanks(Arrays.copyOfRange(path, 72, 200), 128, 128, covered, 30, 1, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{0, 0, 0, 0}), "A query consisting only of the key has no covered flanks.");
		countFlanks(path, 200, 128, Arrays.copyOfRange(path, 72, 200), 0, 0, flankCounts);
		require(Arrays.equals(flankCounts, new int[]{0, 0, 0, 0}), "A read ending at both anchors contributes no flank evidence.");
		checks+=6;
		final Read saved=new Read(path, null, "saved members=0+,1-", 0);
		final Target extracted=fromPath(saved, 200, 128);
		require(extracted.end==200 && extracted.queryId.equals(saved.id) &&
				Arrays.equals(extracted.word, Arrays.copyOfRange(path, 72, 200)),
				"FASTA extraction must preserve the exact full path ID and one-based end convention.");
		boolean rejected=false;
		try{fromPath(saved, 127, 128);}catch(IllegalArgumentException e){rejected=true;}
		require(rejected, "A partial FASTA word must fail instead of padding or shifting the query.");
		checks+=2;
		System.out.println("HOLLOW_ANCHOR_ORIGINS_TEST_PASS checks="+checks);
	}

	/** Keeps scientific input checks enabled without relying on JVM assertion flags. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Retains one small query and primitive aggregate counters. */
	private static final class Target {
		/** Retains an original saved context with no candidate-base designation. */
		Target(final int end_, final byte[] word_){
			this(end_, word_, (byte)'.', null);
		}
		/** Caches both orientations once so the literal scan allocates no per-window keys. */
		Target(final int end_, final byte[] word_, final byte candidate_, final String queryId_){
			assert(word_!=null && word_.length>64) : "Literal matching requires two complete outer 32bp anchors.";
			end=end_;
			candidate=candidate_;
			queryId=queryId_;
			word=word_;
			reverse=Arrays.copyOf(word, word.length);
			AminoAcid.reverseComplementBasesInPlace(reverse);
		}
		final int end;
		final byte candidate;
		final String queryId;
		final byte[] word, reverse;
		byte[] path;
		final long[] counts=new long[9];
	}

	/** Ensures native report threads terminate on validation failures as well as success. */
	private static final class Report implements AutoCloseable {
		Report(final String path){
			writer=new ByteStreamWriter(path, false, false, true);
			writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Anchor report output failed.");}
		final ByteStreamWriter writer;
	}

	private String input, output, reference, summary, label="query";
	private String paths, flanks;
	private int[] ends;
	private boolean test=false;
	private boolean successors=false;
	private boolean queryFasta=false;
	private String scaffold;
	private int queryLength=224;
	private final PrintStream log;
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
