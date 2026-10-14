package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import ukmer.Kmer;

/** Independent exhaustive truth and non-mutating candidate geometry fixtures.
 * @author Fischl */
public final class SameKJoinTest {

	/** Small deterministic tests only; optional directory receives CLI fixtures. */
	public static void main(final String[] args){
		if(args.length==4 && args[0].equals("features")){SameKExportTest.verify(args[1], args[2], Integer.parseInt(args[3])); return;}
		if(args.length==3 && args[0].equals("parity")){parity(args[1], args[2]); return;}
		thirdOccurrence();
		boundariesAndPalindromes();
		exhaustive();
		for(int k : new int[]{3, 31, 32, 33, 62, 64, 95, 96, 124}){geometry(k);}
		invalid();
		pairLabels();
		mergeChains();
		mergeGuards();
		checks+=SameKJoinFeatureTest.run();
		if(args.length>0){fixtures(args[0]);}
		System.err.println("SAME_K_TEST_PASS checks="+checks);
	}

	/** A third repeat copy must not disappear behind a two-hit multiplicity shortcut. */
	private static void thirdOccurrence(){
		final String word=random(64, 91);
		KmerBranchOracle.Result result=oracle(word, "A"+word+"CNNNA"+word+"CNNNT"+word+"G", false);
		check(result.occurrences==3 && result.leftMask==9 && result.rightMask==6,
				"Third-hit conflicting extensions must be retained.");
		check(!result.unambiguous(0, 1), "A present extension is not a unique continuation.");
		result=oracle(word, "A"+word+"CNNNA"+word+"CNNNA"+word+"C", false);
		check(result.occurrences==3 && result.unambiguous(0, 1), "Repeated words with identical flanks are locally unambiguous.");
		final String reversed=rc("A"+word+"C");
		result=oracle(word, reversed, false);
		check(result.occurrences==1 && result.leftMask==1 && result.rightMask==2,
				"Reverse-reference flanks must be transformed into query orientation.");
	}

	/** Unknown flanks, origin wrapping and palindrome orientation are distinct facts. */
	private static void boundariesAndPalindromes(){
		final String word=random(33, 7);
		KmerBranchOracle.Result result=oracle(word, word+"C", false);
		check(result.occurrences==1 && result.unknownFlanks==1 && result.leftMask==0 && result.rightMask==2,
				"Linear boundary is unknown, not a safe missing edge.");
		result=oracle(word, "N"+word+"C", false);
		check(result.unknownFlanks==1 && !result.unambiguous(0, 1), "N flank must block strict local safety.");
		result=oracle("ATAT", "CATATG", false);
		check(result.occurrences==1 && result.unambiguous(1, 2), "Physical palindrome occurrence must not count twice.");
		result=oracle("ATAT", "AATATC", false);
		check(result.occurrences==1 && result.leftMask==5 && result.rightMask==10,
				"Palindrome retains both flank orientations.");
		compare(Arrays.asList(bytes("TAAC"), bytes("ACGT")), Arrays.asList(bytes("ACGTTA")), new boolean[]{true}, 4);
		compare(Arrays.asList(bytes("TAAC")), Arrays.asList(bytes("ACGTTA")), new boolean[]{false}, 4);
		check(oracle("CCCC", "AAAAAA", false).occurrences==0, "Absent word acquired fabricated support.");
	}

	/** Compares the indexed scan with every start/strand on mixed linear/circular references. */
	private static void exhaustive(){
		final Random random=new Random(7731);
		for(int k : new int[]{2, 3, 4, 31, 32, 33, 64, 95, 124}){
			for(int trial=0; trial<8; trial++){
				final String original=random(400, random.nextLong());
				final byte[] first=bytes(original.toLowerCase(java.util.Locale.ROOT));
				first[73]='N';
				final String second=rc(original)+"NN"+original;
				final List<byte[]> references=Arrays.asList(first, bytes(second));
				final boolean[] topology={trial%2==0, trial%3==0};
				final ArrayList<byte[]> queries=new ArrayList<byte[]>();
				for(int i=0; i<12; i++){
					final int start=random.nextInt(original.length()-k+1);
					queries.add(bytes(original.substring(start, start+k)));
				}
				queries.add(bytes(random(k, random.nextLong())));
				compare(queries, references, topology, k);
			}
		}
	}

	/** Builds a scalar oracle using string equality at every position, without seeds or hash tables. */
	private static void compare(final List<byte[]> queries, final List<byte[]> references,
			final boolean[] topology, final int k){
		final KmerBranchOracle.Result[] results=KmerBranchOracle.scan(queries, references, topology, k);
		check(results.length==queries.size(), "Query/result cardinality changed.");
		for(int q=0; q<queries.size(); q++){
			final String word=new String(queries.get(q), StandardCharsets.US_ASCII).toUpperCase(java.util.Locale.ROOT);
			final String reverse=rc(word);
			final long[] expected=new long[5];
			for(int r=0; r<references.size(); r++){
				final String reference=new String(references.get(r), StandardCharsets.US_ASCII).toUpperCase(java.util.Locale.ROOT);
				final String expanded=topology[r] ? reference+reference : reference;
				for(int start=0; start<reference.length(); start++){
					if(start+k>expanded.length()){continue;}
					final String observed=expanded.substring(start, start+k);
					final boolean forward=word.equals(observed), backward=reverse.equals(observed);
					if(!forward && !backward){continue;}
					expected[0]++;
					final int left=(start==0 && !topology[r]) ? -1 : "ACGT".indexOf(reference.charAt((start+reference.length()-1)%reference.length()));
					final int right=(start+k>=reference.length() && !topology[r]) ? -1 : "ACGT".indexOf(reference.charAt((start+k)%reference.length()));
					if(left<0 || right<0){expected[1]++;}
					if(forward){masks(expected, left, right);}
					if(backward){masks(expected, right<0 ? -1 : 3-right, left<0 ? -1 : 3-left);}
				}
			}
			final KmerBranchOracle.Result got=results[q];
			check(Arrays.equals(expected, new long[]{got.occurrences, got.unknownFlanks, got.leftMask, got.rightMask, got.pairMask}),
					"Indexed/exhaustive disagreement for K="+k+" query="+word+" expected="+Arrays.toString(expected));
		}
	}

	/** Scalar truth accumulation independent of Result.add. */
	private static void masks(final long[] result, final int left, final int right){
		assert(result.length==5) : "Scalar fixture layout is count, unknown, left, right and pair masks.";
		if(left>=0){result[2]|=1L<<left;}
		if(right>=0){result[3]|=1L<<right;}
		if(left>=0 && right>=0){result[4]|=1L<<(4*left+right);}
	}

	/** All four endpoint orientations and third competing targets survive without mutation. */
	private static void geometry(final int k){
		final int end=k+40;
		final String genome=random(3*k+150, 1700+k);
		final String a=genome.substring(0, end), b=genome.substring(end-k+1);
		for(int ar=0; ar<2; ar++){
			for(int br=0; br<2; br++){
				final ArrayList<Contig> contigs=new ArrayList<Contig>();
				contigs.add(new Contig(bytes(ar==0 ? a : rc(a)), 71));
				for(int n=0; n<3; n++){contigs.add(new Contig(bytes(br==0 ? b : rc(b)), 90+n));}
				contigs.add(new Contig(bytes("NNNN"), 99));
				final byte[] before=contigs.get(0).bases.clone();
				final Evidence evidence=new Evidence(k);
				evidence.add(genome, 12);
				final boolean[] found=new boolean[3];
				final long count=new SameKJoinCollector(contigs, evidence, k).collect(new SameKJoinCollector.Sink(){
					@Override
					public void accept(final long event, final int se, final int de, final Contig source,
							final boolean sr, final Contig dest, final boolean dr, final int[] ac, final int[] bc){
						check(ac.length==9 && bc.length==9 && ac[0]>0 && bc[0]>0, "Missing endpoint evidence.");
						if(source==contigs.get(0)){
							final String left=new String(SameKJoinCollector.context(source, sr, 0, source.length()), StandardCharsets.US_ASCII);
							final String right=new String(SameKJoinCollector.context(dest, dr, 0, dest.length()), StandardCharsets.US_ASCII);
							if(left.equals(a) && right.equals(b)){found[de/2-1]=true;}
						}
					}
				});
				check(count>=3 && found[0] && found[1] && found[2], "Shared anchor lost alternate partners at K="+k);
				check(Arrays.equals(before, contigs.get(0).bases) && contigs.get(0).id==71 && contigs.get(0).leftEdges==null,
						"Diagnostic collector mutated the assembly.");
				check(new SameKJoinCollector(contigs, new Evidence(k), k).collect(new SameKJoinCollector.Sink(){
					@Override public void accept(long e, int s, int d, Contig x, boolean xr, Contig y, boolean yr, int[] xc, int[] yc){
						throw new AssertionError("Unsupported endpoint became a graph candidate.");
					}
				})==0, "Empty read table must produce no supported candidates.");
			}
		}
	}

	/** Invalid reference metadata and query geometry fail before assigning labels. */
	private static void invalid(){
		for(int mode=0; mode<4; mode++){
			boolean failed=false;
			try{
				if(mode==0){oracle("ACN", "ACGT", false);}
				if(mode==1){KmerBranchOracle.scan(Arrays.asList(bytes("ACGT")), Arrays.asList(bytes("ACG")), new boolean[]{true}, 4);}
				if(mode==2){KmerBranchOracle.scan(Arrays.asList(bytes("ACG")), Arrays.asList(bytes("ACGT")), new boolean[0], 3);}
				if(mode==3){KmerBranchOracle.scan(Arrays.asList(bytes("ACG")), Arrays.asList(bytes("ACGT")), new boolean[]{false}, 4);}
			}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid oracle input was accepted: "+mode);
		}
		check(KmerBranchOracle.scan(new ArrayList<byte[]>(), Arrays.asList(bytes("ACGT")), new boolean[]{false}, 3).length==0,
				"An empty candidate set is a valid census.");
	}

	/** Pair truth is not the local-repeat mask and does not call absent ambiguous context false. */
	private static void pairLabels(){
		final String genome=random(300, 181);
		final FusionJoinDiagnostic.Join good=new FusionJoinDiagnostic.Join(1, 32, 0, 1, 31, 0, 0,
				bytes(genome.substring(20, 140)), bytes(genome.substring(109, 240)));
		final ArrayList<String> refs=new ArrayList<String>(Arrays.asList(genome, rc(genome)));
		final FusionJoinDiagnostic.Hit whole=FusionJoinDiagnostic.find(refs, good.product);
		final FusionJoinDiagnostic.Hit left=FusionJoinDiagnostic.find(refs, good.left), right=FusionJoinDiagnostic.find(refs, good.right);
		check(SameKJoinDiagnostic.pairLabel(good, refs, whole, left, right, false).equals("supported"), "Unique correct pair lost.");
		whole.count=2;
		check(SameKJoinDiagnostic.pairLabel(good, refs, whole, left, right, false).equals("unresolved"), "Repeated whole context was promoted.");
		whole.count=0; right.pos+=50;
		check(SameKJoinDiagnostic.pairLabel(good, refs, whole, left, right, false).equals("contradicted"), "Unique incompatible placements lost.");
		check(SameKJoinDiagnostic.pairLabel(good, refs, whole, left, right, true).equals("unresolved"), "Unimplemented circular pair context became negative.");
		right.count=2;
		check(SameKJoinDiagnostic.pairLabel(good, refs, whole, left, right, false).equals("unresolved"), "Ambiguous arm became a false negative.");
	}

	/** Spelling is invariant to all orientations in a three-unitig chain. */
	private static void mergeChains(){
		final int k=32;
		final String genome=random(320, 701);
		final String[] pieces={genome.substring(0, 130), genome.substring(99, 240), genome.substring(209)};
		for(int mask=0; mask<8; mask++){
			final ArrayList<Contig> contigs=new ArrayList<Contig>();
			for(int i=0; i<3; i++){contigs.add(new Contig(bytes((mask&(1<<i))==0 ? pieces[i] : rc(pieces[i])), i));}
			final SameKJoinMerge.Choices choices=new SameKJoinMerge.Choices(3);
			choices.add(((mask&1)==0 ? 1 : 0), 2+((mask&2)==0 ? 0 : 1));
			choices.add(2+((mask&2)==0 ? 1 : 0), 4+((mask&4)==0 ? 0 : 1));
			final SameKJoinMerge.Result result=SameKJoinMerge.merge(contigs, choices.mates(), k);
			check(result.joins==2 && result.contigs.size()==1 && result.cycleUnitigs==0, "Linear chain lost a join.");
			check(new String(result.contigs.get(0).bases, StandardCharsets.US_ASCII).equals(genome), "Oriented K-1 chain spelling changed.");
			for(int i=0; i<3; i++){
				check(new String(contigs.get(i).bases, StandardCharsets.US_ASCII).equals((mask&(1<<i))==0 ? pieces[i] : rc(pieces[i])),
						"Merge changed source sequence or orientation.");
				check(!contigs.get(i).used() && contigs.get(i).id==i, "Merge changed source ownership or IDs.");
			}
		}
	}

	/** Competitors, same-contig closures and cycles never become arbitrary merges. */
	private static void mergeGuards(){
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(new Contig(bytes("AACCGG"), 0));
		contigs.add(new Contig(bytes("GGAATT"), 1));
		contigs.add(new Contig(bytes("GGTTCC"), 2));
		final SameKJoinMerge.Choices ambiguous=new SameKJoinMerge.Choices(3);
		ambiguous.add(1, 2); ambiguous.add(1, 4);
		SameKJoinMerge.Result result=SameKJoinMerge.merge(contigs, ambiguous.mates(), 3);
		check(result.joins==0 && result.contigs.size()==3, "A competing supported partner was arbitrarily selected.");
		final SameKJoinMerge.Choices closure=new SameKJoinMerge.Choices(3);
		closure.add(0, 1);
		result=SameKJoinMerge.merge(contigs, closure.mates(), 3);
		check(result.joins==0, "Self-contig closure was merged.");
		final SameKJoinMerge.Choices cycle=new SameKJoinMerge.Choices(3);
		cycle.add(0, 3); cycle.add(1, 2);
		result=SameKJoinMerge.merge(contigs, cycle.mates(), 3);
		check(result.joins==0 && result.cycleUnitigs==2 && result.contigs.size()==3, "Cycle was opened or swallowed.");
		check(SameKJoinMerge.merge(new ArrayList<Contig>(), new int[0], 3).contigs.isEmpty(), "Empty census is not a valid no-op.");
		for(int mode=0; mode<3; mode++){
			boolean failed=false;
			try{
				final int[] mates={-1, 2, 1, -1, -1, -1};
				if(mode==0){mates[2]=-1;}
				if(mode==1){mates[1]=6;}
				if(mode==2){mates[1]=-1; mates[2]=-1; mates[0]=4; mates[4]=0;}
				SameKJoinMerge.merge(contigs, mates, 3);
			}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid mates or overlap accepted, mode="+mode);
		}
	}

	/** Writes small ordinary FASTA fixtures for a separate native collect/label smoke run. */
	private static void fixtures(final String directory){
		final String genome=random(8000, 9501);
		write(directory+"/reference.fa", ">truth\n"+genome+"\n");
		final StringBuilder reads=new StringBuilder();
		for(int pos=0, id=0; pos+150<=genome.length(); pos+=10){
			for(int copy=0; copy<2; copy++){
				String read=genome.substring(pos, pos+150);
				if(copy==0 && pos%70==0){read=read.substring(0, 74)+(read.charAt(74)=='A' ? 'C' : 'A')+read.substring(75);}
				reads.append('>').append(id++).append('\n').append(read).append('\n');
			}
		}
		write(directory+"/reads.fa", reads.toString());
	}

	/** Compares sequence multisets, ignoring the deliberately different raw-unitig headers. */
	private static void parity(final String first, final String second){
		final ArrayList<String> a=new ArrayList<String>(), b=new ArrayList<String>();
		for(Contig c : TadpoleGraph.readContigs(first)){a.add(canonical(new String(c.bases, StandardCharsets.US_ASCII)));}
		for(Contig c : TadpoleGraph.readContigs(second)){b.add(canonical(new String(c.bases, StandardCharsets.US_ASCII)));}
		java.util.Collections.sort(a); java.util.Collections.sort(b);
		check(!a.isEmpty() && a.equals(b), "Diagnostic collection changed the matched raw-unitig assembly.");
		System.err.println("SAME_K_ASSEMBLY_PARITY_PASS unitigs="+a.size());
	}

	/** Fixture writes use native output and fail loudly on writer errors. */
	private static void write(final String path, final String text){
		assert(path!=null && text!=null) : "Fixture output path and contents are explicit.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new RuntimeException("Fixture write failed: "+path);}
	}

	/** One-query convenience without hiding the reference topology flag. */
	private static KmerBranchOracle.Result oracle(final String query, final String reference, final boolean circular){
		return KmerBranchOracle.scan(Arrays.asList(bytes(query)), Arrays.asList(bytes(reference)), new boolean[]{circular}, query.length())[0];
	}

	/** Test-only reference count table; no production boxed per-base map. */
	private static final class Evidence implements TadpoleGraph.Counts {
		Evidence(final int k_){k=k_;}
		void add(final String sequence, final int depth){
			assert(depth>0) : "Synthetic evidence adds positive read support.";
			for(int i=0; i+k<=sequence.length(); i++){
				final String word=canonical(sequence.substring(i, i+k));
				map.put(word, depth);
			}
		}
		@Override public int count(final Kmer word){return Math.max(0, map.get(canonical(word.toString())));}
		private final int k;
		private final ObjectIntMap<String> map=new ObjectIntMap<String>(String.class);
	}

	private static String canonical(final String word){final String reverse=rc(word); return word.compareTo(reverse)<=0 ? word : reverse;}
	private static byte[] bytes(final String word){return word.getBytes(StandardCharsets.US_ASCII);}
	private static String rc(final String word){
		final StringBuilder result=new StringBuilder(word.length());
		for(int i=word.length()-1; i>=0; i--){result.append("TGCA".charAt("ACGT".indexOf(word.charAt(i))));}
		return result.toString();
	}
	private static String random(final int length, final long seed){
		assert(length>0) : "Deterministic fixtures use nonempty sequences.";
		final Random random=new Random(seed);
		final StringBuilder result=new StringBuilder(length);
		for(int i=0; i<length; i++){result.append("ACGT".charAt(random.nextInt(4)));}
		return result.toString();
	}
	private static void check(final boolean condition, final String message){checks++; if(!condition){throw new AssertionError(message);}}
	private static int checks;
}
