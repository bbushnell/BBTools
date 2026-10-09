package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import dna.AminoAcid;
import map.ObjectIntMap;
import ukmer.Kmer;

/** Exact geometry, count, orientation and premerge regressions for fusion support.
 * @author Fischl */
public final class FusionKmerSupportTest {

	/** Runs small deterministic tests without reading external data. */
	public static void main(final String[] args){
		final boolean packed=Kmer.PACKED;
		Kmer.PACKED=true;
		try{
			for(int k : new int[]{11, 32, 33, 40, 64, 124}){orientations(k);}
			boundaries();
			branches();
			longFlanks();
			noTrimming();
			selectionAndPremerge();
			parser();
		}finally{Kmer.PACKED=packed;}
		System.out.println("PASS: fusion support orientations, trims, internal missing words, depth, bounds, N, selection, premerge and parser.");
	}

	/** Both strands and both directions inspect the same informative word set. */
	private static void orientations(final int k){
		final int overlap=k, end=300;
		final String genome=random(700, 12+k);
		final String source=genome.substring(0, end), dest=genome.substring(end-overlap);
		final Evidence evidence=new Evidence(k);
		evidence.add(genome, 3);
		for(int orientation=0; orientation<4; orientation++){
			final boolean sr=(orientation&1)!=0, dr=(orientation&2)!=0;
			final Contig a=contig(0, source+"NN", sr), b=contig(1, "NNN"+dest, dr);
			final FusionKmerSupport support=new FusionKmerSupport(evidence, k, 3);
			check(support.supported(a, sr, 2, b, dr, 3, overlap), "Supported trimmed join rejected at K="+k);
			check(support.words==2*k+1, "Wrong word count across full K flanks at K="+k);
			check(support.supported(b, !dr, 3, a, !sr, 2, overlap), "Reciprocal geometry changed support");
			check(!new FusionKmerSupport(evidence, k, 4).supported(a, sr, 2, b, dr, 3, overlap),
					"Count floor ignored");
			final String internal=genome.substring(end-overlap-2, end-overlap-2+k);
			evidence.set(internal, 0);
			check(!support.supported(a, sr, 2, b, dr, 3, overlap), "Internal missing word escaped the scan");
			evidence.set(internal, 3);
		}
	}

	/** The first, interior and last words and both complete flanks are required. */
	private static void boundaries(){
		final String genome=random(96, 99);
		final Contig a=contig(0, genome.substring(0, 64), false);
		final Contig b=contig(1, genome.substring(32), false);
		final Evidence evidence=new Evidence(32);
		evidence.add(genome, 10);
		final FusionKmerSupport support=new FusionKmerSupport(evidence, 32, 1);
		check(support.supported(a, false, 0, b, false, 0, 32), "Exact full window failed");
		check(support.words==65, "K32 exact overlap must inspect65 words");
		for(int offset : new int[]{0, 30, 64}){
			evidence.set(genome.substring(offset, offset+32), 0);
			check(!support.supported(a, false, 0, b, false, 0, 32), "Missing boundary/internal word accepted");
			evidence.set(genome.substring(offset, offset+32), 10);
		}
		check(!support.supported(contig(0, genome.substring(1, 64), false), false, 0, b, false, 0, 32),
				"Missing left flank silently shortened the window");
		check(!support.supported(a, false, 0, contig(1, genome.substring(32, 95), false), false, 0, 32),
				"Missing right flank silently shortened the window");
		a.bases[0]='N';
		check(!support.supported(a, false, 0, b, false, 0, 32), "Unknown flank was skipped");
	}

	/** Rejects both branch directions and hidden dominant alternatives, but tolerates weak noise. */
	private static void branches(){
		final String genome=random(160, 109);
		final Contig a=contig(0, genome.substring(0, 80), false);
		final Contig b=contig(1, genome.substring(48), false);
		for(boolean right : new boolean[]{true, false}){
			final int start=right ? 48 : 47;
			final String original=genome.substring(start, start+32);
			final int pos=right ? 31 : 0;
			final char other=original.charAt(pos)=='A' ? 'C' : 'A';
			final String alternative=original.substring(0, pos)+other+original.substring(pos+1);
			for(int depth : new int[]{1, 5, 10, 100}){
				final Evidence evidence=new Evidence(32);
				evidence.add(genome, 10);
				evidence.set(alternative, depth);
				final boolean accepted=new FusionKmerSupport(evidence, 32, 1).supported(a, false, 0, b, false, 0, 32);
				check(accepted==(depth==1), "Incorrect branch decision: right="+right+", alternative="+depth);
			}
		}
	}

	/** Larger flanks include distant branch evidence without changing the counting K. */
	private static void longFlanks(){
		final String genome=random(1200, 941);
		final int k=32, overlap=55, end=600, flank=128;
		for(int orientation=0; orientation<4; orientation++){
			final boolean sr=(orientation&1)!=0, dr=(orientation&2)!=0;
			final Contig a=contig(0, genome.substring(0, end), sr);
			final Contig b=contig(1, genome.substring(end-overlap), dr);
			final Evidence evidence=new Evidence(k);
			evidence.add(genome, 10);
			final FusionKmerSupport longer=new FusionKmerSupport(evidence, k, 1, flank);
			check(longer.supported(a, sr, 0, b, dr, 0, overlap), "Full longer flanks failed");
			check(longer.words==overlap+2*flank-k+1, "Extended word count is incorrect");
			check(longer.supported(b, !dr, 0, a, !sr, 0, overlap), "Extended reciprocal geometry failed");
			for(int offset : new int[]{end-overlap-flank, end+flank-k}){
				final String word=genome.substring(offset, offset+k);
				evidence.set(word, 0);
				check(!longer.supported(a, sr, 0, b, dr, 0, overlap), "Missing outer word was ignored");
				evidence.set(word, 10);
			}
			final int offset=end-overlap-80;
			final String word=genome.substring(offset, offset+k);
			final String branch=word.substring(0, k-1)+(word.charAt(k-1)=='A' ? 'C' : 'A');
			evidence.set(branch, 10);
			check(new FusionKmerSupport(evidence, k, 1).supported(a, sr, 0, b, dr, 0, overlap),
					"A branch outside auto flanks unexpectedly entered the window");
			check(!longer.supported(a, sr, 0, b, dr, 0, overlap), "A branch inside longer flanks escaped");
			check(new FusionKmerSupport(evidence, k, 1, 1).flank==k, "Minimum flank shortened K context");
			check(!new FusionKmerSupport(evidence, k, 1, end).supported(a, sr, 0, b, dr, 0, overlap),
					"Insufficient full flank was silently clipped");
		}
	}

	/** No-trim eligibility preserves longer exact terminal overlaps while excluding discarded tails. */
	private static void noTrimming(){
		final String genome=random(500, 943);
		for(int trims=0; trims<4; trims++){
			for(boolean allowed : new boolean[]{true, false}){
				final Contig a=contig(0, genome.substring(0, 250)+((trims&1)!=0 ? "NN" : ""), false);
				final Contig b=contig(1, ((trims&2)!=0 ? "NNN" : "")+genome.substring(195), false);
				a.rightBridgeEndpoint=b.leftBridgeEndpoint=true;
				final ArrayList<Contig> contigs=new ArrayList<Contig>();
				contigs.add(a); contigs.add(b);
				final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, 32, 80);
				overlapper.allowTrim=allowed;
				check(overlapper.addEdges()==(allowed || trims==0 ? 1 : 0),
						"No-trim eligibility confused seed offset with discarded bases: "+trims+"/"+allowed);
				if(!allowed && trims==0){
					check(a.rightEdges.get(0).overlap==55, "No-trim lost the longer-than-K terminal overlap");
				}
			}
		}
	}

	/** Discovery vetoes bad joins, and changed read support is rechecked before any merge mutation. */
	private static void selectionAndPremerge(){
		final String genome=random(300, 711);
		final Evidence evidence=new Evidence(25);
		evidence.add(genome, 5);
		final Contig a=contig(0, genome.substring(0, 150), false);
		final Contig b=contig(1, genome.substring(120), false);
		a.rightBridgeEndpoint=b.leftBridgeEndpoint=true;
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(a); contigs.add(b);
		final FusionKmerSupport support=new FusionKmerSupport(evidence, 25, 1);
		final CrossKTipOverlapper rejected=new CrossKTipOverlapper(contigs, 25, 35);
		rejected.support=new FusionKmerSupport(new Evidence(25), 25, 1);
		check(rejected.addEdges()==0, "Unbacked reciprocal overlap created edges");
		final CrossKTipOverlapper accepted=new CrossKTipOverlapper(contigs, 25, 35);
		accepted.support=support;
		check(accepted.addEdges()==1, "Supported reciprocal overlap disappeared");
		final HashMap<Integer, ArrayList<Edge>> incoming=new HashMap<Integer, ArrayList<Edge>>();
		//Incoming indexes own their lists; outbound mutation must not erase them.
		incoming.put(0, new ArrayList<Edge>(b.leftEdges));
		incoming.put(1, new ArrayList<Edge>(a.rightEdges));
		final boolean oldCross=BubblePopper.crossKMerge, oldDirect=BubblePopper.popDirect;
		final FusionKmerSupport oldSupport=BubblePopper.crossKSupport;
		try{
			BubblePopper.crossKMerge=true;
			BubblePopper.popDirect=true;
			BubblePopper.crossKSupport=support;
			evidence.set(genome.substring(119, 144), 0);
			check(new BubblePopper(contigs, incoming, 25).expand(a)==0, "Stale support allowed an actual merge");
			check(!a.used() && !b.used() && a.length()==150 && b.length()==180, "Rejected merge changed sequences");
			evidence.set(genome.substring(119, 144), 5);
			check(new BubblePopper(contigs, incoming, 25).expand(a)==1, "Supported actual merge failed");
			check(new String(a.bases, StandardCharsets.US_ASCII).equals(genome), "Merged product changed");
		}finally{
			BubblePopper.crossKMerge=oldCross;
			BubblePopper.popDirect=oldDirect;
			BubblePopper.crossKSupport=oldSupport;
		}
	}

	/** Fusion-only arguments stay out of the assembly/bridge parsers and defaults remain disabled. */
	private static void parser(){
		final TadpoleMulti.Config defaults=new TadpoleMulti.Config(new String[]{"k=95,64", "out=test.fa"});
		check(!defaults.fusePath && defaults.fusePathFlank==0 && defaults.fuseTrim,
				"Path/trimming defaults changed");
		final TadpoleMulti.Config config=new TadpoleMulti.Config(new String[]{"k=95,64", "out=test.fa",
				"fusepath=t", "fusepathdepth=2", "fusepathflank=128", "fusetrim=f", "prefilter=t", "hashkmers=pair"});
		check(config.fusePath && config.fusePathDepth==2, "Path flags not parsed");
		check(config.fusePathFlank==128 && !config.fuseTrim, "Flank/trimming flags not parsed");
		final TadpoleMulti multi=new TadpoleMulti(config);
		for(String arg : multi.makeArgs(64, false)){
			check(!arg.startsWith("fusepath") && !arg.startsWith("fusetrim"), "Fusion flag leaked into bridge phase");
		}
		final String[] supportArgs=multi.makeFusionSupportArgs(64);
		check(supportArgs[supportArgs.length-2].equals("prefilter=f"), "Support inherits a lossy prefilter");
		check(supportArgs[supportArgs.length-3].equals("hashkmers=explicit"), "Support lacks exact counts");
		check(TadpoleMulti.hasMultipleK(new String[]{"fusepath=t"}), "Dispatcher missed path guard");
		for(String invalid : new String[]{"fusesupportk=124", "fusepathdepth=0", "fusepathflank=-1"}){
			boolean failed=false;
			try{new TadpoleMulti.Config(new String[]{"k=95,64", "out=test.fa", invalid});}
			catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid support option accepted: "+invalid);
		}
	}

	/** Creates one exact-strand fixture with safe default terminal metadata. */
	private static Contig contig(final int id, final String sequence, final boolean reverse){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		if(reverse){AminoAcid.reverseComplementBasesInPlace(bases);}
		final Contig c=new Contig(bases, id);
		c.coverage=5;
		c.leftCode=c.rightCode=Tadpole.DEAD_END;
		return c;
	}

	/** Generates reproducible nonperiodic sequence independently of the count implementation. */
	private static String random(final int length, final long seed){
		final Random random=new Random(seed);
		final StringBuilder out=new StringBuilder(length);
		for(int i=0; i<length; i++){out.append("ACGT".charAt(random.nextInt(4)));}
		return out.toString();
	}

	/** Test failures stay fatal even under accidental -da execution. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Small exact string oracle; allocation is confined to the test, not the production checker. */
	private static final class Evidence implements FusionKmerSupport.Evidence {
		/** Fixes the word length used to populate the independent oracle. */
		Evidence(final int k_){k=k_;}
		/** Adds every reference word at a known depth. */
		void add(final String sequence, final int depth){
			for(int i=0; i+k<=sequence.length(); i++){set(sequence.substring(i, i+k), depth);}
		}
		/** Sets one canonical word's count. */
		void set(final String word, final int depth){counts.put(canonical(word), depth);}
		@Override
		public int count(final Kmer key){return counts.get(canonical(key.toString()));}
		/** A fixed ratio policy makes weak-noise and meaningful-branch fixtures distinct. */
		@Override
		public boolean isJunction(final int max, final int second){return second>0 && second*4>=max;}
		/** Canonicalizes independently of the rolling-key implementation. */
		private static String canonical(final String word){
			final byte[] bases=word.getBytes(StandardCharsets.US_ASCII);
			AminoAcid.reverseComplementBasesInPlace(bases);
			final String reverse=new String(bases, StandardCharsets.US_ASCII);
			return word.compareTo(reverse)<0 ? word : reverse;
		}
		private final int k;
		private final ObjectIntMap<String> counts=new ObjectIntMap<String>(String.class);
	}
}
