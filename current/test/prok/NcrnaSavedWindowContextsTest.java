package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSavedWindowContexts;
import prok.NcrnaSavedWindowContexts.Window;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Native streamed source/context checks, independent of verifier alignment.
 * Adapted from A_SampleBasic2: fixed finite examples, fresh native outputs.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSavedWindowContextsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSavedWindowContextsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				require(parser.parse(arg, a, b), "Unknown context-source fixture parameter: "+arg);
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory without input, overwrite, append or sampling");
			final File root=new File(parser.out1);
			require(root.mkdir(), "The native context-source fixture must use a fresh directory");
			run(root.getAbsolutePath());
			pp.outstream.println("NCRNA_SAVED_WINDOW_CONTEXTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Actual Sources        ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final String root){
		final byte[] first=Prok16sPartialRescueTest.random(100, 101), second=Prok16sPartialRescueTest.random(100, 202);
		first[12]='N';
		second[85]='N';
		final String fasta=root+"/source.fa.gz";
		try(Output out=new Output(fasta)){
			out.writer.print(new ByteBuilder(">shared first contig\n").append(first)
				.append("\n>shared second contig\n").append(second).nl());
		}
		overlaps(fasta, first);
		sidecars(fasta, first);
		oneSource(fasta, first, second);
		mismatches(fasta, first);
	}

	/** New saved context must match the candidate and bind its exact interval;
	 * a genuine N is data, while zero is never allowed in a complete sidecar. */
	private static void sidecars(final String fasta, final byte[] source){
		for(int mode=0; mode<5; mode++){
			final NcrnaSavedWindowContexts contexts=new NcrnaSavedWindowContexts("sidecar", fasta+".must-not-open");
			final Window window=contexts.add("candidate", 0, "shared", 100, 0, 10, 29, 5, slice(source, 10, 29, 0));
			final byte[] bytes=slice(source, 5, 34, 0);
			if(mode==1){bytes[10]=different(bytes[10]);}
			if(mode==2){bytes[0]=0;}
			boolean failed=false;
			try{
				contexts.supply(window, mode==3 ? 6 : 5, mode==4 ? Arrays.copyOf(bytes, 29) : bytes);
				contexts.prepare();
				contexts.sourceFor(window).fill(window.bases, window.origin, 5, 34);
			}catch(IllegalArgumentException expected){failed=true;}
			require(failed==(mode!=0) && contexts.sourceLoads()==0,
				"Complete sidecars must bind original candidate bytes and geometry without source I/O; mode="+mode);
			if(mode==0){check(window, source);}
		}
	}

	/** Reverse-strand saved candidates must supply forward-strand context without
	 * opening an absent source file. This also distinguishes real Ns from holes. */
	private static void overlaps(final String fasta, final byte[] source){
		assert(source.length==100) : "The finite overlap oracle uses explicit100-base original coordinates";
		final NcrnaSavedWindowContexts contexts=new NcrnaSavedWindowContexts("frame", fasta+".must-not-open");
		final Window target=contexts.add("target", 0, "shared", 100, 0, 40, 59, 10, slice(source, 40, 59, 0));
		contexts.add("reverse_supplier", 0, "shared", 100, 1, 30, 69, 0, slice(source, 30, 69, 1));
		contexts.prepare();
		check(target, source);
		contexts.sourceFor(target).fill(target.bases, target.origin, 30, 69);
		require(contexts.sourceLoads()==0 && contexts.sourceBases()==0 && contexts.windows()==2,
			"A completely supplied saved overlap must require zero original-source loads");
	}

	/** Repeated native names do not erase contig ordinals. One streamed file load
	 * fills both strands/all retained contexts, then the second request reuses it. */
	private static void oneSource(final String fasta, final byte[] first, final byte[] second){
		final NcrnaSavedWindowContexts contexts=new NcrnaSavedWindowContexts("exact_frame", fasta);
		final Window a=contexts.add("first", 0, "shared", 100, 0, 10, 29, 5, slice(first, 10, 29, 0));
		final Window b=contexts.add("second", 1, "shared", 100, 1, 60, 79, 5, slice(second, 60, 79, 1));
		contexts.prepare();
		require(a.bases[0]==0 && b.bases[0]==0 && a.bases[12-a.origin]=='N',
			"Unavailable flank bytes must remain zero while genuine genomic Ns are already present");
		contexts.sourceFor(a).fill(a.bases, a.origin, a.origin, a.origin+a.bases.length-1);
		check(a, first);
		check(b, second);
		contexts.sourceFor(b).fill(b.bases, b.origin, b.origin, b.origin+b.bases.length-1);
		require(contexts.sourceLoads()==1 && contexts.sourceBases()==200,
			"Multiple missing-flank requests must stream this exact two-contig source once, without extra alignment or genome retention");
	}

	private static void mismatches(final String fasta, final byte[] first){
		for(int mode=0; mode<4; mode++){
			final NcrnaSavedWindowContexts contexts=new NcrnaSavedWindowContexts("mismatch", fasta);
			final byte[] captured=slice(first, 20, 39, 0);
			if(mode==0){captured[5]=different(captured[5]);}
			final Window window=contexts.add("event", mode==1 ? 2 : 0,
				mode==2 ? "different_name" : "shared", mode==3 ? 101 : 100, 0, 20, 39, 5, captured);
			contexts.prepare();
			boolean failed=false;
			try{contexts.sourceFor(window).fill(window.bases, window.origin, window.origin, window.origin+window.bases.length-1);}
			catch(IllegalArgumentException expected){failed=true;}
			require(failed && contexts.sourceLoads()==1,
				"Wrong source sequence, ordinal, name or length must fail rather than fill apparently compatible flanks: mode="+mode);
		}
		final NcrnaSavedWindowContexts overlap=new NcrnaSavedWindowContexts("overlap_mismatch", fasta+".must-not-open");
		overlap.add("left", 0, "shared", 100, 0, 20, 39, 5, slice(first, 20, 39, 0));
		final byte[] wrong=slice(first, 30, 49, 0);
		wrong[0]=different(wrong[0]);
		overlap.add("right", 0, "shared", 100, 0, 30, 49, 5, wrong);
		boolean failed=false;
		try{overlap.prepare();}catch(IllegalArgumentException expected){failed=true;}
		require(failed && overlap.sourceLoads()==0, "Disagreeing saved overlaps must fail before an original-file read or verification");
	}

	/** Independent original-orientation extraction used only by finite fixtures. */
	private static byte[] slice(final byte[] source, final int start, final int stop, final int strand){
		assert(start>=0 && stop>=start && stop<source.length && (strand==0 || strand==1))
			: "Finite slices use explicit original caller-strand bounds";
		final byte[] oriented=source.clone();
		if(strand==1){AminoAcid.reverseComplementBasesInPlace(oriented);}
		return Arrays.copyOfRange(oriented, start, stop+1);
	}

	private static void check(final Window window, final byte[] source){
		assert(window.bases.length>0) : "Every retained context must include its original nonempty candidate";
		require(Arrays.equals(window.bases, slice(source, window.origin, window.origin+window.bases.length-1, window.strand)),
			"Every retained context byte must match its exact original contig in the caller's strand orientation");
	}
	private static byte different(final byte base){return base=='A' ? (byte)'C' : (byte)'A';}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
}
