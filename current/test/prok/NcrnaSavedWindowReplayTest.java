package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlOrf;
import prok.NcrnaAlignmentCache;
import prok.NcrnaModelAttemptInstrumentSink;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaDdlWindowAnalysis;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;
import shared.Shared;

/** Real ordinary-window versus stored-context verification, including genuine
 * physical edges and interior fragments. Adapted from A_SampleBasic2 and the
 * selected unlimited/physical-edge fixtures. No timing or biological claim.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSavedWindowReplayTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSavedWindowReplayTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				require(parser.parse(arg, a, b), "Unknown saved-window fixture parameter: "+arg);
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory without input, overwrite, append or sampling");
			final File root=new File(parser.out1);
			require(root.mkdir(), "The finite saved-window fixture must use a new output directory");
			shortRnaConstruction();
			for(int mode=0; mode<=4; mode++){
				for(int strand=0; strand<2; strand++){compare(root.getPath(), mode, strand, true);}
			}
			compare(root.getPath(), 4, 0, false);
			pp.outstream.println("NCRNA_SAVED_WINDOW_REPLAY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Actual Verifier        ----------------*/
	/*--------------------------------------------------------------*/

	private static void compare(final String root, final int mode, final int strand, final boolean trim){
		assert(mode>=0 && mode<=4 && (strand==0 || strand==1))
			: "Fixtures cover full, left, right, both physical ends and an interior fragment on each caller strand";
		final byte[] model=Prok16sPartialRescueTest.random(480, 20261115);
		final int from=mode==3 ? 150 : mode==1 || mode==4 ? 240 : 0;
		final int to=mode==3 ? 330 : mode==2 ? 240 : model.length;
		final byte[] observed=Arrays.copyOfRange(model, from, to);
		final int left=mode==1 || mode==3 ? 0 : 300;
		final int right=mode==2 || mode==3 ? 0 : 300;
		final byte[] bases=Prok16sPartialRescueTest.embed(observed, left, right);
		final int first=left+observed.length/2-10, second=first+20;
		final LongHashSet seeds=new LongHashSet(4);
		seeds.add(key(bases, first));
		seeds.add(key(bases, second));
		require(seeds.size()==2, "Two actual distinct seed keys must support every candidate fixture");
		final byte[][] models={model};
		final String[] names={"cobalamin_RF00174_context"};
		final String file=root+"/profiles_"+mode+"_"+strand+"_"+trim+".tsv";
		NcrnaDdlProfileResource.writeSelected(file, new char[1][1024], new int[]{1}, names, models, 1024, new LongHashSet(4));
		final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(file, names, models);
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", models, null, names, seeds,
			17, 109, 300, 9, Integer.MAX_VALUE, true, 2f, 0f, 0f, 0, 0f, 1f, .99f, .99f);
		family.modelPatience=Integer.MAX_VALUE;
		family.seedMinHits=2;
		family.seedDistinct=true;
		family.physicalEdgeMinLen=109;
		family.modelClipRescueId=.85f;
		family.trimAlignmentExtent=trim;
		family.scavengePass2=false;
		final NcrnaScavenger ordinary=NcrnaDdlWindowAnalysis.savedWindowWorker(family, profiles, .99f);
		final int[] window={-1,-1,-1,-1};
		ordinary.setWindowCaptureSink((name, sense, pass, start, stop, source, contig, sequence)->{
			require(window[0]<0 && contig==bases && sense==strand && pass==1,
				"One ordinary candidate must supply the replay's actual original geometry");
			window[0]=start;
			window[1]=stop;
			window[2]=source;
			window[3]=pass;
		});
		final ArrayList<Orf> calls=ordinary.scavenge("context", bases, strand, new ArrayList<int[]>(), new int[]{first,second});
		require(window[0]>=0 && calls.size()==(mode==4 ? 0 : 1),
			"The fixture must really exercise one ordinary accepted full/physical locus or rejected internal fragment: mode="+mode);
		final NcrnaScavenger replay=NcrnaDdlWindowAnalysis.savedWindowWorker(family, profiles, .99f);
		final int flank=replay.savedWindowFlank();
		final int contextStart=Math.max(0, window[0]-flank), contextStop=Math.min(bases.length-1, window[1]+flank);
		final byte[] context=Arrays.copyOfRange(bases, contextStart, contextStop+1);
		final NcrnaDdlOrf result=replay.verifySavedWindow("context", context, contextStart, bases.length,
			strand, window[3], window[0], window[1], window[2]);
		require((result==null)==calls.isEmpty() && replay.alignmentCount()==ordinary.alignmentCount(),
			"Stored real context must preserve actual acceptance and all primary/edge/refinement alignment work");
		if(result!=null){
			final NcrnaDdlOrf expected=(NcrnaDdlOrf)calls.get(0);
			require(result.start+contextStart==expected.start && result.stop+contextStart==expected.stop
				&& result.trnaModel.equals(expected.trnaModel) && result.acceptedIdentity()==expected.acceptedIdentity()
				&& result.ncrnaPartial==expected.ncrnaPartial && result.strand==expected.strand,
				"Context coordinates, original model, identity, physical clipping and strand must reproduce the ordinary proposal");
			if(mode>0){require(result.ncrnaPartial!=0, "The physical-edge fixture must reach clipping, not pass the primary full-model path");}
		}
		if(contextStart>0 && flank>0){
			boolean failed=false;
			try{
				replay.verifySavedWindow("missing", Arrays.copyOfRange(context, 1, context.length), contextStart+1,
					bases.length, strand, window[3], window[0], window[1], window[2]);
			}catch(IllegalArgumentException expected){failed=true;}
			require(failed, "Missing real refinement context must fail before alignment instead of clamping to a false contig edge");
		}
		if(mode==0 && strand==0 && trim){fillOnDemand(family, profiles, bases, left, left+model.length-1, window[2]);}
		if(mode<=3 && strand==0 && trim){checkDiagnostic(family, profiles, bases, window, result, contextStart);}
		if(mode==0 && strand==0 && trim){tinySpanOffEdges(family, profiles, window[2]);}
	}

	/** Literal codon values prove that the RNA-only boundary guard does not alter
	 * available metadata or relax the existing coding-feature construction path. */
	private static void shortRnaConstruction(){
		final byte[] bases={'A','C','G','T','A','C'};
		for(int strand=0; strand<2; strand++){
			for(int length=1; length<=2; length++){
				final Orf left=new NcrnaDdlOrf("left", 0, length-1, strand, bases, ProkObject.RNA);
				final Orf right=new NcrnaDdlOrf("right", bases.length-length, bases.length-1, strand, bases, ProkObject.RNA);
				require(left.start==0 && left.stop==length-1 && left.startCodon==6 && left.stopCodon==-1,
					"Tiny left-edge RNA keeps its real coordinates and available ACG metadata; no stop codon is invented");
				require(right.start==bases.length-length && right.stop==bases.length-1 && right.startCodon==-1 && right.stopCodon==49,
					"Tiny right-edge RNA keeps its real coordinates and available TAC metadata; no start codon is invented");
				final byte[] entire=Arrays.copyOf(bases, length);
				final Orf whole=new NcrnaDdlOrf("whole", 0, length-1, strand, entire, ProkObject.RNA);
				require(whole.startCodon==-1 && whole.stopCodon==-1 && whole.length()==length,
					"A one/two-base contig remains a nonempty RNA candidate without reading outside its real bases");
				boolean rejected=false;
				try{new Orf("coding", 0, length-1, strand, 0, entire, false, ProkObject.CDS);}
				catch(IndexOutOfBoundsException expected){rejected=true;}
				require(rejected, "The RNA metadata guard must not convert a truncated CDS codon into an accepted coding feature");
			}
			final Orf coding=new Orf("coding", 0, 5, strand, 0, bases, false, ProkObject.CDS);
			require(coding.startCodon==6 && coding.stopCodon==49, "Ordinary CDS codon extraction must retain the literal ACG/TAC values");
		}
		boolean rejected=false;
		try{new NcrnaDdlOrf("invalid", -1, 1, 0, bases, ProkObject.RNA);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "RNA codon metadata must not hide invalid candidate coordinates");
	}

	/** The real span-off edge verifier must retain tiny accepted spans rather
	 * than silently reintroducing a three-base or ordinary109-base floor. */
	private static void tinySpanOffEdges(final NcrnaFamily family, final NcrnaDdlProfileResource profiles, final int source){
		final NcrnaModelAttemptInstrumentSink observer=(name, strand, pass, start, stop, rank, model, shared, id, span, accepted)->{};
		final byte[] model=family.library[0];
		for(int length=1; length<=2; length++){
			for(int edge=0; edge<2; edge++){
				// The real failure had a tiny accepted span in a normal window.
				// Entire tiny references instead reach QuantumAligner:204's known
				// long-query prefill limitation before the RNA constructor.
				final byte[] bases=new byte[model.length+64];
				Arrays.fill(bases, (byte)'N');
				final int start=edge==0 ? 0 : bases.length-length;
				System.arraycopy(model, edge==0 ? model.length-length : 0, bases, start, length);
				final NcrnaScavenger diagnostic=NcrnaDdlWindowAnalysis.shuffledWindowWorker(family, profiles, .60f, false);
				final NcrnaDdlOrf actual=diagnostic.verifyShuffledWindow("tiny", bases, 0, bases.length, 0, 1,
					0, bases.length-1, source, null, observer, null);
				require(actual!=null && actual.start==start && actual.stop==start+length-1 && actual.ncrnaPartial!=0,
					"The actual no-span verifier must preserve the complete one/two-base physical-edge match and its clipping flags");
			}
		}
	}

	/** A deliberately empty discovery set rejects ordinary replay before alignment.
	 * The explicit diagnostic must still reproduce real accepted/refined proposals
	 * and physical-edge flags, with and without exact primary-kernel reuse.
	 */
	private static void checkDiagnostic(final NcrnaFamily family, final NcrnaDdlProfileResource profiles,
			final byte[] bases, final int[] window, final NcrnaDdlOrf expected, final int expectedOrigin){
		family.kmerSet.clear();
		final NcrnaScavenger ordinary=NcrnaDdlWindowAnalysis.savedWindowWorker(family, profiles, .99f);
		require(ordinary.verifySavedWindow("seedless", bases, 0, bases.length, 0, 1, window[0], window[1], window[2])==null
			&& ordinary.alignmentCount()==0, "The seed-bypass fixture must actually reject before ordinary alignment");
		final NcrnaAlignmentCache cache=new NcrnaAlignmentCache(family.library, 1, bases.length);
		final long[] attempts={0};
		final NcrnaModelAttemptInstrumentSink observer=(name, strand, pass, start, stop, rank, model, shared, id, span, accepted)->attempts[0]++;
		for(int mode=0; mode<3; mode++){
			final NcrnaScavenger diagnostic=NcrnaDdlWindowAnalysis.shuffledWindowWorker(family, profiles, .99f, true);
			final NcrnaDdlOrf actual=diagnostic.verifyShuffledWindow("seedless", bases, 0, bases.length, 0, 1,
				window[0], window[1], window[2], null, observer, mode==0 ? null : cache);
			require(expected!=null && actual!=null && actual.start==expected.start+expectedOrigin
				&& actual.stop==expected.stop+expectedOrigin && actual.acceptedIdentity()==expected.acceptedIdentity()
				&& actual.ncrnaPartial==expected.ncrnaPartial,
				"Discovery-seed bypass and cached primary kernels must preserve actual refinement and physical-edge acceptance");
		}
		require(attempts[0]>0 && cache.fills()>0 && cache.hits()>0,
			"The diagnostic fixture must execute real primary kernels and then really reuse them");
	}

	/** An exact candidate forces native trimming to read both real flanks. The
	 * on-demand source must supply them inside the same single verification walk. */
	private static void fillOnDemand(final NcrnaFamily family, final NcrnaDdlProfileResource profiles,
			final byte[] bases, final int start, final int stop, final int source){
		assert(start>=10 && stop+10<bases.length) : "The sparse-context fixture needs real interior flanks on both sides";
		final NcrnaScavenger full=NcrnaDdlWindowAnalysis.savedWindowWorker(family, profiles, .99f);
		final NcrnaDdlOrf expected=full.verifySavedWindow("source", bases, 0, bases.length, 0, 1, start, stop, source);
		final int origin=start-10;
		final byte[] context=new byte[stop-start+21];
		System.arraycopy(bases, start, context, 10, stop-start+1);
		final int[] requests={0};
		final NcrnaScavenger sparse=NcrnaDdlWindowAnalysis.savedWindowWorker(family, profiles, .99f);
		final NcrnaDdlOrf actual=sparse.verifySavedWindow("source", context, origin, bases.length,
			0, 1, start, stop, source, (stored, offset, from, to)->{
				require(stored==context && offset==origin && from==start-10 && to==stop+10,
					"The context source must receive the exact native refinement read, not a whole-contig replacement");
				requests[0]++;
				for(int i=from; i<=to; i++){
					require(stored[i-offset]==0 || stored[i-offset]==bases[i], "Already captured bytes must agree with the exact source");
					stored[i-offset]=bases[i];
				}
			});
		require(requests[0]==1 && expected!=null && actual!=null && actual.start+origin==expected.start
			&& actual.stop+origin==expected.stop && actual.acceptedIdentity()==expected.acceptedIdentity()
			&& sparse.alignmentCount()==full.alignmentCount(),
			"Lazy real-flank supply must preserve acceptance, original coordinates and exact single-walk alignment work");
	}

	private static long key(final byte[] bases, final int start){
		assert(start>=0 && start+17<=bases.length) : "The fixture seed must be a complete real forward17mer";
		long key=0;
		for(int i=start; i<start+17; i++){
			final int x=dna.AminoAcid.baseToNumber[bases[i]];
			require(x>=0, "Finite seed examples require unambiguous real bases");
			key=(key<<2)|x;
		}
		return key;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
}
