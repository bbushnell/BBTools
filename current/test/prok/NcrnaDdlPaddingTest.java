package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.AbstractKmerTable;
import map.LongHashSet;
import prok.BuildNcrnaDdlLayoutCounts;
import prok.BuildNcrnaDdlPadding;
import prok.GeneCaller;
import prok.GeneModel;
import prok.NcrnaDdlOracle;
import prok.NcrnaDdlWindowAnalysis;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.NcrnaDdlHoldout;
import prok.NcrnaTextIO;
import prok.Orf;
import stream.Read;

/** Literal padding counts plus a real two-model accepted-envelope fixture.
 * The oracle/verifier are the existing production components, with no ranking
 * profile and no substitute annotation flank.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlPaddingTest{

	public static void run(final String root){
		final byte[] window=bytes("AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCC");
		check(window.length==50, "The literal window consists of16prefix,20occupied and14suffix bases");
		final BuildNcrnaDdlPadding.Counts forward=new BuildNcrnaDdlPadding.Counts();
		BuildNcrnaDdlPadding.addPadding(forward, window, 0, 49, 0, 17, 36, new int[2]);
		final int maskedC=(int)(0x55555&AbstractKmerTable.makeMiddleMask(10, false, 1));
		check(forward.shortCounts[0]==9 && forward.shortCounts[0x5555]==7
			&& forward.longCounts[0]==7 && forward.longCounts[maskedC]==5
			&& sum(forward.shortCounts)==16 && sum(forward.longCounts)==12,
			"Count occurrences only outside the complete accepted envelope, separately by original word span");
		final BuildNcrnaDdlPadding.Counts reverse=new BuildNcrnaDdlPadding.Counts();
		BuildNcrnaDdlPadding.addPadding(reverse, window, 100, 149, 1, 115, 134, new int[2]);
		check(Arrays.equals(forward.shortCounts, reverse.shortCounts) && Arrays.equals(forward.longCounts, reverse.longCounts),
			"One-based forward accepted spans must select the same oriented flanks on the reverse strand");
		final BuildNcrnaDdlPadding.Counts separated=new BuildNcrnaDdlPadding.Counts();
		final byte[] noBridge=bytes("AAAAAAAGGGGGGGGGGAAAAAAA");
		BuildNcrnaDdlPadding.addPadding(separated, noBridge, 0, 23, 0, 8, 17, new int[2]);
		check(sum(separated.shortCounts)==0 && sum(separated.longCounts)==0,
			"Two seven-base flanks may not be concatenated into invented8/10mers across the occupied interior");
		separated.flank(bytes("AAAAANAAAA"), 0, 10);
		check(sum(separated.shortCounts)==0 && sum(separated.longCounts)==0,
			"An undefined base rejects its entire original span even when the middle mask could hide it");
		final long[] ordered=BuildNcrnaDdlLayoutCounts.orderedFrequencies(new int[]{0, 5, 5, 3});
		check(BuildNcrnaDdlLayoutCounts.rawWord(ordered, 0)==1
			&& BuildNcrnaDdlLayoutCounts.rawWord(ordered, 1)==2
			&& BuildNcrnaDdlLayoutCounts.rawWord(ordered, 2)==3,
			"The existing blacklist order is frequency descending and raw-word ascending on ties");
		oracle();
		reuse(root);
	}

	private static void oracle(){
		final byte[] query=Prok16sPartialRescueTest.random(320, 20261120);
		final byte[][] library={Arrays.copyOfRange(query, 40, 200), Arrays.copyOfRange(query, 60, 240)};
		final String[] names={"cobalamin_RF00174_left", "cobalamin_RF00174_right"};
		final LongHashSet seeds=new LongHashSet(4);
		long key=0;
		for(int base=100; base<117; base++){key=(key<<2)|dna.AminoAcid.baseToNumber[query[base]];}
		seeds.add(key);
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", library, null, names, seeds,
			17, 109, query.length, 9, 384, true, 2, 0, 0, 0, 0, 20, .99f, .99f);
		family.modelPatience=8;
		// The320-base window exceeds the ordinary120-base Quantum threshold.
		family.trimAlignmentExtent=false;
		family.scavengePass2=false;
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family);
			final GeneCaller caller=new GeneCaller(60, 10, 10, -1, -1, -1, -1, -1, new GeneModel(false));
			final NcrnaScavenger ordinary=caller.ncrnaScavenger(0);
			final int[] observed={0};
			final BuildNcrnaDdlPadding.Counts counts=new BuildNcrnaDdlPadding.Counts();
			final NcrnaDdlWindowAnalysis analysis=new NcrnaDdlWindowAnalysis("fixture", family, ordinary,
				new NcrnaDdlWindowAnalysis.Sink(){
					@Override
					public void window(final NcrnaDdlWindowAnalysis view, final byte[] window){
						final NcrnaDdlOracle oracle=view.oracle();
						check(oracle.start(0)==41 && oracle.stop(0)==200 && oracle.start(1)==61 && oracle.stop(1)==240,
							"The real verifier must retain both exact models' distinct accepted spans");
						final int[] bounds=new int[2];
						BuildNcrnaDdlPadding.addPadding(counts, window, 0, 319, 0,
							Math.min(oracle.start(0), oracle.start(1)), Math.max(oracle.stop(0), oracle.stop(1)), bounds);
						check(bounds[0]==40 && bounds[1]==240,
							"Padding excludes ALL accepted spans, not merely the ordinary winning model");
						observed[0]++;
					}
					@Override
					public void attempt(final NcrnaDdlWindowAnalysis view, final int rank, final int model,
							final int shared, final float identity, final int span, final boolean passed){}
					@Override
					public void endWindow(final NcrnaDdlWindowAnalysis view){}
				});
			final Read read=new Read(query.clone(), null, "padding_fixture", 0);
			analysis.begin(read);
			try{
				final ArrayList<Orf> calls=ordinary.scavenge(read.id, read.bases, 0, new ArrayList<int[]>());
				check(calls.size()==1 && names[1].equals(calls.get(0).trnaModel), "The longer exact model ranks first by native9mer occurrences; the complete oracle envelope also includes the left model");
			}finally{analysis.end();}
			analysis.finish();
			check(observed[0]==1 && counts.shortSegments==2 && counts.longSegments==2,
				"One real full candidate window must yield exactly its two accepted-envelope padding flanks");
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(saved);
		}
	}

	/** Constructor-only failures precede native caller work and frequency import. */
	private static void reuse(final String root){
		check(new java.io.File(root).mkdir(), "Reuse fixtures require a fresh directory");
		final String header=NcrnaDdlHoldout.PANEL_HEADER+"\n";
		final String pilot="0\tA\ta\t/a.fa\t/a.gff\n";
		final String rest="1\tB\tb\t/b.fa\t/b.gff\n";
		final String sample=root+"/sample.tsv", binding=root+"/binding.tsv", caller=root+"/caller.config";
		final String previous=root+"/pilot";
		write(sample, header+pilot+rest);
		write(binding, "source\tfixture\nbuild\tfixture\n");
		write(caller, "taxonomy=f\n");
		write(previous+".sample.tsv", header+pilot+rest);
		write(previous+".binding.tsv", "source\tfixture\nbuild\tfixture\n");
		write(previous+".inputs.tsv", BuildNcrnaDdlPadding.INPUT_HEADER+"\nA\ta\t/a.fa\t1\t100\t0\t0\t0\t0\t0\n");
		write(previous+".frequencies.tsv.gz", BuildNcrnaDdlPadding.FREQUENCY_HEADER+"\nSHORT8\t0\t1\n");
		construct(caller, sample, binding, previous, root+"/valid", 1, true);
		reject(caller, sample, binding, previous, root+"/overlap", 0, false, "disjoint");
		reject(caller, sample, binding, null, root+"/omitted", 0, true, "every row");
		final String changed=root+"/changed.tsv";
		write(changed, header+pilot+"1\tC\tc\t/c.fa\t/c.gff\n");
		reject(caller, changed, binding, previous, root+"/changed", 1, true, "COMPLETE");
		final String otherBinding=root+"/other.binding.tsv";
		write(otherBinding, "source\tother\nbuild\tfixture\n");
		reject(caller, sample, otherBinding, previous, root+"/different_binding_output", 1, true, "worker-verified");
	}

	private static void construct(final String caller, final String sample, final String binding,
			final String reuse, final String out, final int first, final boolean makeList){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("caller="+caller, "sample="+sample,
			"binding="+binding, "out="+out, "first="+first, "count=1", "makelist="+makeList));
		if(reuse!=null){args.add("reuse="+reuse);}
		new BuildNcrnaDdlPadding(args.toArray(new String[args.size()]), System.err);
	}

	private static void reject(final String caller, final String sample, final String binding,
			final String reuse, final String out, final int first, final boolean makeList, final String reason){
		try{construct(caller, sample, binding, reuse, out, first, makeList);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Reuse must fail for its declared contract: "+reason+"; got "+expected.getMessage());
			return;
		}
		throw new AssertionError("Reuse accepted a broken preparation contract: "+reason);
	}

	private static void write(final String path, final String text){
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(path)){out.writer.print(text);}
	}

	private static long sum(final int[] counts){
		long result=0;
		for(final int count:counts){result+=count;}
		return result;
	}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
}
